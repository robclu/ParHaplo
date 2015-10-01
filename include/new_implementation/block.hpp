// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for a block
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_BLOCK_HPP
#define PARAHAPLO_BLOCK_HPP

#include "operations.hpp"
#include "small_containers.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>
#include <tbb/tbb.h>
#include <tbb/concurrent_unordered_map.h>

#include <iostream>
#include <string>
#include <vector>

// NOTE: All output the the terminal is for debugging and checking at the moment

// TODO : Create unsplittable info and then determine unsplittable blocks

namespace haplo {

// Define some HEX values for the 8 bits comparisons
#define ZERO 0x00
#define ONE  0x01
#define TWO  0x02
#define IH   0x00           // Intricically heterozygous
#define NIH  0x01           // Not intrinsically heterozygous
    
namespace io = boost::iostreams;
using namespace io;

// ----------------------------------------------------------------------------------------------------------
/// @class      Block 
/// @brief      Represents a block of input the for which the haplotypes must be determined
/// @tparam     R       The number of rows in the block
/// @tparam     C       The number of columns in the block
/// @tparam     THI     The number of threads to use for dimension I -- default to 1
/// @tparam     THJ     The number of threads to use for dimension J -- default to 1
// ----------------------------------------------------------------------------------------------------------
template <size_t R, size_t C, size_t THI = 1, size_t THJ = 1>
class Block {
public:
    // ----------------------------------------- TYPES ALIAS'S ----------------------------------------------
    using data_container        = BinaryArray<R * C, 2>;    // R*C Elements, 2 bits per element
    using binary_container_r    = BinaryArray<C>;           // A container of bits C elements long
    using binary_container_c    = BinaryArray<R>;           // A container of bits R elements long
    using atomic_array_r        = std::array<tbb::atomic<int>, R>;
    using atomic_array_c        = std::array<tbb::atomic<int>, C>;
    using atomic_array_2r       = std::array<tbb::atomic<int>, R * 2>;
    using concurrent_umap       = tbb::concurrent_unordered_map<int, byte>;
    
    // ------------------------------------------------------------------------------------------------------
private:
    data_container      _data;                  //!< Container for { '0' | '1' | '-' } data variables
    binary_container_c  _haplo_one;             //!< The binary bits which represent the first haplotype
    binary_container_c  _haplo_two;             //!< The bianry bits which represent the second haplotype
    binary_container_r  _alignment;             //!< The binary bits which represent the alignment of the 
                                                //!< reads to the haplotypes 
    atomic_array_r      _row_mplicities;        //!< The multiplicites of the rows
    atomic_array_c      _col_mplicities;        //!< The multiplicites of the columns
    binary_container_r  _singletons;            //!< If each of the rows is singleton or not
    binary_container_c  _column_types;          //!< The type of a column, IH or NIH
    binary_container_c  _splittable_columns;    //!< If a column is splittable or not
    atomic_array_2r     _row_info;              //!< Start and end positions for the row
    concurrent_umap     _monotones;             //!< Columns that are monotone and do not need to be searched
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to fill the block with data from the input file
    /// @param[in]  data_file       The file to fill the data with
    // ------------------------------------------------------------------------------------------------------
    Block(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns an element from the block at a specified position
    /// @param[in]  row_idx     The index of the row for the element
    /// @param[in]  col_idx     The index of the column for the element
    // ------------------------------------------------------------------------------------------------------
    byte operator()(size_t row_idx, size_t col_idx) const { return _data.get(row_idx * C + col_idx); }
   
    // --------------------- TEMP PRINTING FUNCTIONS ----------------------- 
    void print() const;
    
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the types of each of the columns - if the column is intrisically heterozygous,
    ///             or if it is not intrinsically heterozygous
    // ------------------------------------------------------------------------------------------------------
    void determine_column_types();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the mutplicities of the rows and the columns
    // ------------------------------------------------------------------------------------------------------
    void determine_multiplicities();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data from the input file
    /// @param[in]  data_file   The file to get the data from 
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* data_file);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Finds all duplicates of either a row or a column
    /// @param[in]  checker     The checker function -- for rows or columns
    /// @tparam     CheckerType The type of the checker -- rows or columns 
    // ------------------------------------------------------------------------------------------------------
    template <typename CheckerType>
    void find_duplicates(CheckerType checker);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the start and end indices of each of the rows, if a row is singular, and how
    ///             many ones and zeros are in each column
    // ------------------------------------------------------------------------------------------------------
    void find_params(); 
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_data(size_t row, TP token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes a row, determining all the necessary parameters
    /// @param[in]  row_idx     The row in the data matrix to fill
    // ------------------------------------------------------------------------------------------------------
    void process_row(const size_t row_idx);
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

// ----------------------------------------------- TESTING --------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::print() const 
{
    std::cout << "\n\n|----------STARTING TO PRINT------------|\n";
    std::cout << "|---------------------------------------|";
    
    std::cout << "\n\n|--------------DATA------------|\n\n";
    
    for (size_t r = 0; r < R; ++r) {
        for (size_t c = 0; c < C; ++c) {
            if (_data.get(r * C +c) != 2 ) {
                std::cout << static_cast<unsigned>(_data.get(r * C + c)) << " ";
            } else {
                std::cout << "- ";
            }
        }
        std::cout << "\n";
    }
    
    std::cout << "\n|-----------SINGLETONS----------|";
    
    for (size_t i = 0; i < R; i++) 
        std::cout << static_cast<unsigned>(_singletons.get(i)) << "\n";
    

    std::cout << "\n|------------ROW INFO-----------|";
    
    for (size_t i = 0; i < R; ++i) 
        std::cout << i << " : " << _row_info[2 * i] << " : " << _row_info[2 * i + 1] << "\n";
    

    std::cout << "\n|------------COL TYPES-----------|";
    
    for (auto i = 0; i < C; ++i) 
        std::cout << static_cast<unsigned>(_column_types.get(i)) << " ";
    std::cout << "\n";
    
    std::cout << "\n|------------MONOTONES-----------|";
    
    for (auto i = 0; i < C; ++i) 
        if (_monotones.find(i) != _monotones.end()) std::cout << i << " ";
    std::cout << "\n";

    std::cout << "\n|------------SPLITTABLE-----------|";
    
    for (auto i = 0; i < C; ++i) 
        if (_splittable_columns.get(i) == 1) std::cout << i << " ";
    std::cout << "\n";

    std::cout << "\n|------------HAPLOTYPES-----------|";
    
    std::cout << "h   : ";
    for (auto i = 0; i < C; ++i) std::cout << static_cast<unsigned>(_haplo_one.get(i));
    std::cout << "\nh`  : ";
    for (auto i = 0; i < C; ++i) std::cout << static_cast<unsigned>(_haplo_two.get(i));
    std::cout << "\n|---------------------------------------|\n\n";
        
}


// ------------------------------------------------- PUBLIC -------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
Block<R, C, THI, THJ>::Block(const char* data_file)
: _row_info{0}
{
    fill(data_file);
    find_params();

} 

// ------------------------------------------------- PRIVATE ------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::determine_column_types()
{
    // We can use all the available cores for this 
    constexpr size_t threads_x = (THJ + THI) < C ? (THJ + THI) : C;

    // Over each column in the row
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x)
        {
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, C, threads_x);
     
                // For each column we must traverse downwards and check how many non
                // singular rows there are, and then compare that to the min of the 
                // number of zeros and ones
                for (size_t it_x = 0; it_x < thread_iters_x; ++it_x) {
                    size_t  col_idx             = it_x * threads_x + thread_idx;
                    int     num_elements[2]     = {0, 0};
                    int     num_not_single      = 0;
                    bool    splittable          = true;
                    
                    // For each of the elements in the column (i.e row indices)
                    for (size_t row_idx = 0; row_idx < R; ++row_idx) {
                        auto element_value = _data.get(row_idx * C + col_idx);
                        // Check the element value and update the appropriate
                        // column information container, set start and end values
                        if (element_value == 0) {
                            ++num_elements[0];
                            if (_singletons.get(row_idx) == ZERO) ++num_not_single;
                        } else if (element_value == 1) {
                            ++num_elements[1];
                            if (_singletons.get(row_idx) == ZERO) ++num_not_single;
                        }
                        
                        // Check for splittable condition
                        if (_row_info[row_idx * 2] < col_idx && _row_info[row_idx * 2 + 1] > col_idx)
                            splittable = false;
                    }
                    
                    // Check if the column in monotone, and if the are NIH
                    // if they are none of these, then they are the default (IH)
                    if (num_elements[0] != 0 && num_elements[1] == 0) {
                        _monotones[col_idx] = 0;
                        
                        // Set the values of the haplotypes at these positions
                        _haplo_one.set(col_idx, 0);
                        _haplo_two.set(col_idx, 0);
                        
                        // These columns also fit the IH category
                        _column_types.set(col_idx, NIH);
                    } else if (num_elements[1] != 0 && num_elements[0] == 0) {
                        _monotones[col_idx] = 0;
                        
                        // Set the value of the haplotypes
                        _haplo_one.set(col_idx, 1);
                        _haplo_two.set(col_idx, 1);
                        
                        _column_types.set(col_idx, NIH);
                    } else if (!(std::min(num_elements[0], num_elements[1]) >= (num_not_single / 2))) {
                        _column_types.set(col_idx, NIH);      
                    }
                    
                    // If the column is splittable, then add it to the splittable list
                    if (splittable) _splittable_columns.set(col_idx, 1);
                }
            }
        }
    );
}


template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::fill(const char* data_file) 
{
    // Open file and convert to string (for tokenizer)
    io::mapped_file_source  file(data_file);
    std::string data(file.data(), file.size());

    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> nwline_separator{"\n"};
    
    // Tokenize the data into lines
    tokenizer lines{data, nwline_separator};

    // Counter for the row id 
    size_t row = 0;
   
    // Get the data and store it in the data container 
    for (auto line = lines.begin(); line != lines.end(); ++line) 
        process_data(row++, line);
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::find_params()
{
    // Check that we aren't trying to use more threads than rows or columns
    constexpr size_t threads_y = THI < R ? THI : R;
    
    // Over each of the rows
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y), 
        [&](const tbb::blocked_range<size_t>& thread_ids_y)
        {
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy, R, threads_y);
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t  row_idx = it_y * threads_y + thread_idy;
                    
                    // Prcocess the row, determining all the necessary parameters
                    process_row(row_idx);
                }
            }
        }
    );
    
    // Now we can find which of the columns are IH and which are NIH
    determine_column_types();
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::process_row(const size_t row_idx)
{
    // Determine the number of threads (this function is called by parallel
    // rows so the number of threads each instance can use is limited)
    constexpr size_t threads_x = (THJ / THI) < C ? (THJ / THI) : C;

    // Over each column in the row
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x)
        {
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, C, threads_x);
                int start_idx = -1, end_idx = -1;
     
                // Go through the row elements, checking if each of the element is 
                // a 0 or a 1 so that we can later detrmine if a col is IH or NIH
                for (size_t it_x = 0; it_x < thread_iters_x; ++it_x) {
                    auto col_idx       = it_x * threads_x + thread_idx;
                    auto element_value = _data.get(row_idx * C + col_idx);
                    
                    // If the element has a valid value then we 
                    // can set the start and end index appropriately
                    if (element_value <= 1) {
                        if (start_idx == -1) 
                            start_idx = col_idx;
                        end_idx = col_idx;
                    }
                }
                // Set the row paramters - the start and the end of the row
                // as well as if the row is a singleton or not (has 1 element)
                if (start_idx != end_idx) {
                    _row_info[2 * row_idx]     = start_idx;
                    _row_info[2 * row_idx + 1] = end_idx;
                } else {
                    _singletons.set(row_idx, 1);        // Row is singleton
                }
            }
        }
    );
}

template <size_t R, size_t C, size_t THI, size_t THJ> template <typename TP>
void Block<R, C, THI, THJ>::process_data(size_t row, TP token_pointer) 
{
    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> wspace_separator{" "}; 

    // Create a string from the line token and tokenize it
    std::string line(*token_pointer);
    tokenizer   elements{line, wspace_separator};

    size_t column       = 0;
    size_t row_offset   = row * C;
    
    for (auto& e : elements) {
        // Tokenizer creates a string, but because of the way we tokenized it
        // we know that it only has 1 element, so convert to char
        switch (e[0]) {
            case '0':
                _data.set(row_offset + column, ZERO);
                break;
            case '1':
                _data.set(row_offset + column, ONE);
                break;
            case '-':
                _data.set(row_offset + column, TWO);
                break;
            default:
                std::cerr << "Error reading input data - exiting =(\n";
                exit(1);
        } ++column;
    }
}

 
}           // End namespace haplo
#endif      // PARAHAPLO_BLOCK_HPP
