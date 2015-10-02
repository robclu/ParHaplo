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
#define ZERO    0x00
#define ONE     0x01
#define TWO     0x02
#define THREE   0x03
#define IH      0x00           // Intricically heterozygous
#define NIH     0x01           // Not intrinsically heterozygous
    
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
    using atomic_array_r        = std::array<tbb::atomic<uint>, R>;
    using atomic_array_c        = std::array<tbb::atomic<uint>, C>;
    using atomic_array_2r       = std::array<tbb::atomic<uint>, R * 2>;
    using concurrent_umap       = tbb::concurrent_unordered_map<uint, byte>;
    
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
    atomic_array_2r     _row_info;              //!< Start and end positions for the row
    std::vector<uint>   _splittable_columns;    //!< The indices of the splittable columns
    concurrent_umap     _monotone_columns;      //!< Columns that are monotone and do not need to be searched
    concurrent_umap     _flipped_columns;       //!< Columns that are monotone and do not need to be searched
    
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
    inline byte operator()(size_t row_idx, size_t col_idx) const { return _data.get(row_idx * C + col_idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns a constant reference to the data
    /// @return     A constant reference to the data
    // ------------------------------------------------------------------------------------------------------
    inline const data_container& data() const { return _data; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      If a column is monotone 
    /// @param[in]  col_idx     The index of the column 
    /// @return     If the column at col_idx is monotone or not
    // ------------------------------------------------------------------------------------------------------
    inline bool is_monotone(const uint col_idx) const 
    { 
        return _monotone_columns.find(col_idx) != _monotone_columns.end(); 
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the i'th unsplittable column in the block. For example the 3rd
    ///             unsplittable column may be the 7th column in the block, so this is essentially a mapping
    ///             interface
    /// @param[in]  idx     The index of the unsplittable column that is wanted
    /// @return     The index of the ith unsplittable column's index in the whole block
    // ------------------------------------------------------------------------------------------------------
    inline uint unsplittable_column(const uint idx) const 
    { 
        // Returning 0 if an out of range error occurs (for now)
        return idx < _splittable_columns.size() ? _splittable_columns[idx] : 0; 
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of unsplittable blocks possible from this block
    /// @return     The number of unsplittable blocks from this block
    // ------------------------------------------------------------------------------------------------------
    inline const size_t num_unsplittable_blocks() const { return _splittable_columns.size() - 1; }
   
    // --------------------- TEMP PRINTING FUNCTIONS ----------------------- 
    
    void print() const;
    
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Checks if the column shoul be monotone or not, and then sets it to be monotone or not
    /// @param[in]  col_idx         The index of the column to check if monotone
    /// @param[in]  num_elements    The number of elements in the column which are 0 and which are 1
    /// @param[in]  num_not_single  The number of elements in the columns which are part of a non singular row
    // ------------------------------------------------------------------------------------------------------
    void check_if_monotone(const uint col_idx, const uint* num_elements, const uint num_not_singular); 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Creates a vector of indices of the splittable columns
    /// @param[in]  splittable_info Information about if each of the columns is splittable or not
    // ------------------------------------------------------------------------------------------------------
    void create_splittable_column_vector(const binary_container_c& splittable_info);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the types of each of the columns - if the column is intrisically heterozygous,
    ///             or if it is not intrinsically heterozygous
    // ------------------------------------------------------------------------------------------------------
    void determine_column_types();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data from the input file
    /// @param[in]  data_file   The file to get the data from 
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Goes over each of the rows and processes the row with the row_process function, to
    ///             determine the start and end inidces of the row, and if the row is singular or not
    // ------------------------------------------------------------------------------------------------------
    void find_row_params(); 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Flips all elements of a column if there are more ones than zeros, and records that the
    ///             column has been flipped
    // ------------------------------------------------------------------------------------------------------
    void flip_column_bits(const uint col_idx);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_data(size_t row, TP token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes a row, determining all the necessary parameters, such as the start and end
    ///             indices, and hence if a row is singleton or not
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
        if (_monotone_columns.find(i) != _monotone_columns.end()) std::cout << i << " ";
    std::cout << "\n";

    std::cout << "\n|------------SPLITTABLE-----------|";
    
    for (auto i = 0; i < _splittable_columns.size(); ++i) 
        std::cout << _splittable_columns[i] << " ";
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
    fill(data_file);                    // Get the data from the input file
    determine_column_types();           // Determine which columns are splittable, montone, and IH or NIH
} 

// ------------------------------------------------- PRIVATE ------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::create_splittable_column_vector(const binary_container_c& splittable_info)
{
    bool first_found = false;
    // Go through the splittable columns and make a vector of the start inidices of the splittable columns
    for (size_t col_idx = 0; col_idx < C; ++col_idx) {
        if (splittable_info.get(col_idx) == 1) {
            if (!first_found) {
                if (_monotone_columns.find(col_idx) == _monotone_columns.end()) { // Not a monotone column
                        _splittable_columns.push_back(col_idx);
                        first_found = true;
                }
            } else {
                if (splittable_info.get(col_idx) == 1                           &&  // Column is splittable
                    _monotone_columns.find(col_idx) == _monotone_columns.end()  &&  // Not monotone
                    _column_types.get(col_idx) == IH                            ) { // Is an IH column
                        _splittable_columns.push_back(col_idx);
                }
            } 
            
            // If the last column is NIH, we need to add it
            if (col_idx == C - 1                                            &&      // Last column   
                _column_types.get(col_idx) == NIH                           &&      // NIH column
                _monotone_columns.find(col_idx) == _monotone_columns.end()  &&      // Not monotone
                _splittable_columns.back() != col_idx                       ) {     // Not added
                    _splittable_columns.push_back(col_idx);
            }
        }
    } 
}



template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::check_if_monotone(const uint  col_idx         , 
                                              const uint* num_elements    , 
                                              const uint  num_not_single  )
{
    // Check if the column in monotone, and if the are NIH
    // if they are none of these, then they are the default (IH)
    if (num_elements[0] != 0 && num_elements[1] == 0) {
        _monotone_columns[col_idx] = 0;
        
        // Set the values of the haplotypes at these positions
        _haplo_one.set(col_idx, 0);
        _haplo_two.set(col_idx, 0);
        
        // These columns also fit the IH category
        _column_types.set(col_idx, NIH);
    } else if (num_elements[1] != 0 && num_elements[0] == 0) {
        _monotone_columns[col_idx] = 0;
        
        // Set the value of the haplotypes
        _haplo_one.set(col_idx, 1);
        _haplo_two.set(col_idx, 1);
        
        _column_types.set(col_idx, NIH);
    } else if (!(std::min(num_elements[0], num_elements[1]) >= (num_not_single / 2))) {
        _column_types.set(col_idx, NIH);      
    }    
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::determine_column_types()
{
    // We need the row parameters to be set, so find those
    find_row_params();
    
    // We can use all the available cores for this 
    constexpr size_t threads_x = (THJ + THI) < C ? (THJ + THI) : C;

    // Binary container for if a columns is splittable or not (not by default)
    binary_container_c splittable_info;
    
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
                    uint col_idx            = it_x * threads_x + thread_idx;
                    uint num_elements[2]    = {0, 0};
                    uint num_not_single     = 0;
                    bool splittable         = true;
                    
                    // For each of the elements in the column (i.e row indices)
                    for (size_t row_idx = 0; row_idx < R; ++row_idx) {
                        auto element_value = _data.get(row_idx * C + col_idx);
                        
                        // Add to the (0's | 1's &| non singlton) count
                        if (element_value <= ONE) {
                            element_value == ZERO ? ++num_elements[0] : ++num_elements[1];
                            if (_singletons.get(row_idx) == ZERO) ++num_not_single;
                        }
                        
                        // Check for splittable condition
                        if (_row_info[row_idx * 2]      < col_idx &&         
                            _row_info[row_idx * 2 + 1]  > col_idx )
                            splittable = false;
                    }
                    
                    // Check if the column in monotone
                    check_if_monotone(col_idx, num_elements, num_not_single);
    
                    // If there are more ones than zeros and not monotone, flip all elements
                    if (num_elements[1] > num_elements[0]                           && 
                        _monotone_columns.find(col_idx) == _monotone_columns.end()  ) {
                        flip_column_bits(col_idx);
                    }
                    
                    // If the column is splittable, then add it to the splittable list
                    if (splittable) splittable_info.set(col_idx, 1);
                }
            }
        }
    );
   
    create_splittable_column_vector(splittable_info);
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
void Block<R, C, THI, THJ>::find_row_params()
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
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::flip_column_bits(const uint col_idx)
{
    for (uint row_idx = 0; row_idx < R; ++row_idx) {
        auto element_value = _data.get(row_idx * C + col_idx);
        if (element_value <= ONE) {
            element_value == ZERO 
                ? _data.set(row_idx * C + col_idx, ONE)
                : _data.set(row_idx * C + col_idx, ZERO);
        }
    }
    // Add that this column was flipped
    _flipped_columns[col_idx] = 0;
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
