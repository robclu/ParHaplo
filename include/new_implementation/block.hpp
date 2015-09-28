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

#include <iostream>
#include <string>
#include <vector>

// NOTE: All output the the terminal is for debugging and checking at the moment

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
    using data_container        = BinaryContainer<R * C, 2>;    // R*C Elements, 2 bits per element
    using binary_matrix         = BinaryContainer<R * C, 1>;    // A matrix of bits for each element
    using binary_container_r    = BinaryContainer<C>;           // A container of bits C elements long
    using binary_container_c    = BinaryContainer<R>;           // A container of bits R elements long
    using row_info_container    = std::array<int, R * 2>;       // Start and end of row reads
    using col_info_container    = std::array<int, C * 2>;       // Reads per col and and num 0's
    // ------------------------------------------------------------------------------------------------------
private:
    data_container      _data;                  //!< Container for { '0' | '1' | '-' } data variables
    binary_container_c  _haplo_one;             //!< The binary bits which represent the first haplotype
    binary_container_c  _haplo_two;             //!< The bianry bits which represent the second haplotype
    binary_container_r  _alignment;             //!< The binary bits which represent the alignment of the 
                                                //!< reads to the haplotypes                        
    binary_matrix       _pre_haplotypes;        //!< Potential haplotypes determined by looking at only the 
                                                //!< data from 1 row -- from which the optimal solution is
                                                //!< found
    row_info_container  _row_info;              //!< Information about the start and end positions of the row
    col_info_container  _col_info;
    
    // V2.0
    binary_container_r  _singletons;            //!< If each of the rows is singleton or not
    binary_container_c  _column_types;          //!< 0 if intrinsically heterozygous (IH), 1 if not IH
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
    void print_singletons() const;
    
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Updates the row parameters -- the best value for an element, and the start and end indices
    ///             of the row
    /// @param[in]  row_idx         The index of the row to update
    /// @param[in]  col_idx         The index of the column in the row 
    /// @param[in]  start_idx       The start index of the row (forst non-gap element)
    /// @param[in]  end_idx         The end index of the row (last non-gap element)
    // ------------------------------------------------------------------------------------------------------
    void update_row_params(size_t row_idx, size_t col_idx, int& start_idx, int& end_idx);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data from the input file
    /// @param[in]  data_file   The file to get the data from 
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Filers out the best possible solutions for each of the rows
    // ------------------------------------------------------------------------------------------------------
    void filter(); 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Filters out all the monotone columns and sets the values of the haplotypes at the monotone
    ///             column positions 
    // ------------------------------------------------------------------------------------------------------
    void filter_columns();

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Finds all the singleton rows
    // ------------------------------------------------------------------------------------------------------
    void find_singletons();

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes the data from the file, converting the chars into 2-bit values
    /// @param[in]  row     The row in the data matrix to fill
    /// @param[in]  line    The line which holds the elements as tokens
    /// @tparam     TP      The type of the tokenizer
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_data(size_t row, TP tp);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Filters all columns by finding the monotone columns, as well as counting th
    // ------------------------------------------------------------------------------------------------------
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets all the homozygous variables if all reads in a column had the same value
    // ------------------------------------------------------------------------------------------------------
    void set_homozygous_sites();
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

// ----------------------------------------------- TESTING --------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::print() const 
{
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
    
    // Print column info
    for (size_t c = 0; c < C; ++c)
        std::cout << static_cast<unsigned>(_column_types.get(c)) << " ";
    std::cout << "\n";
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::print_singletons() const 
{
    for (size_t i = 0; i < R; i++) {
        std::cout << static_cast<unsigned>(_singletons.get(i)) << "\n";
    }
    
}

// ------------------------------------------------- PUBLIC -------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
Block<R, C, THI, THJ>::Block(const char* data_file)
: _row_info({0}), _col_info({0})
{
    fill(data_file);
    find_singletons();
    filter_columns();
    //set_homozygous_sites();

    // Debugging ....
    /*
    // Print out the column information
    for (int i = 0; i < 30; ++i) std::cout << "-";
    std::cout << "\nTR :  ";
    for (size_t i = 0; i < 2 * C; i += 2) std::cout << _col_info[i] << " ";
    std::cout << "\nN0 :  ";
    for (size_t i = 1; i < 2 * C; i += 2) std::cout << _col_info[i] << " ";
    std::cout << "\n";    
    for (int i = 0; i < 30; ++i) std::cout << "-";
    std::cout << "\nh :   ";    
    
    // Print the haplotype
    for (size_t i = 0; i < C; ++i) std::cout << static_cast<unsigned>(_haplotype.get(i)) << " ";
    std::cout << "\nh` :  "; 
    for (size_t i = 0; i < C; ++i) {
        if (_homozygous_sites.get(i) == ONE)
            std::cout << static_cast<unsigned>(_haplotype.get(i)) << " ";
        else 
            std::cout << static_cast<unsigned>((_haplotype.get(i) ^ ONE) & ONE) << " ";
    }
    std::cout << "\n";    
    for (int i = 0; i < 30; ++i) std::cout << "-";
    std::cout << "\n";    
    */
} 

// ------------------------------------------------- PRIVATE ------------------------------------------------


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
    for (tokenizer::iterator line = lines.begin(); line != lines.end(); ++line) process_data(row++, line);
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::filter()
{
    // Check that we aren't trying to use more threads than rows
    constexpr size_t threads_y = THI < R ? THI : R;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y), 
        [&](const tbb::blocked_range<size_t>& thread_ids_y)
        {
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy, R, threads_y);
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t  row_id          = it_y * threads_y + thread_idy;
                    int     start_index     = -1;           // First non gap in row
                    int     end_index       =  0;           // Last non gap in row
                    
                    std::cout << "R : " << row_id << " ";   // Debugging
                    
                    // For all of the columns in the row
                    for (size_t col_id = 0; col_id < C; ++col_id) 
                        update_row_params(row_id, col_id, start_index, end_index);
                        
                    // Update the start and end positions of the row
                    _row_info[row_id] = start_index; _row_info[row_id + 1] = end_index;
                    std::cout << " : " << start_index << " " << end_index << "\n";
                }
            }
        }
    );
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::filter_columns()
{
    constexpr size_t threads_x = THJ < C ? THJ : C;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x)
        {
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, C, threads_x);
                
                for (size_t it_x = 0; it_x < thread_iters_x; ++it_x) {
                    size_t col_id               = it_x * threads_x + thread_idx;
                    size_t num_non_singular     = 0;            // Number of 0's and 1's that are non_singluar
                    size_t num_any[2]           = {0, 0};       // Number of 0's and 1's in column
                    
                    for (size_t row_id = 0; row_id < R; ++row_id) {
                        if (_data.get(row_id * C + col_id) == ZERO) {
                            ++num_any[0];
                            if (_singletons.get(row_id) == ZERO) ++num_non_singular;
                        } else if (_data.get(row_id * C + col_id) == ONE) {
                            ++num_any[1];
                            if (_singletons.get(row_id) == ZERO) ++num_non_singular;
                        }
                    }
                    std::cout << col_id << " " << num_any[0] << " " << num_any[1] << " " << num_non_singular << "\n";
                    if (std::min(num_any[0], num_any[1]) < (num_non_singular / 2)) 
                        _column_types.set(col_id, NIH);
                    
                    // TODO: Remove monotone columns and set solution
                }
            }
        }
    );
}


template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::find_singletons() 
{
    constexpr size_t threads_y = THI < R ? THI : R;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y),
        [&](const tbb::blocked_range<size_t>& thread_ids_y)
        {
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy, R, threads_y);  
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t row_id           = it_y * threads_y + thread_idy;
                    size_t col_id           = 0;
                    size_t elements_in_row  = 0;
                    
                    while (elements_in_row <= 1 && col_id < C) {
                        if (_data.get(row_id * C + col_id++) <= 1) 
                            ++elements_in_row;
                    } std::cout << "\n";
                    // Check if only one (or 0 =/) element was found
                    if (elements_in_row <= 1) _singletons.set(row_id, 1);
                }
            }
        }
    );
}

template <size_t R, size_t C, size_t THI, size_t THJ> template <typename TP>
void Block<R, C, THI, THJ>::process_data(size_t row, TP tp) 
{
    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> wspace_separator{" "}; 

    // Create a string from the line token and tokenize it
    std::string line(*tp);
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

/*
template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::update_row_params(size_t row_idx, size_t col_idx, int& start_index, int& end_index)
{
    size_t  element_idx = row_idx * C + col_idx;    
    byte    value       = _data.get(element_idx);
    
    if (value < TWO && start_index == -1) {
        start_index = col_idx;                                  
        
        // Set the value of the aligment and set the haplotype value
        _alignment.set(row_idx, value == ZERO ? ONE : ZERO);
        _pre_haplotypes.set(element_idx, ZERO);
        
        // Increment that this column has a value, and if it is a zero or 1
        ++_col_info[2 *col_idx];
        _col_info[2 * col_idx + 1] = (value == ZERO ?
                                        _col_info[2 * col_idx + 1] + 1  :
                                        _col_info[2 * col_idx + 1]      );
        
        std::cout << static_cast<unsigned>(_pre_haplotypes.get(element_idx)) << " ";
    } else if (value == ZERO && start_index >= 0) {
        end_index = col_idx;
        
        // Update the haplotype based on the aligment and that the read element is a 0
        (_alignment.get(row_idx) == ONE)            ?
            _pre_haplotypes.set(element_idx, ZERO)  :
            _pre_haplotypes.set(element_idx, ONE )  ;
        
        // Increment that the column has a read element for the row, and that it's a 0
        ++_col_info[2 * col_idx];
        ++_col_info[2 * col_idx + 1];
            
        std::cout << static_cast<unsigned>(_pre_haplotypes.get(element_idx)) << " ";
    } else if (value == ONE && start_index >= 0) {
        end_index = col_idx;
        
        // Update the haplotype based on the aligment and that the read element is a 1
        (_alignment.get(row_idx) == ONE)             ?
            _pre_haplotypes.set(element_idx, ONE )   :
            _pre_haplotypes.set(element_idx, ZERO)   ;
        
        // Increment that this column has a read element for the row, but that it's a 1
        ++_col_info[2 * col_idx];
        
        std::cout << static_cast<unsigned>(_pre_haplotypes.get(element_idx)) << " ";
    } else if (value >= TWO) {
       std::cout  << "- "; 
       
    }     
}

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::set_homozygous_sites()
{
    // Check that we aren't trying to use more threads than columns
    constexpr size_t threads_x = THJ < C ? THJ : C;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x), 
        [&](const tbb::blocked_range<size_t>& thread_ids_x)
        {
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, C, threads_x);
                
                // Check the column information to see if this column is homozygous
                for (size_t it_x = 0; it_x < thread_iters_x; ++it_x) {
                    // *2 is because the array holds the number of elements for each column
                    // as well as the number of elements that were a zero for the column
                    size_t col_idx = 2 * it_x * threads_x + thread_idx;
                    
                    if (_col_info[col_idx] == _col_info[col_idx + 1]) {                 // All reads are 0
                        _haplotype.set(col_idx / 2, ZERO);
                        _homozygous_sites.set(col_idx / 2, 1);
                    } else if (_col_info[col_idx] > 0 && _col_info[col_idx +1] == 0) {  // All reads are 1
                        _haplotype.set(col_idx / 2, ONE);
                        _homozygous_sites.set(col_idx / 2, 1);
                    }
                }
            }
        }
    );
}
 */   
}           // End namespace haplo
#endif      // PARAHAPLO_BLOCK_HPP
