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

namespace haplo {

// Define some HEX values for the 8 bits comparisons
#define ZERO 0x00
#define ONE  0x01
#define TWO  0x02
    
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
    using binary_container      = BinaryContainer<C>;           // Just a container of bits
    using row_info_container    = std::array<int, R * 2>;       // Start of row reads
    // ------------------------------------------------------------------------------------------------------
private:
    data_container      _data;                  //!< Container for { '0' | '1' | '-' } data variables
    binary_container    _haplotype;             //!< The binary bits which represent the haplotype
    binary_container    _alignment;             //!< The binary bits which represent the alignment of the 
                                                //!< reads to the haplotype                        
    binary_matrix       _pre_haplotypes;        //!< Potential haplotypes determined by looking at only the 
                                                //!< data from 1 row -- from which the optimal solution is
                                                //!< found
    row_info_container  _row_info;              //!< Information about the start and end positions of the row
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
    /// @brief      Sets the parameters for a potential solutions determined from a single row
    /// @param[in]
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes the data from the file, converting the chars into 2-bit values
    /// @param[in]  row     The row in the data matrix to fill
    /// @param[in]  line    The line which holds the elements as tokens
    /// @tparam     TP      The type of the tokenizer
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_data(size_t row, TP tp);
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

// ------------------------------------------------- PUBLIC -------------------------------------------------

template <size_t R, size_t C, size_t THI, size_t THJ>
Block<R, C, THI, THJ>::Block(const char* data_file)
{
    fill(data_file);
    filter();
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

template <size_t R, size_t C, size_t THI, size_t THJ>
void Block<R, C, THI, THJ>::update_row_params(size_t row_idx, size_t col_idx, int& start_index, int& end_index)
{
    size_t  element_idx = row_idx * C + col_idx;    
    byte    value       = _data.get(element_idx);
    
    if (value < TWO && start_index == -1) {
        start_index = col_idx;
        _alignment.set(row_idx, value == ZERO ? ONE : ZERO);
        _pre_haplotypes.set(element_idx, 0);
        std::cout << static_cast<unsigned>(_pre_haplotypes.get(element_idx)) << " ";
    } else if (value == ZERO && start_index >= 0) {
        end_index = col_idx;
        (_alignment.get(row_idx) == ONE)            ?
            _pre_haplotypes.set(element_idx, ZERO)  :
            _pre_haplotypes.set(element_idx, ONE )  ;
        std::cout << static_cast<unsigned>(_pre_haplotypes.get(element_idx)) << " ";
    } else if (value == ONE && start_index >= 0) {
        end_index = col_idx;
        (_alignment.get(row_idx) == ONE)             ?
            _pre_haplotypes.set(element_idx, ONE )   :
            _pre_haplotypes.set(element_idx, ZERO)   ;
        std::cout << static_cast<unsigned>(_pre_haplotypes.get(element_idx)) << " ";
    } else if (value >= TWO) {
       std::cout  << "- "; 
    }     
}

}           // End namespace haplo
#endif      // PARAHAPLO_BLOCK_HPP
