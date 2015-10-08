// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for a block
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_NEW_BLOCK_HPP
#define PARAHAPLO_NEW_BLOCK_HPP

#include "operations.hpp"
#include "read_info.hpp"
#include "small_containers.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>
#include <tbb/tbb.h>
#include <tbb/concurrent_unordered_map.h>

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>
#include <unordered_map>
#include <utility>

// NOTE: All output the the terminal is for debugging and checking at the moment

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
/// @tparam     Elements    The number of elements in the input data
/// @param      ThreadsX    The threads for the X direction 
/// @param      ThreadsY    The threads for the Y direction 
// ----------------------------------------------------------------------------------------------------------
template <size_t Elements, size_t ThreadsX = 1, size_t ThreadsY = 1>
class Block {
public:
    // ----------------------------------------- TYPES ALIAS'S ----------------------------------------------
    using data_container        = BinaryArray<Elements, 2>;    
    using binary_vector         = BinaryVector<2>;
    using atomic_type           = tbb::atomic<size_t>;
    using atomic_vector         = tbb::concurrent_vector<size_t>;
    using read_info_container   = std::vector<ReadInfo>;
    using col_info_container    = std::vector<size_t>;
    using col_map_type          = std::unordered_map<size_t, std::pair<size_t, size_t>>;
    // ------------------------------------------------------------------------------------------------------
private:
    size_t              _rows;                  //!< The number of reads in the input data
    size_t              _cols;                  //!< The number of SNP sites in the container
    data_container      _data;                  //!< Container for { '0' | '1' | '-' } data variables
    read_info_container _read_info;             //!< Information about each read
    col_map_type        _col_info;             //!< Information about each column (snp site)
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to fill the block with data from the input file
    /// @param[in]  data_file       The file to fill the data with
    // ------------------------------------------------------------------------------------------------------
    Block(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of an element, if it exists, otherwise returns 3
    /// @param[in]  row_idx     The row index of the element
    /// @param[in]  col_idx     The column index of the element
    // ------------------------------------------------------------------------------------------------------
    uint8_t operator()(const size_t row_idx, const size_t col_idx) const;

private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data from the input fill
    /// @param[in]  data_file   The file to get the data from 
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes a line of data
    /// @param      offset          The offset in the data container of the data
    /// @param      line            The data to proces
    /// @tparam     TokenPointer    The token pointer type
    /// @return     The new offset after processing
    // ------------------------------------------------------------------------------------------------------
    template <typename TokenPointer>
    size_t process_data(size_t offset, TokenPointer& line);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the parameters for a column -- the start and end index
    /// @param[in]  col_map     The map used to check if the column starrt has been set
    /// @param[in]  col_idx     The index of the column to set the parameters for
    /// @param[in]  row_idx     The index of the row to update in teh column info
    // ------------------------------------------------------------------------------------------------------
    void set_col_params(const size_t col_idx, const size_t row_idx);
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

// ----------------------------------------------- PUBLIC ---------------------------------------------------

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
Block<Elements, ThreadsX, ThreadsY>::Block(const char* data_file)
: _rows{0}, _cols{0}, _read_info{0}, _col_info{0} 
{
    fill(data_file);                    // Get the data from the input file
} 

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
uint8_t Block<Elements, ThreadsX, ThreadsY>::operator()(const size_t row_idx, const size_t col_idx) const 
{
    // If the element exists
    return _read_info[row_idx].element_exists(col_idx) == true 
        ? _data.get(_read_info[row_idx].offset() + col_idx - _read_info[row_idx].start_index()) : 0x03;
} 

// ------------------------------------------------- PRIVATE ------------------------------------------------

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
void Block<Elements, ThreadsX, ThreadsY>::fill(const char* data_file)
{
    // Open file and convert to string (for tokenizer)
    io::mapped_file_source file(data_file);
    if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
    
    std::string data(file.data(), file.size());

    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> nwline_separator{"\n"};
    
    // Tokenize the data into lines
    tokenizer lines{data, nwline_separator};
    
    // Create a counter for the offset in the data container
    size_t offset = 0;
   
    // Get the data and store it in the data container 
    for (auto& line : lines) {
        offset = process_data(offset, line);
        ++_rows;
    }
    
    if (file.is_open()) file.close();
}

template <size_t Elements, size_t ThreadsX, size_t ThreadsY> template <typename TokenPointer>
size_t Block<Elements, ThreadsX, ThreadsY>::process_data(size_t   offset    ,
                                                         TokenPointer&  line)
{
    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> wspace_separator{" "}; 

    // Tokenize the line
    tokenizer elements{line, wspace_separator};

    std::string read_data;
    size_t start_index = 0, end_index = 0, counter = 0;
    for (auto token : elements) {
        if (counter == 0) 
            start_index = stoul(token);
        else if (counter == 1)
            end_index = stoul(token);
        else
            read_data = token;
        counter++;
    }
    _read_info.emplace_back(_rows, start_index, end_index, offset);

    size_t col_idx = start_index;    
    // Put data into the data vector
    for (const auto& element : read_data) {
        switch (element) {
            case '0':
                _data.set(offset++, ZERO);
                set_col_params(col_idx, _rows);
                break;
            case '1':
                _data.set(offset++, ONE);
                set_col_params(col_idx,_rows);
                break;
            case '-':
                _data.set(offset++, TWO);
                break;
            default:
                std::cerr << "Error reading input data - exiting =(\n";
                exit(1);
        } ++col_idx;
    }
    return offset;
}

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
void Block<Elements, ThreadsX, ThreadsY>::set_col_params(const size_t  col_idx,
                                                         const size_t  row_idx)
{
    if (_col_info.find(col_idx) == _col_info.end()) {
        // Not in map, so set start index to row index
        _col_info[col_idx] = std::make_pair(row_idx, 0);
    } else {
        // In map, so start is set, set end 
        _col_info[col_idx].second = row_idx;
    }
}

}           // End namespace haplo
#endif      // PARAHAPLO_BLOCK_HPP
