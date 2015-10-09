// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for a block
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_BLOCK_HPP
#define PARAHAPLO_BLOCK_HPP

#include "operations.hpp"
#include "read_info.hpp"
#include "snp_info.hpp"
#include "small_containers.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>
#include <tbb/tbb.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/parallel_sort.h>
#include <string>
#include <vector>
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
    using snp_info_container    = tbb::concurrent_unordered_map<size_t, SnpInfo>;
    using concurrent_umap       = tbb::concurrent_unordered_map<size_t, uint8_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    size_t              _rows;                  //!< The number of reads in the input data
    size_t              _cols;                  //!< The number of SNP sites in the container
    size_t              _first_splittable;      //!< 1st nono mono splittable solumn in splittale vector
    data_container      _data;                  //!< Container for { '0' | '1' | '-' } data variables
    read_info_container _read_info;             //!< Information about each read (row)
    snp_info_container  _snp_info;              //!< Information about each snp (col)
    concurrent_umap     _flipped_cols;          //!< Columns which have been flipped
    atomic_vector       _splittable_cols;       //!< A vector of splittable columns
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
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of subblocks in the block
    // ------------------------------------------------------------------------------------------------------
    inline size_t num_subblocks() const { return _splittable_cols.size() - _first_splittable; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the start index of a subblock (or the end index of the previous one) -- returns 0 if
    ///             the given index is out of range
    /// @param[in]  i   The index of the subblock 
    // ------------------------------------------------------------------------------------------------------
    inline size_t subblock(const size_t i) const 
    { 
        return i < _splittable_cols.size() ? _splittable_cols[_first_splittable + i] : 0;
    }
    
    // ------------------------------------------------------------------------------------------------------A
    /// @brief      Returns true if the requested column is monotone, false otherwise -- returns false if the
    ///             index is out of range
    /// @param[in]  i   The index of the column
    // ------------------------------------------------------------------------------------------------------
    inline bool is_monotone(const size_t i) const 
    {
        return i < _cols ? _snp_info.at(i).is_monotone() : false;
    }
    
    // ------------------------------------------------------------------------------------------------------A
    /// @brief      Returns true if the requested column is intrinsically herterozygous, false otherwise -- 
    ///             returns false if the index is out of range
    /// @param[in]  i   The index of the column
    // ------------------------------------------------------------------------------------------------------
    inline bool is_intrin_hetro(const size_t i) const 
    {
        return i < _cols ? (_snp_info.at(i).type() == IH) : false;
    }
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Fills the block with data from the input fill
    /// @param[in]  data_file   The file to get the data from 
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Flips all elements of a column if there are more ones than zeros, and records that the
    ///             column has been flipped
    /// @param[in]  col_idx         The index of the column to flip
    /// @param[in]  col_start_row   The start row of the column
    /// @param[in]  col_end_row     The end row of the column
    // ------------------------------------------------------------------------------------------------------
    void flip_column_bits(const size_t col_idx, const size_t col_start_row, const size_t col_end_row);
    
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
    /// @brief      Processses a snp (column), checking if it is IH or NIH, and if it is montone, or flipping 
    ///             the bits if it has more ones than zeros
    // ------------------------------------------------------------------------------------------------------
    void process_snps(); 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the parameters for a column -- the start and end index
    /// @param[in]  col_idx     The index of the column to set the parameters for
    /// @param[in]  row_idx     The index of the row to update in teh column info
    /// @param[in]  value       The value of the element at row_idx, col_idx
    // ------------------------------------------------------------------------------------------------------
    void set_col_params(const size_t col_idx, const size_t row_idx, const uint8_t value);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sorts the (small in almost all cases) splittable vector, and removes and montone columns 
    ///             from the start of the vector
    // ------------------------------------------------------------------------------------------------------
    void sort_splittable_cols();
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

// ----------------------------------------------- PUBLIC ---------------------------------------------------

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
Block<Elements, ThreadsX, ThreadsY>::Block(const char* data_file)
: _rows{0}, _cols{0}, _first_splittable{0}, _read_info{0}, _splittable_cols{0} 
{
    fill(data_file);                    // Get the data from the input file
    process_snps();                     // Process the SNPs to determine block params
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
    
    // Set the number of columns 
    _cols = _snp_info.size();    
}

template <size_t Elements, size_t ThreadsX, size_t ThreadsY> template <typename TokenPointer>
size_t Block<Elements, ThreadsX, ThreadsY>::process_data(size_t         offset  ,
                                                         TokenPointer&  line    )
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
                set_col_params(col_idx, _rows, ZERO);
                break;
            case '1':
                _data.set(offset++, ONE);
                set_col_params(col_idx,_rows, ONE);
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
                                                         const size_t  row_idx,
                                                         const uint8_t value  )
{
    if (_snp_info.find(col_idx) == _snp_info.end()) {
        // Not in map, so set start index to row index
        _snp_info[col_idx] = SnpInfo(row_idx, 0);
    } else {
        // In map, so start is set, set end 
        _snp_info[col_idx].end_index() = row_idx;
    }
    // Update the value counter
    value == ZERO ? _snp_info[col_idx].zeros()++
                  : _snp_info[col_idx].ones()++;
}

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
void Block<Elements, ThreadsX, ThreadsY>::process_snps()
{
    // We can use all the available cores for this 
    //const size_t threads = (ThreadsX + ThreadsY) < _cols ? (ThreadsX + ThreadsY) : _cols;
    const size_t threads = 1;
    
    // Binary container for if a columns is splittable or not (not by default)
    binary_vector splittable_info(_cols);
    
    // Over each column in the row
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads),
        [&](const tbb::blocked_range<size_t>& thread_ids)
        {
            for (size_t thread_id = thread_ids.begin(); thread_id != thread_ids.end(); ++thread_id) {
                size_t thread_iters = ops::get_thread_iterations(thread_id, _cols, threads);
    
                // For each column we need to "walk" downwards and check how many non singular rows there are 
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t col_idx      = it * threads + thread_id;
                    size_t non_single   = 0;                                // Number of non singular columns
                    bool   splittable   = true;                             // Assume splittable
                    auto&  col_info     = _snp_info[col_idx];
                    
                    // For each of the elements in the column 
                    for (size_t row_idx = col_info.start_index(); row_idx <= col_info.end_index(); ++row_idx) {
                        if (_read_info[row_idx].length() != 1) 
                            ++non_single;
                            
                        // Check for the splittable condition
                        if (_read_info[row_idx].start_index() < col_idx && 
                            _read_info[row_idx].end_index()   > col_idx  )
                                splittable = false;
                    }
                  
                    // If the column fits the non-intrinsically heterozygous criteria, change the type
                    if (std::min(col_info.zeros(), col_info.ones()) >= (non_single / 2))
                        col_info.set_type(NIH);
                        
                    // If there atre more 1's than 0's flip all the bits
                    if (col_info.ones() > col_info.zeros() && !col_info.is_monotone()) 
                        flip_column_bits(col_idx, col_info.start_index(), col_info.end_index());
                    
                    // If the column is splittable, add it to the splittable info 
                    if (splittable && !col_info.is_monotone()) _splittable_cols.push_back(col_idx);
                }
            }
        }
    );
    // Need to sort the splittable columns in ascending order
    sort_splittable_cols();
}

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
void Block<Elements, ThreadsX, ThreadsY>::flip_column_bits(const size_t col_idx       , 
                                                           const size_t col_start_row ,
                                                           const size_t col_end_row   )
{
    for (size_t row_idx = col_start_row; row_idx <= col_end_row; ++row_idx) {
        size_t mem_offset    = _read_info[row_idx].offset() + col_idx - _read_info[row_idx].start_index(); 
        auto   element_value = _data.get(mem_offset);
        
        // Check that the element is not a gap, then set it
        if (element_value <= ONE) {                             
            element_value == ZERO 
                ? _data.set(mem_offset, ONE)
                : _data.set(mem_offset, ZERO);
        }
    }
    // Add that this column was flipped
    _flipped_cols[col_idx] = 0;
}

template <size_t Elements, size_t ThreadsX, size_t ThreadsY>
void Block<Elements, ThreadsX, ThreadsY>::sort_splittable_cols()
{
    // Sort the splittable vector -- this is still NlgN, so not really parallel
    tbb::parallel_sort(_splittable_cols.begin(), _splittable_cols.end(), std::less<size_t>());
    
    // Set the start index to be the first non-monotone column
    // tbb doesn't have erase and it'll be slow to erase from the front
    while (_snp_info[_splittable_cols[_first_splittable]].is_monotone()) ++_first_splittable;
    
    // Check that the last column is in the vector (just some error checking incase)
    if (_splittable_cols[_splittable_cols.size() - 1] != _cols -1) 
        _splittable_cols.push_back(_cols - 1);
}


}           // End namespace haplo
#endif      // PARAHAPLO_BLOCK_HPP
