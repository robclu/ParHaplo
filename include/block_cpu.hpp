// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library -- cpu implementation
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_BLOCK_IMPLEMENTATION_CPU_HPP
#define PARAHAPLO_BLOCK_IMPLEMENTATION_CPU_HPP

#include "block_interface.hpp"
#include "operations.hpp"
#include "read.hpp"
#include "subblock_info.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>
#include <tbb/tbb.h>

#include <bitset>
#include <iostream>
#include <vector>

namespace io = boost::iostreams;

namespace haplo {

// Specialization for CPU block implementation    
template <size_t Rows, size_t Cols, size_t Cores>
class BlockImplementation<Rows, Cols, Cores, Device::CPU> {
public:
    // ----------------------------------- TYPEDEFS ---------------------------------------------------------
    using data_type         = haplo::Data;
    using read_type         = haplo::Read;
    using subinfo_type      = haplo::SubBlockInfo;
    using data_container    = std::array<data_type, Rows * Cols>;
    using read_container    = std::array<read_type, Rows>;
    using subinfo_container = std::vector<subinfo_type>;
    using colinfo_container = std::bitset<Cols>;
    using reference_type    = data_type&;
    using value_type        = data_type;
    // ------------------------------------------------------------------------------------------------------
private:
    data_container      _data;              //!< The data for the block
    read_container      _read_info;         //!< Information for all the reads
    subinfo_container   _subblock_info;     //!< Information for sub-blocks
    colinfo_container   _column_info;       //!< Information for if each column is splittable or not
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- uses the given filename to get the data. This constructor must be used to 
    ///             ensure that the block implementation has valid data
    /// @param[in]  data_file    The name of the data file to get the input from
    // ------------------------------------------------------------------------------------------------------
    explicit BlockImplementation(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of cores available to the block
    /// @return     The number of cores available to the block 
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t num_cores() { return Cores; }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of rows in the block 
    /// @return     The number of rows in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t rows() { return Rows; }
    
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of columns in the block
    /// @return     The number of columns in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t cols() { return Cols; }         
  
    //-------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size (numbre of elements) in the block
    /// @return     The number of elements in the block
    //-------------------------------------------------------------------------------------------------------
    static constexpr size_t size() { return Rows * Cols; } 

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the information for the sub block i of this block
    /// @param[in]  i   The index of the sub-block for which the information must be given
    /// @return     The subblock information for the sub-block with index i
    // ------------------------------------------------------------------------------------------------------
    const subinfo_type& subblock_info(const size_t i) const { return _subblock_info.at(i); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a refernce to the element of the block at row and column -- does not implement bound
    ///             checking as this is optimozed for performance, .at() can be implemented at a later stage 
    ///             if bound checking access is required
    /// @param[in]  row     The row of the element in the block
    /// @param[in]  col     The column of the element in the block 
    // ------------------------------------------------------------------------------------------------------
    reference_type operator()(const size_t row, const size_t col) { return _data[row * Cols + col]; }
 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the read at row i of the block -- nor error checking for performance, again .at can
    ///             be implemented at a later stage
    /// @param[in]  The index of the read to get from the block
    /// @return     The read at row i in the block
    // ------------------------------------------------------------------------------------------------------
    const read_type& operator[](const size_t i) const { return _read_info[i]; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Loads data into the Block from the given file. The function will look for as many elements 
    ///             as are specified by the Block dimensions, so for a Block with 3 rows and 6 columns, the
    ///             function will look for 3 lines each with 6 elements.                                     \n
    ///                                                                                                      \n
    ///             The file is loaded as a memory mapped file since the input data files will likely be     \n
    ///             huge, which should save have a significant performance incease.                          \n
    ///                                                                                                      \n
    ///             See boost memory mapped files for reference:                                             \n
    ///                 http://www.boost.org/doc/libs/1_38_0/libs/iostreams/doc/index.html                   \n
    ///                                                                                                      \n
    ///             Or if you download the boost libraries the source for the memory mapped files is at:     \n
    ///                 boost/iostreams/device/mapped_file.hpp                                               \n
    /// @param[in]  filename        The name of the file to get the input data from.
    // ------------------------------------------------------------------------------------------------------
    void fill(const char* data_file);
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Creates a list of column indices which are unsplittable, which can then be used by
    ///             different threads to create filtered (unnecessary data removed) sub-blocks which can be 
    ///             sovled  with ILP using a branch and bound implementation
    // ------------------------------------------------------------------------------------------------------
    void determine_subblock_info(); 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The column information stores all columns as splittable by default since this is the 
    ///             default initialization of the bitset struct which is being used, so this functions changes
    ///             all the non-splittable column information in the _column_info struct
    // ------------------------------------------------------------------------------------------------------
    void find_unsplittable_columns(); 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the start (first 0 | 1) and end (last 0 | 1)  positions of each of the reads 
    // ------------------------------------------------------------------------------------------------------
    void get_read_info(); 
};

// ------------------------------------------------ IMPLEMENTATIONS -----------------------------------------

// ---------------------------------------------------- PUBLIC ----------------------------------------------

template <size_t Rows, size_t Cols, size_t Cores>
BlockImplementation<Rows, Cols, Cores, Device::CPU>::BlockImplementation(const char* data_file) 
{
    fill(data_file);                    // Fill the data file with 
    get_read_info();                    // Get the read information for the block
    find_unsplittable_columns();        // Find all the columns which cannot be split
    determine_subblock_info();          // Determine the information for sub-blocks 
}

template <size_t Rows, size_t Cols, size_t Cores>
void BlockImplementation<Rows, Cols, Cores, Device::CPU>::fill(const char* data_file)
{
    using namespace io;
    
    // Create a readonly memory mapped file to get the data
    io::mapped_file_source mapped_input(data_file);
    const char* input_data = mapped_input.data();

    // Check that we arent't using more threads than rows
    constexpr size_t threads_to_use = Rows < Cores ? Rows : Cores;
    
    // Parallel tasks to get the input data from the
    // memory mapped file into the class _data container
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices) 
        {
            // This first loop is actually in parallel
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                // Determine the number of iterations this thread must  perform
                size_t thread_iters = ops::get_thread_iterations(idx, Rows, threads_to_use);
                
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t row_offset = threads_to_use * it + idx;      // Offset due to row in data
                    size_t map_offset = row_offset * (Cols * 2);        // Add whitespace and \n char
                    
                    for (size_t elem_idx = map_offset; elem_idx < map_offset + (Cols * 2); elem_idx += 2) {
                        _data[elem_idx / 2] = *(input_data + elem_idx);
                    }
                }
            }
        }
    );
    mapped_input.close();   
}

// ---------------------------------------------------- PRIVATE ---------------------------------------------

template <size_t Rows, size_t Cols, size_t Cores>
void BlockImplementation<Rows, Cols, Cores, Device::CPU>::determine_subblock_info()
{
    // Incase the container is not empty
    if (_subblock_info.size() != 0) _subblock_info.clear();
    
    for (int i = 1, start = 0; i < _column_info.size(); ++i) {
        if (_column_info[i] == 0) {                 // We have found the end, add an element to unsplittable
            _subblock_info.emplace_back(start, i);
            start = i;                              // Update the new start
        }
    }
}

template <size_t Rows, size_t Cols, size_t Cores>
void BlockImplementation<Rows, Cols, Cores, Device::CPU>::find_unsplittable_columns()
{
    // Check that we aren't trying to use more threads than there are columns
    size_t threads_to_use = Cores > Cols ? Cols : Cores;
    
    // Create some threads where each one determines if a column is splittable
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices) 
        {
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                size_t thread_iters = ops::get_thread_iterations(idx, Cols, threads_to_use);
                
                for (size_t it = 0; it < thread_iters; ++it) {
                    // Go through each of the reads and check if idx lies between the 
                    // start and the end position of the read, which makes the col not
                    // splittable. Since bitset default intializes to all 0's, a 0 in the 
                    // _column_indo array means that a columns is not splittable
                    size_t  index      = threads_to_use * it + idx;
                    bool    splittable = true;                                      // Default as splittable
                    size_t  i          = 0;
                    while (splittable && i < Rows) {
                        if (index > _read_info[i].start() && index < _read_info[i].end()) {
                            _column_info[index] = 1;
                            splittable          = false;                            // Break while
                        } ++i;
                    }
                }
            }
        }
    );
}

template <size_t Rows, size_t Cols, size_t Cores>
void BlockImplementation<Rows, Cols, Cores, Device::CPU>::get_read_info() 
{
    // Check that we aren't trying to use more threads than there are rows
    size_t threads_to_use = Cores < Rows ? Cores : Rows;
    
    // Create some threads where each one will get the information of a row
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices)
        {
            // As above, this is a loop is parallelized by tbb
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                // Determine the number of iterations for each thread
                size_t thread_iters = ops::get_thread_iterations(idx, Rows, threads_to_use);

                for (size_t it = 0; it < thread_iters; ++it) {
                    int read_start = -1, read_end = -1, counter = 0;
                    
                    size_t start_offset = (threads_to_use * it + idx) * Cols;
                    size_t end_offset   = start_offset + Cols;
                    
                    for (auto elem = _data.begin() + start_offset; elem != _data.begin() + end_offset; ++elem) {
                        if (elem->value() != 2 ) {
                            if (read_start != -1)           // If a start position has been found for the read
                                read_end = counter;
                            else                            // If a start position has not been found 
                                read_start = counter;
                        }
                        ++counter;
                    }
                    if (read_end == -1) read_end = read_start;     // Case for a read with only a single element
                    _read_info[it * threads_to_use + idx] =  Read(read_start, read_end);
                }
            }
        }
    );
}

}           // End namespace haplo

#endif      // PARAHAPLO_BLOCK_CPU_HPP
