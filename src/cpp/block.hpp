// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BLOCK_HPP
#define PARAHAPLO_CPP_BLOCK_HPP

#include "block_expressions.hpp"
#include "block_exceptions.hpp"
#include "operations.hpp"
#include "read.hpp"
#include "subblock_info.hpp"

#include <tbb/tbb.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>

#include <bitset>
#include <cassert>
#include <string>
#include <iostream>

namespace io = boost::iostreams;

namespace haplo {
 
// ----------------------------------------------------------------------------------------------------------
/// @class  Block 
/// @brief  Defines a block - which is an input to the paplotype assemebly problem. Operations can be 
///         performed on the block, such as decomposition, row, column and duplicate removal to simplify 
///         the block making the problem easier to compute. Block inherits BlockExpression so that we can
///         do operations on the block that look as their mathematical expressions would.
/// @tparam Rows    The number of rows in the block
/// @tparam Cols    The number of columns in the block
// ----------------------------------------------------------------------------------------------------------
template <size_t Rows, size_t Cols>
class Block : public BlockExpression<Block<Rows, Cols>> {
public:
    // --------------------------------------- Typedefs -----------------------------------------------------
    using container_type    = typename BlockExpression<Block<Rows, Cols>>::container_type;
    using size_type         = typename BlockExpression<Block<Rows, Cols>>::size_type;
    using value_type        = typename BlockExpression<Block<Rows, Cols>>::value_type;
    using reference         = typename BlockExpression<Block<Rows, Cols>>::reference;
    // ------------------------------------------------------------------------------------------------------
    using container_subinfo     = std::vector<SubBlockInfo>;
    using container_read        = std::vector<Read>;

    // Incease we want to use a constant expression version    
    static constexpr size_t num_rows = Rows;
    static constexpr size_t num_cols = Cols;
private:
    container_type      _data;                  //!< Container for the data
    container_subinfo   _subblock_info;         //!< Information for the unsplittable sub-blocks
public:
    std::bitset<Cols>   _column_info;           //!< Information for if the block are/aren't splittable
    container_read      _read_info;             //!< Vector for holding the information for each read (row)
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Defult constructor, intitlizes the data to empty
    // ------------------------------------------------------------------------------------------------------
    Block() : _data(Rows * Cols) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructs a Block from an input data file - calls the fill() function.
    /// @param[in]  filename        The name of the file which has the input data
    /// @param[in]  num_threads     The number of threads to use for loading the data
    // ------------------------------------------------------------------------------------------------------
    Block(const std::string filename, const size_t num_threads = 1);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructs a Block from a block expression
    /// @param[in]  expression  The expression to create the Block from
    /// @tparam     Expression  The type of expression 
    // ------------------------------------------------------------------------------------------------------
    template <typename Expression>
    Block(BlockExpression<Expression> const& expression) 
    {
        Expression const& expr = expression;
        _data.resize(expr.size());
        for (size_type i = 0; i < expr.size(); ++i) {
            _data[i] = expr[i];                                     // Put data from expr into this block
        }
    }
    
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
    ///                                                                                                      \n
    ///             The function also uses threads to improve performance - the number of threads to use needs
    ///             to be specified if the default of 1 is not used.
    /// @param[in]  filename        The name of the file to get the input data from.
    /// @param[in]  num_threads     The number of threads to use 
    // ------------------------------------------------------------------------------------------------------
    void fill(const std::string fileanme, const size_t num_threads = 1);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the block (total number of elements)
    /// @return     The total number of elements in the block
    // ------------------------------------------------------------------------------------------------------
    inline size_type size() const { return _data.size(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of rows in the Block
    /// @return     The number of rows in the Block 
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_type rows() { return Rows; }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of columns in the Block
    /// @return     The number of columns in the Block 
    // ------------------------------------------------------------------------------------------------------
    static constexpr size_type columns() { return Cols; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Assigns a pointer to data to the vector so that we can use the vector as data rather than 
    ///             the raw pointer. There shouldn't be any performance loss, but we will have to test that
    /// @param[in]  data_source     A reference to the data array
    /// @tparam     N               The number of elements in the data array
    template <size_t N>
    inline void assign_data(const haplo::Data (&data_source)[N])
    {
        static_assert( N == (Rows * Cols), "Invalid data source size for block" );
        _data.assign(data_source, data_source + N);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the Block's data
    /// @return     A reference to the Block's data
    // ------------------------------------------------------------------------------------------------------
    inline container_type const& get_data() const { return _data; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Print the block data (for debugging and verifying operation results quickly)
    // ------------------------------------------------------------------------------------------------------
    void print() const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Interface for getting the read information of the block 
    /// @return     A constant reference to the read info 
    // ------------------------------------------------------------------------------------------------------
    container_read const& read_info() const { return _read_info; }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Interface for gettting the unsplittable sub block information 
    // ------------------------------------------------------------------------------------------------------
    container_subinfo const& subblock_info() const { return _subblock_info; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the element at position i
    /// @param[in]  i   The index of the element in the data container
    /// @return     The value of the element at position i in the data container
    // ------------------------------------------------------------------------------------------------------
    inline value_type operator[](size_type i) const { return _data[i]; }
private: 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the start (forst 0 | 1) and end (last 0 | 1)  positions of each of the reads 
    ///             (rows)
    /// @param[in]  num_threads     The number of thereads to use for getting the row information
    // ------------------------------------------------------------------------------------------------------
    void get_read_info(const size_t num_threads);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      The column information stores all columns as splittable by default since this is the 
    ///             default initialization of the bitset struct which is being used, so this functions changes
    ///             all the non-splittable column information in the _column_info struct
    /// @param[in]  num_threads     The number of threads to use
    // ------------------------------------------------------------------------------------------------------
    void find_unsplittable_columns(const size_t num_threads = 1);
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Creates a list of column infices which are unsplittable, which can then be used by
    ///             different threads to create filtered (unnecessary data removed) sub-blocks which can be 
    ///             sovled  with ILP
    /// @return     A vector of SubBlockInfo (an interface to make the code more readable) which is the
    ///             information for the splittable column.
    // ------------------------------------------------------------------------------------------------------
    void get_subblock_info();
};

// ---------------------------------------- IMPLEMENTATIONS -------------------------------------------------

// -------------------------------------------- PUBLIC ------------------------------------------------------

template <size_t Rows, size_t Cols>
Block<Rows, Cols>::Block(const std::string filename, size_t num_threads) 
{
    fill(filename, num_threads);                                    // Fill the block with information
    get_read_info(num_threads);                                     // Get the read information for the block
    find_unsplittable_columns(num_threads);                         // Find all the unsplittable columns 
    get_subblock_info();                                            // Get all the sub-block information
}

template <size_t Rows, size_t Cols>
void Block<Rows, Cols>::fill(const std::string filename, const size_t num_threads) 
{
    using namespace io;
    
    // Create a readonly memory mapped file to get the data
    io::mapped_file_source mapped_input(filename.c_str());
    const char* input_data = mapped_input.data();

    // Make sure that we have enough space for all the data. We need to use 
    // resize since the vector must hold data that gets modified. If we just 
    // allocate, then a thread may try to insert at pos 100 when the vector 
    // only has 40 elements due to the random order of thread execution, which 
    // will cause undefined behaviour
    _data.resize(Rows * Cols);

    // Check that we arent't using more threads than rows
    size_t threads_to_use = Rows < num_threads ? Rows : num_threads;
    
    // Parallel tasks to get the input data from the memory
    // mapped file into the class _data vector 
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

template <size_t Rows, size_t Cols>
void Block<Rows, Cols>::print() const 
{
    for (int r = 0; r < Rows; ++r) {
        for (int c = 0; c < Cols; ++c) {
            std::cout << static_cast<unsigned int>(_data[r * Cols + c].value()) << " ";
        }   
        std::cout << "\n";
    }
}

// --------------------------------------------- PRIVATE ----------------------------------------------------

template <size_t Rows, size_t Cols>
void Block<Rows, Cols>::get_read_info(const size_t num_threads) 
{
    // Check that we aren't trying to use more threads than there are rows
    size_t threads_to_use = (num_threads < Rows ? num_threads : Rows);
    
    _read_info.resize(Rows);                                    // Allocate some space in the vector
    
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

template <size_t Rows, size_t Cols>
void Block<Rows, Cols>::find_unsplittable_columns(const size_t num_threads)
{
    // Check that we aren't trying to use more threads than there are columns
    size_t threads_to_use = (num_threads > Cols ? Cols : num_threads);
    
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
                        }
                        ++i;
                    }
                }
            }
        }
    );
}

template <size_t Rows, size_t Cols>
void Block<Rows, Cols>::get_subblock_info()
{
    // The class doesn't intitialize the size of the unsplittable_info
    // container so we assume 0 (we could clear it here, but that's runtime 
    // overhead that we don't want, so if this function is only called once
    // then the behaviour will be as expected)
    for (int i = 1, start = 0; i < _column_info.size(); ++i) {
        if (_column_info[i] == 0) {                 // We have found the end, add an element to unsplittable
            _subblock_info.emplace_back(start, i);
            start = i;                              // Update the new start
        }
    }
}

}           // End namespace haplo

#endif      // PARAHAPLO_CPP_BLOCK_HPP
