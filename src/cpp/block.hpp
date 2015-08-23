// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BLOCK_HPP
#define PARAHAPLO_CPP_BLOCK_HPP

#include "block_expressions.hpp"
#include "block_exceptions.hpp"
#include "read.hpp"

#include <tbb/tbb.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>

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
private:
    container_type      _data;          //!< Container for the data - just std::vector<T>(Rows * Cols)
public:
    std::vector<Read>   _read_info;     //!< Vector for holding the information for each read (row)
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
    /// @brief  Gets the size of the block (total number of elements)
    /// @return The total number of elements in the block
    // ------------------------------------------------------------------------------------------------------
    inline size_type size() const { return _data.size(); }
    
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
    /// @brief      Gets a data element from the Block.
    /// @param[in]  r       The row in the block of the data element to get
    /// @param[in]  c       The column in the block of the data element to get
    /// @return     The value of the element at position [r, c] in the block
    // ------------------------------------------------------------------------------------------------------
    inline unsigned int operator()(const int r, const int c) const { return _data[r * Cols + c].value(); } 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the start (forst 0 | 1) and end (last 0 | 1)  positions of each of the reads 
    ///             (rows)
    /// @param[in]  num_threads     The number of thereads to use for getting the row information
    // ------------------------------------------------------------------------------------------------------
    void get_read_info(const size_t num_threads);

private:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines how many iterations each thread must perform in a given situation. Say for
    ///             example there is a block with 9 rows and we have 4 threads, and each thread does some 
    ///             operations on a row, then each thread will need to peform 2 iterations, except for one of
    ///             the threads, which will need to perform 3. So the iteration mapping for the threads would
    ///             be:                                                                                     \n\n
    ///             Thread Id  | Rows to operate on | Iterations                                            \n\n
    ///             0          | 0, 4, 8            | 3
    ///             1          | 1, 5               | 2
    ///             2          | 2, 6               | 2
    ///             3          | 3, 7               | 2
    /// @param[in]  thread_id       The thread number 
    /// @param[in]  total_ops       The total number of operations (9 in the above example)
    /// @param[in]  num_threads     The number of threads being used
    /// @return     The total number of iterations for the thread
    // ------------------------------------------------------------------------------------------------------
    size_t get_thread_iterations(const size_t thread_id, const size_t total_ops, const size_t num_threads) const;
};

// ---------------------------------------- IMPLEMENTATIONS -------------------------------------------------

// -------------------------------------------- PUBLIC ------------------------------------------------------

template <size_t Rows, size_t Cols>
Block<Rows, Cols>::Block(const std::string filename, size_t num_threads) 
{
    fill(filename, num_threads);
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
                size_t thread_iters = get_thread_iterations(idx, Rows, threads_to_use);
                
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
            std::cout << _data[r * Cols + c].value() << " ";
        }   
        std::cout << "\n";
    }
}

// --------------------------------------------- PRIVATE ----------------------------------------------------

template <size_t Rows, size_t Cols>
size_t Block<Rows, Cols>::get_thread_iterations(const size_t thread_id  , 
                                                const size_t total_ops  , 
                                                const size_t num_threads) const 
{
    return (total_ops / num_threads) + 
            (thread_id < (total_ops % num_threads) ? 1 : 0);
}

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
                size_t thread_iters = get_thread_iterations(idx, Rows, threads_to_use);

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

#endif      // PARAHAPLO_CPP_BLOCK_HPP
