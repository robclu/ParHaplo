// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BLOCK_HPP
#define PARAHAPLO_CPP_BLOCK_HPP

#include "block_expressions.hpp"
#include "block_exceptions.hpp"
#include "../general/input_parser.hpp"

#include <omp.h>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>

#include <cassert>
#include <string>
#include <iostream>

namespace sp = boost::spirit;
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
template <std::size_t Rows, std::size_t Cols>
class Block : public BlockExpression<Block<Rows, Cols>> {
public:
    // --------------------------------------- Typedefs -----------------------------------------------------
    using container_type    = typename BlockExpression<Block<Rows, Cols>>::container_type;
    using size_type         = typename BlockExpression<Block<Rows, Cols>>::size_type;
    using value_type        = typename BlockExpression<Block<Rows, Cols>>::value_type;
    using reference         = typename BlockExpression<Block<Rows, Cols>>::reference;
    // ------------------------------------------------------------------------------------------------------
private:
    container_type  _data;          //!< Container for the data - just std::vector<T>(Rows * Cols)
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Defult constructor, intitlizes the data to empty
    // ------------------------------------------------------------------------------------------------------
    Block() : _data(Rows * Cols) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructs a Block from an input data file - calls the fill() function.
    /// @param[in]  filename        The name of the file which has the input data
    /// @param[in]  num_threads     The number of threads to use for loading the data
    // ------------------------------------------------------------------------------------------------------
    Block(const std::string filename, const std::size_t num_threads = 1);
    
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
    void fill(const std::string fileanme, const std::size_t num_threads = 1);
    
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
    template <std::size_t N>
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
};

// ---------------------------------------- Implementations -------------------------------------------------

template <std::size_t Rows, std::size_t Cols>
Block<Rows, Cols>::Block(const std::string filename, std::size_t num_threads) 
{
    fill(filename, num_threads);
}

template <std::size_t Rows, std::size_t Cols>
void Block<Rows, Cols>::fill(const std::string filename, const std::size_t num_threads) 
{
    using namespace sp;
    using namespace io;
    
    // Create a readonly memory mapped file to get the data
    io::mapped_file_source mapped_input(filename.c_str());
    const char* input_data = mapped_input.data();

    int counter = 0;
    for (auto& elem : mapped_input) {
        if (counter++ % 2 == 0 ) {
            std::cout  <<  elem;   
        }
    } 
    std::cout << std::endl;
    _data.reserve(Rows * Cols);                            // Make sure we have enough space for all the data
 
#pragma omp parallel num_threads( num_threads < Rows ? num_threads : Rows )
{
    int thread_idx          = omp_get_thread_num();
    int num_used_threads    = omp_get_num_threads(); 
    
    // Determine total number of iterations for the thread
    int iters = (Rows / num_used_threads)    +                      // Each thred gets this many iterations 
                 (thread_idx < (Rows % num_used_threads)    ?       // If row / thread is not an integer
                            static_cast<int>(1)             :       // add iteration if tidx is less then rem
                            static_cast<int>(0)             );      // otherwise don't add another iteration

    for (int it = 0; it < iters; ++it) {
        int row_offset  = num_used_threads * it + thread_idx; 
        int map_offset  = row_offset * (Cols * 2);                  // Add whitespace and newline character
        
        // Put elements from the line into _data vecor
        for (int index = map_offset; index < map_offset + (Cols * 2); index += 2) {
            _data.insert(_data.begin() + (index / 2), *(input_data + index));
            if (thread_idx == 1 ) {
                std::cout << index << " : " << *(input_data + index) << " ";
            }
        }       
    }
#pragma omp barrier         // Wait for all threads to complete
}                           // End omp parallel
    // Debugging - print
    print();
}

template <std::size_t Rows, std::size_t Cols>
void Block<Rows, Cols>::print() const 
{
    for (int r = 0; r < Rows; ++r) {
        for (int c = 0; c < Cols; ++c) {
            std::cout << _data[r * Cols + c].value() << " ";
        }   
        std::cout << "\n";
    }
}

}           // End namespace haplo

#endif      // PARAHAPLO_CPP_BLOCK_HPP
