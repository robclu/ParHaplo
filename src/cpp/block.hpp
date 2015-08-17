// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BLOCK_HPP
#define PARAHAPLO_CPP_BLOCK_HPP

#include "block_expressions.hpp"
#include "block_exceptions.hpp"
#include "../general/input_parser.hpp"

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
/// @tparam T       The type of data used by the block
/// @tparam Rows    The number of rows in the block
/// @tparam Cols    The number of columns in the block
// ----------------------------------------------------------------------------------------------------------
template <typename T, std::size_t Rows, std::size_t Cols>
class Block : public BlockExpression<T, Block<T, Rows, Cols>> {
public:
    // --------------------------------------- Typedefs -----------------------------------------------------
    using container_type    = typename BlockExpression<T, Block<T, Rows, Cols>>::container_type;
    using size_type         = typename BlockExpression<T, Block<T, Rows, Cols>>::size_type;
    using value_type        = typename BlockExpression<T, Block<T, Rows, Cols>>::value_type;
    using reference         = typename BlockExpression<T, Block<T, Rows, Cols>>::reference;
    // ------------------------------------------------------------------------------------------------------
private:
    container_type  _data;          //!< Container for the data - just std::vector<T>(Rows * Cols)
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Defult constructor, intitlizes the data to empty
    // ------------------------------------------------------------------------------------------------------
    Block() : _data(Rows * Cols) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructs a Block from a string which gives the location of the input data. The file is 
    ///             loaded as a memory mapped file - since the input data files will likely be huge, this 
    ///             should save have a significant performance incease.                                      \n
    ///                                                                                                      \n
    ///             See boost memory mapped files for reference:                                             \n
    ///                 http://www.boost.org/doc/libs/1_38_0/libs/iostreams/doc/index.html                   \n
    ///                                                                                                      \n
    ///             Or if you download the boost libraries the source for the memory mapped files is at:     \n
    ///                 boost/iostreams/device/mapped_file.hpp
    /// @param      filename        The name of the file which has the input data
    /// @param      num_elements    The number of elements to read from the file
    // ------------------------------------------------------------------------------------------------------
    Block(const std::string filename, const std::size_t num_elements);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructs a Block from a block expression
    /// @param[in]  expression  The expression to create the Block from
    /// @tparam     Expression  The type of expression 
    // ------------------------------------------------------------------------------------------------------
    template <typename Expression>
    Block(BlockExpression<T, Expression> const& expression) 
    {
        Expression const& expr = expression;
        _data.resize(expr.size());
        for (size_type i = 0; i < expr.size(); ++i) {
            _data[i] = expr[i];                                     // Put data from expr into this block
        }
    }
    
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
    inline void assign_data(const T (&data_source)[N])
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
    /// @brief      Gets input data for a block from a file
    /// @param[in]  filename    The file to get the input data from
    // ------------------------------------------------------------------------------------------------------
};

// ---------------------------------------- Implementations -------------------------------------------------

template <typename T, std::size_t Rows, std::size_t Cols>
Block<T, Rows, Cols>::Block(const std::string filename, std::size_t num_elements) 
{
    using namespace sp;
    using namespace io;
    
    // Create a readonly memory mapped file to get the data
    io::mapped_file_source mapped_input_data(filename.c_str());
    
    _data.reserve(num_elements);                            // Make sure we have enough space for all the data
    
    try {
        if (!qi::parse( mapped_input_data.begin()           ,
                        mapped_input_data.end()             ,
                        sp::ascii::char_ >> sp::qi::space   ,
                        _data                               ) ) {
            throw BlockInputParseError(filename);
        }
    } catch (BlockInputParseError& e) {
        std::cout << e.what() << std::endl;
    }        
}

}           // End namespace haplo

#endif      // PARAHAPLO_CPP_BLOCK_HPP
