// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file for the block class for the parahaplo library
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_BLOCK
#define PARAHAPLO_CPP_BLOCK

#include "block_expressions.hpp"

#include <cassert>

namespace phap {
 
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


}           // End namespace phap

#endif      // PARAHAPLO_CPP_BLOCK
