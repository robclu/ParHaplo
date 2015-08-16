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
    /// @brief  Gets the size of the block (total number of elements)
    /// @return The total number of elements in the block
    // ------------------------------------------------------------------------------------------------------
    size_type size() const { return _data.size(); }
};

}           // End namespace phap

#endif      // PARAHAPLO_CPP_BLOCK
