// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for UnsplittableBlock
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP
#define PARHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP

#include "block.hpp"

#include <tbb/atomic.h>

#include <bitset>
#include <numeric>

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class  UnsplittableBlock
/// @brief  A block which cannot be split
/// @tparam Expression  The expression (Block for example) which is the base of the UnsplittableBlock
// ----------------------------------------------------------------------------------------------------------
template <typename Expression>
class UnsplittableBlock : public BlockExpression<UnsplittableBlock<Expression>> {
public:
    // -------------------------------------- Typedefs ------------------------------------------------------
    using container_type    = std::vector<haplo::Data>;
    using size_type         = typename container_type::size_type;
    using value_type        = typename container_type::value_type;           
    using reference         = typename container_type::reference;
    // ------------------------------------------------------------------------------------------------------
    // Define a type alias for a container that keeps info about the singleton rows
    using singleton_container   = std::array<bool, Expression::num_rows>;
    using data_container        = std::vector<uint8_t>;
private:
    Expression const&   _expression;        //!< The expression that's the base of this class
    singleton_container _singleton_info;    //!< Array of bools for which rows of _expression are singletons
                                            //!< 1 = not singleton (so we can count the number of not singles)
    data_container      _data;              //!< Data for the unsplittable block 
    size_type           _size;              //!< The size of the unsplittable block (number of elements)
    size_type           _cols;              //!< The number of columns
    size_type           _rows;              //!< The number of rows 
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the expression, and the data to 0
    /// @param[in]  expression  The base expression from which this class derives -- for example a Block, this
    ///             Unsplittable block is then a sub-block of that block.
    /// @param[in]  index       The index of the unsplittable-block from the base expression -- for example,
    ///             a Block may be splittable into 3 unsplittable (smaller) blocks, so there will be 3 indices
    ///             (0, 1, 2)
    /// @param[in]  num_threads The number of threads to use for all functions for the class
    // ------------------------------------------------------------------------------------------------------
    UnsplittableBlock(const Expression& expression      , 
                      const size_t      index       = 0 ,
                      const size_t      num_threads = 1 );
    
    // ------------------------------------------------------------------------------------------------------
    //! @brief      Returns the size of the Expression
    //! @return     The size of the Expression
    // ------------------------------------------------------------------------------------------------------
    size_type size() const { return _size; }

    // ------------------------------------------------------------------------------------------------------
    //! @brief      Gets an element from the data
    //! @param[in]  i   The element to fetch.
    //! @return     The value of the element at position i of the data.
    // ------------------------------------------------------------------------------------------------------
    value_type operator[](size_type i) const { _data[i]; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the number of rows in the UnsplittableBlock 
    /// @return     The number of rows in the UnsplittableBlock
    // ------------------------------------------------------------------------------------------------------
    size_type rows() const { return _rows; } 
   // ------------------------------------------------------------------------------------------------------ 
   /// @brief      Gets the number of rows in the UnsplittableBlock                                          // ------------------------------------------------------------------------------------------------------
   /// @return     The number of rows in the UnsplittableBlock                                               /// @brief      Prints the UnsplittableBlock
   // ------------------------------------------------------------------------------------------------------ // ------------------------------------------------------------------------------------------------------
    void print() const;
    
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      First removes all the singleton rows (rows with only 1 value) from the Unsplittable block
    ///             from which the number of rows in the unsplittable block, and hence the total size of the
    ///             block, can be determined.
    /// @param[in]  index           The index of the unsplittable block inforamtion from Expression to get the
    ///             data from
    /// @param[in]  num_threads     The number of threads to use for function 
    // ------------------------------------------------------------------------------------------------------
    void determine_params(const size_t index = 0, const size_t num_threads = 1);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Uses the determined parameters to fill the data for the sub-block from the _expression
    /// @param[in]  index       The index of the unsplittable block in the Expression
    /// @param[in]  num_threads The number of threads to use for filling the data
    // ------------------------------------------------------------------------------------------------------
    void fill(const size_t index = 0, const size_t num_threads = 1);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Maps rows from the unsplittable block to the Expression. Say for example the Expression
    ///             has 10 rows, and 5 are singular (2, 3, 6, 8, 9) and must be removed, so the _singletons
    ///             array will look like :                                                                  \n\n
    ///             [1, 1, 0, 0, 1, 1, 0, 1, 0, 0]                                                          \n\n
    ///             So then the mapping is                                                                  \n\n
    /// 
    ///             | UnsplittableRow | Expression Row  | (0 indexing)                                      \n
    ///             |       0         |         0       |                                                   \n
    ///             |       1         |         1       |                                                   \n
    ///             |       2         |         4       |                                                   \n
    ///             |       3         |         5       |                                                   \n
    ///             |       4         |         7       |                                                   \n 
    /// 
    ///             The functions maps an UnSplittable row to an expression row
    /// @param[in]  u_row   The unspittable row
    /// @return     The expresison row
    // ------------------------------------------------------------------------------------------------------
    int row_map(const size_t urow);
};

// ---------------------------------------- IMPLEMENTATIONS -------------------------------------------------

// ------------------------------------------- PUBLIC -------------------------------------------------------

template <typename Expression>
UnsplittableBlock<Expression>::UnsplittableBlock(const Expression&  expression  , 
                                                 const size_t       index       ,
                                                 const size_t       num_threads )
: _expression(expression)
{
   determine_params(index, num_threads);                    
   fill(index, num_threads);
}

template <typename Expression>
void UnsplittableBlock<Expression>::print() const 
{
    for (int r = 0; r < _rows; ++r) {
        for (int c = 0; c < _cols; ++c) 
            std::cout << static_cast<unsigned int>(_data[r * _cols + c]) << " ";
        std::cout << "\n";
    }
}
// ------------------------------------------ PRIVATE -------------------------------------------------------

template <typename Expression>
void UnsplittableBlock<Expression>::determine_params(const size_t index, const size_t num_threads)
{
    // Number of rows in the Expression 
    constexpr size_t rows = Expression::num_rows;
    
    // Create threads where each checks for singletons
    const size_t threads_to_use = num_threads > rows ? rows : num_threads;

    // Information for the sub-block
    const haplo::SubBlockInfo& info = _expression.subblock_info()[index];
        
    // Check which rows are singletons
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices) 
        {
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                size_t thread_iters = util::get_thread_iterations(idx, rows, threads_to_use);
                
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t num_elements = 0;
                    // Now we need to go through all elements in the row 
                    // of Expression and check of there is only 1 element
                    for (int col = info.start(); col <= info.end() && num_elements < 2; ++col) {
                        if (_expression(it * threads_to_use + idx, col) != 2) 
                            ++num_elements;         // Not a gap, so increment num_elements
                    }
                    // If we have found more than a single element then the row needs to be added
                    // to the UnsplittableBlock, so set the relevant value 1 (not singletons)
                    _singleton_info[it * threads_to_use + idx] = num_elements > 1 ? 1 : 0;
                }
            }
        }
    );
    
    // This acclally counts the number of rows which aren't singletons
    _rows = std::accumulate(_singleton_info.begin(), _singleton_info.end(), 0);
    _cols = info.columns();
    _size = _rows * info.columns();
    _data.resize(_size, 0);    
}

template <typename Expression>
void UnsplittableBlock<Expression>::fill(const size_t index, const size_t num_threads)
{
    // Check that we aren't using too many threads
    const size_t threads_to_use = num_threads > _rows ? _rows : num_threads;
 
    // Reference to the sub-block info for this unsplittable block
    const haplo::SubBlockInfo& info = _expression.subblock_info()[index];
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices)
        {
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                size_t thread_iters = util::get_thread_iterations(idx, _rows, threads_to_use);
                
                for (size_t it = 0; it < thread_iters; ++it) {
                    int expression_row = row_map(it * threads_to_use + idx);
                    for (size_t col = info.start(); col <= info.end(); ++col) {
                        _data[(it * threads_to_use + idx) * _cols   +               // Row offset
                              (col - info.start())                  ]               // Column offset
                        = _expression(expression_row, col);                         // Value
                    }                  
                }               
            }
        }
    );
}

template <typename Expression>
int UnsplittableBlock<Expression>::row_map(const size_t unsplittable_row)
{
    int counter = -1, index = 0;
    while (counter != unsplittable_row) {
        if (_singleton_info[index++] == 1) ++counter;
    }
    return --index;
}

}       // End namespace haplo

#endif  // PARAHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP
