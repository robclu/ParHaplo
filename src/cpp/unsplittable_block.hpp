// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for UnsplittableBlock
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP
#define PARHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP

#include "block.hpp"

#include <tbb/atomic.h>

#include <bitset>

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
private:
    Expression const&   _expression;    //!< The expression that's the base of this class
    container_type      _data;          //!< Data for the unsplittable block 
    size_type           _size;          //!< The size of the unsplittable block (number of elements)
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
    
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      First removes all the singleton rows (rows with only 1 value) from the Unsplittable block
    ///             and then adds the remaining data to the _data container 
    /// @param[in]  index           The index of the unsplittable block inforamtion from Expression to get the
    ///             data from
    /// @param[in]  num_threads     The number of threads to use for function 
    // ------------------------------------------------------------------------------------------------------
    void fill(const size_t index = 0, const size_t num_threads = 1);
    
};

// ---------------------------------------- IMPLEMENTATIONS -------------------------------------------------

template <typename Expression>
UnsplittableBlock<Expression>::UnsplittableBlock(const Expression&  expression  , 
                                                 const size_t       index       ,
                                                 const size_t       num_threads )
: _expression(expression), _data(0)
{
   fill(index, num_threads);                   // Fill the data using the Expression
}

template <typename Expression>
void UnsplittableBlock<Expression>::fill(const size_t index, const size_t num_threads)
{
    // Number of rows in the Expression 
    constexpr size_t rows = Expression::num_rows;
    
    // Create an atomic variable to represent which of the 
    // rows from expression are not singleton and must 
    // therefore be added to the unsplittable element
    std::bitset<rows> unsplittable_rows;
    
    // Create threads where each checks for singletons
    const size_t threads_to_use = num_threads > rows ? rows : num_threads;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices) 
        {
            // Get the information for this unsplittable (sub-)block
            const haplo::SubBlockInfo& info = _expression.subblock_info()[index];
            
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
                    // to the UnsplittableBlock, so set the relevant value in unsplittable_rows
                    if (num_elements > 1) unsplittable_rows[idx + threads_to_use + it] = 1;
                }
            }
        }
    );
    for (int i = 0; i < unsplittable_rows.size(); ++i) std::cout << unsplittable_rows[i] << " : ";
    std::cout << "\n";
}

}       // End namespace haplo

#endif  // PARAHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP
