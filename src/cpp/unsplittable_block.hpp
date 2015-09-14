// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for UnsplittableBlock
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP
#define PARHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP

#include "block.hpp"
#include "equality_checker.hpp"

#include <tbb/concurrent_unordered_map.h>

#include <bitset>
#include <numeric>
#include <iomanip>
#include <ios>

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
    using container_type    = typename Expression::container_type;
    using size_type         = typename container_type::size_type;
    using value_type        = typename container_type::value_type;           
    using reference         = typename container_type::reference;
    // ------------------------------------------------------------------------------------------------------
    // Define a type alias for a container that keeps info about the singleton rows
    using singleton_container   = std::array<bool, Expression::num_rows>;
    // Define a type alias for a container for a concurrent unordered_map for 
    // storing duplicate rows and columns, and the multiplicities
    using concurrent_map        = tbb::concurrent_unordered_map<uint16_t, uint16_t>;
    
private:
    Expression const&   _expression;        //!< The expression that's the base of this class
    singleton_container _singleton_info;    //!< Array of bools for which rows of _expression are singletons
                                            //!< 1 = not singleton (so we can count the number of not singles)
    concurrent_map      _row_mplicity;      //!< Multiplicities of each of the rows
    concurrent_map      _col_mplicity;      //!< Multiplicities of each of the columns
    container_type      _data;              //!< Data for the unsplittable block 
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
    /// @brief      Gets the number of columns in the UnsplittableBlock 
    /// @return     The number of columns in the UnsplittableBlock
    // ------------------------------------------------------------------------------------------------------
    size_type columns() const { return _cols; }  

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the multiplicity of row i, if it has a multiplicity, otherwise returns 0
    /// @return     The multiplicity of row i
    // ------------------------------------------------------------------------------------------------------
    const uint16_t& row_multiplicity(const uint16_t i) const 
    { 
        return _row_mplicity.find(i) != _row_mplicity.end() ? _row_mplicity.at(i) : 0; 
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the multiplicity of column i, if it has a multiplicity, otherwise returns 0
    /// @return     The multiplicity of column i
    // ------------------------------------------------------------------------------------------------------
    const uint16_t& col_multiplicity(const uint16_t i) const 
    { 
        return _col_mplicity.find(i) != _col_mplicity.end() ? _col_mplicity.at(i) : 0; 
    }
    
    // -------------------------------------- DEBUGGING FUNCTIONS -------------------------------------------
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Prints the UnsplittableBlock
    // ------------------------------------------------------------------------------------------------------
    void print() const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Prints the multiplicities of the non-duplicate rows in the unsplittable block
    // ------------------------------------------------------------------------------------------------------
    void print_row_multiplicities() const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Prints the multiplicities of the non-duplicate columns in the unsplittable block
    // ------------------------------------------------------------------------------------------------------
    void print_col_multiplicities() const;
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Removes all duplicate columns or rows
    /// @param[in]  num_threads     The number of threads to use
    /// @tparam     Checker         The fuctor which checks either row or column equivalence
    // ------------------------------------------------------------------------------------------------------
    template <typename Checker>
    void remove_duplicates(const size_t num_threads, Checker checker);
    
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
    // Create row and column equivalence checker functors
    EqualityChecker<check::ROWS>    row_checker;
    EqualityChecker<check::COLUMNS> col_checker;
    
    determine_params(index, num_threads);           // Set the internal parameters             
    fill(index, num_threads);                       // Get the data
    remove_duplicates(num_threads, row_checker);    // Remove all duplicate rows
    remove_duplicates(num_threads, col_checker);    // Remove all duplicate columns
}

template <typename Expression>
void UnsplittableBlock<Expression>::print() const 
{
    for (int r = 0; r < _rows; ++r) {
        for (int c = 0; c < _cols; ++c) 
            std::cout << static_cast<unsigned int>(_data[r * _cols + c].value()) << " ";
        std::cout << "\n";
    }
}

template <typename Expression>
void UnsplittableBlock<Expression>::print_row_multiplicities() const 
{
    std::cout << "\nPrinting row multiplicities for unsplittable block:\n";
    std::cout << std::left << std::setw(12) << "Row ID" << std::setw(12) << "Multiplicity\n";
    std::cout << "------------------------\n";
    for (size_t i = 0; i < _rows; ++i) {
        if (_row_mplicity.find(i) != _row_mplicity.end())
            std::cout   << std::left        <<  std::setw(12)                                   << i  
                        << std::setw(12)    <<  static_cast<unsigned int>(_row_mplicity.at(i))  << "\n";
    }
}

template <typename Expression>
void UnsplittableBlock<Expression>::print_col_multiplicities() const 
{
    std::cout << "\nPrinting col multiplicities for unsplittable block:\n";
    std::cout << std::left << std::setw(12) << "Col ID" << std::setw(12) << "Multiplicity\n";
    std::cout << "------------------------\n";
    for (size_t i = 0; i < _cols; ++i) {
        if (_col_mplicity.find(i) != _col_mplicity.end())
            std::cout   << std::left        <<  std::setw(12)                                   << i  
                        << std::setw(12)    <<  static_cast<unsigned int>(_col_mplicity.at(i))  << "\n";
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
                size_t thread_iters = ops::get_thread_iterations(idx, rows, threads_to_use);
                
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t num_elements = 0;
                    // Now we need to go through all elements in the row 
                    // of Expression and check of there is only 1 element
                    for (int col = info.start(); col <= info.end() && num_elements < 2; ++col) {
                        if (_expression(it * threads_to_use + idx, col).value() != 2) 
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
    _data.resize(_size, static_cast<uint8_t>(0));    
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
                size_t thread_iters = ops::get_thread_iterations(idx, _rows, threads_to_use);
                
                for (size_t it = 0; it < thread_iters; ++it) {
                    int expression_row = row_map(it * threads_to_use + idx);
                    for (size_t col = info.start(); col <= info.end(); ++col) {
                        _data[(it * threads_to_use + idx) * _cols   +               // Row offset
                              (col - info.start())                  ]               // Column offset
                        = _expression(expression_row, col).value();                 // Value
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
    
template <typename Expression> template <typename Checker>
void UnsplittableBlock<Expression>::remove_duplicates(const size_t num_threads, Checker checker)
{
    // This can check either row or column equivalence, the following 
    // variables are defined:
    // elements : This is the total number of elements for comparison, so
    //            if rows are beinf checked then elements = number of columns
    //            since for 2 rows, each corresponding column must be check, 
    //            and vice versa for columns
    // indexes  : This is the total number of rows or columns (an index) to check,
    //            if columns are being checked then indexes = _cols and if rows
    //            are being checked the  indexes = _rows

    // Create an unordered map and use the index of a duplicate column
    // as the key, and the column of which it is a duplicate as the value
    concurrent_map duplicates;
   
    // Set the number of elements for each comparison, the total number of comparisons
    // (indexes), stride and the total number of threads to use depending on if column
    // or row duplicates are being checked
    const size_t elements       = Checker::type == check::ROWS ? _cols : _rows;
    const size_t indexes        = Checker::type == check::ROWS ? _rows : _cols;
    const size_t stride         = _cols;    
    const size_t threads_to_use = num_threads > indexes ? indexes : num_threads;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_to_use),
        [&](const tbb::blocked_range<size_t>& thread_indices)
        {
            for (size_t idx = thread_indices.begin(); idx != thread_indices.end(); ++idx) {
                size_t thread_iters = ops::get_thread_iterations(idx, indexes, threads_to_use);
            
                for (uint16_t it = 0; it < thread_iters; ++it) {
                    uint16_t current_idx                = ops::thread_map(idx, threads_to_use, it);
                    uint16_t current_idx_multiplicity   = 1;
                    // Check if this thread is a duplicate (already been found by another thread)
                    if (duplicates.find(current_idx) == duplicates.end()) {
                        // Not a duplicate so we must look for duplicates for each of the remaining rows/cols
                        for (uint8_t other_idx = current_idx + 1; other_idx < indexes; ++other_idx) {
                            // Check if the rows/columns are equal
                            if (checker(&_data[0], current_idx, other_idx, elements, stride)) {
                                // Add that other_idx (a row or column ahead of current_idx) 
                                // is a duplicate of current_idx, and use other_idx as a key
                                // so that when it comes to that idx, no work is done
                                duplicates[other_idx] = current_idx;
                                // Add to the multiplicity of the current_idx as duplicate found.
                                // Use the temp variable since it's almost definitely faster to write
                                // to the concurrent container once rather than on each iteration
                                ++current_idx_multiplicity;
                            }
                        }
                        // Write the multiplicity of the idx to the concurrent row or column container
                        Checker::type == check::ROWS 
                                       ? _row_mplicity[current_idx] = current_idx_multiplicity
                                       : _col_mplicity[current_idx] = current_idx_multiplicity;
                    }
                }
            }
        }
    );
}

}       // End namespace haplo

#endif  // PARAHAPLO_CPP_UNSPLITTABLE_BLOCK_HPP
