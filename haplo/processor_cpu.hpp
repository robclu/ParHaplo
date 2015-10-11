// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo processor cpu functionality, processes rows/columns 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_PROCESSOR_CPU_HPP
#define PARAHAPLO_PROCESSOR_CPU_HPP

#include "devices.hpp"
#include "operations.hpp"
#include "tree_cpu.hpp"
#include "processor.hpp"

#include <tbb/tbb.h>
#include <iostream>

namespace haplo {

// ------------------------------------------------- ROWS : DUPLICATES  -------------------------------------

template <typename FriendType>
class Processor<FriendType, proc::row_dups, devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using friend_type       = FriendType;
    // ------------------------------------------------------------------------------------------------------
private:
    friend_type& _friend;           //!< The friend class this class has access to to process
    
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the friend class to operate on
    /// @param[in]  friend_class    The class to do the processing for
    // ------------------------------------------------------------------------------------------------------
    Processor(friend_type& friend_class) : _friend(friend_class) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Operator to invoke the processing on the friend class, the processing looks through all
    ///             elements of two rows to check for equivalence
    /// @param      row_idx     The index of the column in the friend class to process
    // ------------------------------------------------------------------------------------------------------
    void operator()(const size_t row_idx);
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Compares two columns, to check if they are equal and if they result in any node links
    /// @param[in]  row_idx_top     The index of the top row in the comparison
    /// @param[in]  row_idx_bot     The index of the bottom row in the comparison
    /// @param[in]  threads_x       The number of threads that can be used to across the rows
    /// @return     If the columns are equal
    // ------------------------------------------------------------------------------------------------------
    bool compare_rows(const size_t   row_idx_top     , 
                      const size_t   row_idx_bot     , 
                      const size_t   threads_x       ); 
};
  
// --------------------------------------- IMPLEMENTATION ---------------------------------------------------

template <typename FriendType>
void Processor<FriendType, proc::row_dups, devices::cpu>::operator()(const size_t row_idx) 
{
    constexpr size_t THX = friend_type::THREADS_X; constexpr size_t THY = friend_type::THREADS_Y;
    
    const size_t threads_y = THY < (_friend._rows - row_idx - 1) 
                           ? THY : (_friend._rows - row_idx - 1);
    
    tbb::atomic<size_t> multiplicity{1};        // The number of rows equal to this row
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y),
        [&](const tbb::blocked_range<size_t>& thread_ids_y) 
        {   
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy                   , 
                                                                   _friend._rows - row_idx - 1 , 
                                                                   threads_y                    ); 
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t row_idx_bot = row_idx + it_y * threads_y + thread_idy + 1;
                    
                    // If the row below this is not a duplicate
                    if (_friend._duplicate_rows.find(row_idx_bot) == _friend._duplicate_rows.end()) {
                        // The number of threads that we can use to process the two columns
                        const size_t threads_x = (THX / threads_y + (THX % threads_y)) < _friend._cols 
                                               ? (THX / threads_y + (THX % threads_y)) : _friend._cols;
                                                    
                        // Check if the columns are duplicates (also determines link values)
                        if (compare_rows(row_idx, row_idx_bot, threads_x) == true) {
                            // Bot is a duplicate if this row
                            _friend._duplicate_rows[row_idx_bot] = row_idx;
                            ++multiplicity;
                        }
                    }
                }
            }
        }
    );
    // Set the multiplicity of the row
    _friend._row_multiplicities[row_idx] = multiplicity;  
}

template <typename FriendType>
bool Processor<FriendType, proc::row_dups, devices::cpu>::compare_rows(const size_t row_idx_top  ,
                                                                       const size_t row_idx_bot  ,
                                                                       const size_t threads_x    )
{
    tbb::atomic<bool> rows_equal{true};
    
    size_t start_col, end_col;
    
    // Look for early exit
    if (_friend._read_info[row_idx_top].start_index() != _friend._read_info[row_idx_bot].start_index()) {
        return false;
    } else { 
        start_col = std::min(_friend._read_info[row_idx_top].start_index(),
                             _friend._read_info[row_idx_bot].start_index());
    }
    if (_friend._read_info[row_idx_top].end_index() != _friend._read_info[row_idx_bot].end_index()) {
        return false;
    } else { 
        end_col = std::max(_friend._read_info[row_idx_top].end_index(),
                           _friend._read_info[row_idx_bot].end_index());
    }    
    const size_t total_cols = end_col - start_col + 1;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, total_cols, threads_x); 
                
                for (size_t it_x = 0; it_x < thread_iters_x && rows_equal; ++it_x) {
                    size_t col_idx = it_x * threads_x + thread_idx + start_col;
                   
                    // If the columns aren't duplicates
                    if (_friend._duplicate_cols.find(col_idx) == _friend._duplicate_cols.end()) {
                        // If the values at these two positions are equivalent
                        auto value_top = _friend(row_idx_top, col_idx);
                        auto value_bot = _friend(row_idx_bot, col_idx);
                        
                        // Check that that all threads terminate here and that it is not a source of error
                        if (value_top != value_bot) {
                            rows_equal = false;
                            break;
                        }
                    }
                }
            }
        }
    );
    return static_cast<bool>(rows_equal); 
}

// ------------------------------- COLUMNS : DUPLICATES AND NODE LINKS  -------------------------------------

template <typename FriendType>
class Processor<FriendType, proc::col_dups_links, devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using friend_type   = FriendType;
    using tree_type     = Tree<friend_type, devices::cpu>;
    // ------------------------------------------------------------------------------------------------------
private:
    friend_type&     _friend;           //!< The friend class this class has access to to process
    tree_type&       _tree;             //!< The tree to initialize 
    
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Contructor -- sets the friend class to operate on
    /// @param[in]  friend_class    The class to do the processing for
    /// @param[in]  tree            The tree to initialize
    // ------------------------------------------------------------------------------------------------------
    Processor(friend_type& friend_class, tree_type& tree) : _friend(friend_class), _tree(tree) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Operator to invoke the processing on the friend class, the processing looks through all
    ///             columns to the right of col_idx and determines equivalence and if the columns are linked
    /// @param[in]  col_idx     The index of the column in the friend class to process
    // ------------------------------------------------------------------------------------------------------
    void operator()(const size_t col_idx);
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Compares two columns, to check if they are equal and if they result in any node links
    /// @param[in]  col_idx_left    The index of the left column
    /// @param[in]  col_idx_right   The index of the right column 
    /// @param[in]  threads_y       The number of threads that can be used to look down the columns
    /// @return     If the columns are equal
    // ------------------------------------------------------------------------------------------------------
    bool compare_columns(const size_t   col_idx_left    , 
                         const size_t   col_idx_right   , 
                         const size_t   threads_y       );      
};

// ---------------------------------------- IMPEMENTATION ---------------------------------------------------

template <typename FriendType>
void Processor<FriendType, proc::col_dups_links, devices::cpu>::operator()(const size_t col_idx)
{
    constexpr size_t THX = friend_type::THREADS_X; constexpr size_t THY = friend_type::THREADS_Y;
    
    const size_t threads_x = THX < (_friend._cols - col_idx - 1) 
                           ? THX : (_friend._cols - col_idx - 1);
    tbb::atomic<size_t> multiplicity{1};        // The number of columns equal tothis column
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx                   , 
                                                                   _friend._cols - col_idx - 1 , 
                                                                   threads_x                    ); 
                
                for (size_t it_x = 0; it_x < thread_iters_x; ++it_x) {
                    size_t col_idx_right = col_idx + it_x * threads_x + thread_idx + 1;
                    
                    // If the column to the right is not a duplicate
                    if (_friend._duplicate_cols.find(col_idx_right) == _friend._duplicate_cols.end()) {
                        // The number of threads that we can use to process the two columns
                        const size_t threads_y = (THY / threads_x + (THY % threads_x)) < _friend._rows 
                                               ? (THY / threads_x + (THY % threads_x)) : _friend._rows;
                        // Check if the columns are duplicates (also determines link values)
                        if (compare_columns(col_idx, col_idx_right, threads_y) == true) {
                            // Right is a duplicate of col_idx
                            _friend._duplicate_cols[col_idx_right] = col_idx;
                            ++multiplicity;
                        }
                    }
                }
            }
        }
    );
    // Set the weight of the node in the tree to be the multiplicity
    _tree.node(col_idx).weight() = multiplicity;
}

template <typename FriendType>
bool Processor<FriendType, proc::col_dups_links, devices::cpu>::compare_columns(const size_t col_idx_left   ,   
                                                                                const size_t col_idx_right  ,
                                                                                const size_t threads_y      )
{
    tbb::atomic<bool> cols_equal{true}; tbb::atomic<bool> link_created{false};
    
    // Find the start row for the comparison 
    size_t start_row = std::min(_friend._snp_info[col_idx_left].start_index(),
                                _friend._snp_info[col_idx_right].start_index());
    
    // Find the end row for the comparison
    size_t end_row = std::max(_friend._snp_info[col_idx_left].end_index(),
                              _friend._snp_info[col_idx_right].end_index());

    // We could have exited early, but then we would not be able to determine the links
    const size_t total_rows = end_row - start_row + 1;
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y),
        [&](const tbb::blocked_range<size_t>& thread_ids_y) 
        {   
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy, total_rows, threads_y); 
              
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t row_idx = it_y * threads_y + thread_idy + start_row;
                    
                    // If the rows aren't duplicates
                    if (_friend._duplicate_rows.find(row_idx) == _friend._duplicate_rows.end()) {
                        // If the values at these two positions are equivalent
                        auto value_left   = _friend(row_idx, col_idx_left);
                        auto value_right  = _friend(row_idx, col_idx_right);
                        
                        if (value_left == value_right) {
                            // Could be a duplicate, check if the value are valid for a link
                            if (value_left <= 1 && value_right <= 1) {
                                if (!link_created) {
                                    _tree.create_link(col_idx_left, col_idx_right);
                                    link_created = true;
                                }
                                // Optimal if the values are the same
                                _tree.link(col_idx_left, col_idx_right).homo_weight()
                                += _friend._row_multiplicities[row_idx];
                            }
                        } else if (value_left != value_right) {
                            //  Definitely can't be a duplicate 
                            cols_equal = false;
                            // Check if the values are valid for a link
                            if (value_left <= 1 && value_right <= 1) {
                                if (!link_created) {
                                    _tree.create_link(col_idx_left, col_idx_right);
                                    link_created = true;
                                }
                                // Optimal if the values are opposite
                                _tree.link(col_idx_left, col_idx_right).hetro_weight()
                                += _friend._row_multiplicities[row_idx];
                            }                            
                        }
                    }
                }
            }
        }
    );
    if (!cols_equal) {
        // Determine the contribution to the worst case value for these columns
        const size_t worst_case_value = _tree.link_max(col_idx_left, col_idx_right) -
                                        _tree.link_min(col_idx_left, col_idx_right);
        
        // Add to their worst case values
        _tree.node(col_idx_left).worst_case_value().fetch_and_add(worst_case_value * 
                                                                  _tree.node(col_idx_right).weight());
        _tree.node(col_idx_right).worst_case_value().fetch_and_add(worst_case_value);
        
        // Check if the right node is the global worst case
        if (_tree.node(col_idx_right).worst_case_value() * _tree.node(col_idx_right).weight() > _tree.max_worst_case()) {
            _tree.max_worst_case()  = _tree.node(col_idx_right).worst_case_value() * _tree.node(col_idx_right).weight();
            _tree.start_node()      = col_idx_right;
        }
    }
    return static_cast<bool>(cols_equal);
}

}               // End namespace haplo
#endif          // PARAHAPLO_PROCESSOR_CPU_HPP
