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
    constexpr size_t THI = friend_type::THREADS_I; constexpr size_t THJ = friend_type::THREADS_J;
    
    const size_t threads_y = THI < (_friend._rows - row_idx - 1) 
                           ? THI : (_friend._rows - row_idx - 1);
    
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
                        const size_t threads_x = (THJ / threads_y + (THJ % threads_y)) < _friend._cols 
                                               ? (THJ / threads_y + (THJ % threads_y)) : _friend._cols;
                                                    
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
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, _friend._cols, threads_x); 
                
                for (size_t it_x = 0; it_x < thread_iters_x && rows_equal; ++it_x) {
                    size_t col_idx = it_x * threads_x + thread_idx;
                   
                    // If the columns aren't duplicates
                    if (_friend._duplicate_cols.find(col_idx) == _friend._duplicate_cols.end()) {
                        // If the values at these two positions are equivalent
                        const size_t offset_top = row_idx_top * _friend._cols + col_idx;
                        const size_t offset_bot = row_idx_bot * _friend._cols + col_idx;
                        
                        if (_friend._data.get(offset_top) != _friend._data.get(offset_bot)) {
                            // Not duplicate rows, we can return
                            rows_equal = false;
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
    using tree_type     = Tree<devices::cpu>;
    using friend_type   = FriendType;
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
    constexpr size_t THI = friend_type::THREADS_I; constexpr size_t THJ = friend_type::THREADS_J;
    
    const size_t threads_x = THJ < (_friend._cols - col_idx - 1) 
                           ? THJ : (_friend._cols - col_idx - 1);
    
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
                        const size_t threads_y = (THI / threads_x + (THI % threads_x)) < _friend._rows 
                                               ? (THI / threads_x + (THI % threads_x)) : _friend._rows;
                                                    
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
    _tree.node_weight(col_idx) = multiplicity;
    // Set the multiplicity of the column
    _friend._col_multiplicities[col_idx] = multiplicity; 
}

template <typename FriendType>
bool Processor<FriendType, proc::col_dups_links, devices::cpu>::compare_columns(const size_t col_idx_left   ,   
                                                                                const size_t col_idx_right  ,
                                                                                const size_t threads_y      )
{
    tbb::atomic<bool>   cols_equal{true};
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y),
        [&](const tbb::blocked_range<size_t>& thread_ids_y) 
        {   
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy, _friend._rows, threads_y); 
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t row_idx = it_y * threads_y + thread_idy;
                    
                    // If the rows aren't duplicates
                    if (_friend._duplicate_rows.find(row_idx) == _friend._duplicate_rows.end()) {
                       // If the values at these two positions are equivalent
                        const size_t offset_left  = row_idx * _friend._cols + col_idx_left;
                        const size_t offset_right = row_idx * _friend._cols + col_idx_right;
                        
                        if (_friend._data.get(offset_left) == _friend._data.get(offset_right)) {
                            // Could be a duplicate, check if the value are valid for a link
                            if (_friend._data.get(offset_left) <= 1 && _friend._data.get(offset_right) <= 1) {
                                // Optimal if the values are the same
                                _tree.link<links::homo>(col_idx_left, col_idx_right)
                                    += _friend._row_multiplicities[row_idx];
                            }
                        } else if (_friend._data.get(offset_left) != _friend._data.get(offset_right)) {
                            //  Definitely can't be a duplicate 
                            cols_equal = false;
                            // Check if the values are valid for a link
                            if (_friend._data.get(offset_left) <= 1 && _friend._data.get(offset_right) <= 1) {
                                // Optimal if the values are opposite
                                _tree.link<links::hetro>(col_idx_left, col_idx_right)
                                    += _friend._row_multiplicities[row_idx];
                            }                            
                        }
                    }
                }
            }
        }
    );
    // Determine the contribution to the worst case value for these columns
    const size_t worst_case_value = _tree.link_max(col_idx_left, col_idx_right) -       // Max of the links
        std::min(_tree.link<links::homo>(col_idx_left, col_idx_right) ,                 // Min of the links
                 _tree.link<links::hetro>(col_idx_left, col_idx_right));
    
    // Add to their worst case values
    _tree.node_worst_case(col_idx_left)  += worst_case_value;
    _tree.node_worst_case(col_idx_right) += worst_case_value;

    return static_cast<bool>(cols_equal);
}

// ------------------------------------- COLUMNS : REMOVE MONTONES ------------------------------------------

template <typename FriendType>
class Processor<FriendType, proc::col_rem_mono, devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using tree_type     = Tree<devices::cpu>;
    using friend_type   = FriendType;
    // ------------------------------------------------------------------------------------------------------
private:
    friend_type&     _friend;           //!< The friend class this class has access to to process
    
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Contructor -- sets the friend class to operate on
    /// @param[in]  friend_class    The class to do the processing for
    // ------------------------------------------------------------------------------------------------------
    Processor(friend_type& friend_class) : _friend(friend_class) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Operator to invoke the processing on the friend class, the processing looks through all
    ///             the monotone columns and removes them
    /// @param[in]  col_idx     The index of the column in the friend class to process
    // ------------------------------------------------------------------------------------------------------
    void operator()();
};

// ------------------------------------------ IMPLEMENTATION ------------------------------------------------

template <typename FriendType>
void Processor<FriendType, proc::col_rem_mono, devices::cpu>::operator()() 
{
    const size_t cols_before            = _friend._cols;
    size_t       elements               = _friend._cols * _friend._rows;
    size_t       cols_removed           = 0;
    size_t       col_idx_with_removal   = 0;
    
    // Unroll the column deteltion process
    while (col_idx_with_removal < elements) {
        size_t col_idx_without_removal = (col_idx_with_removal + cols_removed) % cols_before;
        
        if (_friend._monotone_cols.find(col_idx_without_removal) != _friend._monotone_cols.end() ) {
            // Montone column so remove it 
            _friend._data.remove_element(col_idx_with_removal);
            ++cols_removed;
            --elements;
        } else {
            ++col_idx_with_removal;     // Just move to the next element
        }
    }
    _friend._cols -= _friend._monotone_cols.size();
}

}               // End namespace haplo
#endif          // PARAHAPLO_PROCESSOR_CPU_HPP
