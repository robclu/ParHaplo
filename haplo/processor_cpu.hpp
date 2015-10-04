// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo processor cpu functionality, processes rows/columns 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_PROCESSOR_CPU_HPP
#define PARAHAPLO_PROCESSOR_CPU_HPP

#include "devices.hpp"
#include "node_container_cpu.hpp"
#include "operations.hpp"
#include "processor.hpp"

#include <tbb/tbb.h>

#include <iostream>

namespace haplo {

// ------------------------------------------------- ROWS : DUPLICATES  -------------------------------------

template <typename FriendType>
class Processor<FriendType, proc::ROW_DUP, devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using friend_type       = FriendType;
    // ------------------------------------------------------------------------------------------------------
private:
    friend_type* _friend;           //!< The friend class this class has access to to process
    
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the friend class to operate on
    /// @param[in]  friend_class    The class to do the processing for
    // ------------------------------------------------------------------------------------------------------
    Processor(FriendType* friend_class) : _friend(friend_class) {}
    
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
void Processor<FriendType, proc::ROW_DUP, devices::cpu>::operator()(const size_t row_idx) 
{
    constexpr size_t THI = friend_type::THREADS_I; constexpr size_t THJ = friend_type::THREADS_J;
    
    const size_t threads_y = THI < (_friend->_rows - row_idx - 1) 
                           ? THI : (_friend->_rows - row_idx - 1);
    
    tbb::atomic<size_t> multiplicity{1};        // The number of rows equal to this row
    
    std::cout << threads_y << " R : " << row_idx << "\n";
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y),
        [&](const tbb::blocked_range<size_t>& thread_ids_y) 
        {   
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy                   , 
                                                                   _friend->_rows - row_idx - 1 , 
                                                                   threads_y                    ); 
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t row_idx_bot = row_idx + it_y * threads_y + thread_idy + 1;
                    
                    // If the row below this is not a duplicate
                    if (_friend->_duplicate_rows.find(row_idx_bot) == _friend->_duplicate_rows.end()) {
                        // The number of threads that we can use to process the two columns
                        const size_t threads_x = (THJ / threads_y + (THJ % threads_y)) < _friend->_cols 
                                               ? (THJ / threads_y + (THJ % threads_y)) : _friend->_cols;
                                                    
                        // Check if the columns are duplicates (also determines link values)
                        if (compare_rows(row_idx, row_idx_bot, threads_x) == true) {
                            // Bot is a duplicate if this row
                            _friend->_duplicate_rows[row_idx_bot] = row_idx;
                            ++multiplicity;
                        }
                    }
                }
            }
        }
    );
    // Set the multiplicity of the row
    _friend->_row_multiplicities[row_idx] = multiplicity;  
}

template <typename FriendType>
bool Processor<FriendType, proc::ROW_DUP, devices::cpu>::compare_rows(const size_t row_idx_top  ,
                                                                      const size_t row_idx_bot  ,
                                                                      const size_t threads_x    )
{
    tbb::atomic<bool> rows_equal{true};
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx, _friend->_cols, threads_x); 
                
                for (size_t it_x = 0; it_x < thread_iters_x && rows_equal; ++it_x) {
                    size_t col_idx = it_x * threads_x + thread_idx;
                   
                    // If the columns aren't duplicates
                    if (_friend->_duplicate_cols.find(col_idx) == _friend->_duplicate_cols.end()) {
                        // If the values at these two positions are equivalent
                        const size_t offset_top = row_idx_top * _friend->_cols + col_idx;
                        const size_t offset_bot = row_idx_bot * _friend->_cols + col_idx;
                        
                        if (_friend->_data.get(offset_top) != _friend->_data.get(offset_bot)) {
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
class Processor<FriendType, proc::COL_DUP_LINKS, devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using node_container    = NodeContainer<devices::cpu>;
    using friend_type       = FriendType;
    // ------------------------------------------------------------------------------------------------------
private:
    friend_type* _friend;           //!< The friend class this class has access to to process
    
public:    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Contructor -- sets the friend class to operate on
    /// @param[in]  friend_class    The class to do the processing for
    // ------------------------------------------------------------------------------------------------------
    Processor(FriendType* friend_class) : _friend(friend_class) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Operator to invoke the processing on the friend class, the processing looks through all
    ///             columns to the right of col_idx and determines equivalence and if the columns are linked
    /// @param      col_idx     The index of the column in the friend class to process
    /// @param      nodes       The nodes to set the parameters for
    // ------------------------------------------------------------------------------------------------------
    void operator()(const size_t col_idx, node_container& nodes);
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Compares two columns, to check if they are equal and if they result in any node links
    /// @param[in]  col_idx_left    The index of the left column
    /// @param[in]  col_idx_right   The index of the right column 
    /// @param[in]  threads_y       The number of threads that can be used to look down the columns
    /// @param[in]  nodes           The nodes to set the links for, if any valid links
    /// @return     If the columns are equal
    // ------------------------------------------------------------------------------------------------------
    bool compare_columns(const size_t       col_idx_left    , 
                         const size_t       col_idx_right   , 
                         const size_t       threads_y       ,
                         node_container&    nodes           ); 
};

// ---------------------------------------- IMPEMENTATION ---------------------------------------------------

template <typename FriendType>
void Processor<FriendType, proc::COL_DUP_LINKS, devices::cpu>::operator()(const size_t      col_idx     ,
                                                                          node_container&   nodes       )
{
    constexpr size_t THI = friend_type::THREADS_I; constexpr size_t THJ = friend_type::THREADS_J;
    
    const size_t threads_x = THJ < (_friend->_cols - col_idx - 1) 
                           ? THJ : (_friend->_cols - col_idx - 1);
    
    tbb::atomic<size_t> multiplicity{1};        // The number of columns equal tothis column
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_x),
        [&](const tbb::blocked_range<size_t>& thread_ids_x) 
        {   
            for (size_t thread_idx = thread_ids_x.begin(); thread_idx != thread_ids_x.end(); ++thread_idx) {
                size_t thread_iters_x = ops::get_thread_iterations(thread_idx                   , 
                                                                   _friend->_cols - col_idx - 1 , 
                                                                   threads_x                    ); 
                
                for (size_t it_x = 0; it_x < thread_iters_x; ++it_x) {
                    size_t col_idx_right = col_idx + it_x * threads_x + thread_idx + 1;
                    
                    // If the column to the right is not a duplicate
                    if (_friend->_duplicate_cols.find(col_idx_right) == _friend->_duplicate_cols.end()) {
                        // The number of threads that we can use to process the two columns
                        const size_t threads_y = (THI / threads_x + (THI % threads_x)) < _friend->_rows 
                                               ? (THI / threads_x + (THI % threads_x)) : _friend->_rows;
                                                    
                        // Check if the columns are duplicates (also determines link values)
                        if (compare_columns(col_idx, col_idx_right, threads_y, nodes) == true) {
                            // Right is a duplicate of col_idx
                            _friend->_duplicate_cols[col_idx_right] = col_idx;
                            ++multiplicity;
                        }
                    }
                }
            }
        }
    );
    // Set the multiplicity of the column
    _friend->_col_multiplicities[col_idx] = multiplicity; 
}

template <typename FriendType>
bool Processor<FriendType, proc::COL_DUP_LINKS, devices::cpu>::compare_columns(const size_t     col_idx_left    ,   
                                                                               const size_t     col_idx_right   ,
                                                                               const size_t     threads_y       , 
                                                                               node_container&  nodes           )
{
    tbb::atomic<bool> cols_equal{true};
    
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads_y),
        [&](const tbb::blocked_range<size_t>& thread_ids_y) 
        {   
            for (size_t thread_idy = thread_ids_y.begin(); thread_idy != thread_ids_y.end(); ++thread_idy) {
                size_t thread_iters_y = ops::get_thread_iterations(thread_idy, _friend->_rows, threads_y); 
                
                for (size_t it_y = 0; it_y < thread_iters_y; ++it_y) {
                    size_t row_idx = it_y * threads_y + thread_idy;
                    
                    // If the rows aren't duplicates
                    if (_friend->_duplicate_rows.find(row_idx) == _friend->_duplicate_rows.end()) {
                       // If the values at these two positions are equivalent
                        const size_t offset_left  = row_idx * _friend->_cols + col_idx_left;
                        const size_t offset_right = row_idx * _friend->_cols + col_idx_right;
                        
                        if (_friend->_data.get(offset_left) == _friend->_data.get(offset_right)) {
                            // Could be a duplicate, check if the value are valid for a link
                            if (_friend->_data.get(offset_left) <= 1 && _friend->_data.get(offset_right) <= 1) {
                                // Optimal if the values are the same
                                nodes.link(col_idx_left, col_idx_right)._homo_weight += 
                                                                    _friend->_row_multiplicities[row_idx];
                            }
                        } else if (_friend->_data.get(offset_left) != _friend->_data.get(offset_right)) {
                            //  Definitely can't be a duplicate 
                            cols_equal = false;
                            // Check if the values are valid for a link
                            if (_friend->_data.get(offset_left) <= 1 && _friend->_data.get(offset_right) <= 1) {
                                // Optimal if the values are opposite
                                nodes.link(col_idx_left, col_idx_right)._hetro_weight += 
                                                                    _friend->_row_multiplicities[row_idx];
                            }                            
                        }
                    }
                }
            }
        }
    );
    return static_cast<bool>(cols_equal);
}
}               // End namespace haplo
#endif          // PARAHAPLO_PROCESSOR_CPU_HPP
