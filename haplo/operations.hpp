// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo operations -- general utilities
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_OPERATIONS_HPP
#define PARAHAPLO_OPERATIONS_HPP

namespace haplo {
namespace ops   {

// ----------------------------------------------------------------------------------------------------------
/// @brief      Determines how many iterations each thread must perform in a given situation. Say for
///             example there is a block with 9 rows and we have 4 threads, and each thread does some 
///             operations on a row, then each thread will need to peform 2 iterations, except for one of
///             the threads, which will need to perform 3. So the iteration mapping for the threads would
///             be:                                                                                     \n\n
///             Thread Id  | Rows to operate on | Iterations                                            \n\n
///             0          | 0, 4, 8            | 3
///             1          | 1, 5               | 2
///             2          | 2, 6               | 2
///             3          | 3, 7               | 2
/// @param[in]  thread_id       The thread number 
/// @param[in]  total_ops       The total number of operations (9 in the above example)
/// @param[in]  num_threads     The number of threads being used
/// @return     The total number of iterations for the thread
// ----------------------------------------------------------------------------------------------------------
inline size_t get_thread_iterations(const size_t thread_id  , 
                                    const size_t total_ops  , 
                                    const size_t num_threads)
{
    return (total_ops / num_threads) + 
            (thread_id < (total_ops % num_threads) ? 1 : 0);
}

// ----------------------------------------------------------------------------------------------------------
/// @brief      Determines a mapping when assigning threads to rows or columns -- assuming that all threads in 
///             the group are assigned a row/col of the data and then if the thread performs multiple 
///             iterations it's assigned a thead to operator on data which has an index the number of threads 
///             in the thread group ahead of its previous data. 
/// @param[in]  thread_idx  The index of the thread in the thread group
/// @param[in]  num_threads The total number of threads in the thread group
/// @param[in]  thread_it   The iteration of the thread
/// @return     The row/col index for the given thread id based on its iteration and the number of threads 
///             used
// ----------------------------------------------------------------------------------------------------------
inline size_t thread_map(const size_t thread_idx, const size_t num_threads, const size_t thread_it) 
{
    return thread_it * num_threads + thread_idx;
}

}           // End namespace ops
}           // End namespace haplo

#endif      // PARAHAPLO_BLOCK_OPERATIONS_HPP
