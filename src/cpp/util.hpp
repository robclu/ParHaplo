// ----------------------------------------------------------------------------------------------------------
/// @file   util.hpp
/// @brief  Header file for utility functions for parahaplo
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_UTIL_HPP
#define PARAHAPLO_CPP_UTIL_HPP

namespace haplo {
namespace util  {

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
/// @brief      Determines the row when assigning threads to rows -- assuming that all threads in the group
///             are assigned a row of the dat and then if the thread performs multiple iterations its assigned
///             a thead to operator on data which is the number of threads in the group ahead of its previous
///             data. 
/// @param[in]  thread_idx  The index of the thread in the thread group
/// @param[in]  num_threds  The total number of threads in the thread group
/// @param[in]  thread_it   The iteration of the thread
/// @return     The rows index for the given thread id based on its iteration and the number of threads used
// ----------------------------------------------------------------------------------------------------------
inline size_t thread_row(const size_t thread_idx, const size_t num_threads, const size_t thread_it) 
{
    return thread_it * num_threads + thread_idx;
}

}           // End namespace util
}           // End namespace haplo

#endif      // PARAHAPLO_CPP_UTIL_HPP
