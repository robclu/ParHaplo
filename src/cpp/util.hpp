// ----------------------------------------------------------------------------------------------------------
/// @file   util.hpp
/// @brief  Header file for utility functions for parahaplo
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_UTIL_HPP
#define PARAHAPLO_CPP_UTIL_HPP

namespace haplo {
namespace util  {

// ------------------------------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------------------------------
size_t get_thread_iterations(const size_t thread_id  , 
                             const size_t total_ops  , 
                             const size_t num_threads)
{
    return (total_ops / num_threads) + 
            (thread_id < (total_ops % num_threads) ? 1 : 0);
}

}           // End namespace util
}           // End namespace haplo

#endif      // PARAHAPLO_CPP_UTIL_HPP
