// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo parallel sorter class cpu implementation 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_SORTER_CPU_HPP
#define PARHAPLO_SORTER_CPU_HPP

#include "buffer.hpp"
#include "comparators.hpp"
#include "sorter.hpp"

#include <tbb/parallel_invoke.h>

#include <algorithm>
#include <iostream>

namespace haplo {
namespace sorts {
    
static constexpr int outplace    = 0;
static constexpr int inplace     = 1;
static constexpr int start_case  = 2;

}

// Specialization for CPU using TBB 
template <>
class Sorter<devices::cpu> {
public:
private:
    const size_t _sort_term;
    const size_t _merge_term;
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor
    // ------------------------------------------------------------------------------------------------------
    Sorter() noexcept :_sort_term(200), _merge_term(500) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Operator to invoke the sorting
    /// @param[in]
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void operator()(DataType* input_start, DataType* input_end, Comparator& comparator) const;
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Parallel sort function
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void para_sort(DataType*    input_start  , DataType*  input_end , 
                   DataType*    output_start , int        sort_type ,
                   Comparator&  comparator                          ) const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Parallel merge function 
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void para_merge(DataType*   first_start  , DataType*  first_end    , 
                    DataType*   second_start , DataType*  second_end   ,
                    DataType*   output_start , bool       clean_mem    ,
                    Comparator& comparator                             ) const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Serial sort
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void serial_sort(DataType*   input_start , DataType*  input_end   , 
                     DataType*   output_start, int        sort_type   ,
                     Comparator& comparator                           ) const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Serial merge function 
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void serial_merge(DataType* first_start  , DataType*   first_end    , 
                      DataType* second_start , DataType*   second_end   ,
                      DataType* output_start , Comparator& comparator   ) const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Serial clean up function, to remove extra buffers
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType>
    void free_memory(DataType* start, DataType* end) const 
    { 
        while (end != start) { --end; (*end).~DataType(); } 
    }
    
};
    
// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::operator()(DataType*   input_start, DataType* input_end  , 
                                      Comparator& comparator                        ) const 
{
    // If we can allocate enough memory for aanother array
    if (Buffer<DataType> output_array = Buffer<DataType>(sizeof(DataType) * (input_end - input_start))) {
        // Sort the array not-ipplce
        para_sort(input_start, input_end, output_array.fetch_data(), sorts::start_case, comparator);
    } else {
        // Have to use serial sort
        std::stable_sort(input_start, input_end, comparator);
    }
}

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::para_sort(DataType*   input_start , DataType*   input_end, 
                                     DataType*   output_start, int         sort_type,
                                     Comparator& comparator                         ) const 
{
    if (input_end - input_start <= _sort_term) {                                            // Serial sort
        serial_sort(input_start, input_end, output_start, sort_type, comparator);
    } else {                                                                                // Parallel sort
        // Do the parallel sorting
        DataType* input_mid  = input_start  + ((input_end - input_start) / 2);
        DataType* output_mid = output_start + (input_mid - input_start);
        DataType* output_end = output_start + (input_end - input_start);
            
        // Call the parallel function on each hald of the input arrays
        tbb::parallel_invoke([=]{ para_sort(input_start, input_mid, output_start, !sort_type, comparator); } ,
                             [=]{ para_sort(input_mid  , input_end, output_mid  , !sort_type, comparator); } );
        
        const bool is_start_case = (sort_type == sorts::start_case);
        
        // Do the merging
        if (sort_type != sorts::outplace)
            para_merge(output_start, output_mid, output_mid, output_end, input_start, is_start_case, comparator);
        else
            para_merge(input_start, input_mid, input_mid, input_end, output_start, sorts::outplace, comparator);
    }
}

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::para_merge(DataType*   first_start  , DataType*  first_end , 
                                      DataType*   second_start , DataType*  second_end, 
                                      DataType*   output_start , bool       clean_mem ,
                                      Comparator& comparator                          ) const
{
    auto first_size  = first_end  - first_start;
    auto second_size = second_end - second_start;
    
    if (first_size + second_size <= _merge_term) {                                  // Serial merge
        serial_merge(first_start, first_end, second_start, second_end, output_start, comparator);
        if (clean_mem) {                                                            // If we must clean buffers
            free_memory(first_start, first_end); free_memory(second_start, second_end);
        }
    } else {                                                                        // Parallel merge
        // Case for large arrays -- find the middle of the arrays
        DataType* first_mid, *second_mid;
        
        if (first_size < second_size) {
            second_mid = second_start + ((second_end - second_start) / 2);
            first_mid  = std::upper_bound(first_start, first_end, *second_mid, comparator);
        } else {
            first_mid  = first_start + ((first_end - first_start) / 2);
            second_mid = std::lower_bound(second_start, second_end, *first_mid, comparator);
        }
        
        DataType* output_mid = output_start + ((first_mid - first_start) + (second_mid - second_start));
        
        // Call the function in parallel
        tbb::parallel_invoke(
            [=]{para_merge(first_start, first_mid, second_start, second_mid, output_start, clean_mem, comparator);},
            [=]{para_merge(first_mid  , first_end, second_mid  , second_end, output_mid  , clean_mem, comparator);} 
        );
    }
}

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::serial_sort(DataType*    input_start   , DataType*    input_end     , 
                                       DataType*    output_start  , int          sort_type     ,
                                       Comparator&  comparator                                  ) const 
{
    // Change to std::sort if more performace is required
    std::stable_sort(input_start, input_end, comparator);
    if (sort_type != sorts::start_case) {                            // If we haven't cleaned the memory
        DataType* output_end = output_start + (input_end - input_start);
        if (sort_type != sorts::outplace) {
            // Initialize the temp buffer
            for (; output_start < output_end; ++output_start) new(&*output_start) DataType;
        } else {
            // Initialize the temp buffer and move the data into it
            for (; output_start < output_end; input_start++, ++output_start) 
                new(&*output_start) DataType(std::move(*input_start));
        }
    }
}

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::serial_merge(DataType*  first_start  , DataType*   first_end       , 
                                        DataType*  second_start , DataType*   second_end      , 
                                        DataType*  output_start , Comparator& comparator      ) const
{
    if (first_start != first_end) {                                 // While not at the end of first
        if (second_start != second_end) {                           // While not at the end of second
            while (true) {                                          // Keep going
                if (comparator(*second_start, *first_start)) {      // If second is more important that first
                    *output_start = std::move(*second_start);       // Move second into the output
                    ++output_start; 
                    if (++second_start == second_end) break;        // Break if we are at the end of second
                } else {                                            // First ise more important
                    *output_start = std::move(*first_start);        // Move first into the result
                    ++output_start; 
                    if (++first_start == first_end) 
                        goto move_second;                           // Move second to the output
                }
            }
        }
        // Second becomes first
        second_start = first_start;                 
        second_end   = first_end;                                   
    }
move_second:
    std::move(second_start, second_end, output_start);              // Move second into the output
}

}               // End namespace haplo
#endif          // PARHAPLO_SORTER_CPU_HPP

