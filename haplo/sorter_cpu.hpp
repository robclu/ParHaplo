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

namespace haplo {
    
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
    void operator()(DataType* input_start, DataType* input_end, Comparator comparator);
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Parallel sort function
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void para_sort(DataType* input_start, DataType* input_end, DataType* output_start, Comparator comparator);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Parallel merge function 
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void para_merge(DataType* first_start  , DataType*  first_end    , 
                    DataType* second_start , DataType*  second_end   ,
                    DataType* output_start , Comparator comparator   );
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Serial sort
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void serial_sort(DataType*  input_start , 
                     DataType*  input_end   , 
                     DataType*  output_start, 
                     Comparator comparator  );
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Serial merge function 
    // ------------------------------------------------------------------------------------------------------
    template <typename DataType, typename Comparator>
    void serial_merge(DataType* first_start  , DataType*  first_end    , 
                    DataType* second_start , DataType*  second_end   ,
                    DataType* output_start , Comparator comparator   );
};
    
// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::operator()(DataType* input_start, DataType* input_end, Comparator comparator)
{
    // If we can allocate enough memory for aanother array
    if (Buffer<DataType> output_array = Buffer<DataType>(sizeof(DataType) * (input_end - input_start))) {
        // Sort the array not-ipplce
        para_sort(input_start, input_end, output_array.fetch_data(), comparator);
    } else {
        // Failed to get the memory -- create an inplace implementation if necessary
    }
}

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::para_sort(DataType*   input_start   , 
                                     DataType*   input_end     , 
                                     DataType*   output_start  , 
                                     Comparator  comparator    )
{
    if (input_end - input_start < sort_term) {                      // Serial sort
    } else {                                                        // Parallel sort
        // Do the parallel sorting
        DataType* input_mid  = input_start  + ((input_end - input_start) / 2);
        DataType* output_mid = output_start + (input_mid - input_start);
            
        // Call the parallel function on each hald of the input arrays
        tbb::parallel_invoke([=]{ para_sort(input_start, input_mid, output_start, comparator); } ,
                             [=]{ para_sort(input_mid  , input_end, output_mid  , comparator); } );
        
        // Merge the arrays
        para_merge(input_start, input_mid, input_mid, input_end, output_start, comparator);
    }
}

template <typename DataType, typename Comparator>
void Sorter<devices::cpu>::para_merge(DataType* first_start  , DataType*  first_end      , 
                                      DataType* second_start , DataType*  second_end     , 
                                      DataType* output_start , Comparator comparator     )
{
    auto first_size  = first_end  - first_start;
    auto second_size = second_end - second_start;
    
    if (first_size +second_size < merge_term) {                 // Serial merge
    } else {                                                    // Parallel merge
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
                [=]{para_merge(first_start, first_mid, second_start, second_mid, output_start, comparator);},
                [=]{para_merge(first_mid  , first_end, second_mid  , second_end, output_mid  , comparator);} 
        );
    }
}

}               // End namespace haplo
#endif          // PARHAPLO_SORTER_CPU_HPP

