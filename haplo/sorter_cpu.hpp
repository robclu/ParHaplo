// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo parallel sorter class cpu implementation 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_SORTER_CPU_HPP
#define PARHAPLO_SORTER_CPU_HPP

#include "comparators.hpp"
#include "sorter.hpp"

#include <tbb/parallel_invoke.hpp>

#include <algorithm>

namespace haplo {
           
// Specialization for CPU using TBB 
template <>
class Sorter<devices::cpu> {
public:
private:
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor
    // ------------------------------------------------------------------------------------------------------
    Sorter() noexcept {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Operator to invoke the sorting
    /// @param[in]
    // ------------------------------------------------------------------------------------------------------
    template <typename Iterator, typename Comparator>
    void operator()(Iterator input_start, Iterator input_end, Comparator comparator);
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Parallel merge function 
    // ------------------------------------------------------------------------------------------------------
    template <typename Iterator, typename Comparator>
    void merge(Iterator fh_start , Iterator fh_end    , 
               Iterator sh_start , Iterator sh_end    ,
               Iterator res_start, Comparator comparator   );
    
}
    
// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

template <typename Iterator, typename Comparator>
void Sorter<devices::cpu>::operator()(Iterator input_start, Iterator input_end, Comparator comparator)
{
    // We can check here is we can get enough memory for the non-inplace (faster) sorting method
    
    // Create the resultant array
    typename Iterator::value_type output_array[end - start];

    Iterator input_mid  = input_start  + ((input_end - input_start) / 2);
    Iterator output_mid = output_array + (input_mid - input_start);
    Iterator output_end = output_array + (input_end - input_start);
    
    // Call the parallel function on each hald of the input arrays
    tbb::parallel_invoke([=]{operator()(input_start, input_mid, output_start, comparator);},
                         [=]{operator()(input_mid  , input_end, output_mid  , comparator);} );
    
    // Merge the arrays
    merge(input_start, input_mid, input_mid, input_end, output_start, comparator);
}

template <typename Iterator, typename Comparator>
void Sorter<devices::cpu>::merge(Iterator first_start  , Iterator first_end   , 
                                 Iterator second_start , Iterator second_end  , 
                                 Iterator output_start , Comparator comparator     )
{
    // Check for small array and do serial sorting 
    
    // Case for large arrays -- find the middle of the arrays
    Iterator first_mid, second_mid;
   
    auto first_size  = first_end  - first_start;
    auto second_size = second_end - second_start;
    
    if (first_size < second_size) {
        second_mid = second_start + ((second_end - second_start) / 2);
        first_mid  = std::upper_bound(first_start, first_end, *second_mid, comparator);
    } else {
        first_mid  = first_start + ((first_end - first_start) / 2);
        second_mid = std::lower_bound(second_start, second_end, *first_mid, comparator);
    }
    
    Iterator output_mid = output_start + ((first_mid - first_start) + (second_mid - second_start));
    
    // Call the function in parallel
    tbb::parallel_invoke([=]{merge(first_start, first_mid, second_start, second_mid, output_start, comparator);},
                         [=]{merge(first_mid  , first_end, second_mid  , second_end, output_start, comparator);} );
}

}               // End namespace haplo
#endif          // PARHAPLO_SORTER_CPU_HPP

