// ----------------------------------------------------------------------------------------------------------
/// @file Header file for the parahaplo fragment class which stores the contribution of a fragment to
///       the MEC score 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_FRAGMENT_H
#define PARAHAPLO_FRAGMENT_H

#include "cuda_defs.h"

struct ALIGN(16) Fragment {
    uint64_t score : 62;
    uint64_t set   : 2 ;
    size_t   index ;
    
    CUDA_HD 
    Fragment() noexcept : score{0}, set{0}, index{0} {}
};

#endif          // PARAHAPLO_FRAGMENT_H

