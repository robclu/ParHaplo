// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo snp info class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_EDGE_H
#define PARHAPLO_EDGE_H

#include "cuda_defs.h"
#include <stdint.h>

namespace haplo {
    
struct ALIGN(16) Edge {
    float   distance;       // Distance between the fragments
    int32_t f1;             // The first fragment index
    int32_t f2;             // The second fragment index
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor to set the variables
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    Edge() noexcept : distance{0.0f}, f1{0}, f2{0} {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Assignment operator
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    void operator=(const Edge& rhs) 
    {
        distance = rhs.distance;
        f1       = rhs.f1;
        f2       = rhs.f2;
    }

};

}               // End namespace haplo
#endif          // PARAHPLO_EDGE_H
           
