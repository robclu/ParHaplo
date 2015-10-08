// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo bound calculator class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_BOUNDER_HPP
#define PARHAPLO_BOUNDER_HPP

#include "bounder.hpp"

namespace haplo {

//-----------------------------------------------------------------------------------------------------------
/// @struct     Bounds  
/// @brief      Simple wrapper struct which allows the bounding operator to return both bounds at the same
///             time
//-----------------------------------------------------------------------------------------------------------
struct Bounds {
    
    size_t lower;       //!< The lower bound
    size_t upper;       //!< The upper bound

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for the bound 
    /// @param[in]  lower   The lower bound
    /// @param[in]  upper   The upper bound
    // ------------------------------------------------------------------------------------------------------
    explicit Bounds(const size_t lower_bound, const size_t upper_bound)
    : lower(lower_bound), upper(upper_bound) {}
};

// ----------------------------------------------------------------------------------------------------------
/// @class      Bounder 
/// @brief      Provides the bounding calculations for a node
/// @tparam     DeviceType  The type of device used for the computation of the bound
// ----------------------------------------------------------------------------------------------------------
template <uint8_t DeviceType>
class Bounder;


}               // End namespace haplo
#endif          // PARHAPLO_BOUNDER_HPP
