// ----------------------------------------------------------------------------------------------------------
/// @file   usplit_block_interface.hpp
/// @brief  Header file for the unsplittable block interface class for the parahaplo library -- this just 
///         defines an 'empty' 
///         base which allows the specific block implementation to be selected at compile time 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_USPLIT_BLOCK_INTERFACE_HPP
#define PARAHAPLO_USPLIT_BLOCK_INTERFACE_HPP

#include "devices.hpp"

namespace haplo{

// ----------------------------------------------------------------------------------------------------------
/// @class      UnsplittableBlockInterface 
/// @brief      Empty base class for unsplittable blocks, which provides an interface for selecting the
///             implementation method -- CPU, GPU, or PHI
/// @tparam     Implementation  The type of implementation to use -- CPU, GPU, or PHI
// ----------------------------------------------------------------------------------------------------------
template <typename Implementation>
class UnsplittableBlockInterface : Implementation {
    
};

// ----------------------------------------------------------------------------------------------------------
/// @class      UnsplittableBlockImplementation 
/// @brief      UnsplittableBlock implementation class which can be specialized to provide CPU, GPU and PHI
///             implementations
/// @tparam     DeviceType      The type of device which is used -- CPU, GPU or PHI
/// @tparam     BaseBlock       The block which is the base of the unsplittable block 
// ----------------------------------------------------------------------------------------------------------
template <Device DeviceType, typename BaseBlock>
class UnsplittableBlockImplementation;
        
}               // End namespace haplo

#endif          // PARAHAPLO_USPLIT_BLOCK_INTERFACE
