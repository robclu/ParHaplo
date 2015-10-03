/// @file   unsplittables_block.hpp
/// @brief  Header file for the unsplittable block class definition
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_UNSPLITTABLE_BLOCK_HPP
#define PARAHAPLO_UNSPLITTABLE_BLOCK_HPP

#include "block.hpp"

namespace haplo   {
namespace devices {
    
static constexpr byte cpu = 0;
static constexpr byte gpu = 1;
static constexpr byte phi = 2;

}           // End namespace devices

// ----------------------------------------------------------------------------------------------------------
/// @class      UnsplittableBlock   
/// @brief      General class for creating sub blocks from an entire block, which can then be solved in
///             parallel
/// @tparam     THI         The number of threads to use in the I dimension (for rows...)
/// @tparam     THJ         The number of threads to use in the J dimension (for columns...)
/// @tparam     DeviceType  The type of device to use 
// ----------------------------------------------------------------------------------------------------------
template <typename Block, size_t THI, size_t THJ, byte DeviceType>
class UnsplittableBlock;

// Specializations are in their respective files
    
}               // End namespace haplo
#endif          // PARAHAPLO_UNSPLITTABLE_BLOCK_HPP
