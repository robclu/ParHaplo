// ----------------------------------------------------------------------------------------------------------
/// @file   usplit_block.hpp
/// @brief  Header file to include all unsplittable block implementations and define the unsplittable types
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_USPLIT_BLOCK_HPP
#define PARAHAPLO_USPLIT_BLOCK_HPP

#include "usplit_block_cpu.hpp"

namespace haplo {

class Void {};

// Put and type aliases for unsplittable blocks here
template <Device DeviceType, typename BaseBlock = Void>
using UnsplittableBlock = UnsplittableBlockImplementation<DeviceType, BaseBlock>;

}

#endif      // PARAHAPLO_USPLIT_BLOCK_HPP
