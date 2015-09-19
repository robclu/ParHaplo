// ----------------------------------------------------------------------------------------------------------
/// @file   usplit_block.hpp
/// @brief  Header file to include all unsplittable block implementations and define the unsplittable types
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_USPLIT_BLOCK_HPP
#define PARAHAPLO_USPLIT_BLOCK_HPP

#include "usplit_block_cpu.hpp"

namespace haplo {

// Put and type aliases for unsplittable blocks here
template <typename BaseBlock, Device DeviceType = BaseBlock::device_type, size_t Cores = BaseBlock::num_cores()>
using UnsplittableBlock = UnsplittableBlockImplementation<BaseBlock, DeviceType, Cores>;


}

#endif      // PARAHAPLO_USPLIT_BLOCK_HPP
