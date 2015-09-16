// ----------------------------------------------------------------------------------------------------------
/// @file   block.hpp
/// @brief  Header file to include all block implementations and define the types
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_BLOCK_HPP
#define PARAHAPLO_BLOCK_HPP

#include "block_cpu.hpp"

namespace haplo {
    
// Alias for Block
template <size_t Rows, size_t Cols, size_t Cores, haplo::Device DeviceType>
using Block = BlockInterface<BlockImplementation<Rows, Cols, Cores, DeviceType>>;

}

#endif      // PARAHAPLO_BLOCK_HPP