// ----------------------------------------------------------------------------------------------------------
/// @file   subblock_block.hpp
/// @brief  Header file for the subblock class definition
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_SUB_BLOCK_HPP
#define PARAHAPLO_SUB_BLOCK_HPP

#include "block.hpp"

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      SubBlock   
/// @brief      General class for creating sub blocks from an entire block, which can then be solved in
///             parallel
/// @tparam     ThreadsX    The number of threads to use in the X dimesion (columns)
/// @tparam     ThreadsY    The number of threads to use in the Y dimension (rows)
/// @tparam     DeviceType  The type of device to use 
// ----------------------------------------------------------------------------------------------------------
template <typename Block, size_t ThreadsX, size_t ThreadsY, uint8_t DeviceType>
class SubBlock;

// Specializations are in their respective files
    
}               // End namespace haplo
#endif          // PARAHAPLO_SUB_BLOCK_HPP
