// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo processor functionality, processes rows/columns 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_PROCESSOR_HPP
#define PARAHAPLO_PROCESSOR_HPP

#include "devices.hpp"

namespace haplo {
namespace proc {
    
static constexpr uint8_t ROW_DUP           = 0x00;
static constexpr uint8_t COL_DUP_LINKS     = 0x02;

}       // End namespace proc

template <typename FriendType, uint8_t ProcessType, uint8_t DeviceType>
class Processor;

}               // End namespace haplo
#endif          // PARAHAPLO_PROCESSOR_HPP
