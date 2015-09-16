// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo devices enum -- defines the available computational devices
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_DEVICES_HPP
#define PARHAPLO_DEVICES_HPP

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @enum       Device
/// @brief      The available devices 
// ----------------------------------------------------------------------------------------------------------
enum class Device : uint8_t {
 
CPU,            //!< For multi-threaded CPU implementations  using Intel TBB with multi-core CPUs
GPU,            //!< For multi-threaded GPU implementations using CUDA with Nvidia GPUs
PHI             //!< For multi-threaded PHI implementations using Intel TBB with Intel Xeon Phi
    
};

}           // End namespace haplo

#endif      // PARAHAPLO_DEVICES_HPP
