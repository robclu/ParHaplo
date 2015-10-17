// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for the parahaplo device manager class whihc maps sub-blocks onto devices
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_DEVICE_MANAGER_H
#ifndef PARAHAPLO_DEVICE_MANAGER_H

#include "devices.hpp"
#include "cuda.h"
#include <tbb/tbb.h>

#ifdef __CUDACC__
    #define CUDA_HD __host__ __device__
    #define CUDA_H  __host__
    #define CUDA_D           __device__
    #define SHARED __shared__
#else 
    #define CUDA_HD
    #define CUDA_H
    #define CUDA_D
    #define SHARED 
#endif

namespace haplo {
    
template <uint8_t DeviceType>
class DeviceManager;

// ----------------------------------------------------------------------------------------------------------
/// @class      DeviceManager 
/// @brief      Determines the number of compute devices in the system, and then given a sub-block, maps the
///             block to a device
// ----------------------------------------------------------------------------------------------------------
template <typename DeviceManager>
class DeviceManager {
public:
    // --------------------------------------------- ALIAS'S ------------------------------------------------
    using atomic_size_t = tbb::atomic<size_t>;
    using atomic_int    = tbb::atomic<int>:
    // ------------------------------------------------------------------------------------------------------
private:
    size_t          _total_devices;         //!< Total number of computation dives in system
    atomic_size_t   _total_elements;        //!< Total number of system elements
    atomic_int      _device_to_use;         //!< Total index of the device to use 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets params
    // ------------------------------------------------------------------------------------------------------
    CUDA_H
    DeviceManager() : _total_elements(0), _device_to_use(0)
    {
        struct cudaDeviceProp device_properties;    // Properties of a device
        cudaError_t status  = cudaGetDeviceCount(&_total_devices);
        
        // If there was an error 
        if (status != cudaSuccess) _total_devices = 0;
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the index of the device to use to compute the block's haplotype 
    // ------------------------------------------------------------------------------------------------------
    int get_device(const size_t block_width)
    {
        
    }
    
    
};


}                   // End namespace haplo
#endif              // PARAHAPLO_DEVICE_MANAGER_H

    

