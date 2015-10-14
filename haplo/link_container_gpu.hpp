// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo gpu link contianer  
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_LINK_CONTAINER_GPU_HPP
#define PARHAPLO_LINK_CONTAINER_GPU_HPP

#include "linkv2.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace haplo {
    
template <>
class LinkContainer<devices::gpu> {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using link_container_host      = thrust::host_vector<LinkV2>; 
    using link_container_device    = thrust::device_vector<LinkV2>; 
    // ------------------------------------------------------------------------------------------------------
private:
    link_container_host         _host_links;        //!< Links for the CPU
    link_container_device       _device_links;      //!< Links for the GPU
    size_t                      _nodes;             //!< The number of nodes got which there may be links
    
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief Constructor 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    LinkContainer(const size_t num_nodes) noexcept {}
    

    // ------------------------------------------------------------------------------------------------------
    /// @brief  Inserts a link into the container, with thread safety
    // ------------------------------------------------------------------------------------------------------
    __host__ template <typename MutexType>
    inline void push_back(const size_t idx_one, const size_t idx_two, MutexType& mutex)
    {
        if (!link_exists(idx_one, idx_two)) {
            auto key_value = key(idx_one, idx_two);
            std::lock_guard<MutexType>  lock(mutex);
            _host_links.push_back(
        
    }
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Determines the key for an index pair
    // ------------------------------------------------------------------------------------------------------
    __host__
    inline size_t key(const size_t idx_one, const size_t idx_two) 
    {
        size_t mem_offset = 0, idx_low = 0, idx_high = 0;
        #ifdef __CUDAACC__ 
            idx_low  = min(idx_one, idx_two);
            idx_high = max(idx_one, idx_two);
        #else
            idx_low  = std::min(idx_one, idx_two);
            idx_high = std::max(idx_one, idx_two);            
        #endif
            mem_offset = 
    }
    

    
    
};



};              // End namespace haplo  
#endif          // PARAHAPLO_LINK_CONTAINER_GPU_HPP

