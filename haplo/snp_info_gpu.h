// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo snp info class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_SNP_INFO_GPU_H
#define PARHAPLO_SNP_INFO_GPU_H

#ifdef __CUDACC__
    #define CUDA_HD __host__ __device__
    #define CUDA_H  __host__
    #define CUDA_D  __device__
    #define SHARED  __shared__
#else
    #define CUDA_HD
    #define CUDA_H
    #define CUDA_D
    #define SHARED
#endif

#include "snp_info.hpp"
#include <stdint.h>

namespace haplo {
           
// ----------------------------------------------------------------------------------------------------------
/// @class      SnpInfoGpu
/// @brief      Stores some information about a SNP for the GPU
// ----------------------------------------------------------------------------------------------------------
class SnpInfoGpu {
private:
    size_t      _start_idx;
    size_t      _end_idx;
    uint8_t     _type;          //!< IH or NIH
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    SnpInfoGpu()
    : _start_idx{0}, _end_idx{0}, _type{0} {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor -- sets the values
    /// @param[in]  snp_info        The cpu side snp info to create this snp from
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    SnpInfoGpu(const SnpInfo& snp_info_cpu) noexcept
    : _start_idx(snp_info_cpu.start_index()), 
      _end_idx(snp_info_cpu.end_index())    , 
      _type(snp_info_cpu.type())            {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the staet index of the read 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t start_index() const { return _start_idx; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the staet index of the read 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t& start_index() { return _start_idx; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the end index of the read
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t end_index() const { return _end_idx; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the end index of the read
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t& end_index() { return _end_idx; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the type of the snp
    /// @param[in]  value   The value to set the type to
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline void set_type(const uint8_t value) { _type = value & 0x03; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the type of the snp
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t type() const { return _type; }    

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the length of the read 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t length() const { return _end_idx - _start_idx + 1; }
};

}           // End namespace haplo
#endif      // PARAHAPLO_SNP_INFO_GPU_H
