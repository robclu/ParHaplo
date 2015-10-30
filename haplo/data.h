// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for a data struct for the GPU
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_DATA_H
#define PARHAPLO_DATA_H

#include "cuda_defs.h"
#include "read_info.h"
#include "snp_info_gpu.h"
#include <thrust/host_vector.h>

namespace haplo {
 
struct Data {    
public:
    //------------------------------------- ALIAS's ---------------------------------------------------------
    using read_info_type                = ReadInfo;
    using snp_info_type                 = SnpInfoGpu;
    using small_type                    = uint8_t;
    using small_container               = thrust::host_vector<small_type>;
    //------------------------------------------------------------------------------------------------------- 

    size_t          snps;               // Number of snps to serach
    size_t          reads;              // Number of reads for the sub-block
    small_type*     data;               // The actual data 
    read_info_type* read_info;          // The information for each read for fast access
    snp_info_type*  snp_info;           // The information for each snp for fast access

    //------------------------------------------------------------------------------------------------------- 
    /// @brief  Constructor to initialize avriables
    //------------------------------------------------------------------------------------------------------- 
    CUDA_HD
    Data(const size_t num_snps, const size_t num_reads)
    : snps(num_snps), reads(num_reads) {}
};

}
#endif          // PARAHAPLO_DATA_H
