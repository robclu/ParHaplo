#include "cuda.h"
#include <stdint.h>
#include "read_info_gpu.cuh"
#include "snp_info_gpu.h"

namespace haplo {
    
__global__ void search_tree(uint8_t* data       , size_t data_size, 
                            ReadInfo* read_info , size_t read_size, 
                            SnpInfoGpu* snp_info, size_t snps_size, 
                            size_t * result)
{
    // Check that this is working 
    for (size_t i = 0; i < read_size; ++i) {
        if (read_info[i].start_index() != 0) 
            *result = read_info[i].start_index();
    }
}

}               // End namespace haplo
