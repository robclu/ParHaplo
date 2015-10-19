// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree node for the gpu search tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_NODE_H
#define PARHAPLO_TREE_NODE_H

#include <stdint.h>

#ifdef __CUDAACC__ 
    #define CUDA_HD __host__ __device__
#else 
    #define CUDA_HD
#endif

namespace haplo {
    
struct TreeNode {
    
    size_t          haplo_idx   ;
    size_t          root_idx    ; 
    size_t          node_idx    ;
    size_t          alignments  ;
    unsigned int    ubound      ;
    unsigned int    lbound      ;
    uint8_t         value       ;
    size_t*         indices     ;
    size_t*         read_ids    ;
    uint8_t*        read_values ;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Assignment operator
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    void operator=(const TreeNode& other) const
    {
        other.haplo_idx  = haplo_idx    ;
        other.root_idx   = root_idx     ;
        other.alignments = alignements  ;
        other.ubound     = ubound       ;
        other.lbound     = lbound       ;
        other.value      = value        ;
        
        // Pointer copy
        if (other.alignments != alignments) {
            free(other.indices); free(other.read_ids); free(other.read_values); 
          
            other.indices       = new size_t[alignments];
            other.read_ids      = new size_t[alignments];
            other.read_values   = new uint8_t[alignments];
        }
        
        for (size_t i = 0; i < alignments; ++i) {
            other.indices[i]     = indices[i];
            other.read_ids[i]    = read_ids[i];
            other.read_values[i] = read_values[i];
        }
    }
};

}

#endif 
