// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for internal graph struct which holds the graph variables
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_GRAPH_INTERNAL_H
#define PARHAPLO_GRAPH_INTERNAL_H

#include "edge.h"
#include "fragment.h"

namespace haplo     {
namespace internal  {
 
struct ALIGN(16) Graph {    
public:
    Edge*           edges;                  // The edges for the graph
    size_t*         set_one;                // The first partition
    size_t*         set_two;                // The second partition
    size_t*         mec_score;              // The mec score for the solution
    size_t*         snp_scores_one;         // Contribution of each snp (best case)
    size_t*         snp_scores_two;         // Contribution of each snp (worst case)
    uint8_t*        haplo_one;              // The first haplotype -- for set 1
    uint8_t*        haplo_two;              // The second haplotype -- for set 2
    uint8_t*        haplo_one_temp;         // A temporary haplotype
    uint8_t*        haplo_two_temp;         // A temporary haplotype
    Fragment*       fragments;              // The fragments for the partitions
    size_t          set_one_size;           // Number of fragments in p1
    size_t          set_two_size;           // Number of fragments in p2 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Initializes the class variables 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    Graph() noexcept : set_one_size{0}, set_two_size{0} {} 
};

}
}
#endif          // PARAHAPLO_GRAPH_INTERNAL_H
