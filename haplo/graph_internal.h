// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for internal graph struct which holds the graph variables
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_GRAPH_INTERNAL_H
#define PARHAPLO_GRAPH_INTERNAL_H

#include "edge.h"

namespace haplo     {
namespace internal  {
 
struct ALIGN(16) Graph {    
public:
    Edge*           edges;                  // The edges for the graph
    size_t*         set_one;                // The first partition
    size_t*         set_two;                // The second partition
    size_t          set_one_size;            // Number of fragments in p1
    size_t          set_two_size;            // Number of fragments in p2 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Initializes the class variables 
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    Graph() noexcept : set_one_size{0}, set_two_size{0} {} 
};

}
}
#endif          // PARAHAPLO_GRAPH_INTERNAL_H
