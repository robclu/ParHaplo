// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for internal graph struct which holds the graph variables
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_GRAPH_INTERNAL_H
#define PARHAPLO_GRAPH_INTERNAL_H

#include "edge.h"

namespace haplo     {
namespace internal  {
 
struct Graph {    
public:
    Edge*           edges;              // The edges for the graph
    
    // Need tp still add partitions
};

}
}
#endif          // PARAHAPLO_GRAPH_INTERNAL_H
