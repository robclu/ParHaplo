// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo graph class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_GRAPH_H
#define PARHAPLO_GRAPH_H

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Graph
/// @brief      Graph class which is used to search the space for the potential haplotypes
/// @tparam     DeviceType  The type of device to use the node on -- so that we can optimize functions for the
///             different implementations and so that each one can be excluded from compilation if necessary
/// @tparam     SubBlockType    The type of the sublock from which this tree is derived
// ----------------------------------------------------------------------------------------------------------
template <typename SubBlockType, uint8_t DeviceType>
class Graph;

}           // End namespace haplo
#endif      // PARAHAPLO_GRAPH_H

