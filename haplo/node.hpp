// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo nade class, which hold properties of a bit of a potention solution
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_HPP
#define PARHAPLO_NODE_HPP

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Node    
/// @brief      Represents a possible solution at a SNP site, and holds properties which define its
///             contribution to the final solution and how it is related to other nodes
/// @tparam     DeviceType  The type of device to use the node on -- so that we can optimize functions for the
///             different implementations and so that each one can be excluded from compilation if necessary
// ----------------------------------------------------------------------------------------------------------
template <byte DeviceType>
class Node;

}           // End namespace haplo
#endif      // PARAHAPLO_NODE_HPP

