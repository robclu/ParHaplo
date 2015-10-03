// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo node container class and the structs which it uses
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_CONTAINER_HPP
#define PARHAPLO_NODE_CONTAINER_HPP

#include "devices.hpp"

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @struct     Node
/// @brief      Each node has a weight and and index, the index represents the position which the node models
///             in the haplotype and the weight is the significance of the variable
// ----------------------------------------------------------------------------------------------------------
struct Node {
    
    size_t _weight;     //!< The weight of the node (how important it is)
    size_t _index;      //!< The index of the node (which position in the haplotype it represents)
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    Node() : _weight(1), _index(0) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the node so that it cant be used bya sorting function
    /// @return     The value (weight of the node)
    // ------------------------------------------------------------------------------------------------------
    inline size_t value() const { return _weight; }
};

// ----------------------------------------------------------------------------------------------------------
/// @struct     Link
/// @brief      A link between two nodes, there is a homozygous component -- how strongly correlated the nodes
///             are (that they should have the same value) -- and a heterozygous component -- how stongly they
///             should be different.
// ----------------------------------------------------------------------------------------------------------
struct Link {
  
    size_t _homo_weight;        //!< Weight of the link if the nodes have the same ideal values
    size_t _hetro_weight;       //!< Weight of the link if the nodes have different ideal values

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    Link() : _homo_weight(0), _hetro_weight(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the link so that it cant be used bya sorting function
    /// @return     The value (maximum weight of the node)
    // ------------------------------------------------------------------------------------------------------
    inline size_t value() const { return std::max(_homo_weight, _hetro_weight); }
};

// ----------------------------------------------------------------------------------------------------------
/// @class      NodeContainer 
/// @brief      Holds nodes for the tree that needs to be searched. The container is structured as follows:
///             The information for each of the nodes is stored first
///
///             [{weight0,index0}, {weight1,index1}, ..., {weightN,indexN}]
///
///             Then the link weights for the connections between the nodes is stored:
///
///             [{how01,hew01}, {how02,hew02},...,{how0N,hew0N},{how12,hew12},....{how(N-1)N,hew(N-1)N}]
///
///             where:
///
///             howAB = homozygous weight between node A and B
///             hewAB = heterozygous weight between node A and B
/// @tparam     DeviceType      The type of device being used
// ----------------------------------------------------------------------------------------------------------
template <uint8_t DeviceType>
class NodeContainer;

}           // End namespace haplo
#endif      // PARAHAPLO_NODE_CONTAINER_HPP

