// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo node container class and the structs which it uses
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_CONTAINER_HPP
#define PARHAPLO_NODE_CONTAINER_HPP

#include "devices.hpp"

#include <tbb/tbb.h>

#include <memory>

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @class     Link
/// @brief      A link between two nodes, there is a homozygous component -- how strongly correlated the nodes
///             are (that they should have the same value) -- and a heterozygous component -- how stongly they
///             should be different.
// ----------------------------------------------------------------------------------------------------------
class Link {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using atomic_type = tbb::atomic<size_t>;
private:  
    atomic_type     _homo_weight;        //!< Weight of the link if the nodes have the same ideal values
    atomic_type     _hetro_weight;       //!< Weight of the link if the nodes have different ideal values
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    Link() noexcept : _homo_weight{0}, _hetro_weight{0} {}
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for link class
    // ------------------------------------------------------------------------------------------------------
    ~Link() noexcept {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the homozygous weight 
    /// @return     A reference to the homozygous weight
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& homo_weight() { return _homo_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the heteroygous weight 
    /// @return     A reference to the heteroygous weight
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& hetro_weight() { return _hetro_weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Const ccessor for the homozygous weight 
    /// @return     A tosnt reference to the homozygous weight
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& homo_weight() const { return _homo_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Const ccessor for the heteroygous weight 
    /// @return     A const reference to the heteroygous weight
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& hetro_weight() const { return _hetro_weight; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the link so that it cant be used bya sorting function
    /// @return     The value (maximum weight of the node)
    // ------------------------------------------------------------------------------------------------------
    inline size_t value() const { return std::max(_homo_weight, _hetro_weight); }
};

// ----------------------------------------------------------------------------------------------------------
/// @class     Node
/// @brief      Each node has a weight and and index, the index represents the position which the node models
///             in the haplotype and the weight is the significance of the variable
// ----------------------------------------------------------------------------------------------------------
class Node {
public:
    // -------------------------------------- ALIAS'S -------------------------------------------------------
    using atomic_type   = tbb::atomic<size_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    atomic_type     _weight;            //!< The weight of the node (how important it is)
    atomic_type     _worst_case;        //!< The worst case contribution to the score
    atomic_type     _haplo_pos;         //!< The position in the haplotype the node represents
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    Node() noexcept : _weight{1}, _worst_case{0}, _haplo_pos{0} {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for node class
    // ------------------------------------------------------------------------------------------------------
    ~Node() noexcept {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the weight
    /// @return     A reference to the weight
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& weight() { return _weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the weight
    /// @return     A constant reference to the weight
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& weight() const { return _weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the haplo position
    /// @return     A reference to haplo position
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& position() { return _haplo_pos; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the haplo position
    /// @return     A reference to haplo position
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& position() const { return _haplo_pos; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the node so that it cant be used by sorting function
    /// @return     The value (weight of the node)
    // ------------------------------------------------------------------------------------------------------
    inline size_t value() const { return _weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      The worst case value of the node (const reference)
    /// @return     A const reference to the worst case value of the node
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& worst_case_value() const { return _worst_case; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The worst case value of the node
    /// @return     A reference to the worst case value of the node
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& worst_case_value() { return _worst_case; };
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

