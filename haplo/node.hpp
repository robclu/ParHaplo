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
/// @class      Node
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
    size_t          _num_elements;      //!< The number of elements in the col the node represents
    uint8_t         _haplo_value;       //!< The value of the haplotype for this position
    uint8_t         _is_intrin;         //!< If the node is intrinsically heterozygoud
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for initialization
    // ------------------------------------------------------------------------------------------------------
    Node() noexcept : 
    _weight{1}, _worst_case{0}, _haplo_pos{0}, _num_elements(0),_haplo_value{0}, _is_intrin{1} {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for node class
    // ------------------------------------------------------------------------------------------------------
    ~Node() noexcept {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Copy constructor
    /// @param[in]  other       The other node to copy from
    // ------------------------------------------------------------------------------------------------------
    Node(const Node& other) noexcept 
    : _weight(other._weight)            , _worst_case(other._worst_case)   , _haplo_pos(other._haplo_pos)  , 
      _num_elements(other._num_elements),  _haplo_value(other._haplo_value), _is_intrin(other._is_intrin) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Move constructor
    /// @param[in]  other       The other node to copy from
    // ------------------------------------------------------------------------------------------------------
    Node(Node&& other) noexcept 
    : _weight(std::move(other._weight))             ,    
      _worst_case(std::move(other._worst_case))     , 
      _haplo_pos(std::move(other._haplo_pos))       ,
      _num_elements(std::move(other._num_elements)) ,
      _haplo_value(std::move(other._haplo_value))   ,
      _is_intrin(std::move(other._is_intrin)) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Copy assigment operator
    /// @param[in]  other       The other node to copy from
    // ------------------------------------------------------------------------------------------------------
    void operator=(const Node& other) 
    {
        _weight = other.weight()            ; _worst_case = other.worst_case_value(); 
        _haplo_pos = other.position()       ; _num_elements = other.elements();
        _haplo_value = other.haplo_value()  ; _is_intrin = other.type();
    }
    
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
    /// @brief      Accessor for the haplo value
    /// @return     The value of the haplotype at this position
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t haplo_value() const { return _haplo_value; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Accessor for the haplo value
    /// @return     A reference to haplo value
    // ------------------------------------------------------------------------------------------------------
    inline void set_haplo_value(const size_t value) { _haplo_value = value & 0x01; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the value of the node so that it cant be used by sorting function
    /// @return     The value (weight of the node)
    // ------------------------------------------------------------------------------------------------------
    inline size_t value() const { return _weight; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      The worst case value of the node (const reference)
    /// @return     A const reference to the worst case value of the node
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type worst_case_value() const { return _worst_case; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The worst case value of the node
    /// @return     A reference to the worst case value of the node
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& worst_case_value() { return _worst_case; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The type of the column -- IH or not IH
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t type() const { return _is_intrin; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The type of the column -- IH or not IH
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t& type() { return _is_intrin; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The number of elements in the column the node represents
    // ------------------------------------------------------------------------------------------------------
    inline size_t elements() const { return _num_elements; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The number of elements in the column the node represents
    // ------------------------------------------------------------------------------------------------------
    inline size_t& elements() { return _num_elements; };
};

// ----------------------------------------------------------------------------------------------------------
/// @class      NodeContainer 
/// @brief      Holds nodes for the tree that needs to be searched.
/// @tparam     DeviceType      The type of device being used
// ----------------------------------------------------------------------------------------------------------
template <uint8_t DeviceType>
class NodeContainer;

}           // End namespace haplo
#endif      // PARAHAPLO_NODE_CONTAINER_HPP

