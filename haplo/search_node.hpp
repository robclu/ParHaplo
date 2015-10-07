// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo search node class -- defines a lightweight node for a search tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_SEARCH_NODE_HPP
#define PARHAPLO_SEARCH_NODE_HPP

#include <cstdint>

namespace haplo {

class SearchNode {
private:
    uint16_t    _index : 15;    //!< The index the node represents
    uint16_t    _value : 1;     //!< If the value is a 1 or a 0
    uint16_t    _upper_bound;   //!< Node upper bound
    uint16_t    _lower_bound;   //!< Node lower bound
    size_t      _left;          //!< The index of the left node
    size_t      _right;         //!< The index of the right node
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    explicit SearchNode() noexcept 
    : _index(0), _value(0), _upper_bound(0), _lower_bound(0), _left(0), _right(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for when the parameters are known
    // ------------------------------------------------------------------------------------------------------
    explicit SearchNode(const uint16_t index       , const uint8_t  value       ,  
                        const uint16_t upper_bound , const uint16_t lower_bound )
    : _index(index & 0x7FFF)    , _value(value & 0x01)      , 
      _upper_bound(upper_bound) , _lower_bound(lower_bound) ,
      _left(0)                  , _right(0)                 {}

    inline void set_value(const uint16_t value) { _value = value & 0x01; }
    
    inline void set_index(const uint16_t index) { _index = index & 0x7FFF; }
    
    inline uint8_t value() const { return _value & 0x01; }
    
    inline uint16_t index() const { return _index & 0x7FFF; }
    
    inline uint16_t& upper_bound() { return _upper_bound; }
    
    inline uint16_t& lower_bound() { return _lower_bound; }
    
    inline uint16_t upper_bound() const { return _upper_bound; }
    
    inline uint16_t lower_bound() const { return _lower_bound; }
    
    inline size_t& left() { return _left; }
    
    inline size_t& right() { return _right; }
    
    inline size_t left() const { return _left; }
    
    inline size_t right() const { return _right; }
};

}                   // End namespace haplo
#endif              // PARHAPLO_SEARCH_NODE_HPP
