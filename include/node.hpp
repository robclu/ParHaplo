// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo nade class for a tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_HPP
#define PARHAPLO_NODE_HPP

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Node
/// @brief      Defines a node of a tree, which is used to calculate a bound for branches stemming from the
///             tree, which can be used to either prune the tree or explore the branch further
// ----------------------------------------------------------------------------------------------------------
class Node {
private:
    short       _x_index;           //!< The index of the variable x for which the value should be determined
    short       _y_index;           //!< The index of the variable y for which the value should be determined   
    uint8_t     _x_value : 1;       //!< The value of x for the given x_index
    uint8_t     _y_value : 1;       //!< The value of the y for the given y_index
    uint8_t     _score   : 6;       //!< The score of this node (0 or 1, but use 6 bits for padding)          
public:
    Node(short x_index, short y_index, uint8_t x_value, uint8_t y_value)
    : _x_index(x_index), _y_index(y_index), _x_value(x_value & 1), _y_value(y_value & 1) {}
    
    inline uint8_t x_value() const { return _x_value; }
    
    inline uint8_t y_value() const { return _y_value; }
};

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_HPP
