// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo nade class for a tree
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_HPP
#define PARHAPLO_NODE_HPP

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Node
/// @brief      Defines a node of a tree, which is used to calculate a bound for branches stemming from the
///             tree, which can be used to either prune the tree or explore the branch further. The nodes do
///             not have links to the next node as the tree is designed to be represented as a vector with a
///             "guess" size -- since for this problem we can make a pretty accurate guess of the size -- once
///             or two reallloations (and the cache improvement of hvaing all tree nodes right next to each 
///             other) will be better performance than having the nodes randomly scattered through memory
// ----------------------------------------------------------------------------------------------------------
class Node {
private:
    short       _x_index;           //!< The index of the variable x for which the value should be determined
    short       _y_index;           //!< The index of the variable y for which the value should be determined   
    uint8_t     _x_value : 1;       //!< The value of x for the given x_index
    uint8_t     _y_value : 1;       //!< The value of the y for the given y_index
    uint8_t     _t_value : 1;       //!< The value of the constraint parameter
    uint8_t     _score   : 6;       //!< The score of this node (0 or 1, but use 6 bits for padding)          
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor for a node
    // ------------------------------------------------------------------------------------------------------
    Node() : _x_index(0), _y_index(0), _x_value(0), _y_value(0), _t_value(0), _score(0) {};
 
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for when the y variable index and value is know
    /// @param[in]  y_index     The index of the y variable
    /// @param[in]  y_value     The value of the y variable
    // ------------------------------------------------------------------------------------------------------
    Node (const short y_index, const uint8_t y_value)
    : _y_index(y_index), _y_value(y_value & 1), _x_index(0), _x_value(0), _t_value(0), _score(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for when the both the x and y variables and indices are known
    /// @param[in]  x_index     The index of the x variable
    /// @param[in]  x_value     The value of the x variable
    /// @param[in]  y_index     The index of the y variable
    /// @param[in]  y_value     The value of the y variable
    // ------------------------------------------------------------------------------------------------------
    Node(short x_index, short y_index, uint8_t x_value, uint8_t y_value)
    : _x_index(x_index), _y_index(y_index), _x_value(x_value & 1), _y_value(y_value & 1), _t_value(0) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the x variable 
    /// @param[in]  value       The value to set the variable x to
    // ------------------------------------------------------------------------------------------------------
    inline void set_x_value(uint8_t value) { _x_value = value & 1; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the value of the y variable 
    /// @param[in]  value       The value to set the variable y to
    // ------------------------------------------------------------------------------------------------------
    inline void set_y_value(uint8_t value) { _y_value = value & 1; } 
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the index of the x variable 
    /// @param[in]  value       The value to set the x index to
    // ------------------------------------------------------------------------------------------------------
    inline void set_x_index(short index) { _x_index = index; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the index of the y variable 
    /// @param[in]  value       The value to set the y index to
    // ------------------------------------------------------------------------------------------------------
    inline void set_y_index(short index) { _y_index = index; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the x variable
    /// @return     The value of the x variable
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t x_value() const { return _x_value; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the value of the y variable
    /// @return     The value of the y variable
    // ------------------------------------------------------------------------------------------------------
    inline uint8_t y_value() const { return _y_value; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Computes the score of the node based on the x and y value using the first computation 
    ///                 score = (1 - x_j - y_i + 2t_i,j)
    ///             where:
    ///                 x_j     = value of the variable x_j = _x_value
    ///                 y_i     = value of the variable y_i = _y_value
    ///                 t_i,j   = constrait
    /// @return     The score of the node based on the parameters
    // ------------------------------------------------------------------------------------------------------
    uint8_t score_zero();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Computes the score of the node based on the x and y value using the second computation 
    ///                 score = (y_i + x_j - 2t_i,j)
    ///             where:
    ///                 x_j     = value of the variable x_j = _x_value
    ///                 y_i     = value of the variable y_i = _y_value
    ///                 t_i,j   = constrait
    /// @return     The score of the node based on the parameters
    // ------------------------------------------------------------------------------------------------------
    uint8_t score_one();
private:
    void determine_constraint();
};

// ---------------------------------------------- IMPLEMENTATIONS -------------------------------------------

void Node::determine_constraint()
{
    // If one of the values is 0, then the constraint is 0
    if (!(_x_value & 0x01) || !(_y_value & 0x01)) _t_value = 0x00;
    else if ((_x_value & 0x01) && (_y_value & 0x01)) _t_value = 0x01;
}

uint8_t Node::score_zero()
{
    determine_constraint();
    return 1 - _x_value - _y_value + (2 *_t_value);
}

uint8_t Node::score_one()
{
    determine_constraint();
    return _y_value + _x_value - (2 *_t_value);
}

}           // End namespace haplo
#endif      // PARAHAPLO_NODE_HPP
