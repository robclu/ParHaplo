// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class - this does the branching and bounding when solving the
///         haplotype 
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_CPU_HPP
#define PARHAPLO_TREE_CPU_HPP

#include "devices.hpp"
#include "node.hpp"

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Tree 
/// @brief      Defines a tree which is searched for the optimal haplotype of an unsplittable block --
///             the implementation method (hardware device which is used) is selected at compile time. The
///             tree uses a vector of nodes rather than a linked list of nodes since the number of nodes can
///             be easily estimated and have the nodes next to each other in memory has significantly better 
///             performance than if the nodes were randomly distributed in memory -- as in the linked list 
///             case
/// @tparam     DeviceType  The type fo device to use -- CPU, GPU, PHI
// ----------------------------------------------------------------------------------------------------------
template <haplo::Device DeviceType>
class Tree;

// Specialize for CPUs
template <> class Tree<haplo::Device::CPU> {
public:
    // --------------------------------------------- ALIAS'S ------------------------------------------------
    using node_container = std::vector<Node>;
    // ------------------------------------------------------------------------------------------------------
    node_container  _nodes;                 //!< Nodes for this tree
    short           _score;                 //!< The score of the tree -- how good of a solution it is
    short           _root_index;            //!< The index of the root variable for the tree
    uint8_t         _root_value;           //!< The value of the root -- a 0 or 1
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Defualt constructor for a tree
    // ------------------------------------------------------------------------------------------------------
    Tree() : _root_index(0), _root_value(0), _score(0), _nodes(0) {};
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for the tree, given a default index and the number of nodes which which it
    ///             should determine the optimal variable values
    /// @param[in]  root_index      The index of variable which is the root of the tree 
    /// @param[in]  index_value     The value of the variable the root index represents (a 0 or a 1)
    Tree(short root_index, uint8_t root_value, size_t num_nodes) 
    : _root_index(root_index)                       , 
      _root_value(root_value)                       , 
      _score(0)                                     , 
      _nodes(num_nodes, Node(root_index, root_value))
    {
        std::cout << "CPU tree\n";
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Searches a tree, trying to find the optimal configuration of the possible nodes given the
    ///             root of the tree -- the nodes with the best score
    /// @param[in]  block       The block for which the tree was created to find the optimal solution for
    /// @tparam     BlockType   The type of the block
    // ------------------------------------------------------------------------------------------------------
    template <typename BlockType>
    void solve(const BlockType& block);

    void print() const;    
};

// --------------------------------------------- IMPLEMENTATIONS --------------------------------------------

template <typename BlockType>
void Tree<Device::CPU>::solve(const BlockType& block)
{
    std::cout << "solving\n";
    
    // Create 2 comprison nodes 
    Node node_0(_root_index, _root_value);
    Node node_1(_root_index, _root_value);
    
    // We need to determine the best possible x variables, 
    // so go through all the possible inidices
    for (short x_index = 0; x_index < _nodes.size(); ++x_index) {
         // Set the two comparison nodes to have different x values
         node_0.set_x_index(x_index); node_1.set_x_index(x_index);
         node_0.set_x_value(0); node_1.set_x_value(1);
         
         // Set the index of the corresponding node in the node container
         _nodes[x_index].set_x_index(x_index);
         
         // Get the value of the block element -- 0 or 1 -- since the comptation function 
         // for the node is different depending on what the block element value is
         uint8_t block_element_value = block(_root_index, x_index).value();
         
         // Determine the better node based on the score for the two node
         uint8_t best_value = node_0.score(block_element_value) < node_1.score(block_element_value) ? 0 : 1;
         
         _nodes[x_index].set_x_value(best_value);
    }
}

void Tree<Device::CPU>::print() const 
{
    for (const auto& node : _nodes) {
        std::cout << static_cast<unsigned>(node.y_index()) << " " << static_cast<unsigned>(node.y_value()) << " "
                  << static_cast<unsigned>(node.x_index()) << " " << static_cast<unsigned>(node.x_value()) << "\n";
    }
}

}           // End namespace haplo
    
#endif      // PARAHAPLO_TREE_HPP
