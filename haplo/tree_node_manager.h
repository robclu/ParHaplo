// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree node manager for the gpu implementation -- managers where the nodes
///         are placed in memory
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_NODE_MANAGER_H
#define PARHAPLO_TREE_NODE_MANAGER_H

#include "tree_node.h"

#ifdef __CUDACC__
    #define CUDA_HD __host__ __device__
    #define CUDA_H  __host__
    #define CUDA_D  __device__
    #define SHARED  __shared__
#else
    #define CUDA_HD
    #define CUDA_H
    #define CUDA_D
    #define SHARED
#endif

namespace haplo {

// ----------------------------------------------------------------------------------------------------------
/// @brief      Manages the nodes in the GPU search tree
// ----------------------------------------------------------------------------------------------------------
class NodeManagerGpu {
private:
    size_t      _num_nodes;             //!< The number of nodes in the tree
    size_t      _next_node_index;       //!< Index of the next node
public:
    TreeNode*   nodes;                 //!< The nodes to manage
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor -- min 4 nodes
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    explicit NodeManagerGpu(const size_t num_nodes) noexcept
    : _num_nodes(num_nodes), _next_node_index(1) {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the node manager -- number of nodes
    /// @return     The number of nodes being managed
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline size_t& num_nodes() { return _num_nodes; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Resizes the node container 
    /// @param[in]  num_nodes   The number of nodes to resize to
    // ------------------------------------------------------------------------------------------------------
    //void resize(const size_t num_nodes) 
    //{
    //    size_t elements = _nodes.size();
    //    _nodes.reserve(num_nodes);
    //    while (elements++ < num_nodes) _nodes.push_back(SearchNode());
    //}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the next node, adn makes space for another after it -- the left and
    ///             subnodes of the tree
    /// @return     The index of the next available subnode
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline TreeNode& get_next_node() 
    { 
        // Very bare -- need to improve
        #ifdef __CUDAACC__
            return nodes[atomicAdd(&_next_node_index, 2)];
        #else 
            size_t before = _next_node_index;
            _next_node_index += 2;
            return nodes[before];
        #endif
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the next node, adn makes space for another after it -- the left and
    ///             subnodes of the tree
    /// @return     The index of the next available subnode
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline TreeNode& node(const size_t i) { return nodes[i]; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the index of the next node, adn makes space for another after it -- the left and
    ///             subnodes of the tree
    /// @return     The index of the next available subnode
    // ------------------------------------------------------------------------------------------------------
    CUDA_HD
    inline TreeNode* node_ptr(const size_t i) { return &nodes[i]; }
    
    CUDA_HD
    inline const TreeNode* node_ptr(const size_t i) const { return &nodes[i]; }

};

}                   // End namespace haplo
#endif              // PARHAPLO_SEARCH_NODE_HPP
