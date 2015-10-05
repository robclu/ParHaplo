// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo node container cpu class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_NODE_CONTAINER_CPU_HPP
#define PARHAPLO_NODE_CONTAINER_CPU_HPP

#include "node_container.hpp"

#include <vector>

namespace haplo {

// Specialization for the CPU 
template <>
class NodeContainer<devices::cpu> {
public:
    // ------------------------------------------ ALIAS'S ---------------------------------------------------
    using info_container    = std::vector<Node>;
    using link_container    = std::vector<Link>;
    // ------------------------------------------------------------------------------------------------------
private:
    size_t              _nodes;                 //!< The number of nodes
    info_container      _node_info;             //!< Information for each of the nodes
    link_container      _node_links;            //!< Data for each of the nodes
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for creating a node container
    /// @param[in]  nodes   The number of nodes in the container
    // ------------------------------------------------------------------------------------------------------
    NodeContainer(const size_t nodes) noexcept
    : _nodes(nodes), _node_info(nodes), _node_links((nodes - 1) * nodes / 2)
    {
        size_t position = 0;
        for (auto& info : _node_info) info.position() = position++;
    }
   
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Destructor for node container class
    // ------------------------------------------------------------------------------------------------------
    ~NodeContainer() noexcept {}

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Move constructor
    /// @param[in]  other  The other node container to assign this container to
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t ContainerInstance>
    NodeContainer(NodeContainer<ContainerInstance>&& other) noexcept
    : _nodes(other._nodes), _node_info(std::move(other._node_info)), _node_links(std::move(other._node_links))    
    {
        other._nodes   = 0;
    } 

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Returns the number of nodes in the container 
    /// @return     The number of nodes in the node container 
    // ------------------------------------------------------------------------------------------------------
    inline size_t num_nodes() const { return _nodes; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the node information (nodes) for the container 
    /// @return     The nodes for the container
    // ------------------------------------------------------------------------------------------------------
    inline const info_container& nodes() const { return _node_info; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the node links for the container 
    /// @return     The node links for the container
    // ------------------------------------------------------------------------------------------------------
    inline const link_container& links() const { return _node_links; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Resets the number of nodes in the container (doesn't do any memory removal, just changes
    ///             the number of nodes for the index mapping
    /// @param[in]  new_num_nodes   The new number of nodes in the container
    // ------------------------------------------------------------------------------------------------------
    inline void set_num_nodes(const size_t new_num_nodes) { _nodes = new_num_nodes; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the position in the haplotype a node represents
    /// @param[in]  index   The index of the node
    // ------------------------------------------------------------------------------------------------------
    inline size_t& haplo_pos(const size_t index) { return _node_info[index].position(); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a node (it's information which can be compared to links and other nodes
    /// @param[in]  index   The index of the node to get
    /// @return     The node at the given index
    // ------------------------------------------------------------------------------------------------------
    const Node& operator[](size_t index) const { return _node_info[index]; };
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of the node at container position index (const refernce)
    /// @param[in]  index       The index of the node to get the weight for
    /// @return     The weight of the node at index i
    // ------------------------------------------------------------------------------------------------------
    inline const size_t& weight(const size_t index) const { return _node_info[index].weight(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of the node at container position index (non const refernce)
    /// @param[in]  index       The index of the node to get the weight for
    /// @return     The weight of the node at index i
    // ------------------------------------------------------------------------------------------------------
    inline size_t& weight(const size_t index) { return _node_info[index].weight(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes (the link between nodes 0 and 1 is the same as the link 
    ///             between nodes 1 and 0 -- always use the smaller index as the first parameter
    /// @param[in]  node_idx_a    The index of the first node (this must always be < node_idx_b)
    /// @param[in]  node_idx_b    The index of the second node (this must always be > node_idx_a)
    /// @return     The link between the nodes
    // ------------------------------------------------------------------------------------------------------
    inline const Link& link(const size_t node_idx_a, const size_t node_idx_b) const
    {
        const size_t link_idx = _node_links.size() - ((_nodes - node_idx_a) * (_nodes - node_idx_a - 1) / 2) 
                                + node_idx_b - node_idx_a - 1;
        return _node_links[link_idx];
        
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a non-cost reference to link between two nodes (the link between nodes 0 and 1 is 
    ///             the same as the link between nodes 1 and 0 -- always use the smaller index as the first 
    ///             parameter
    /// @param[in]  node_idx_a    The index of the first node (this must always be < node_idx_b)
    /// @param[in]  node_idx_b    The index of the second node (this must always be > node_idx_a)
    /// @return     The link between the nodes
    // ------------------------------------------------------------------------------------------------------
    inline Link& link(const size_t node_idx_a, const size_t node_idx_b)
    {
        const size_t link_idx = _node_links.size() - ((_nodes - node_idx_a) * (_nodes - node_idx_a - 1) / 2) 
                                 + node_idx_b - node_idx_a - 1;
        return _node_links[link_idx]; 
    }
};

}           // End namespace haplo
#endif      // PARAHAPLO_NODE_CONTAINER_CPU_HPP

