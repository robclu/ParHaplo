// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo tree class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_TREE_CPU_HPP
#define PARHAPLO_TREE_CPU_HPP

#include "bounder_cpu.hpp"
#include "node_manager_cpu.hpp"
#include "node_selector_cpu.hpp"
#include "tree.hpp"

#include <iostream>
#include <limits>

namespace haplo {
namespace links {
    
static constexpr uint8_t homo   = 0x00;
static constexpr uint8_t hetro  = 0x01;

}               // End namespace links

// Update atomic varibale to min
template <typename T1, typename T2>
void atomic_min_update(tbb::atomic<T1>& atomic_var, T2 value)
{
    T1 state;
    do {
        state = atomic_var;         // Capture state
        if (state <= value) break;  // Exit earlt
    } while (atomic_var.compare_and_swap(value, state) != state);
}

// ----------------------------------------------------------------------------------------------------------
/// @class      Tree    
/// @brief      Holds nodes which can then be searched to find the optimal haplotypes
/// @tparam     DeviceType  The type of device to use the node on -- so that we can optimize functions for the
///             different implementations and so that each one can be excluded from compilation if necessary
// ----------------------------------------------------------------------------------------------------------
template <>
class Tree<devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using node_container    = NodeContainer<devices::cpu>;              // Container for the nodes
    using link_container    = LinkContainer<devices::cpu>;              // Container for the links
    using manager_type      = NodeManager<devices::cpu>;                // Manager for the search nodes
    using bounder_type      = Bounder<devices::cpu>;                    // Bound calculator type
    using selector_type     = NodeSelector<devices::cpu>;               // Node selector type
    using atomic_type       = tbb::atomic<size_t>;
    using atomic_binvec     = tbb::concurrent_vector<bool>;             // Atomic binary vector
    // ------------------------------------------------------------------------------------------------------
private:
    atomic_type         _num_nih;                   //!< Number of NIH columns
    atomic_type         _start_node;                //!< The node at which to start the search
    atomic_type         _start_node_worst_case;     //!< The worst case value of the start node
    node_container      _nodes;                     //!< The nodes in the tree
    link_container      _links;                     //!< Links between the nodes of the tree
    atomic_binvec       _haplo_one;                 //!< The one haplotpye (h)
    atomic_binvec       _haplo_two;                 //!< The other haplotype (h')
    atomic_binvec       _alignments;                //!< The alignemnts of the reads
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Default constructor
    // ------------------------------------------------------------------------------------------------------
    Tree() noexcept 
    : _num_nih(0), _start_node(0), _start_node_worst_case(0), _nodes(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for a tree
    /// @param[in]  nodes   The number of nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    Tree(const size_t nodes) noexcept 
    : _num_nih(0), _start_node(0), _start_node_worst_case(0), _nodes(nodes) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Desctructor
    // ------------------------------------------------------------------------------------------------------
    ~Tree() noexcept {}
  
    // ------------------------------------------------------------------------------------------------------
    /// @brief      The mazimum worst case value for the tree
    /// @return     A reference to the maximim worst case value 
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& max_worst_case() { return _start_node_worst_case; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      The index of the start node
    /// @return     A reference to the start node index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& start_node() { return _start_node; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the size of the tree (the number of nodes
    /// @return     The size of of the tree
    // ------------------------------------------------------------------------------------------------------
    inline size_t size() const { return _nodes.num_nodes(); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Resizes the tree to a certain number of nodes
    /// @param[in]  num_nodes   The number of nodes to create for the tree
    // ------------------------------------------------------------------------------------------------------
    inline void resize(const size_t num_nodes) 
    {
        if (_nodes.num_nodes() != num_nodes) _nodes.resize(num_nodes);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the nodes of the tree
    /// @return     The nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    inline const node_container& nodes() const { return _nodes; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Creates a link for the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)    
    // ------------------------------------------------------------------------------------------------------
    inline void create_link(const size_t node_idx_lower, const size_t node_idx_upper)
    {
        _links.insert(node_idx_lower, node_idx_upper);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the links of the tree
    /// @return     The links for  the tree
    // ------------------------------------------------------------------------------------------------------
    inline const link_container& links() const { return _links; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    /// @tparam     LinkType          The type of the link to get
    // ------------------------------------------------------------------------------------------------------
    template <uint8_t LinkType>
    inline atomic_type& link(const size_t node_idx_lower, const size_t node_idx_upper);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    // ------------------------------------------------------------------------------------------------------    
    inline Link& link(const size_t node_idx_lower, const size_t node_idx_upper) 
    {
        return _links.at(node_idx_lower, node_idx_upper);
    }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the max of a link between two nodes of the tree
    /// @param[in]  node_idx_one    The index of the one of the nodes
    /// @param[in]  node_idx_two    The index of another node
    // ------------------------------------------------------------------------------------------------------
    inline size_t link_max(const size_t node_idx_one, const size_t node_idx_two)
    {
        return _links.exists(node_idx_one, node_idx_two) ? _links.at(node_idx_one, node_idx_two).max() : 0;
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the min of a link between two nodes of the tree
    /// @param[in]  node_idx_one    The index of the one of the nodes
    /// @param[in]  node_idx_two    The index of another node
    // ------------------------------------------------------------------------------------------------------
    inline size_t link_min(const size_t node_idx_one, const size_t node_idx_two)
    {
        return _links.exists(node_idx_one, node_idx_two) ? _links.at(node_idx_one, node_idx_two).min() : 0;
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a node from the tree
    /// @param[in]  idx     The index of the node
    // ------------------------------------------------------------------------------------------------------
    inline Node& node(const size_t idx) { return _nodes[idx]; }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets a node from the tree
    /// @param[in]  idx     The index of the node
    // ------------------------------------------------------------------------------------------------------
    inline const Node& node(const size_t idx) const { return _nodes[idx]; }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of a node 
    /// @param[in]  idx     The index of the node
    /// @return     The weight of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& node_weight(const size_t idx) { return _nodes.weight(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the worst case value of a node
    /// @param[in]  idx     The index of the node
    /// @return     The worst case value of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& node_worst_case(const size_t idx) { return _nodes.worst_case_value(idx); }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the worst case value of a node
    /// @param[in]  idx     The index of the node
    /// @return     The worst case value of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type node_worst_case(const size_t idx) const { return _nodes.worst_case_value(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the haplotype position of a node -- the position in the haplotype a node represents
    /// @param[in]  node_idx    The index of the node to get the haplotype position of
    /// @return     The position the node represents in the haplotype
    // ------------------------------------------------------------------------------------------------------
    inline atomic_type& node_haplo_pos(const size_t node_idx) { return _nodes.haplo_pos(node_idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the weight of a node 
    /// @param[in]  idx     The index of the node
    /// @return     The weight of the node at the index
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& node_weight(const size_t idx) const { return _nodes.weight(idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the haplotype position of a node -- the position in the haplotype a node represents
    /// @param[in]  node_idx    The index of the node to get the haplotype position of
    /// @return     The position the node represents in the haplotype
    // ------------------------------------------------------------------------------------------------------
    inline const atomic_type& node_haplo_pos(const size_t node_idx) const { return _nodes.haplo_pos(node_idx); }
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Searches the tree for the best solution 
    /// @param[in]  min_upper_bound The min (or starting) upper bound
    /// @param[in]  num_nih         The number of NIH columns (that don't need ot be searched)
    /// @tparam     BranchCores     The number of cores available for parallel brach search
    /// @tparam     OpCores         The number of cores available for the operations
    // ------------------------------------------------------------------------------------------------------
    template <size_t BranchCores, size_t OpCores>
    void explore(const size_t max_upped_bound, const size_t num_nih);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts the solution into a binary vector from the unordered nodes
    // ------------------------------------------------------------------------------------------------------
    void format_haplotypes();
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Moves down the sub-nodes of the current root node of a subtree tree
    /// @param[in]  node_manager    The manager of the nodes
    /// @param[in]  node_selector   The selector for the nodes
    /// @param[in]  bounder         The bound calculator object
    /// @param[in]  min_upper_bound The lowest upper bound so far
    /// @param[in]  start_index     The index of the start node in the search nodes
    /// @param[in]  num_subnodes    The number of subnodes to search
    /// @tparam     BranchCores     The number of cores available for parallel brach search
    /// @tparam     OpCores         The number of cores available for the operations
    /// @return     The index of the optimal node from the previous iteration 
    // ------------------------------------------------------------------------------------------------------
    template <size_t BranchCores, size_t OpCores>
    size_t search_subnodes(manager_type&  node_manager, selector_type& node_selector  ,
                           bounder_type&  bounder     , atomic_type&   min_upper_bound,
                           const size_t   start_index , const size_t   num_subnodes );
};

// -------------------------------------- IMPLEMENTATIONS ---------------------------------------------------

template <>
inline tbb::atomic<size_t>& Tree<devices::cpu>::link<links::homo>(const size_t node_idx_lower, 
                                                                  const size_t node_idx_upper)
{
    return _links.at(node_idx_lower, node_idx_upper).homo_weight();
}

template <>
inline tbb::atomic<size_t>& Tree<devices::cpu>::link<links::hetro>(const size_t node_idx_lower, 
                                                                   const size_t node_idx_upper)
{
    return _links.at(node_idx_lower, node_idx_upper).hetro_weight();
}

template <size_t BranchCores, size_t OpCores>
void Tree<devices::cpu>::explore(const size_t min_upper_bound, const size_t num_nih) 
{  
    manager_type    node_manager(_nodes.num_nodes());                           // Create a node manager
    selector_type   node_selector(_nodes, _links, _start_node);                 // Create a node selector
    bounder_type    bound_calculator(_nodes, _links, node_manager.nodes());     // Create a bound calculator
    atomic_type     min_ubound{0};                                              // Starting min upper bound
    
    // Set the upper bound and the number of NIH columns (nodes)
    min_ubound = min_upper_bound; _num_nih = num_nih;
    
    // DEBUGGING for the moment
    std::cout << " - - - - - - - EXPLORING TREE - - - - - - -\n";
   
    // For the first node in the tree                    
    auto& root_node = node_manager.node(0);             // Get the root 
    root_node.set_index(_start_node);                   // Set the index of the root node
    root_node.set_value(0);                             // Setting the value to 0
    root_node.left()  = 1; root_node.right() = 2;
  
    std::cout << "SN : " << _start_node << "\n";
    
    // Start node's upper bound is the total number of elements 
    root_node.upper_bound() = min_upper_bound - _nodes[0].elements(); root_node.lower_bound() = 0;
  
    std::cout << "SUB : " << root_node.upper_bound() << "\n";
    
    // Pass the upper bounds to the subnodes
    auto& left_node  = node_manager.node(1);
    auto& right_node = node_manager.node(2);
    
    // Need to do max upper found calculation
    left_node.upper_bound()  = root_node.upper_bound(); left_node.lower_bound()  = 0;
    right_node.upper_bound() = root_node.upper_bound(); right_node.lower_bound() = 0;
    
    // Make left and right point back to root so that we can go backwards out of the recursion
    left_node.root() = 0; right_node.root() = 0;
    left_node.set_type(0); right_node.set_type(1);
    
    // Make sure there is enough memory for another level of nodes
    node_manager.add_node_level(2);
    
    // Search the sutrees, start with 2 subtrees -- this runs until the solution is found
    search_subnodes<BranchCores, OpCores>(node_manager, node_selector, bound_calculator, min_ubound, 1, 2);
}

template <size_t BranchCores, size_t OpCores>
size_t Tree<devices::cpu>::search_subnodes(manager_type&  node_manager     , selector_type& node_selector   , 
                                           bounder_type&  bound_calculator , atomic_type&   min_ubound      ,
                                           const size_t   start_index      , const size_t   num_subnodes    )
{
    // Check how many branch cores we need
    const size_t branch_cores = BranchCores > num_subnodes ? num_subnodes : BranchCores;
    atomic_type  num_branches{0};                                               // Branches to search
    atomic_type  min_lbound{0};                                                 // Best lower bound
    atomic_type  best_index{0};                                                 // Index of best node
    atomic_type  best_value{0};                                                 // Index of best node
    const size_t search_idx   = node_selector.select_node();                    // Index in node array
    const size_t haplo_idx    = _nodes[search_idx].position();                  // Haplo var index
  
    std::cout << search_idx << "\n";
    
    min_lbound = std::numeric_limits<size_t>::max();                            // Set LB 
    
    // Get the index of the 
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, branch_cores),
        [&](const tbb::blocked_range<size_t>& threads)
        {
            for (size_t thread_id = threads.begin(); thread_id != threads.end(); ++thread_id) {
                size_t thread_iters = ops::get_thread_iterations(thread_id, num_subnodes, branch_cores);
                for (size_t it = 0; it < thread_iters; ++it) {
                    const size_t node_idx   = start_index + it * branch_cores + thread_id;
                    auto& node              = node_manager.node(node_idx);              // Node to search
                    
                    constexpr size_t bound_threads = OpCores / BranchCores == 0 
                                                ? 1 : OpCores / BranchCores;
                                                
                    // Get the bounds for the node and update them            
                    auto bounds = bound_calculator.calculate<bound_threads>(node, haplo_idx, search_idx);
                    node.upper_bound() -= (_nodes[search_idx].elements() - bounds.lower);
                    node.lower_bound() += bounds.lower;
                    
                    std::cout << node_idx << " " << node.lower_bound() << " " << node.upper_bound() << 
                         " " << _nodes[search_idx].elements() << " BNDS\n";
                    
                    // If the node is not going to be printed, then create children
                    if (node.lower_bound() <= min_ubound && search_idx != node_selector.last_search_index() - 1) {
                        size_t left_child_idx = node_manager.get_next_node();
                        auto& left_child  = node_manager.node(left_child_idx);
                        auto& right_child = node_manager.node(left_child_idx + 1);
                        
                        // Set the start bounds of the left child node
                        left_child.set_bounds(node.bounds());
                        right_child.set_bounds(node.bounds());
                  
                        // Make the children point back to this node
                        left_child.root() = node_idx; right_child.root() = node_idx;
                        left_child.set_type(0); right_child.set_type(1);
                   
                        num_branches.fetch_and_add(2);                  // 2 more branches next it
                    }
                    atomic_min_update(min_ubound, node.upper_bound());
                    atomic_min_update(min_lbound, node.lower_bound());
                        
                    if (node.lower_bound() == min_lbound) {
                        best_index = node_idx; best_value = node.type();
                    }
                }
            }
        }
    );
   
    std::cout << "UB : " << min_ubound << "\n";
    std::cout << "LB : " << min_lbound << "\n";
    std::cout << "LI : " << node_selector.last_search_index() << "\n";
    
    // If we do not have a terminating case, then we must recurse
    if (num_branches > 2 && search_idx != node_selector.last_search_index() - _num_nih - 1) {
        // Make sure that we have enough space for the next level of nodes
        node_manager.add_node_level(num_branches);
        
        best_index = search_subnodes<BranchCores, OpCores>(node_manager                 ,       
                                                           node_selector                , 
                                                           bound_calculator             ,  
                                                           min_ubound                   , 
                                                           start_index + num_subnodes   , 
                                                           num_branches                 );
    } 
    // Otherwise set the best value in the haplo node
    _nodes[search_idx].set_haplo_value(node_manager.node(best_index).type());
   
    std::cout << "Index      : "  << haplo_idx << "\n";
    std::cout << "Num NIH    : "  << _num_nih  << "\n";
    std::cout << "Best Value : "  << static_cast<unsigned>(node_manager.node(best_index).type()) << "\n";
    std::cout << "\n";
    return node_manager.node(best_index).root();
}

void Tree<devices::cpu>::format_haplotypes() 
{
    // Go through the nodes and put the correct position into the binary vector 
    _haplo_one.resize(_nodes.num_nodes()); _haplo_two.resize(_nodes.num_nodes());
    
    for (auto& node : _nodes.nodes()) {
        if (node.type() == 1) {                                             // Node is IH
            _haplo_one[node.position()] = node.haplo_value();
            _haplo_two[node.position()] = !node.haplo_value();
        } else {                                                            // Node is NIH
            _haplo_one[node.position()] = node.haplo_value();
            _haplo_two[node.position()] = node.haplo_value();            
        }
    }
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_CPU_HPP

