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
/// @tparam     SubBlockType    The type of the sublock from which this tree is derived
// ----------------------------------------------------------------------------------------------------------
template <typename SubBlockType>
class Tree<SubBlockType, devices::cpu> {
public:
    // ----------------------------------------------- ALIAS'S ----------------------------------------------
    using sub_block_type    = SubBlockType;                             // Type of the subblock
    using node_container    = NodeContainer<devices::cpu>;              // Container for the nodes
    using link_container    = LinkContainer<devices::cpu>;              // Container for the links
    using manager_type      = NodeManager<devices::cpu>;                // Manager for the search nodes
    using bounder_type      = Bounder<devices::cpu>;                    // Bound calculator type
    using selector_type     = NodeSelector<devices::cpu>;               // Node selector type
    using atomic_type       = tbb::atomic<size_t>;
    // ------------------------------------------------------------------------------------------------------
private:
    sub_block_type&     _sub_block;                 //!< The sub-block from which the tree is derived
    atomic_type         _start_node;                //!< The node at which to start the search
    atomic_type         _start_node_worst_case;     //!< The worst case value of the start node
    node_container      _nodes;                     //!< The nodes in the tree
    link_container      _links;                     //!< Links between the nodes of the tree
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Basic constructor
    /// @param[in]  sub_block   The sub_block for which this tree solves the haplotypes
    // ------------------------------------------------------------------------------------------------------
    Tree(sub_block_type& sub_block) noexcept 
    : _sub_block(sub_block), _start_node(0), _start_node_worst_case(0), _nodes(0) {}
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor for a tree
    /// @param[in]  sub_block   The sub block for which this tree solves teh haplotype
    /// @param[in]  nodes       The number of nodes in the tree
    // ------------------------------------------------------------------------------------------------------
    Tree(sub_block_type& sub_block, const size_t nodes) noexcept 
    : _sub_block(sub_block), _start_node(0), _start_node_worst_case(0), _nodes(nodes) {}
    
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
    /// @brief      Gets the links of the tree
    /// @return     The links for  the tree
    // ------------------------------------------------------------------------------------------------------
    inline const link_container& links() const { return _links; }
    
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
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    // ------------------------------------------------------------------------------------------------------    
    inline Link& link(const size_t node_idx_lower, const size_t node_idx_upper) 
    {
        return _links.at(node_idx_lower, node_idx_upper);
    }

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the link between two nodes of the tree
    /// @param[in]  node_idx_lower    The index of the lower node (index with a lower value)
    /// @param[in]  node_idx_upper    The index of the upper node (index with a higher value)
    // ------------------------------------------------------------------------------------------------------    
    inline const Link& link(const size_t node_idx_lower, const size_t node_idx_upper) const
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
    /// @brief      Searches the tree for the best solution 
    /// @tparam     BranchCores     The number of cores available for parallel brach search
    /// @tparam     OpCores         The number of cores available for the operations
    // ------------------------------------------------------------------------------------------------------
    template <size_t BranchCores, size_t OpCores>
    void explore();
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
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the aligments after setting the haplotypes
    /// @tparam     BranchCores     The number of cores available for parallel brach search
    /// @tparam     OpCores         The number of cores available for the operations
    // ------------------------------------------------------------------------------------------------------
    template <size_t BranchCores, size_t OpCores>
    void set_alignments();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the values of the NIH columns
    /// @tparam     BranchCores     The number of cores available for parallel brach search
    /// @tparam     OpCores         The number of cores available for the operations
    // ------------------------------------------------------------------------------------------------------
    template <size_t BranchCores, size_t OpCores>
    void determine_nih_nodes();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Updates the counts for which NIH values are best
    /// @param[in]  counts      The counts of the respective options
    /// @param[in]  row_idx     The row index to check against
    /// @param[in]  elem_value  The valu of the element in the sub block at the haplo index and row index
    // ------------------------------------------------------------------------------------------------------
    void update_nih_count(size_t* counts, const size_t row_idx, const uint8_t elem_value);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sets the values of the haplotype at the NIH columns 
    /// @param[in]  haplo_idx   The index of the haplotype to set
    /// @param[in]  value       The value to set the haplotype to 
    // ------------------------------------------------------------------------------------------------------
    void set_haplo_nih_result(const size_t haplo_idx, const size_t value);
};

// -------------------------------------- IMPLEMENTATIONS ---------------------------------------------------

template <typename SubBlockType> template <size_t BranchCores, size_t OpCores>
void Tree<SubBlockType, devices::cpu>::explore()
{  
    manager_type    node_manager(_nodes.num_nodes());                           // Create a node manager
    selector_type   node_selector(_nodes, _links, _start_node);                 // Create a node selector
    bounder_type    bound_calculator(_nodes, _links, node_manager.nodes());     // Create a bound calculator
    atomic_type     min_ubound{0};                                              // Starting min upper bound
    
    // Set the starting upper bound 
    min_ubound = _sub_block._elements;
    
    // DEBUGGING for the moment
    //std::cout << " - - - - - - - EXPLORING TREE - - - - - - -\n";
   
    // For the first node in the tree                    
    auto& root_node = node_manager.node(0);             // Get the root 
    root_node.set_index(_start_node);                   // Set the index of the root node
    root_node.set_type(0);                              // Setting the type to be 0 (a "left" node)
    root_node.left()  = 1; root_node.right() = 2;
  
    //std::cout << "SN : " << _start_node << "\n";
    
    // Start node's upper bound is the total number of elements 
    root_node.upper_bound() = min_ubound - _nodes[0].elements(); root_node.lower_bound() = 0;
  
    //std::cout << "SUB : " << root_node.upper_bound() << "\n";
    
    // Pass the upper bounds to the subnodes
    auto& left_node  = node_manager.node(1);
    auto& right_node = node_manager.node(2);
    
    // Need to do max upper found calculation
    left_node.upper_bound()  = root_node.upper_bound(); left_node.lower_bound()  = 0;
    right_node.upper_bound() = root_node.upper_bound(); right_node.lower_bound() = 0;
    
    // Make left and right point back to root so that we can go backwards out of the recursion
    left_node.root() = 0 ; right_node.root() = 0;
    left_node.set_type(0); right_node.set_type(1);
    
    // Make sure there is enough memory for another level of nodes
    node_manager.add_node_level(2);
    
    // Check that the haplotypes have enough memory
    if (_sub_block._haplo_one.size() < _sub_block._cols) _sub_block._haplo_one.resize(_sub_block._cols);
    if (_sub_block._haplo_two.size() < _sub_block._cols) _sub_block._haplo_two.resize(_sub_block._cols);
    
    // Search the sutrees, start with 2 subtrees -- this runs until the solution is found
    search_subnodes<BranchCores, OpCores>(node_manager, node_selector, bound_calculator, min_ubound, 1, 2);
    
    // Set the value in the haplotype
    _sub_block._haplo_one.set(_nodes[0].position(), 0); _sub_block._haplo_two.set(_nodes[0].position(), 1);
    
    // Set all the alignments 
    set_alignments<BranchCores, OpCores>();
    
    // Determine the values of all the NIH columns
    if (_sub_block._num_nih > 0) determine_nih_nodes<BranchCores, OpCores>();
}

template <typename SubBlockType> template <size_t BranchCores, size_t OpCores>
size_t Tree<SubBlockType, devices::cpu>::search_subnodes(
                                           manager_type&  node_manager     , selector_type& node_selector   , 
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
  
    min_lbound = std::numeric_limits<size_t>::max();                            // Set LB 
    
    std::cout << "Searching -- " << search_idx << " -- " << num_subnodes << "\n";
    
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
                    
                    const size_t last_idx = node_selector.last_search_index() - _sub_block._num_nih - 1;
                    
                    // If the node is not going to be printed, then create children
                    if (node.lower_bound() <= min_ubound && search_idx < last_idx) {
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
    
    // If we do not have a terminating case, then we must recurse
    if (num_branches > 2 && search_idx != node_selector.last_search_index() - _sub_block._num_nih - 1) {
        // Make sure that we have enough space for the next level of nodes
        node_manager.add_node_level(num_branches);
        
        best_index = search_subnodes<BranchCores, OpCores>(node_manager                 ,       
                                                           node_selector                , 
                                                           bound_calculator             ,  
                                                           min_ubound                   , 
                                                           start_index + num_subnodes   , 
                                                           num_branches                 );
    } 

    // Set the haplotypes 
    _sub_block._haplo_one.set(haplo_idx, node_manager.node(best_index).type());
    _sub_block._haplo_two.set(haplo_idx, !node_manager.node(best_index).type());
   
    return node_manager.node(best_index).root();
}

template <typename SubBlockType> template <size_t BranchCores, size_t OpCores>
void Tree<SubBlockType, devices::cpu>::set_alignments()
{
    tbb::concurrent_vector<bool> alignment_set(_sub_block._rows);
    
    // The aligments are set based on the order in which the haplotype position values were determined,
    // so we go through the haplotype nodes and then set the aligments accordingly
    for (size_t node_idx = 0; node_idx <= _nodes.num_nodes() - _sub_block._num_nih; ++node_idx) {
        
        // Get the index of the haplotype, it's value and the number of threads to use
        const size_t  haplo_idx   = _nodes[node_idx].position();
        const uint8_t haplo_value = _sub_block._haplo_one.get(haplo_idx);
        const size_t  elements    = _sub_block._snp_info[haplo_idx].length();
        const size_t  row_start   = _sub_block._snp_info[haplo_idx].start_index();
        const size_t  threads     = BranchCores + OpCores > elements ? elements : BranchCores + OpCores;
        
        tbb::parallel_for( 
            tbb::blocked_range<size_t>(0, threads),
            [&](const tbb::blocked_range<size_t>& thread_ids) 
            {
                for (size_t thread_id = thread_ids.begin(); thread_id != thread_ids.end(); ++thread_id) {
                size_t thread_iters = ops::get_thread_iterations(thread_id, elements, threads); 
              
                    for (size_t it = 0; it < thread_iters; ++it) {
                        size_t row_idx = it * threads + thread_id + row_start;  
                        
                        // Check if the data has the same value as the haplotype
                        if (haplo_value == _sub_block(row_idx, haplo_idx) && 
                            _sub_block(row_idx, haplo_idx) <= 1           &&
                            alignment_set[row_idx - row_start] == false    )
                        {
                            _sub_block._alignments.set(row_idx, 1); alignment_set[row_idx] = true;
                        } else if (haplo_value != _sub_block(row_idx, haplo_idx) &&
                                   _sub_block(row_idx, haplo_idx) <= 1           &&
                                   alignment_set[row_idx - row_start] == false   )
                        {
                            _sub_block._alignments.set(row_idx, 0); alignment_set[row_idx] = true;
                        }
                    }
                }
            }
        ); 
    }
    //for (size_t i = 0; i < _sub_block._alignments.size(); ++i)
    //    std::cout << static_cast<unsigned>(_sub_block._alignments.get(i)) << "\n";
    
}

template <typename SubBlockType> template <size_t BranchCores, size_t OpCores>
void Tree<SubBlockType, devices::cpu>::determine_nih_nodes()
{
    const size_t threads   = BranchCores + OpCores > _sub_block._num_nih 
                                ? _sub_block._num_nih : BranchCores + OpCores;
    const size_t start_idx = _nodes.num_nodes() - _sub_block._num_nih - _sub_block._duplicate_cols.size(); 
   
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, threads),
        [&](const tbb::blocked_range<size_t>& thread_ids)
        {
            for (size_t thread_id = thread_ids.begin(); thread_id != thread_ids.end(); ++thread_id) {
                size_t thread_iters = ops::get_thread_iterations(thread_id, _sub_block._num_nih, threads); 
              
                for (size_t it = 0; it < thread_iters; ++it) {
                    size_t node_idx  = it * threads + thread_id + start_idx;
                    size_t haplo_idx = _nodes[node_idx].position();
                    auto&  snp_info  = _sub_block._snp_info[haplo_idx];
                    
                    // Make sure this column can be modified
                    if (snp_info.type() == 1) {
                        size_t counts[3] = {0, 0, 0};
                        size_t start = snp_info.start_index(), end = snp_info.end_index();
                    
                        for (size_t row_idx = start; row_idx <= end; ++row_idx) {
                            update_nih_count(counts, row_idx, _sub_block(row_idx, haplo_idx)); 
                        }
                        
                        // Determine the max element
                        size_t best_result = counts[0] >= counts[1] ? 0 : 1;
                        best_result = counts[best_result] >= counts[2] ? best_result : 2;
                        
                        // Set the results 
                        set_haplo_nih_result(haplo_idx, best_result);
                    }
                }
            }
        }
    );
}
   
template <typename SubBlockType>
void Tree<SubBlockType, devices::cpu>::update_nih_count(size_t*       counts    , 
                                                        const size_t  row_idx   , 
                                                        const uint8_t elem_value)
{
    switch (elem_value) {
        case 0:
            if (_sub_block._alignments.get(row_idx) == 1) {
                ++counts[0]; ++counts[2];
            } else {
               ++counts[0];
            }
            break;
        case 1:
            if (_sub_block._alignments.get(row_idx) == 0) {
                ++counts[1]; ++counts[2];
            } else {
                ++counts[1];
            }
            break;
        default: break;
    }
}

template <typename SubBlockType>
void Tree<SubBlockType, devices::cpu>::set_haplo_nih_result(const size_t haplo_idx, const size_t value)
{
    switch(value) {
        case 0:
            _sub_block._haplo_one.set(haplo_idx, 0);
            _sub_block._haplo_two.set(haplo_idx, 0);
            break;
        case 1:
            _sub_block._haplo_one.set(haplo_idx, 1);
            _sub_block._haplo_two.set(haplo_idx, 0);            
            break;
        case 2:
            _sub_block._haplo_one.set(haplo_idx, 0);
            _sub_block._haplo_two.set(haplo_idx, 1);                
            break;
        default: break;
    }
}

}           // End namespace haplo
#endif      // PARAHAPLO_TREE_CPU_HPP

