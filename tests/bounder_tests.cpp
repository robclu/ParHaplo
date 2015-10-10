// ----------------------------------------------------------------------------------------------------------
/// @file   bounder_tests.cpp
/// @brief  Test suite for parahaplo bound calculation functor tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BounderTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/bounder_cpu.hpp"
#include "../haplo/node_selector_cpu.hpp"
#include "../haplo/node_manager_cpu.hpp"

BOOST_AUTO_TEST_SUITE( BounderSuite )

BOOST_AUTO_TEST_CASE( canCorrectlyDetermineBounds )
{
    const size_t num_nodes = 3;
    
    using link_container = haplo::LinkContainer<haplo::devices::cpu>;
    link_container links;
    
    // Create a node container with 3 nodes
    haplo::NodeContainer<haplo::devices::cpu> nodes(num_nodes);
    
    // Create a node manager for search nodes
    haplo::NodeManager<haplo::devices::cpu> node_manager(6);
    
    // Set the weights
    for (size_t i = 0; i < num_nodes; ++i) nodes.weight(i) = 1;
  
    // Set some links manually 
    links.insert(0, 1, haplo::Link(1, 2));  // 01 : 1 2
    links.insert(0, 2, haplo::Link(1, 0));  // 02 : 1 0
    links.insert(1, 2, haplo::Link(1, 0));  // 12 : 1 0
    
    // This should order the nodes as (choosing 1 to be the start)
    //  1 - 2 - 0

    // Set some parameters of the search nodes in the node manager
    node_manager.node(0).left()  = 1;
    node_manager.node(0).right() = 2;
    
    node_manager.node(1).root() = 0;
    size_t next_index = node_manager.get_next_node();
    
    node_manager.node(1).left()  = next_index;
    node_manager.node(1).right() = next_index + 1;
    
    node_manager.node(3).root() = 1;
    node_manager.node(4).root() = 1;
    node_manager.node(4).set_type(1);
    
    // Create the selector -- start at node 1
    haplo::NodeSelector<haplo::devices::cpu> selector(nodes, links, 1);

    // Create the bounder object 
    haplo::Bounder<haplo::devices::cpu> bound_calculator(nodes, links, node_manager.nodes());
    
    auto idx        = selector.select_node();
    auto node       = nodes[idx];
    auto end_idx    = selector.last_selected_index();
    
    auto bounds = bound_calculator.calculate<1>(node_manager.node(1), node.position(), end_idx);

    BOOST_CHECK( bounds.lower == 0 );
    BOOST_CHECK( bounds.upper == 1 );
    
    // Select a new node
    idx     = selector.select_node();
    node    = nodes[idx];
    end_idx = selector.last_selected_index();
    
    bounds = bound_calculator.calculate<4>(node_manager.node(3), node.position(), end_idx);
    
    BOOST_CHECK( bounds.lower == 2 );
    BOOST_CHECK( bounds.upper == 1 );
    
    bounds = bound_calculator.calculate<4>(node_manager.node(4), node.position(), end_idx);
    
    BOOST_CHECK( bounds.lower == 1 );
    BOOST_CHECK( bounds.upper == 2 );
}

BOOST_AUTO_TEST_SUITE_END()
