// ----------------------------------------------------------------------------------------------------------
/// @file   node_manager_tests.cpp
/// @brief  Test suite for parahaplo node manager tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE NodeManagerTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/node_manager_cpu.hpp"

BOOST_AUTO_TEST_SUITE( NodeManagerSuite )
    
BOOST_AUTO_TEST_CASE( canCreateANodeManager )
{
    // Node manager for 10 nodes
    haplo::NodeManager<haplo::devices::cpu> node_manager(10);
    
    BOOST_CHECK( node_manager.num_nodes() == 10 );
}

BOOST_AUTO_TEST_CASE( canCreateDefaultManagerAndResize )
{
    haplo::NodeManager<haplo::devices::cpu> node_manager;
    node_manager.resize(20);
    
    BOOST_CHECK( node_manager.num_nodes() == 20 ); 
}

BOOST_AUTO_TEST_CASE( canCreateLevelsOfNodesCorectly )
{
    haplo::NodeManager<haplo::devices::cpu> node_manager(4);   
    
    // Add a level of nodes with 2 nodes -- in this case this does nothing
    node_manager.add_node_level(2);
    
    // Default to 3 nodes min, so next node in node 4
    auto next_index  = node_manager.get_next_node();
    auto& left_node  = node_manager.node(next_index);
    auto& right_node = node_manager.node(next_index + 1);
    
    left_node.set_index(4);
    right_node.set_index(5);
        
    BOOST_CHECK( node_manager.node(3).index() == 4 );
    BOOST_CHECK( node_manager.node(4).index() == 5 );
    
    // Ask for a new level with 4 nodes (so 7 nodes total)
    // only 4 allocated initially so manager should re-allocate
    node_manager.add_node_level(4);
    
    next_index = node_manager.get_next_node();
    auto& left_node_1  = node_manager.node(next_index);
    auto& right_node_1 = node_manager.node(next_index + 1);    

    left_node_1.set_index(9);
    right_node_1.set_index(7);
    
    BOOST_CHECK( node_manager.node(5).index()  == 9 );
    BOOST_CHECK( node_manager.node(6).index()  == 7 );
    BOOST_CHECK( node_manager.num_nodes()      == 16  );
}


BOOST_AUTO_TEST_SUITE_END()
