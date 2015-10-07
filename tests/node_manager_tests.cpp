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

BOOST_AUTO_TEST_CASE( canGetAndSetNodes )
{
    haplo::NodeManager<haplo::devices::cpu> node_manager;
    node_manager.resize(20);    
    
    // Will be a pointer to the first node
    auto node = node_manager.get_new_node();
    
    node->set_index(4);
    
    BOOST_CHECK( node_manager.node(0).index() == 4 );
}


BOOST_AUTO_TEST_SUITE_END()
