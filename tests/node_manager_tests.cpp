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
    // Search node manager for 10 nodes
    haplo::NodeManager<haplo::devices::cpu> node_manager(10);
    
    // Get a new node -- not dynamic allocation
    haplo::SearchNode* new_node = node_manager.get_new_node();
    
}

BOOST_AUTO_TEST_SUITE_END()
