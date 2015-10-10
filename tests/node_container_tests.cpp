// ----------------------------------------------------------------------------------------------------------
/// @file   node_container_tests.cpp
/// @brief  Test suite for parahaplo node container tets
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE NodeContainerTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/node_container_cpu.hpp"

BOOST_AUTO_TEST_SUITE( NodeContainerSuite )

BOOST_AUTO_TEST_CASE( canCopyConstructNode )
{
    haplo::Node node_one;
    node_one.weight()   = 3;
    node_one.position() = 4;
    
    haplo::Node node_two(node_one);
    
    BOOST_CHECK( node_two.weight()   == 3 );
    BOOST_CHECK( node_two.position() == 4 );
}

BOOST_AUTO_TEST_CASE( canMoveConstructNode )
{
    haplo::Node node_one;
    node_one.weight()   = 3;
    node_one.position() = 4;
    
    haplo::Node node_two(std::move(node_one));
    
    BOOST_CHECK( node_two.weight()   == 3 );
    BOOST_CHECK( node_two.position() == 4 );
}

BOOST_AUTO_TEST_CASE( canCreateNodeContainer )
{
    haplo::NodeContainer<haplo::devices::cpu> nodes(12);
    
    BOOST_CHECK( nodes.num_nodes() == 12 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetNodePosition )
{
    haplo::NodeContainer<haplo::devices::cpu> nodes(12); 
    
    BOOST_CHECK( nodes.haplo_pos(2) == 2 );
    
    nodes.haplo_pos(2) = 4;
    
    BOOST_CHECK( nodes.haplo_pos(2) == 4 );
}

BOOST_AUTO_TEST_CASE( canSetAndGetNodeContainerWeights )
{
    haplo::NodeContainer<haplo::devices::cpu> nodes(12);
 
    nodes.weight(0)  = 2;
    nodes.weight(3)  = 4;
    nodes.weight(11) = 17;
    
    BOOST_CHECK( nodes.weight(0)  == 2 );
    BOOST_CHECK( nodes.weight(1)  == 1 );
    BOOST_CHECK( nodes.weight(2)  == 1 );
    BOOST_CHECK( nodes.weight(3)  == 4 );
    BOOST_CHECK( nodes.weight(4)  == 1 );
    BOOST_CHECK( nodes.weight(5)  == 1 );
    BOOST_CHECK( nodes.weight(6)  == 1 );
    BOOST_CHECK( nodes.weight(7)  == 1 );
    BOOST_CHECK( nodes.weight(8)  == 1 );
    BOOST_CHECK( nodes.weight(9)  == 1 );
    BOOST_CHECK( nodes.weight(10) == 1 );
    BOOST_CHECK( nodes.weight(11) == 17 );
}


BOOST_AUTO_TEST_CASE( canMoveNodeContainer )
{
    using node_container = haplo::NodeContainer<haplo::devices::cpu>;
    node_container nodes(4);
    
    nodes.weight(0)  = 2;
    nodes.weight(2)  = 7;    
    
    node_container new_nodes(std::move(nodes));
    
    BOOST_CHECK( nodes.num_nodes()     == 0 );
    BOOST_CHECK( new_nodes.weight(0)   == 2 );
    BOOST_CHECK( new_nodes.weight(2)   == 7 );
}

BOOST_AUTO_TEST_CASE( canGetWorstCaseValueOfANode )
{
    using node_container = haplo::NodeContainer<haplo::devices::cpu>;
    node_container nodes(4);
       
    nodes.worst_case_value(2) = 12;
    
    BOOST_CHECK( nodes.worst_case_value(2) == 12 );
}

BOOST_AUTO_TEST_SUITE_END()
