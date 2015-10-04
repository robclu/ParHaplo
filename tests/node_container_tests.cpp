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

BOOST_AUTO_TEST_CASE( canSetAndGetNodeContainerLinks )
{
    haplo::NodeContainer<haplo::devices::cpu> nodes(4);
 
    nodes.link(0, 1).homo_weight()   = 1;
    nodes.link(0, 1).hetro_weight()  = 2;
    nodes.link(2, 3).homo_weight()   = 3;
    nodes.link(2, 3).hetro_weight()  = 4;
    
    BOOST_CHECK( nodes.link(0, 1).homo_weight()  == 1 );
    BOOST_CHECK( nodes.link(0, 1).hetro_weight() == 2 );
    BOOST_CHECK( nodes.link(0, 2).homo_weight()  == 0 );
    BOOST_CHECK( nodes.link(0, 2).hetro_weight() == 0 );
    BOOST_CHECK( nodes.link(0, 3).homo_weight()  == 0 );
    BOOST_CHECK( nodes.link(0, 3).hetro_weight() == 0 );
    BOOST_CHECK( nodes.link(1, 2).homo_weight()  == 0 );
    BOOST_CHECK( nodes.link(1, 2).hetro_weight() == 0 );
    BOOST_CHECK( nodes.link(1, 3).homo_weight()  == 0 );
    BOOST_CHECK( nodes.link(1, 3).hetro_weight() == 0 );
    BOOST_CHECK( nodes.link(2, 3).homo_weight()  == 3 );
    BOOST_CHECK( nodes.link(2, 3).hetro_weight() == 4 );
}

BOOST_AUTO_TEST_CASE( canCompareNodesAndLinks )
{
    haplo::NodeContainer<haplo::devices::cpu> nodes(4);
       
    nodes.link(0, 1).homo_weight()   = 1;
    nodes.link(0, 1).hetro_weight()  = 2;
    
    nodes.weight(0)  = 2;
    nodes.weight(2)  = 7;
    
    BOOST_CHECK( nodes[0].value() == nodes.link(0, 1).value() );
    BOOST_CHECK( nodes[2].value() >  nodes.link(0, 1).value() );
}

BOOST_AUTO_TEST_CASE( canMoveNodeContainerWithAssigmentOperator )
{
    using node_container = haplo::NodeContainer<haplo::devices::cpu>;
    node_container nodes(4);
       
    nodes.link(0, 1).homo_weight()   = 1;
    nodes.link(0, 1).hetro_weight()  = 2;
    
    nodes.weight(0)  = 2;
    nodes.weight(2)  = 7;    
    
    node_container new_nodes(std::move(nodes));
    
    BOOST_CHECK( nodes.num_nodes()                  == 0 );
    BOOST_CHECK( new_nodes[0].weight()              == 2 );
    BOOST_CHECK( new_nodes.link(0, 1).homo_weight() == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
