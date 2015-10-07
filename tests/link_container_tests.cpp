// ----------------------------------------------------------------------------------------------------------
/// @file   link_container_tests.cpp
/// @brief  Test suite for parahaplo link container tets
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE LinkContainerTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/link_container_cpu.hpp"

BOOST_AUTO_TEST_SUITE( LinkContainerSuite )

BOOST_AUTO_TEST_CASE( canCopyConstructLink )
{
    haplo::Link link_one;
    link_one.homo_weight()  = 3;
    link_one.hetro_weight() = 4;
    
    haplo::Link link_two(link_one);
    
    BOOST_CHECK( link_two.homo_weight()  == 3 );
    BOOST_CHECK( link_two.hetro_weight() == 4 );
}


BOOST_AUTO_TEST_CASE( canUseLinkContainer )
{
    using link_container = haplo::LinkContainer<haplo::devices::cpu>;
    link_container links;
    
    haplo::Link new_link;
    new_link.homo_weight() = 12;
    
    links.insert(1, 1, new_link);
    
    BOOST_CHECK( links.exists(1, 1) == true  );
    BOOST_CHECK( links.exists(0, 1) == false );
}

BOOST_AUTO_TEST_SUITE_END()
