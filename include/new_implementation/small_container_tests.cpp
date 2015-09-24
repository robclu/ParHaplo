// ----------------------------------------------------------------------------------------------------------
/// @file   small_container_tests.cpp
/// @brief  Test suite for parahaplo small container tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SmallContainerTests
#endif
#include <boost/test/unit_test.hpp>

#include "small_containers.hpp"


BOOST_AUTO_TEST_SUITE( BinaryContainerSuite )
    
BOOST_AUTO_TEST_CASE( canCreateBinaryContainer ) 
{
    haplo::BinaryContainer<12> bits;
    
    // Set some bits
    bits.set(1, 1);
    bits.set(2, 1);
    bits.set(7, 1);
    bits.set(9, 1);
    
    BOOST_CHECK( bits.get(1) == 1 );
    BOOST_CHECK( bits.get(2) == 1 );
    BOOST_CHECK( bits.get(3) == 0 );
    BOOST_CHECK( bits.get(7) == 1 );
    BOOST_CHECK( bits.get(8) == 0 );
    BOOST_CHECK( bits.get(9) == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
