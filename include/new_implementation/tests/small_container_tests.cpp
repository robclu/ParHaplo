// ----------------------------------------------------------------------------------------------------------
/// @file   small_container_tests.cpp
/// @brief  Test suite for parahaplo small container tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SmallContainerTests
#endif
#include <boost/test/unit_test.hpp>

#include "../small_containers.hpp"


BOOST_AUTO_TEST_SUITE( BinaryContainerSuite )
    
BOOST_AUTO_TEST_CASE( canCreateBinaryContainerWith1BitPerElement ) 
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<12> elements;
    
    // Set some elements
    elements.set(1, 1);
    elements.set(2, 1);
    elements.set(7, 1);
    elements.set(9, 1);    
    elements.set(3, 0);
    
    BOOST_CHECK( elements.get(1) == 1 );
    BOOST_CHECK( elements.get(2) == 1 );
    BOOST_CHECK( elements.get(3) == 0 );
    BOOST_CHECK( elements.get(7) == 1 );
    BOOST_CHECK( elements.get(8) == 0 );
    BOOST_CHECK( elements.get(9) == 1 );
    BOOST_CHECK( sizeof(elements) == 2);
}

BOOST_AUTO_TEST_CASE( canCreateBinaryContainerWith2BitsPerElement ) 
{
    // Define the conatiner to use 2 bits per element and 
    // 18 elements so 36 bits = 5 bytes
    haplo::BinaryContainer<18, 2> elements;

    // Set some elements
    elements.set(1 , 1);
    elements.set(3 , 2);
    elements.set(7 , 1);
    elements.set(14, 3);
    elements.set(8 , 0);
    
    BOOST_CHECK( elements.get(1)  == 1 );
    BOOST_CHECK( elements.get(2)  == 0 );
    BOOST_CHECK( elements.get(3)  == 2 );
    BOOST_CHECK( elements.get(7)  == 1 );
    BOOST_CHECK( elements.get(8)  == 0 );
    BOOST_CHECK( elements.get(14) == 3 );
    BOOST_CHECK( sizeof(elements) == 5);
}

BOOST_AUTO_TEST_SUITE_END()
