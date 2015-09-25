// ----------------------------------------------------------------------------------------------------------
/// @file   block_tests.cpp
/// @brief  Test suite for parahaplo block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "block.hpp"

constexpr char* input_file = "input.txt";

BOOST_AUTO_TEST_SUITE( BlockSuite )
    
BOOST_AUTO_TEST_CASE( canCreateABlock )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::Block<10, 7> block(input_file);
    
    BOOST_CHECK( sizeof(block)  == 18 );
    BOOST_CHECK( block(0, 0)    == 0 );
    BOOST_CHECK( block(0, 1)    == 2 );
    BOOST_CHECK( block(1, 0)    == 1 );
    BOOST_CHECK( block(1, 1)    == 0 );
    BOOST_CHECK( block(2, 2)    == 0 );
    BOOST_CHECK( block(5, 3)    == 0 );
    BOOST_CHECK( block(8, 5)    == 2 );  
}

BOOST_AUTO_TEST_SUITE_END()
