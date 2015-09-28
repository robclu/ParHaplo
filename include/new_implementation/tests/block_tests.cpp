// ----------------------------------------------------------------------------------------------------------
/// @file   block_tests.cpp
/// @brief  Test suite for parahaplo block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../block.hpp"

constexpr char* input_1 = "input_unfiltered_1.txt";
constexpr char* input_2 = "input_unfiltered_2.txt";
constexpr char* input_3 = "input_singleton_rows.txt";

BOOST_AUTO_TEST_SUITE( BlockSuite )
    
BOOST_AUTO_TEST_CASE( canCreateABlock1 )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::Block<10, 12> block(input_1);
    
    block.print();
    
    BOOST_CHECK( block(0, 0 ) == 1 );
    BOOST_CHECK( block(0, 1 ) == 1 );
    BOOST_CHECK( block(1, 0 ) == 2 );
    BOOST_CHECK( block(1, 1 ) == 0 );
    BOOST_CHECK( block(2, 2 ) == 0 );
    BOOST_CHECK( block(5, 3 ) == 1 );
    BOOST_CHECK( block(8, 4 ) == 0 );  
    BOOST_CHECK( block(8, 5 ) == 2 );  
    BOOST_CHECK( block(9, 8 ) == 0 );
    BOOST_CHECK( block(9, 11) == 1 );
}

BOOST_AUTO_TEST_CASE( canDetermineSingletonColumns )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::Block<6, 8> block(input_3);
    
    block.print();
}

BOOST_AUTO_TEST_SUITE_END()
