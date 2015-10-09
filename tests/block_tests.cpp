// ----------------------------------------------------------------------------------------------------------
/// @file   block_tests.cpp
/// @brief  Test suite for parahaplo block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/block.hpp"

static constexpr const char* input_1 = "input_files/input_zero.txt";

BOOST_AUTO_TEST_SUITE( BlockSuite )
    
BOOST_AUTO_TEST_CASE( canCreateABlockAndGetData )
{
    // Define for 28 elements with 1 core for each dimension
    using block_type = haplo::Block<28>; 
    
    block_type block(input_1);
        
    BOOST_CHECK( block(0, 0 ) == 1 );
    BOOST_CHECK( block(0, 1 ) == 0 );
    BOOST_CHECK( block(1, 0 ) == 3 );
    BOOST_CHECK( block(1, 1 ) == 1 );
    BOOST_CHECK( block(1, 2 ) == 0 );
    BOOST_CHECK( block(2, 1 ) == 0 );
    BOOST_CHECK( block(2, 2 ) == 0 );
    BOOST_CHECK( block(2, 3 ) == 0 );
    BOOST_CHECK( block(3, 1 ) == 0 );
    BOOST_CHECK( block(3, 2 ) == 1 );
    BOOST_CHECK( block(5, 3 ) == 1 );
    BOOST_CHECK( block(8, 4 ) == 0 );  
    BOOST_CHECK( block(8, 5 ) == 2 );  
    BOOST_CHECK( block(8, 10) == 2 );
    BOOST_CHECK( block(9, 9 ) == 2 );
    BOOST_CHECK( block(9, 11) == 0 );
}

BOOST_AUTO_TEST_CASE( canDetermineMonotoneColumns )
{
    // Define for 28 elements with 4 cores for each dimension
    using block_type = haplo::Block<28, 4, 4>; 
    
    block_type block(input_1);    
    
    BOOST_CHECK( block.is_monotone(0)  == true  ); 
    BOOST_CHECK( block.is_monotone(1)  == false ); 
    BOOST_CHECK( block.is_monotone(2)  == false ); 
    BOOST_CHECK( block.is_monotone(3)  == false ); 
    BOOST_CHECK( block.is_monotone(4)  == false ); 
    BOOST_CHECK( block.is_monotone(5)  == true  ); 
    BOOST_CHECK( block.is_monotone(6)  == true  ); 
    BOOST_CHECK( block.is_monotone(7)  == false ); 
    BOOST_CHECK( block.is_monotone(8)  == true  ); 
    BOOST_CHECK( block.is_monotone(9)  == true  ); 
    BOOST_CHECK( block.is_monotone(10) == false ); 
    BOOST_CHECK( block.is_monotone(11) == false ); 
}

BOOST_AUTO_TEST_CASE( canDetermineSplittableColumns )
{
    using block_type = haplo::Block<28, 4, 4>;
    
    block_type block(input_1);
    
    BOOST_CHECK( block.num_subblocks() == 4  );
    BOOST_CHECK( block.subblock(0)     == 1  );
    BOOST_CHECK( block.subblock(1)     == 3  );
    BOOST_CHECK( block.subblock(2)     == 4  );
    BOOST_CHECK( block.subblock(3)     == 11 );
}

BOOST_AUTO_TEST_SUITE_END()
