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

static constexpr const char* input_1 = "input_files/input_unfiltered_1.txt";
static constexpr const char* input_2 = "input_files/input_unfiltered_2.txt";
static constexpr const char* input_3 = "input_files/input_singleton_rows.txt";

BOOST_AUTO_TEST_SUITE( BlockSuite )
    
BOOST_AUTO_TEST_CASE( canCreateABlockAndGetData )
{
    // Define a block type with 10 rows, 12 columns, 1 I-dim thread (default) and 1 J-dim thread (default)
    using block_type = haplo::Block<10, 12>; 
    
    block_type block(input_1);
    
    // Debugging 
    //block.print();
    
    // NOTE: Any column with more ones than 0's, the bits get flipped
    
    BOOST_CHECK( block(0, 0 ) == 1 );
    BOOST_CHECK( block(0, 1 ) == 0 );
    BOOST_CHECK( block(1, 0 ) == 2 );
    BOOST_CHECK( block(1, 1 ) == 1 );
    BOOST_CHECK( block(2, 2 ) == 0 );
    BOOST_CHECK( block(5, 3 ) == 1 );
    BOOST_CHECK( block(8, 4 ) == 0 );  
    BOOST_CHECK( block(8, 5 ) == 2 );  
    BOOST_CHECK( block(9, 8 ) == 0 );
    BOOST_CHECK( block(9, 11) == 0 );
}

BOOST_AUTO_TEST_CASE( canDetermineMonotoneColumns )
{
    // Define a block type with 10 rows, 12 columns, 2 I-dim threads and 2 J-dim threads
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    block_type block(input_1);
   
    // Debugging 
    //block.print();
    
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
    // Define a block type with 10 rows, 12 columns, 4 I-dim threads and 4 J-dim threads
    using block_type = haplo::Block<10, 12, 4, 4>;
    
    block_type block(input_1);
    
    block.print();
    
    BOOST_CHECK( block.num_unsplittable_blocks() == 3  );
    BOOST_CHECK( block.unsplittable_column(0)    == 1  );
    BOOST_CHECK( block.unsplittable_column(1)    == 3  );
    BOOST_CHECK( block.unsplittable_column(2)    == 4  );
    BOOST_CHECK( block.unsplittable_column(3)    == 11 );
}

BOOST_AUTO_TEST_SUITE_END()
