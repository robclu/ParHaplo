// ----------------------------------------------------------------------------------------------------------
/// @file   unsplittable_block_tests.cpp
/// @brief  Test suite for parahaplo unsplittable block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/unsplittable_block_cpu.hpp"

static constexpr const char* input_1 = "input_files/input_unfiltered_1.txt";
static constexpr const char* input_2 = "input_files/input_unfiltered_2.txt";
static constexpr const char* input_3 = "input_files/input_singleton_rows.txt";
static constexpr const char* input_4 = "input_files/input_duplicate_rows.txt";
static constexpr const char* input_5 = "input_files/input_duplicate_cols.txt";

BOOST_AUTO_TEST_SUITE( UnsplittableBlockSuite )
   
// NOTE: This test doesn't actually test anything, but the outputs of 
//       the test should show that an out of range exception was thrown
BOOST_AUTO_TEST_CASE( errorIsThrownForOutOfRangeUnsplittableBlock  )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> ublock(block, 7);
}

BOOST_AUTO_TEST_CASE( canCreateUnsplittableBlockCorrectlyAndGetData1 )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> ublock(block, 0);

    ublock.print();
    

    BOOST_CHECK( ublock(0, 0) == 1 );
    BOOST_CHECK( ublock(0, 1) == 0 );
    BOOST_CHECK( ublock(0, 2) == 2 );
    BOOST_CHECK( ublock(1, 0) == 0 );
    BOOST_CHECK( ublock(1, 1) == 0 );
    BOOST_CHECK( ublock(1, 2) == 0 );
    BOOST_CHECK( ublock(2, 0) == 0 );
    BOOST_CHECK( ublock(2, 1) == 1 );
    BOOST_CHECK( ublock(2, 2) == 2 );
}


BOOST_AUTO_TEST_CASE( canCreateUnsplittableBlockCorrectlyAndGetData2 )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> ublock(block, 1);

    ublock.print();
    

    BOOST_CHECK( ublock(0, 0) == 0 );
    BOOST_CHECK( ublock(0, 1) == 1 );
    BOOST_CHECK( ublock(1, 0) == 1 );
    BOOST_CHECK( ublock(1, 1) == 0 );
}

BOOST_AUTO_TEST_CASE( canCreateUnsplittableBlockCorrectlyAndGetData3 )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 2, 2, haplo::devices::cpu> ublock(block, 2);

    ublock.print();
    
    BOOST_CHECK( ublock(0, 0) == 2 );
    BOOST_CHECK( ublock(0, 1) == 2 );
    BOOST_CHECK( ublock(0, 2) == 2 );
    BOOST_CHECK( ublock(0, 3) == 2 );
    BOOST_CHECK( ublock(0, 4) == 2 );
    BOOST_CHECK( ublock(0, 5) == 2 );
    BOOST_CHECK( ublock(0, 6) == 1 );
    BOOST_CHECK( ublock(0, 7) == 1 );
    BOOST_CHECK( ublock(1, 0) == 2 );
    BOOST_CHECK( ublock(1, 1) == 2 );
    BOOST_CHECK( ublock(1, 2) == 2 );
    BOOST_CHECK( ublock(1, 3) == 2 );
    BOOST_CHECK( ublock(1, 4) == 2 );
    BOOST_CHECK( ublock(1, 5) == 2 );
    BOOST_CHECK( ublock(1, 6) == 0 );
    BOOST_CHECK( ublock(1, 7) == 0 ); 
    BOOST_CHECK( ublock(2, 0) == 0 );
    BOOST_CHECK( ublock(2, 1) == 2 );
    BOOST_CHECK( ublock(2, 2) == 2 );
    BOOST_CHECK( ublock(2, 3) == 1 );
    BOOST_CHECK( ublock(2, 4) == 0 );
    BOOST_CHECK( ublock(2, 5) == 1 );
    BOOST_CHECK( ublock(2, 6) == 2 );
    BOOST_CHECK( ublock(2, 7) == 0 ); 
    BOOST_CHECK( ublock(3, 0) == 2 );
    BOOST_CHECK( ublock(3, 1) == 1 );
    BOOST_CHECK( ublock(3, 2) == 0 );
    BOOST_CHECK( ublock(3, 3) == 0 );
    BOOST_CHECK( ublock(3, 4) == 0 );
    BOOST_CHECK( ublock(3, 5) == 2 );
    BOOST_CHECK( ublock(3, 6) == 0 );
    BOOST_CHECK( ublock(3, 7) == 0 );  
}


BOOST_AUTO_TEST_CASE( canDetermineDuplicateRowsAndMultiplicites )
{
    // Define a 12x12 block with a 4x4 CPU core grid
    using block_type = haplo::Block<12, 12, 4, 4>;
    
    // First create the block (using duplicate row input)
    block_type block(input_4);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 4, 4, haplo::devices::cpu> ublock(block, 0);

    ublock.print();
}

BOOST_AUTO_TEST_CASE( canDetermineDuplicateColumnsAndMultiplicities )
{
    // Define a 10x14 block with a 4x4 CPU core grid
    using block_type = haplo::Block<10, 14, 4, 4>;
    
    // First create the block (using duplicate column input)
    block_type block(input_5);
    
    // Create an unsplittable block from the block out of range -- just to illustrate out of range error
    haplo::UnsplittableBlock<block_type, 4, 4, haplo::devices::cpu> ublock(block, 0);

    ublock.print();
}

BOOST_AUTO_TEST_SUITE_END()
