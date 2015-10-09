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
static constexpr const char* input_6 = "input_files/input_simulated_converted_1.txt";


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

    std::cout << "HERE";

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

    BOOST_CHECK( ublock.size() == 4 );
    BOOST_CHECK( ublock(0, 0)  == 0 );
    BOOST_CHECK( ublock(0, 1)  == 1 );
    BOOST_CHECK( ublock(1, 0)  == 1 );
    BOOST_CHECK( ublock(1, 1)  == 0 );
}

BOOST_AUTO_TEST_CASE( canRemoveMonotoneColumns )
{
    // Define a 10x12 block with a 2x2 CPU core grid
    using block_type = haplo::Block<10, 12, 2, 2>;
    
    // First create the block
    block_type block(input_1);
   
    // Use 4 cores for x and 4 for y 
    haplo::UnsplittableBlock<block_type, 4, 4, haplo::devices::cpu> ublock(block, 2);
    
    BOOST_CHECK( ublock.size() == 16 );
    BOOST_CHECK( ublock(0, 0)  == 2 );
    BOOST_CHECK( ublock(0, 1)  == 2 );
    BOOST_CHECK( ublock(0, 2)  == 1 );
    BOOST_CHECK( ublock(0, 3)  == 1 );
    BOOST_CHECK( ublock(1, 0)  == 2 );
    BOOST_CHECK( ublock(1, 1)  == 2 );
    BOOST_CHECK( ublock(1, 2)  == 0 );
    BOOST_CHECK( ublock(1, 3)  == 0 );
    BOOST_CHECK( ublock(2, 0)  == 0 );
    BOOST_CHECK( ublock(2, 1)  == 1 );
    BOOST_CHECK( ublock(2, 2)  == 2 );
    BOOST_CHECK( ublock(2, 3)  == 0 );
    BOOST_CHECK( ublock(3, 0)  == 2 );
    BOOST_CHECK( ublock(3, 1)  == 0 );
    BOOST_CHECK( ublock(3, 2)  == 0 );
    BOOST_CHECK( ublock(3, 3)  == 0 );
}


BOOST_AUTO_TEST_CASE( canInitializeTreeWhenDuplicateRowsInInput )
{
    // Define a 12x12 block with a 4x4 CPU core grid
    using block_type = haplo::Block<12, 12, 4, 4>;
    
    // First create the block (using duplicate row input)
    block_type block(input_4);
    
    haplo::UnsplittableBlock<block_type, 4, 4, haplo::devices::cpu> ublock(block, 0);
    
    // Get the tree which the block created
    auto tree = ublock.tree();
   
    // The node weights are the column multiplicities
    BOOST_CHECK( tree.node_weight(0) == 1 );
    BOOST_CHECK( tree.node_weight(1) == 1 );
    BOOST_CHECK( tree.node_weight(2) == 1 );
    BOOST_CHECK( tree.node_weight(3) == 1 );
   
    // The worst case values are the links between the nodes
    BOOST_CHECK( tree.node_worst_case(0) == 3 );
    BOOST_CHECK( tree.node_worst_case(1) == 3 );
    BOOST_CHECK( tree.node_worst_case(2) == 6 );
    BOOST_CHECK( tree.node_worst_case(3) == 2 );
}

BOOST_AUTO_TEST_CASE( canInitializeTreeWhenDuplicateColumnsInInput )
{
    // Define a 10x14 block with a 4x4 CPU core grid
    using block_type = haplo::Block<10, 14, 4, 4>;
    
    // First create the block (using duplicate column input)
    block_type block(input_5);

    haplo::UnsplittableBlock<block_type, 4, 4, haplo::devices::cpu> ublock(block, 0);
   
    auto tree = ublock.tree();

    // The node weights are the column multiplicities
    BOOST_CHECK( tree.node_weight(0) == 2 );
    BOOST_CHECK( tree.node_weight(1) == 1 );
    BOOST_CHECK( tree.node_weight(2) == 2 );
    BOOST_CHECK( tree.node_weight(3) == 1 );
    BOOST_CHECK( tree.node_weight(4) == 1 );
    
    // The worst case values are the links between the nodes
    // 1 and 3 are duplicates - we don't care about their values
    BOOST_CHECK( tree.node_worst_case(0) == 3 );
    BOOST_CHECK( tree.node_worst_case(2) == 3 );
    BOOST_CHECK( tree.node_worst_case(4) == 4 );
}

BOOST_AUTO_TEST_CASE( canFindTreeStartNodeForSearch ) 
{
    // Define a 10x14 block with a 4x4 CPU core grid
    using block_type = haplo::Block<10, 14, 4, 4>;
    
    // First create the block (using duplicate column input)
    block_type block(input_5);

    haplo::UnsplittableBlock<block_type, 4, 4, haplo::devices::cpu> ublock(block, 0);
   
    auto tree = ublock.tree();    

    BOOST_CHECK( tree.max_worst_case()  == 6 );
    BOOST_CHECK( tree.start_node()      == 2 );
}
BOOST_AUTO_TEST_SUITE_END()
