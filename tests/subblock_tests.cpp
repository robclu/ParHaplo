// ----------------------------------------------------------------------------------------------------------
/// @file   subblock_tests.cpp
/// @brief  Test suite for parahaplo sub-block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE SubBlockTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/subblock_cpu.hpp"

static constexpr const char* input_1 = "input_files/input_zero.txt";

BOOST_AUTO_TEST_SUITE( SubBlockSuite )
   
// NOTE: This test doesn't actually test anything, but the outputs of 
//       the test should show that an out of range exception was thrown
BOOST_AUTO_TEST_CASE( errorIsThrownForOutOfRangeSubBlock  )
{
    // Define a block for a with 4 CPU cores
    using block_type = haplo::Block<28, 2, 2>;
    
    // First create the block
    block_type block(input_1);
    
    // Create a sub-block from the block out of range -- just to illustrate out of range error
    haplo::SubBlock<block_type, 2, 2, haplo::devices::cpu> sub_block(block, 7);
}

BOOST_AUTO_TEST_CASE( canCreateSubBlockCorrectlyAndGetData1 )
{
    // Define a block for a with 4 CPU cores
    using block_type = haplo::Block<28, 2, 2>;
    
    // First create the block
    block_type block(input_1);

    haplo::SubBlock<block_type, 2, 2, haplo::devices::cpu> subblock(block, 0);

    BOOST_CHECK( subblock(0, 0) == 1 );
    BOOST_CHECK( subblock(0, 1) == 0 );
    BOOST_CHECK( subblock(0, 2) == 3 );
    BOOST_CHECK( subblock(1, 0) == 0 );
    BOOST_CHECK( subblock(1, 1) == 0 );
    BOOST_CHECK( subblock(1, 2) == 0 );
    BOOST_CHECK( subblock(2, 0) == 0 );
    BOOST_CHECK( subblock(2, 1) == 1 );
    BOOST_CHECK( subblock(2, 2) == 3 );
}

BOOST_AUTO_TEST_CASE( canRemoveMonotoneColumns )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_1);
    subblock_type   sub_block(input_1, 2);
               
    BOOST_CHECK( sub_block(0, 0)  == 3 );
    BOOST_CHECK( sub_block(0, 1)  == 3 );
    BOOST_CHECK( sub_block(0, 2)  == 1 );
    BOOST_CHECK( sub_block(0, 3)  == 1 );
    BOOST_CHECK( sub_block(1, 0)  == 3 );
    BOOST_CHECK( sub_block(1, 1)  == 3 );
    BOOST_CHECK( sub_block(1, 2)  == 0 );
    BOOST_CHECK( sub_block(1, 3)  == 0 );
    BOOST_CHECK( sub_block(2, 0)  == 0 );
    BOOST_CHECK( sub_block(2, 1)  == 1 );
    BOOST_CHECK( sub_block(2, 2)  == 2 );
    BOOST_CHECK( sub_block(2, 3)  == 0 );
    BOOST_CHECK( sub_block(3, 0)  == 3 );
    BOOST_CHECK( sub_block(3, 1)  == 0 );
    BOOST_CHECK( sub_block(3, 2)  == 0 );
    BOOST_CHECK( sub_block(3, 3)  == 0 );
}

/*

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
*/
BOOST_AUTO_TEST_SUITE_END()
