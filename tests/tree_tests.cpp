// ----------------------------------------------------------------------------------------------------------
/// @file   tree_tests.cpp
/// @brief  Test suite for parahaplo tree tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE TreeTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/block.hpp"
#include "../haplo/subblock_cpu.hpp"
#include "../haplo/tree_cpu.hpp"

static constexpr const char* input_zero = "input_files/input_zero.txt";

BOOST_AUTO_TEST_SUITE( TreeSuite )
    
BOOST_AUTO_TEST_CASE( canCreateATreeAndGetAndSetLinks )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_zero);
    subblock_type   sub_block(block, 2); 
    
    // Create a tree with 12 nodes
    haplo::Tree<subblock_type, haplo::devices::cpu> tree(sub_block, 12);

    tree.create_link(0, 1);
    tree.create_link(0, 2);
    tree.create_link(0, 3);
    
    tree.link(0, 1).homo_weight()  = 5;
    tree.link(0, 1).hetro_weight() = 4;
    
    tree.link(0, 3).hetro_weight() += 4;
    
    BOOST_CHECK( tree.link(0, 1).homo_weight()  == 5 );
    BOOST_CHECK( tree.link(1, 0).homo_weight()  == 5 );
    BOOST_CHECK( tree.link(0, 2).homo_weight()  == 0 );
    BOOST_CHECK( tree.link(0, 1).hetro_weight() == 4 );
    BOOST_CHECK( tree.link(0, 2).hetro_weight() == 0 );
    BOOST_CHECK( tree.link(0, 3).hetro_weight() == 4 );
    BOOST_CHECK( tree.link_max(0, 1)            == 5 );
    BOOST_CHECK( tree.link_min(0, 1)            == 4 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetNodeHaploTypePosition )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_zero);
    subblock_type   sub_block(block, 2); 
    
    // Create a tree with 12 nodes
    haplo::Tree<subblock_type, haplo::devices::cpu> tree(sub_block, 12);
    
    BOOST_CHECK( tree.node(4).position() == 4 );
    
    tree.node(4).position() = 12;
    
    BOOST_CHECK( tree.node(4).position() == 12 );
}


BOOST_AUTO_TEST_CASE( canGetAndSetNodeWeight )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_zero);
    subblock_type   sub_block(block, 2); 
    
    // Create a tree with 12 nodes
    haplo::Tree<subblock_type, haplo::devices::cpu> tree(sub_block, 12);
    
    BOOST_CHECK( tree.node(4).weight() == 1 );
    
    tree.node(4).weight() = 12;
    
    BOOST_CHECK( tree.node(4).weight() == 12 );
}

BOOST_AUTO_TEST_CASE( canGetAndSetNodeWorstCaseValue )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_zero);
    subblock_type   sub_block(block, 2); 
    
    // Create a tree with 12 nodes
    haplo::Tree<subblock_type, haplo::devices::cpu> tree(sub_block, 12); 
    
    tree.node(2).worst_case_value() = 12;
    
    BOOST_CHECK( tree.node(2).worst_case_value() == 12 );    
    BOOST_CHECK( tree.node(4).worst_case_value() == 0  );    
}

BOOST_AUTO_TEST_CASE( canResizeTree )
{
    using block_type    = haplo::Block<28, 4, 4>;
    using subblock_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>;
    
    block_type      block(input_zero);
    subblock_type   sub_block(block, 2); 
    
    // Create a tree with 12 nodes
    haplo::Tree<subblock_type, haplo::devices::cpu> tree(sub_block);  
    
    BOOST_CHECK( tree.size() == 0 );    
    
    tree.resize(12);
    
    BOOST_CHECK( tree.size() == 12 );    
}

BOOST_AUTO_TEST_SUITE_END()
