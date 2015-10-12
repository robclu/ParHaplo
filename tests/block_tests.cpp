// ----------------------------------------------------------------------------------------------------------
/// @file   block_tests.cpp
/// @brief  Test suite for parahaplo block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE BlockTests
#endif
#include <boost/test/unit_test.hpp>
#include <chrono>

//#include "../haplo/block.hpp"
#include "../haplo/subblock_cpu.hpp"

using namespace std::chrono;

static constexpr const char* input_1 = "input_files/input_zero.txt";
static constexpr const char* input_6 = "input_files/input_six.txt";
static constexpr const char* input_7 = "tests_files/output_7.txt";

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

BOOST_AUTO_TEST_CASE( canCreateSubBlocks )
{
    using block_type     = haplo::Block<28, 4, 4>;
    using sub_block_type = haplo::SubBlock<block_type, 4, 4, haplo::devices::cpu>; 
    
    block_type block(input_1);
    
    sub_block_type sub_1(block, 0);
    sub_block_type sub_2(block, 1);
    sub_block_type sub_3(block, 2);
    
    // Time the selection
    high_resolution_clock::time_point start = high_resolution_clock::now();     
    
    // Solve each of them     
    sub_1.find_haplotypes();
    sub_2.find_haplotypes();
    sub_3.find_haplotypes();
    
    high_resolution_clock::time_point solve_end = high_resolution_clock::now();
    duration<double> solve_time = duration_cast<duration<double>>(solve_end - start);
    
    std::cout << "SOLVE TIME : " << solve_time.count() << " seconds\n";
    
    // reconstruct the haplotypes
    block.merge_haplotype(sub_1); 
    block.merge_haplotype(sub_2); 
    block.merge_haplotype(sub_3); 
    
    high_resolution_clock::time_point reconstruct_end = high_resolution_clock::now();
    duration<double> reconstruction_time = duration_cast<duration<double>>(reconstruct_end - solve_end);
    
    std::cout << "MERGE TIME : " << reconstruction_time.count() << " seconds\n";
    
    block.print_haplotypes();
}

BOOST_AUTO_TEST_CASE( canSolve )
{
    using block_type     = haplo::Block<5603, 4, 4>;
    using sub_block_type = haplo::SubBlock<block_type, 64, 64, haplo::devices::cpu>; 
    
    block_type block(input_6);

    //block.print();
    std::cout << "\n\n";
    
    std::cout << block.num_subblocks() << "\n";
    sub_block_type sub_1(block, 1);
    sub_1.print();
    
    std::cout <<"\n\n";
    //sub_block_type sub_2(block, 1);
    
    // Time the selection
    high_resolution_clock::time_point start = high_resolution_clock::now();     
    
    // Solve each of them     
    sub_1.find_haplotypes();
  //  sub_2.find_haplotypes();
    
    high_resolution_clock::time_point solve_end = high_resolution_clock::now();
    duration<double> solve_time = duration_cast<duration<double>>(solve_end - start);
    
    std::cout << "SOLVE TIME : " << solve_time.count() << " seconds\n";
    
    // reconstruct the haplotypes
    block.merge_haplotype(sub_1); 
//    block.merge_haplotype(sub_2); 
    
    high_resolution_clock::time_point reconstruct_end = high_resolution_clock::now();
    duration<double> reconstruction_time = duration_cast<duration<double>>(reconstruct_end - solve_end);
    
    std::cout << "MERGE TIME : " << reconstruction_time.count() << " seconds\n";
    
    block.print_haplotypes();
    block.print_col_types();
}

BOOST_AUTO_TEST_SUITE_END()
