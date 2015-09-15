// ----------------------------------------------------------------------------------------------------------
/// @file   parahaplo_cpp_tests.cpp
/// @brief  Test suites for the c++ implementation of the parahaplo library using Bost.Unit.
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE       ParahaploCppTests
#include <boost/test/unit_test.hpp>

#include "cpp/unsplittable_block.hpp"

#include <iostream>

// ------------------------------------------ DATA TESTS ----------------------------------------------------

BOOST_AUTO_TEST_SUITE( DataTestSuite )

BOOST_AUTO_TEST_CASE( canCreateDataTypeAndGetValueForCorrectInput )
{
    haplo::Data data_zero('0');
    haplo::Data data_one( '1');
    haplo::Data data_dash('-');
    
    BOOST_CHECK( data_zero.value() == 0 );
    BOOST_CHECK( data_one.value()  == 1 );
    BOOST_CHECK( data_dash.value() == 2 );
}

BOOST_AUTO_TEST_CASE( canCreateDataTypeAndGetValueForIncorrectInput )
{
    haplo::Data data_wrong_1('b');
    haplo::Data data_wrong_2('z');
    
    BOOST_CHECK( data_wrong_1.value() == 3 );
    BOOST_CHECK( data_wrong_2.value() == 3 );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( BlockOperations )

BOOST_AUTO_TEST_CASE( canDetermineIfTwoDataRowsAreEqual )
{
    // Create a vectors (rows) of data 
    std::vector<haplo::Data> rows{ '1', '0', '1', '-', '-',
                                   '1', '0', '1', '-', '-' };

    size_t stride       = 5,                // 5 elements between the start of each row
           row_length   = 5;                // The number of elements in each row
    
    // Create an instance of the checking function for rows
    haplo::EqualityChecker<haplo::check::ROWS> checker;
    
    bool rows_equal = checker(&rows[0], 0, 1, row_length, stride);
    
    BOOST_CHECK( rows_equal == true );
}

BOOST_AUTO_TEST_CASE( canDetermineIfTwoDataColumnsAreEqual )
{
    // Create data 4 rows and 2 columns
    std::vector<haplo::Data> data{ '1', '1', 
                                   '1', '1', 
                                   '-', '-',
                                   '0', '0'};
    
    size_t stride     = 2,          // Number of columns in data 
           col_length = 4;          // Number of elements in the columns
   
    haplo::EqualityChecker<haplo::check::COLUMNS> checker;

    bool columns_equal = checker(&data[0], 0, 1, col_length, stride);
    
    BOOST_CHECK( columns_equal == true );
}

BOOST_AUTO_TEST_SUITE_END()

// ------------------------------------ BINARY CONTAINER TESTS ----------------------------------------------

BOOST_AUTO_TEST_SUITE( BinaryContainerTestSuite )
    
BOOST_AUTO_TEST_CASE( canCreateBinaryContainerAndModifyValues )
{
    haplo::BinaryContainer container;
    
    BOOST_CHECK( container.get(0) == 0 );
    
    // Modify the value of the second bit
    container.set(1, 1);
    
    BOOST_CHECK( container.get(1) == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
   
// ---------------------------------------- BNB VARIABLE TESTS ----------------------------------------------

BOOST_AUTO_TEST_SUITE( BnbVariableTestSuite )
    
BOOST_AUTO_TEST_CASE( canCreateAndModify1DVariableArray )
{
    // 10 elements 
    haplo::BnbVariable<1> x(10);
    
    x.set(4, 0);
    
    BOOST_CHECK( x.get(4) == 0 );
    
    // Modify the value of the second bit
    x.set(4, 1);
    
    BOOST_CHECK( x.get(4) == 1 );
}

BOOST_AUTO_TEST_CASE( canCreateAndModify2DVariableArray )
{
    // 10 elements by 10 elements
    haplo::BnbVariable<2> t(10, 10);
    
    t.set(2, 1, 0);
    
    BOOST_CHECK( t.get(2, 1) == 0 );
    
    // Modify the value of the second bit
    t.set(2, 1, 1);
    
    BOOST_CHECK( t.get(2, 1) == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
    
// ---------------------------------------- BLOCK TESTS -----------------------------------------------------

BOOST_AUTO_TEST_SUITE( BlockTestSuite )
    
BOOST_AUTO_TEST_CASE( canCreateABlockOfAnyTypeAnSize )
{
    haplo::Block<3, 4>    block_34;
    haplo::Block<9, 9>    block_99;
    
    BOOST_CHECK( block_34.size() == 12 );
    BOOST_CHECK( block_99.size() == 81 );
}

BOOST_AUTO_TEST_CASE( canGetNumberOfRows )
{
    haplo::Block<3, 4> block_34;
    
    BOOST_CHECK( block_34.rows() == 3 );
}

BOOST_AUTO_TEST_CASE( canGetNumberOfColumns )
{
    haplo::Block<3, 4> block_34;
    
    BOOST_CHECK( block_34.columns() == 4 );
}

BOOST_AUTO_TEST_CASE( canAssignAndGetBlockData )
{
    haplo::Block<2, 2>  block_22;
    haplo::Data         data[4] = {'0', '-', '1', '1'};
    
    // Assign data - static assert for dimensions mismatch
    block_22.assign_data(data);
    
    // Get a reference to the data
    const std::vector<haplo::Data>& block_22_data = block_22.get_data();
    
    BOOST_CHECK( block_22_data[0] == 0 );
    BOOST_CHECK( block_22_data[1] == 2 );
    BOOST_CHECK( block_22_data[2] == 1 );
    BOOST_CHECK( block_22_data[3] == 1 );
}

BOOST_AUTO_TEST_CASE( canCreateBlockFromInputFile ) 
{ 
    const std::string input_file = "data/test_input_file.txt";
    haplo::Block<10, 7> block_10_7(input_file, 8);         // Use 8 threads
    
    BOOST_CHECK( block_10_7(0, 0) == 0 );
    BOOST_CHECK( block_10_7(0, 1) == 2 );
    BOOST_CHECK( block_10_7(1, 0) == 1 );
    BOOST_CHECK( block_10_7(1, 1) == 0 );
    BOOST_CHECK( block_10_7(2, 2) == 0 );
    BOOST_CHECK( block_10_7(5, 3) == 0 );
    BOOST_CHECK( block_10_7(8, 5) == 2 );
}

BOOST_AUTO_TEST_CASE( canGetReadInfoCorrectly )
{
    const std::string input_file = "data/test_input_file.txt";
    haplo::Block<10, 7> block_10_7(input_file, 8);         // Use 8 threads

    // Constructor determines read info 
    // when the input file is given 
    
    BOOST_CHECK( block_10_7.read_info()[0].start()    == 0 );
    BOOST_CHECK( block_10_7.read_info()[0].end()      == 0 );
    BOOST_CHECK( block_10_7.read_info()[0].length()   == 1 );
    BOOST_CHECK( block_10_7.read_info()[4].start()    == 2 );
    BOOST_CHECK( block_10_7.read_info()[4].end()      == 3 );
    BOOST_CHECK( block_10_7.read_info()[4].length()   == 2 ); 
    BOOST_CHECK( block_10_7.read_info()[8].start()    == 3 );
    BOOST_CHECK( block_10_7.read_info()[8].end()      == 6 );
    BOOST_CHECK( block_10_7.read_info()[8].length()   == 4 );
}

BOOST_AUTO_TEST_CASE( canGetUnsplittableSubBlockInfo )
{
    const std::string input_file = "data/test_input_file.txt";
    haplo::Block<10, 7> block_10_7(input_file, 8);         // Use 8 threads

    // Constructor determines subblock info 
    // when the input file is given  

    BOOST_CHECK( block_10_7.subblock_info()[0].start() == 0 );    
    BOOST_CHECK( block_10_7.subblock_info()[0].end()   == 2 );    
    BOOST_CHECK( block_10_7.subblock_info()[1].start() == 2 );    
    BOOST_CHECK( block_10_7.subblock_info()[1].end()   == 3 );    
    BOOST_CHECK( block_10_7.subblock_info()[2].start() == 3 );    
    BOOST_CHECK( block_10_7.subblock_info()[2].end()   == 6 );    
}

BOOST_AUTO_TEST_SUITE_END()
    
// ----------------------------------- UNSPLITTABLE BLOCK TESTS ---------------------------------------------

BOOST_AUTO_TEST_SUITE( UnsplittableBlockSuite )

BOOST_AUTO_TEST_CASE( canCreateAnUnsplittableBlockAndDetermineRowMultiplicities )
{
    using block_type = haplo::Block<13, 7>;
    
    const std::string input_file = "data/test_input_file_duplicate_rows.txt";
    block_type block(input_file, 8);                        // Use 8 threads 
    
    // Create an Unsplittable block 
    haplo::UnsplittableBlock<block_type> usb(block);
    
    BOOST_CHECK( usb.row_multiplicity(0) == 2 );
    BOOST_CHECK( usb.row_multiplicity(2) == 2 );
    BOOST_CHECK( usb.row_multiplicity(4) == 1 );
}

BOOST_AUTO_TEST_CASE( canCreateAnUnsplittableBlockAndDetermineColumnMultiplicities )
{
    using block_type = haplo::Block<13, 9>;
    
    const std::string input_file = "data/test_input_file_duplicate_cols.txt";
    block_type block(input_file, 8);                        // Use 8 threads 
    
    // Create an Unsplittable block 
    haplo::UnsplittableBlock<block_type> usb(block);
    
    BOOST_CHECK( usb.col_multiplicity(0) == 2 );
    BOOST_CHECK( usb.col_multiplicity(2) == 2 );
    BOOST_CHECK( usb.col_multiplicity(4) == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
    
