// ----------------------------------------------------------------------------------------------------------
/// @file   equality_checker_tests.cpp
/// @brief  Test suite for parahaplo block tests
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE EqualityChckerTests
#endif
#include <boost/test/unit_test.hpp>

#include "../equality_checker.hpp"
#include "../small_containers.hpp"

BOOST_AUTO_TEST_SUITE( EqualityCheckerSuite )
    
BOOST_AUTO_TEST_CASE( canCheckEquivalentRows )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<8, 1> row1;
   
    row1.set(1, 1);
    row1.set(5, 1);
    
    // Row1 looks like : 01000100, which can be 2 rows with 4 elements each
    //  : 0100
    //  : 0100
   
    haplo::EqualityChecker<haplo::check::ROWS> row_checker;
    
    BOOST_CHECK( row_checker(row1, 0, 1, 4, 4) == true );
}

BOOST_AUTO_TEST_CASE( canCheckNonEquivalentRows )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<8, 1> row1;
   
    row1.set(1, 1);
    row1.set(4, 1);
    
    // Row1 looks like : 01000100, which can be 2 rows with 4 elements each
    //  : 0100
    //  : 1000
   
    haplo::EqualityChecker<haplo::check::ROWS> row_checker;
    
    BOOST_CHECK( row_checker(row1, 0, 1, 4, 4) == false );
}

BOOST_AUTO_TEST_CASE( canCheckEquivalentColumns )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<16, 1> cols;
   
    cols.set(1 , 1);
    cols.set(5 , 1);
    cols.set(9 , 1);
    cols.set(13, 1);
    
    // Row1 looks like : 0100010001000100, which can be 2 rows with 4 elements each
    //  : 0100
    //  : 0100
    //  : 0100
    //  : 0100
   
    haplo::EqualityChecker<haplo::check::COLUMNS> col_checker;
    
    BOOST_CHECK( col_checker(cols, 0, 2, 4, 4) == true );
    BOOST_CHECK( col_checker(cols, 0, 3, 4, 4) == true );
    BOOST_CHECK( col_checker(cols, 2, 3, 4, 4) == true );
}

BOOST_AUTO_TEST_CASE( canCheckNonEquivalentColumns )
{
    // Define the binary container to use 1 bit per element (default setting)
    haplo::BinaryContainer<16, 1> cols;
   
    cols.set(1 , 1);
    cols.set(5 , 1);
    cols.set(9 , 1);
    cols.set(13, 1);
    
    // Row1 looks like : 0100010001000100, which can be 2 rows with 4 elements each
    //  : 0100
    //  : 0100
    //  : 0100
    //  : 0100
   
    haplo::EqualityChecker<haplo::check::COLUMNS> col_checker;
    
    BOOST_CHECK( col_checker(cols, 0, 1, 4, 4) == false );
}

BOOST_AUTO_TEST_SUITE_END()
