// ----------------------------------------------------------------------------------------------------------
/// @file   data_checker_tests.cpp
/// @brief  Test suite for parahaplo data converter tets
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE DataConverterTests
#endif
#include <boost/test/unit_test.hpp>

#include "../haplo/data_converter.hpp"

static constexpr const char* input_1 = "input_files/input_unfiltered_1.txt";

BOOST_AUTO_TEST_SUITE( DataConverterSuite )
    
BOOST_AUTO_TEST_CASE( canCreateDataConverter )
{
    haplo::DataConverter converter(input_1);
    
    converter.print();
    // rest of test ...
 }


BOOST_AUTO_TEST_SUITE_END()
