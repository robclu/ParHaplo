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

static constexpr const char* input_1    = "input_files/input_simulated_2.txt";
static constexpr const char* input_2    = "input_files/input_simulated_1.txt";
static constexpr const char* output_1   = "output_files/output_simulated_2.txt"; 
BOOST_AUTO_TEST_SUITE( DataConverterSuite )
    
BOOST_AUTO_TEST_CASE( canCreateDataConverter )
{
    haplo::DataConverter converter(input_1);
    
    converter.print();
    // rest of test ...

    converter.write_data_to_file(output_1);
 }


BOOST_AUTO_TEST_SUITE_END()
