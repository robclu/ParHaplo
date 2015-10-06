// ----------------------------------------------------------------------------------------------------------
/// @file   data_checker_tests.cpp
/// @brief  Test suite for parahaplo data converter tets
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE DataConverterTests
#endif
#include <boost/test/unit_test.hpp>
#include <vector>

#include "../haplo/data_converter.hpp"

#define ZERO    0x00
#define ONE     0x01
#define TWO     0x02
#define THREE   0x03

static constexpr const char* input_1    = "input_files/input_simulated_2.txt";
static constexpr const char* input_2    = "input_files/input_simulated_1.txt";
static constexpr const char* input_3    = "input_files/input_simulated_3.txt";
static constexpr const char* input_4    = "input_files/input_dataset_1.txt";
static constexpr const char* output_1   = "output_files/output_simulated_2.txt";

BOOST_AUTO_TEST_SUITE( DataConverterSuite )
    
/*BOOST_AUTO_TEST_CASE( canCreateDataConverter )
{
    haplo::DataConverter converter(input_1);
    
    //converter.print();

    converter.write_data_to_file(output_1);
}

BOOST_AUTO_TEST_CASE( canConvertDataFromBinary )
{
    
    haplo::DataConverter converter(input_1);
    //converter.print();
    haplo::BinaryArray<30, 2> input;
    std::vector<char> output;
    for(size_t i = 0; i < 30; ++i){
        if(i % 2 == 0){
            input.set(i, ONE);
        }
        else {
            input.set(i, ZERO);
        }
    
    }
    
    output = converter.convert_data_from_binary(input);
    BOOST_CHECK(output.at(0) == 't');
    BOOST_CHECK(output.at(1) == 'a');
    BOOST_CHECK(output.at(2) == 'c');
    BOOST_CHECK(output.at(3) == 'c');
    BOOST_CHECK(output.at(4) == 'a');
    BOOST_CHECK(output.at(5) == 'a');
    BOOST_CHECK(output.at(6) == 'a');
    BOOST_CHECK(output.at(7) == 'g');
    BOOST_CHECK(output.at(8) == 'a');
    BOOST_CHECK(output.at(9) == 'a');
    
    //debugging
    //std::cout << std::endl << "Before:" << std::endl;
    //for (auto i = 0; i < 30; ++i) std::cout << static_cast<unsigned>(input.get(i));
    
    //std::cout << std::endl << "After:" << std::endl;
    //for (auto i = 0; i < 30; ++i) std::cout << output.at(i);
}*/

/*BOOST_AUTO_TEST_CASE( canMapBinaryToChar )
{
    haplo::DataConverter converter(input_1);
    char test_0;
    char test_1;
    char test_2;
    char test_3;
    test_0 = converter.convert_byte_to_char(ZERO);
    test_1 = converter.convert_byte_to_char(ONE);
    test_2 = converter.convert_byte_to_char(TWO);
    test_3 = converter.convert_byte_to_char(THREE);
    
    BOOST_CHECK(test_0 == 'a');
    BOOST_CHECK(test_1 == 'c');
    BOOST_CHECK(test_2 == 't');
    BOOST_CHECK(test_3 == 'g');
}

BOOST_AUTO_TEST_CASE( canMapCharToBinary )
{
    haplo::DataConverter converter(input_1);
    haplo::byte test_0;
    haplo::byte test_1;
    haplo::byte test_2;
    haplo::byte test_3;
    test_0 = converter.convert_char_to_byte('a');
    test_1 = converter.convert_char_to_byte('c');
    test_2 = converter.convert_char_to_byte('t');
    test_3 = converter.convert_char_to_byte('g');

    BOOST_CHECK(test_0 == ZERO);
    BOOST_CHECK(test_1 == ONE);
    BOOST_CHECK(test_2 == TWO);
    BOOST_CHECK(test_3 == THREE);
}*/

BOOST_AUTO_TEST_CASE( canConvertDataset )
{
    haplo::DataConverter converter(input_4);
    converter.printMap();
}


BOOST_AUTO_TEST_SUITE_END()
