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
#include <fstream>

#include "../haplo/data_converter.hpp"

#define ZERO    0x00
#define ONE     0x01
#define TWO     0x02
#define THREE   0x03

static constexpr const char* input_1    = "input_files/input_simulated_1.txt";
static constexpr const char* input_2    = "input_files/input_simulated_2.txt";
static constexpr const char* input_3    = "input_files/input_simulated_3.txt";
static constexpr const char* input_4    = "input_files/input_dataset_1.txt";
static constexpr const char* input_5    = "input_files/input_dataset_2.txt";
static constexpr const char* output_1   = "output_files/output_simulated_1.txt";
static constexpr const char* output_2   = "output_files/output_simulated_2.txt";
static constexpr const char* output_3   = "output_files/output_simulated_3.txt";
static constexpr const char* output_4   = "output_files/output_dataset_1.txt";

BOOST_AUTO_TEST_SUITE( DataConverterSuite )

// Do tests on three different inputs for the simulated data
BOOST_AUTO_TEST_CASE( canCreateDataConverter )
{
    
    std::vector<const char*> inputs = {input_1, input_2, input_3};
    std::vector<const char*> outputs = {output_1, output_2, output_3};

    for(size_t i = 0; i < inputs.size(); ++i){
        haplo::DataConverter converter(inputs.at(i), outputs.at(i));
        
        //if want to see output in command window
        //converter.print_simulated();
        
        size_t number_of_lines_input = 0;
        size_t number_of_lines_output = 0;
        std::string line;
        std::ifstream input;
        input.open(input_1);
        
        while (std::getline(input, line))
            ++number_of_lines_input;
        
        std::ifstream output(output_1);
        
        while (std::getline(output, line))
            ++number_of_lines_output;
        
        BOOST_CHECK(number_of_lines_input == number_of_lines_output);
        
    }
    
}

/*// Tests conversion from binary haplotype to actg haplotype
BOOST_AUTO_TEST_CASE( canConvertDataFromBinary )
{
    
    haplo::DataConverter converter(input_1, output_1);
    //converter.print_simulated();
    haplo::BinaryArray<10, 2> input;
    std::vector<char> output;
    // create test haplotype
    for(size_t i = 0; i < 10; ++i){
        if(i % 2 == 0){
            input.set(i, ONE);
        }
        else {
            input.set(i, ZERO);
        }
    
    }
    
    output = converter.convert_data_from_binary(input);
    BOOST_CHECK(output.at(0) == 'g');
    BOOST_CHECK(output.at(1) == 'c');
    BOOST_CHECK(output.at(2) == 'g');
    BOOST_CHECK(output.at(3) == 'c');
    BOOST_CHECK(output.at(4) == 'g');
    BOOST_CHECK(output.at(5) == 'c');
    BOOST_CHECK(output.at(6) == 'c');
    BOOST_CHECK(output.at(7) == 'g');
    BOOST_CHECK(output.at(8) == 't');
    BOOST_CHECK(output.at(9) == 'c');
    
    //debugging
    //std::cout << std::endl << "Before:" << std::endl;
    //for (auto i = 0; i < 10; ++i) std::cout << static_cast<unsigned>(input.get(i));
    
    //std::cout << std::endl << "After:" << std::endl;
    //for (auto i = 0; i < 10; ++i) std::cout << output.at(i);
}

// Tests conversion from atcg haplotype to binary haplotype
BOOST_AUTO_TEST_CASE( canConvertDataToBinary )
{
    haplo::DataConverter converter(input_2, output_2);
    //converter.print_simulated();
    std::vector<char> input = {'a', 't', 't', 't', 'g', 'a', 'c', 'g', 'a', 't'};
    std::vector<size_t> output;
    
    output = converter.convert_data_to_binary(input);
    BOOST_CHECK(output.at(0) == 0);
    BOOST_CHECK(output.at(1) == 1);
    BOOST_CHECK(output.at(2) == 0);
    BOOST_CHECK(output.at(3) == 0);
    BOOST_CHECK(output.at(4) == 0);
    BOOST_CHECK(output.at(5) == 1);
    BOOST_CHECK(output.at(6) == 0);
    BOOST_CHECK(output.at(7) == 0);
    BOOST_CHECK(output.at(8) == 1);
    BOOST_CHECK(output.at(9) == 0);
    
    
}

BOOST_AUTO_TEST_CASE( canMapBinaryToChar )
{
    haplo::DataConverter converter(input_1, output_1);
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
    haplo::DataConverter converter(input_1, output_1);
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
}

BOOST_AUTO_TEST_CASE( canConvertDataset )
{
    haplo::DataConverter converter(input_4, input_5, output_4);
    //converter.print_dataset();
}*/


BOOST_AUTO_TEST_SUITE_END()
