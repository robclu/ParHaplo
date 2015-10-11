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
static constexpr const char* input_4    = "input_files/input_simulated_4.txt";
static constexpr const char* input_5    = "input_files/input_simulated_5.txt";
static constexpr const char* input_6    = "input_files/input_dataset_1.txt";
static constexpr const char* input_7    = "input_files/input_dataset_2.txt";
static constexpr const char* output_1   = "output_files/output_simulated_1.txt";
static constexpr const char* output_2   = "output_files/output_simulated_2.txt";
static constexpr const char* output_3   = "output_files/output_simulated_3.txt";
static constexpr const char* output_4   = "output_files/output_simulated_4.txt";
static constexpr const char* output_5   = "output_files/output_simulated_5.txt";
static constexpr const char* output_6   = "output_files/output_dataset_1.txt";
static constexpr const char* answer_letters_4    = "input_files/input_simulated_4_answer.txt";
static constexpr const char* answer_letters_5    = "input_files/input_simulated_5_answer.txt";
static constexpr const char* answer_binary_4    = "output_files/output_simulated_4_answer.txt";
static constexpr const char* answer_binary_5    = "output_files/ouput_simulated_5_answer.txt";

BOOST_AUTO_TEST_SUITE( DataConverterSuite )

// Do tests on three different inputs for the simulated data
BOOST_AUTO_TEST_CASE( canCreateDataConverter )
{
    
    std::vector<const char*> inputs = {input_1, input_2, input_3, input_4, input_5};
    std::vector<const char*> outputs = {output_1, output_2, output_3, output_4, output_5};
    std::vector<const char*> answers_letters = {answer_letters_4, answer_letters_5};
    std::vector<const char*> answers_binary = {answer_binary_4, answer_binary_5};


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
        
        BOOST_CHECK(number_of_lines_input + 1 == number_of_lines_output);
        
        std::vector<char> input_string;
        std::vector<size_t> output_value;
        
        if(i == 3 || i == 4){
            
            std::ifstream infile(answers_letters.at(i-3));
            std::ofstream outfile(answers_binary.at(i-3));
            std::string line;
            while(infile >> line){
                for(int len = 0; len < line.length(); ++len){
                    input_string.push_back(line[len]);
                }
                
                output_value = converter.convert_data_to_binary(input_string);
                
                for(int vec = 0; vec < output_value.size(); ++vec){
                    outfile << output_value.at(vec);
                }
                outfile << std::endl;
                    
                input_string.clear();
                output_value.clear();
                
            
            }
        
          }
    

        
    }
    
}

// Tests conversion from binary haplotype to actg haplotype
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
    haplo::DataConverter converter(input_6, input_7, output_6);
    //converter.print_dataset();
}

// Check all cases of cigar value are processed properly
BOOST_AUTO_TEST_CASE( canProcessCigarValue )
{
    haplo::DataConverter converter(input_6, input_7, output_6);
    // example 1 : mid padding, end in S/H, start and mid M
    size_t start_position = 2;
    size_t end_position = 10;
    std::string cigar_value = "2M2P3M4S";
    std::string sequence = "ACTGACTGA";
    converter.process_cigar_value(start_position, end_position, cigar_value, sequence);
    BOOST_CHECK(sequence == "AC22TGA");
    BOOST_CHECK(start_position == 2);
    BOOST_CHECK(end_position == 8);
    BOOST_CHECK(sequence.length() == (end_position - start_position + 1));
    
    // example 2 : mid S/H, beginning padding, mid M
    start_position = 4;
    end_position = 9;
    cigar_value = "2N4M2M3H";
    sequence = "ACTGAC";
    converter.process_cigar_value(start_position, end_position, cigar_value, sequence);
    BOOST_CHECK(sequence == "22ACTGAC");
    BOOST_CHECK(start_position == 4);
    BOOST_CHECK(end_position == 11);
    BOOST_CHECK(sequence.length() == (end_position - start_position + 1));
    
    // example 3 : beginning S/H, mid padding, mid and end M
    start_position = 5;
    end_position = 11;
    cigar_value = "2S3M2N2M";
    sequence = "ACTGACT";
    converter.process_cigar_value(start_position, end_position, cigar_value, sequence);
    BOOST_CHECK(sequence == "TGA22CT");
    BOOST_CHECK(start_position == 5);
    BOOST_CHECK(end_position == 11);
    BOOST_CHECK(sequence.length() == (end_position - start_position + 1));
    
    // example 4 : check I has same action as M
    start_position = 5;
    end_position = 11;
    cigar_value = "2M3I2N2I";
    sequence = "ACTGACT";
    converter.process_cigar_value(start_position, end_position, cigar_value, sequence);
    BOOST_CHECK(sequence == "ACTGA22CT");
    BOOST_CHECK(start_position == 5);
    BOOST_CHECK(end_position == 13);
    BOOST_CHECK(sequence.length() == (end_position - start_position + 1));
    

}


BOOST_AUTO_TEST_SUITE_END()
