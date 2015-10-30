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

static constexpr const char* input_1    = "new_inputs/geraci_0.1/100_3_0.1_0.4/input_1.txt";
static constexpr const char* input_2    = "new_inputs/geraci_0.1/100_3_0.1_0.4/input_2.txt";
static constexpr const char* input_3    = "new_inputs/geraci_0.1/100_3_0.1_0.4/input_3.txt";

static constexpr const char* input_4    = "new_inputs/geraci_0.1/100_5_0.1_0.4/input_1.txt";
static constexpr const char* input_5    = "new_inputs/geraci_0.1/100_5_0.1_0.4/input_2.txt";
static constexpr const char* input_6    = "new_inputs/geraci_0.1/100_5_0.1_0.4/input_3.txt";

static constexpr const char* input_7    = "new_inputs/geraci_0.1/100_8_0.1_0.4/input_1.txt";
static constexpr const char* input_8    = "new_inputs/geraci_0.1/100_8_0.1_0.4/input_2.txt";
static constexpr const char* input_9    = "new_inputs/geraci_0.1/100_8_0.1_0.4/input_3.txt";

static constexpr const char* input_10    = "new_inputs/geraci_0.1/100_10_0.1_0.4/input_1.txt";
static constexpr const char* input_11    = "new_inputs/geraci_0.1/100_10_0.1_0.4/input_2.txt";
static constexpr const char* input_12    = "new_inputs/geraci_0.1/100_10_0.1_0.4/input_3.txt";

static constexpr const char* input_13    = "new_inputs/geraci_0.1/350_3_0.1_0.4/input_1.txt";
static constexpr const char* input_14    = "new_inputs/geraci_0.1/350_3_0.1_0.4/input_2.txt";
static constexpr const char* input_15    = "new_inputs/geraci_0.1/350_3_0.1_0.4/input_3.txt";

static constexpr const char* input_16    = "new_inputs/geraci_0.1/350_5_0.1_0.4/input_1.txt";
static constexpr const char* input_17    = "new_inputs/geraci_0.1/350_5_0.1_0.4/input_2.txt";
static constexpr const char* input_18    = "new_inputs/geraci_0.1/350_5_0.1_0.4/input_3.txt";

static constexpr const char* input_19    = "new_inputs/geraci_0.1/350_8_0.1_0.4/input_1.txt";
static constexpr const char* input_20    = "new_inputs/geraci_0.1/350_8_0.1_0.4/input_2.txt";
static constexpr const char* input_21    = "new_inputs/geraci_0.1/350_8_0.1_0.4/input_3.txt";

static constexpr const char* input_22    = "new_inputs/geraci_0.1/350_10_0.1_0.4/input_1.txt";
static constexpr const char* input_23    = "new_inputs/geraci_0.1/350_10_0.1_0.4/input_2.txt";
static constexpr const char* input_24    = "new_inputs/geraci_0.1/350_10_0.1_0.4/input_3.txt";

static constexpr const char* input_25    = "new_inputs/geraci_0.1/700_3_0.1_0.4/input_1.txt";
static constexpr const char* input_26    = "new_inputs/geraci_0.1/700_3_0.1_0.4/input_2.txt";
static constexpr const char* input_27    = "new_inputs/geraci_0.1/700_3_0.1_0.4/input_3.txt";

static constexpr const char* input_28    = "new_inputs/geraci_0.1/700_5_0.1_0.4/input_1.txt";
static constexpr const char* input_29    = "new_inputs/geraci_0.1/700_5_0.1_0.4/input_2.txt";
static constexpr const char* input_30    = "new_inputs/geraci_0.1/700_5_0.1_0.4/input_3.txt";

static constexpr const char* input_31    = "new_inputs/geraci_0.1/700_8_0.1_0.4/input_1.txt";
static constexpr const char* input_32    = "new_inputs/geraci_0.1/700_8_0.1_0.4/input_2.txt";
static constexpr const char* input_33    = "new_inputs/geraci_0.1/700_8_0.1_0.4/input_3.txt";

static constexpr const char* input_34    = "new_inputs/geraci_0.1/700_10_0.1_0.4/input_1.txt";
static constexpr const char* input_35    = "new_inputs/geraci_0.1/700_10_0.1_0.4/input_2.txt";
static constexpr const char* input_36    = "new_inputs/geraci_0.1/700_10_0.1_0.4/input_3.txt";

static constexpr const char* output_1    = "new_outputs/geraci_0.1/100_3_0.1_0.4/output_1";
static constexpr const char* output_2    = "new_outputs/geraci_0.1/100_3_0.1_0.4/output_2";
static constexpr const char* output_3    = "new_outputs/geraci_0.1/100_3_0.1_0.4/output_3";

static constexpr const char* output_4    = "new_outputs/geraci_0.1/100_5_0.1_0.4/output_1";
static constexpr const char* output_5    = "new_outputs/geraci_0.1/100_5_0.1_0.4/output_2";
static constexpr const char* output_6    = "new_outputs/geraci_0.1/100_5_0.1_0.4/output_3";

static constexpr const char* output_7    = "new_outputs/geraci_0.1/100_8_0.1_0.4/output_1";
static constexpr const char* output_8    = "new_outputs/geraci_0.1/100_8_0.1_0.4/output_2";
static constexpr const char* output_9    = "new_outputs/geraci_0.1/100_8_0.1_0.4/output_3";

static constexpr const char* output_10    = "new_outputs/geraci_0.1/100_10_0.1_0.4/output_1";
static constexpr const char* output_11    = "new_outputs/geraci_0.1/100_10_0.1_0.4/output_2";
static constexpr const char* output_12    = "new_outputs/geraci_0.1/100_10_0.1_0.4/output_3";

static constexpr const char* output_13    = "new_outputs/geraci_0.1/350_3_0.1_0.4/output_1";
static constexpr const char* output_14    = "new_outputs/geraci_0.1/350_3_0.1_0.4/output_2";
static constexpr const char* output_15    = "new_outputs/geraci_0.1/350_3_0.1_0.4/output_3";

static constexpr const char* output_16    = "new_outputs/geraci_0.1/350_5_0.1_0.4/output_1";
static constexpr const char* output_17    = "new_outputs/geraci_0.1/350_5_0.1_0.4/output_2";
static constexpr const char* output_18    = "new_outputs/geraci_0.1/350_5_0.1_0.4/output_3";

static constexpr const char* output_19    = "new_outputs/geraci_0.1/350_8_0.1_0.4/output_1";
static constexpr const char* output_20    = "new_outputs/geraci_0.1/350_8_0.1_0.4/output_2";
static constexpr const char* output_21    = "new_outputs/geraci_0.1/350_8_0.1_0.4/output_3";

static constexpr const char* output_22    = "new_outputs/geraci_0.1/350_10_0.1_0.4/output_1";
static constexpr const char* output_23    = "new_outputs/geraci_0.1/350_10_0.1_0.4/output_2";
static constexpr const char* output_24    = "new_outputs/geraci_0.1/350_10_0.1_0.4/output_3";

static constexpr const char* output_25    = "new_outputs/geraci_0.1/700_3_0.1_0.4/output_1";
static constexpr const char* output_26    = "new_outputs/geraci_0.1/700_3_0.1_0.4/output_2";
static constexpr const char* output_27    = "new_outputs/geraci_0.1/700_3_0.1_0.4/output_3";

static constexpr const char* output_28    = "new_outputs/geraci_0.1/700_5_0.1_0.4/output_1";
static constexpr const char* output_29    = "new_outputs/geraci_0.1/700_5_0.1_0.4/output_2";
static constexpr const char* output_30    = "new_outputs/geraci_0.1/700_5_0.1_0.4/output_3";

static constexpr const char* output_31    = "new_outputs/geraci_0.1/700_8_0.1_0.4/output_1";
static constexpr const char* output_32    = "new_outputs/geraci_0.1/700_8_0.1_0.4/output_2";
static constexpr const char* output_33    = "new_outputs/geraci_0.1/700_8_0.1_0.4/output_3";

static constexpr const char* output_34    = "new_outputs/geraci_0.1/700_10_0.1_0.4/output_1";
static constexpr const char* output_35    = "new_outputs/geraci_0.1/700_10_0.1_0.4/output_2";
static constexpr const char* output_36    = "new_outputs/geraci_0.1/700_10_0.1_0.4/output_3";

static constexpr const char* answer_letters_1    = "new_inputs/geraci_0.1/100_3_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_2    = "new_inputs/geraci_0.1/100_3_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_3    = "new_inputs/geraci_0.1/100_3_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_4    = "new_inputs/geraci_0.1/100_5_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_5    = "new_inputs/geraci_0.1/100_5_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_6    = "new_inputs/geraci_0.1/100_5_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_7    = "new_inputs/geraci_0.1/100_8_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_8    = "new_inputs/geraci_0.1/100_8_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_9    = "new_inputs/geraci_0.1/100_8_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_10    = "new_inputs/geraci_0.1/100_10_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_11    = "new_inputs/geraci_0.1/100_10_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_12    = "new_inputs/geraci_0.1/100_10_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_13    = "new_inputs/geraci_0.1/350_3_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_14    = "new_inputs/geraci_0.1/350_3_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_15    = "new_inputs/geraci_0.1/350_3_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_16    = "new_inputs/geraci_0.1/350_5_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_17    = "new_inputs/geraci_0.1/350_5_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_18    = "new_inputs/geraci_0.1/350_5_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_19    = "new_inputs/geraci_0.1/350_8_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_20    = "new_inputs/geraci_0.1/350_8_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_21    = "new_inputs/geraci_0.1/350_8_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_22    = "new_inputs/geraci_0.1/350_10_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_23    = "new_inputs/geraci_0.1/350_10_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_24    = "new_inputs/geraci_0.1/350_10_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_25    = "new_inputs/geraci_0.1/700_3_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_26    = "new_inputs/geraci_0.1/700_3_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_27    = "new_inputs/geraci_0.1/700_3_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_28    = "new_inputs/geraci_0.1/700_5_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_29    = "new_inputs/geraci_0.1/700_5_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_30    = "new_inputs/geraci_0.1/700_5_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_31    = "new_inputs/geraci_0.1/700_8_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_32    = "new_inputs/geraci_0.1/700_8_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_33    = "new_inputs/geraci_0.1/700_8_0.1_0.4/input_3_soln.txt";

static constexpr const char* answer_letters_34    = "new_inputs/geraci_0.1/700_10_0.1_0.4/input_1_soln.txt";
static constexpr const char* answer_letters_35    = "new_inputs/geraci_0.1/700_10_0.1_0.4/input_2_soln.txt";
static constexpr const char* answer_letters_36    = "new_inputs/geraci_0.1/700_10_0.1_0.4/input_3_soln.txt";


static constexpr const char* answer_binary_1    = "new_outputs/geraci_0.1/100_3_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_2    = "new_outputs/geraci_0.1/100_3_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_3    = "new_outputs/geraci_0.1/100_3_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_4    = "new_outputs/geraci_0.1/100_5_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_5    = "new_outputs/geraci_0.1/100_5_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_6    = "new_outputs/geraci_0.1/100_5_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_7    = "new_outputs/geraci_0.1/100_8_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_8    = "new_outputs/geraci_0.1/100_8_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_9    = "new_outputs/geraci_0.1/100_8_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_10    = "new_outputs/geraci_0.1/100_10_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_11    = "new_outputs/geraci_0.1/100_10_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_12    = "new_outputs/geraci_0.1/100_10_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_13    = "new_outputs/geraci_0.1/350_3_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_14    = "new_outputs/geraci_0.1/350_3_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_15    = "new_outputs/geraci_0.1/350_3_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_16    = "new_outputs/geraci_0.1/350_5_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_17    = "new_outputs/geraci_0.1/350_5_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_18    = "new_outputs/geraci_0.1/350_5_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_19    = "new_outputs/geraci_0.1/350_8_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_20    = "new_outputs/geraci_0.1/350_8_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_21    = "new_outputs/geraci_0.1/350_8_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_22    = "new_outputs/geraci_0.1/350_10_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_23    = "new_outputs/geraci_0.1/350_10_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_24    = "new_outputs/geraci_0.1/350_10_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_25    = "new_outputs/geraci_0.1/700_3_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_26    = "new_outputs/geraci_0.1/700_3_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_27    = "new_outputs/geraci_0.1/700_3_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_28    = "new_outputs/geraci_0.1/700_5_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_29    = "new_outputs/geraci_0.1/700_5_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_30    = "new_outputs/geraci_0.1/700_5_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_31    = "new_outputs/geraci_0.1/700_8_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_32    = "new_outputs/geraci_0.1/700_8_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_33    = "new_outputs/geraci_0.1/700_8_0.1_0.4/output_3_soln.txt";

static constexpr const char* answer_binary_34    = "new_outputs/geraci_0.1/700_10_0.1_0.4/output_1_soln.txt";
static constexpr const char* answer_binary_35    = "new_outputs/geraci_0.1/700_10_0.1_0.4/output_2_soln.txt";
static constexpr const char* answer_binary_36    = "new_outputs/geraci_0.1/700_10_0.1_0.4/output_3_soln.txt";

// These test files correspond to bigger bam files and are not included due to their size
/*static constexpr const char* input_7    = "input_files/input_dataset_1.txt";
static constexpr const char* input_8    = "input_files/input_dataset_2.txt";
static constexpr const char* output_7   = "output_files/chr20_output.txt";*/




BOOST_AUTO_TEST_SUITE( DataConverterSuite )

// Do tests on three different inputs for the simulated data
BOOST_AUTO_TEST_CASE( canCreateDataConverter )
{
    
    std::vector<const char*> inputs = {input_1, input_2, input_3, input_4, input_5, input_6, input_7, input_8, input_9, input_10, input_11, input_12, input_13, input_14, input_15, input_16, input_17, input_18, input_19, input_20, input_21, input_22, input_23, input_24, input_25, input_26, input_27, input_28, input_29, input_30, input_31, input_32, input_33, input_34, input_35, input_36};
    
    std::vector<const char*> outputs = {output_1, output_2, output_3, output_4, output_5, output_6, output_7, output_8, output_9, output_10, output_11, output_12, output_13, output_14, output_15, output_16, output_17, output_18, output_19, output_20, output_21, output_22, output_23, output_24, output_25, output_26, output_27, output_28, output_29, output_30, output_31, output_32, output_33, output_34, output_35, output_36};

    std::vector<const char*> answers_letters = {answer_letters_1, answer_letters_2, answer_letters_3, answer_letters_4, answer_letters_5, answer_letters_6, answer_letters_7, answer_letters_8, answer_letters_9, answer_letters_10, answer_letters_11, answer_letters_12, answer_letters_13, answer_letters_14, answer_letters_15, answer_letters_16, answer_letters_17, answer_letters_18, answer_letters_19, answer_letters_20, answer_letters_21, answer_letters_22, answer_letters_23, answer_letters_24, answer_letters_25, answer_letters_26, answer_letters_27, answer_letters_28, answer_letters_29, answer_letters_30, answer_letters_31, answer_letters_32, answer_letters_33, answer_letters_34, answer_letters_35, answer_letters_36};

    std::vector<const char*> answers_binary = {answer_binary_1, answer_binary_2, answer_binary_3, answer_binary_4, answer_binary_5, answer_binary_6, answer_binary_7, answer_binary_8, answer_binary_9, answer_binary_10, answer_binary_11, answer_binary_12, answer_binary_13, answer_binary_14, answer_binary_15, answer_binary_16, answer_binary_17, answer_binary_18, answer_binary_19, answer_binary_20, answer_binary_21, answer_binary_22, answer_binary_23, answer_binary_24, answer_binary_25, answer_binary_26, answer_binary_27, answer_binary_28, answer_binary_29, answer_binary_30, answer_binary_31, answer_binary_32, answer_binary_33, answer_binary_34, answer_binary_35, answer_binary_36};


    int total_elements_in_file = 0;
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
        
        std::vector<char> input_string;
        std::vector<size_t> output_value;
        
        if(i >= 0){
            std::ifstream infile(answers_letters.at(i));
            std::ofstream outfile(answers_binary.at(i));
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
                    

// used for bigger datasets
/*BOOST_AUTO_TEST_CASE( canConvertDataset )
{
    haplo::DataConverter converter(input_7, input_8, output_7);
    //converter.print_dataset();
}*/



BOOST_AUTO_TEST_SUITE_END()
