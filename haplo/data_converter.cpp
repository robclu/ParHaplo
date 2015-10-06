#include "data_converter.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

namespace haplo {

namespace io = boost::iostreams;
using namespace io;
    
#define ZERO    0x00
#define ONE     0x01
#define TWO     0x02
#define THREE   0x03

DataConverter::DataConverter(const char* data_file):
_data(0), _rows(0), _columns(0), _aBase(0), _cBase(0), _tBase(0), _gBase(0)
{
    convert_dataset_to_binary(data_file);
}
    
void DataConverter::printMap() const
{
    int counter = 0;
    for (auto& x: _chr1_ref_and_alt_seq){
        counter++;
        std::cout << std::endl << counter << ": " <<  x.first << std::endl;
        std::cout << ": ";
        (x.second).print();
        std::cout << std::endl;
        

    }
    
    std::cout << "chromosome 2" << std::endl;
    
    for (auto& x: _chr2_ref_and_alt_seq){
        counter++;
        std::cout << counter << ": " <<  x.first << ": ";
        (x.second).print();
        std::cout << std::endl;
        
        
    }
    
    
}


void DataConverter::print() const 
{
    //for (const auto& element : _data) std::cout << element;
    
}

void DataConverter::convert_simulated_data_to_binary(const char* data_file)
{
    // Open file and convert to string (for tokenizer)
    io::mapped_file_source file(data_file);
    if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
    
    std::string data(file.data(), file.size());
    
    // Create a tokenizer to tokenize by newline character
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> nwline_separator{"\n"};
    
    // Tokenize the data into lines and create a 
    tokenizer lines{data, nwline_separator};

    _columns = lines.begin()->length();
    
    _aBase.resize(_columns);
    _cBase.resize(_columns);
    _tBase.resize(_columns);
    _gBase.resize(_columns);
    
    // Get the data and store it in the data container
    for (const auto& line : lines) {
       // std::cout << line << std::endl;
        find_base_occurrance(line);
        _rows++;
    }
    
    determine_simulated_ref_sequence();
    
    for (const auto& line : lines) {
        process_line(line);
    }
    
    if (file.is_open()) file.close();
}
    
void DataConverter::convert_dataset_to_binary(const char* data_file)
{
        // Open file and convert to string (for tokenizer)
        io::mapped_file_source file(data_file);
        if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
        
        std::string data(file.data(), file.size());
        
        // Create a tokenizer to tokenize by newline character
        using tokenizer = boost::tokenizer<boost::char_separator<char>>;
        boost::char_separator<char> nwline_separator{"\n"};
        
        // Tokenize the data into lines and create a
        tokenizer lines{data, nwline_separator};
        
        //_columns = lines.begin()->length();
        size_t header_length = 5;
    
        for (const auto& line : lines) {
            //std::cout << "line: " << line << std::endl;
            if(header_length == 0) {
                //std::cout << header_length << std::endl;
                determine_dataset_ref_sequence(line);
                _rows++;
            }
            else {
                  header_length--;
            }
          
            // std::cout << "i am here 1" << std::endl;
        }
    
        //std::cout << "i am here 3" << std::endl;
        /*for (const auto& line : lines) {
            process_line(line);
            //row++;
        }*/
        
        if (file.is_open()) file.close();
}

    
template <typename TP>    
void DataConverter::find_base_occurrance(const TP& line)
{
    // count occurances of the bases (a,c,t,g) in each column of a line
    for(size_t i = 0; i < _columns; ++i) {
        if(line[i] == 'a')
            _aBase.push_back(_aBase.at(i)++);
        else if(line[i] == 'c')
            _cBase.push_back(_cBase.at(i)++);
        else if(line[i] == 't')
            _tBase.push_back(_tBase.at(i)++);
        else if(line[i] == 'g')
            _gBase.push_back(_gBase.at(i)++);
    }
}

void DataConverter::determine_simulated_ref_sequence()
{
    size_t base_occurance[4];
    // spaces between for and the bracket
    // ++i not i++
    // not i=0, i = 0
    // use size_t, int wont be big enough
    for (size_t i = 0 ; i < _columns; i++) {
        base_occurance[0] = {_aBase.at(i)};
        base_occurance[1] = {_cBase.at(i)};
        base_occurance[2] = {_tBase.at(i)};
        base_occurance[3] = {_gBase.at(i)};
        
        std::sort(base_occurance, base_occurance + 4);
        
        if(_aBase.at(i) == base_occurance[3])
            _refSeq.push_back('a');
        else if(_cBase.at(i) == base_occurance[3])
            _refSeq.push_back('c');
        else if(_tBase.at(i) == base_occurance[3])
            _refSeq.push_back('t');
        else if(_gBase.at(i) == base_occurance[3])
            _refSeq.push_back('g');
        
        if(_aBase.at(i) == base_occurance[2])
            _altSeq.push_back('a');
        else if(_cBase.at(i) == base_occurance[2])
            _altSeq.push_back('c');
        else if(_tBase.at(i) == base_occurance[2])
            _altSeq.push_back('t');
        else if(_gBase.at(i) == base_occurance[2])
            _altSeq.push_back('g');
    }

    
}
    
template <typename TP>
void DataConverter::determine_dataset_ref_sequence(const TP& token_pointer)
{
        //Create a tokenizer to tokenize by newline character and another by whitespace
        using tokenizer = boost::tokenizer<boost::char_separator<char>>;
        boost::char_separator<char> wspace_separator{" "};

        tokenizer   elements{token_pointer, wspace_separator};
        size_t column_counter = 0;
        size_t position = 0;
        char ref_base;
        char alt_base;
        size_t haplo_one;
        size_t haplo_two;
    
        for (auto& e : elements) {
            
            if(column_counter == 0){
                
                if(e == "chr1")
                    _chromosome = 1;
                else if(e == "chr2")
                    _chromosome = 2;
                else if(e == "chr3")
                    _chromosome = 3;
                else if(e == "chr4")
                    _chromosome = 4;
                else if(e == "chr5")
                    _chromosome = 5;
                else if(e == "chr6")
                    _chromosome = 6;
                else if(e == "chr7")
                    _chromosome = 7;
                else if(e == "chr8")
                    _chromosome = 8;
                else if(e == "chr9")
                    _chromosome = 9;
                else if(e == "chr10")
                    _chromosome = 10;
                else if(e == "chr11")
                    _chromosome = 11;
                else if(e == "chr12")
                    _chromosome = 12;
                else if(e == "chr13")
                    _chromosome = 13;
                else if(e == "chr14")
                    _chromosome = 14;
                else if(e == "chr15")
                    _chromosome = 15;
                else if(e == "chr16")
                    _chromosome = 16;
                else if(e == "chr17")
                    _chromosome = 17;
                else if(e == "chr18")
                    _chromosome = 18;
                else if(e == "chr19")
                    _chromosome = 19;
                else if(e == "chr20")
                    _chromosome = 20;
                else if(e == "chr21")
                    _chromosome = 21;
                else if(e == "chr22")
                    _chromosome = 22;
                else
                    std::cout << "Invalid chromosome" << std::endl;
            }
            else if(_chromosome >=1 && column_counter > 0){
                
                //position column
                if(column_counter == 1){
                    position = std::stoul(e);
                }
                //ref sequence column (first character)
                else if(column_counter == 3){
                    ref_base = e[0];
                }
                //alt sequence column (first character)
                else if(column_counter == 4){
                    alt_base = e[0];
                }
                //ground truth column 1|1
                else if(column_counter == 9){
                    if(e[1] == '|') {
                        haplo_one = e[0] - 48;
                        haplo_two = e[2] - 48;
                    }
                   
                   
                }
 
            }
            else {
                        
                std::cerr << "Error reading input data - exiting =(\n";
                exit(1);
                    
            }
            
            column_counter++;
        }
    
        storeBaseData(_chromosome, position, ref_base, alt_base, true, haplo_one, haplo_two);

}
    
void DataConverter::storeBaseData(size_t chromosome, size_t position, char ref_base, char alt_base, bool real, size_t haplo_one, size_t haplo_two)
{
    if(chromosome == 1)
        _chr1_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 2)
        _chr2_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 3)
        _chr3_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 4)
        _chr4_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 5)
        _chr5_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 6)
        _chr6_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 7)
        _chr7_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 8)
        _chr8_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 9)
        _chr9_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 10)
        _chr10_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 11)
        _chr11_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 12)
        _chr12_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 13)
        _chr13_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 14)
        _chr14_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 15)
        _chr15_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 16)
        _chr16_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 17)
        _chr17_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 18)
        _chr18_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 19)
        _chr19_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 20)
        _chr20_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 21)
        _chr21_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else if(chromosome == 22)
        _chr22_ref_and_alt_seq[position] = Base(ref_base, alt_base, real, haplo_one, haplo_two);
    else
        std::cout << "Invaid chromosomoe" << std::endl;

 }

    
byte DataConverter::convert_char_to_byte(char input)
{
    switch(input){
            case 'A':
                return ZERO;
            case 'C':
                return ONE;
            case 'T':
                return TWO;
            case 'G':
                return THREE;
    }
    
}
                    
char DataConverter::convert_byte_to_char(byte input)
{
    switch(input){
        case ZERO:
            return 'a';
        case ONE:
            return 'c';
        case TWO:
            return 't';
        case THREE:
            return 'g';
    }
    
}
    
template <typename TP>
void DataConverter::process_line(const TP& line) 
{
    //test examples
    //std::string referenceAllele = "ttcaaaataatgcactgtgaccaacccttt";
    //std::string alternateAllele = "cagctgtgctgcacaaacataagtaatcaa";

    // { brackets go on the same line inside functions 
    // only use { on a new line for the start of a function
    for (size_t i = 0; i < _columns; ++i) {
        if(line[i] == _refSeq.at(i)) {
            _data.push_back('1');
        }
        else if(line[i] == '-') {
            _data.push_back('-');
        }
        else {
            _data.push_back('0');
        }
        // Add in the space
        _data.push_back(' ');
    }
    // Add in the newline character
    _data.push_back('\n');
}

void DataConverter::write_data_to_file(const char* filename)
{
    // Create the parameters for the output file
    io::mapped_file_params file_params(filename);
    
    // Resize the file
    file_params.new_file_size = sizeof(char) * _data.size();
    
    // Open file and check that it opened
    io::mapped_file_sink output_file(file_params);
    if (!output_file.is_open()) throw std::runtime_error("Could not open output file =(!\n");
    
    // Copy the data to the file
    std::copy(_data.begin(), _data.end(), output_file.data());
    
    // Close the file
    if (output_file.is_open()) output_file.close();
}
    
}// End namespcea haplo
