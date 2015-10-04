#include "data_converter.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

namespace haplo {

namespace io = boost::iostreams;
using namespace io;

DataConverter::DataConverter(const char* data_file):
_data(0), _rows(0), _columns(0), _aBase(0), _cBase(0), _tBase(0), _gBase(0)
{
    convert_data(data_file);
    //write_data_to_file();
}

void DataConverter::print() const 
{
    for (const auto& element : _data) std::cout << element;
}

void DataConverter::convert_data(const char* data_file)
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
        find_base_occurrance(line);
        _rows++;
    }
    
    determine_ref_sequence();
    
    for (const auto& line : lines) {
        process_line(line);
    }
    
    if (file.is_open()) file.close();
}
    
template <typename TP>    
void DataConverter::find_base_occurrance(const TP& line)
{
    // Add some comments about what's happening
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

void DataConverter::determine_ref_sequence()
{
    // Underscores, not camelcase -- base_occurance
    size_t baseOccurance[4];
    // spaces between for and the bracket
    // ++i not i++
    // not i=0, i = 0
    // use size_t, int wont be big enough
    for (size_t i = 0 ; i < _columns; i++) {
        baseOccurance[0] = {_aBase.at(i)};
        baseOccurance[1] = {_cBase.at(i)};
        baseOccurance[2] = {_tBase.at(i)};
        baseOccurance[3] = {_gBase.at(i)};
        
        std::sort(baseOccurance, baseOccurance + 4);
        
        if(_aBase.at(i) == baseOccurance[3])
            _refSeq.push_back('a');
        else if(_cBase.at(i) == baseOccurance[3])
            _refSeq.push_back('c');
        else if(_tBase.at(i) == baseOccurance[3])
            _refSeq.push_back('t');
        else if(_gBase.at(i) == baseOccurance[3])
            _refSeq.push_back('g');
        
        if(_aBase.at(i) == baseOccurance[2])
            _altSeq.push_back('a');
        else if(_cBase.at(i) == baseOccurance[2])
            _altSeq.push_back('c');
        else if(_tBase.at(i) == baseOccurance[2])
            _altSeq.push_back('t');
        else if(_gBase.at(i) == baseOccurance[2])
            _altSeq.push_back('g');
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
        if(line.at(i) == _refSeq.at(i)) {
            _data.push_back('1');
        } else if(line[i] == '-') {
            _data.push_back('-');
        } else {
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
