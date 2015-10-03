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
    write_data_to_file();
}

void DataConverter::print() const 
{
    int counter = 0;
    
    for(int index = 0; index < _data.size(); index++)
    {
        std::cout << _data.at(index) << " ";
        counter++;
        if(counter == _columns)
        {
           std::cout << std::endl;
            counter = 0;
        }
    }
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
    
    // Tokenize the data into lines
    tokenizer lines{data, nwline_separator};

    std::string line(*(lines.begin()));
    
    _columns = line.length();
    
    _aBase.resize(_columns);
    _cBase.resize(_columns);
    _tBase.resize(_columns);
    _gBase.resize(_columns);
    
    
    // Get the data and store it in the data container
    for (auto line = lines.begin(); line != lines.end(); ++line)
    {
        _rows++;
        find_base_occurrance(_rows, line);
    }
    
    determine_ref_sequence();
    
    size_t rows = 0;
    
    for (auto line = lines.begin(); line != lines.end(); ++line)
    {
        process_line(rows++, line);
    }
    
    if (file.is_open()) file.close();


}
    
template <typename TP>    
void DataConverter::find_base_occurrance(size_t line_number, TP token_pointer)
{
    std::string line(*token_pointer);
    
        for(int i = 0; i < _columns; i++)
        {
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
    int baseOccurance[4];
    for(int i=0 ; i < _columns; i++)
    {
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
    
    /*for(int j= 0; j < _refSeq.size(); j++)
    {
        std::cout << _refSeq.at(j);
        
    }
    std::cout << std::endl;
    for(int j= 0; j < _refSeq.size(); j++)
    {
        std::cout << _altSeq.at(j);
        
    }
    std::cout << std::endl;*/
    
}
    
    
template <typename TP>
void DataConverter::process_line(size_t line_number, TP token_pointer) 
{
    //test examples
    //std::string referenceAllele = "ttcaaaataatgcactgtgaccaacccttt";
    //std::string alternateAllele = "cagctgtgctgcacaaacataagtaatcaa";
    
    std::string line(*token_pointer);


    for (int i = 0; i < _columns; i++) {

        if(line.at(i) == _refSeq.at(i))
        {
            _data.push_back('1');
        }
        else if(line[i] == '-')
        {
            _data.push_back('-');
        }
        else
        {
            _data.push_back('0');
            
        }
        
    }
}
    
void DataConverter::write_data_to_file()
{
    std::ofstream myfile;
    myfile.open("output_files/input_simulated_2_ouput.txt");
    
    int counter = 0;
    
    for(int index = 0; index < _data.size(); index++)
    {
        myfile << _data.at(index) << " ";
        counter++;
        if(counter == _columns)
        {
            myfile << std::endl;
            counter = 0;
        }
    }

    myfile.close();

}
    
}// End namespcea haplo
