#include "data_converter.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>

#include <iostream>
#include <vector>
#include <fstream>

namespace haplo {

namespace io = boost::iostreams;
using namespace io;

DataConverter::DataConverter(const char* data_file):
_data(0)
{
    convert_data(data_file);   
}

void DataConverter::print() const 
{
    int counter = 0;
    
    //size of the file
    for(int rows = 0; rows < 27*30; rows++)
    {
        std::cout << _data.at(rows) << " ";
        counter++;
        if(counter == 30)
        {
           std::cout << std::endl;
            counter = 0;
        }
    }
    // Print funcitonns ...
    // 
}

void DataConverter::convert_data(const char* data_file)
{
    // Open file and convert to string (for tokenizer)
    io::mapped_file_source file(data_file);
    if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
    
    std::string data(file.data(), file.size());
    
    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> nwline_separator{"\n"};
    
    // Tokenize the data into lines
    tokenizer lines{data, nwline_separator};
    
    // Counter for the row id
    size_t row = 0;
    
    // Get the data and store it in the data container
    for (auto line = lines.begin(); line != lines.end(); ++line)
        process_line(row++, line);
    
    if (file.is_open()) file.close();


}

template <typename TP>
void DataConverter::process_line(size_t line_number, TP token_pointer) 
{
    std::string referenceAllele = "ttcaaaataatgcactgtgaccaacccttt";
    std::string alternateAllele = "cagctgtgctgcacaaacataagtaatcaa";
    
    std::string line(*token_pointer);

    
    for (int i = 0; i < line.length(); i++) {

        if(line.at(i) == referenceAllele[i])
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

}               // End namespcea haplo
