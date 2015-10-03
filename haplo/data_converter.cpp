#include "data_converter.hpp"

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/tokenizer.hpp>

#include <iostream>

namespace haplo {

namespace io = boost::iostreams;
using namespace io;

DataConverter::DataConverter(const char* data_file)
: _data(0)
{
    convert_data(data_file);   
}

void DataConverter::print() const 
{
    // Print funcitonns ...
    // 
}

void DataConverter::convert_data(const char* data_file)
{
    // Change this too
    
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
    size_t line_number = 0;
   
    // Get the data and store it in the data container 
    for (auto line = lines.begin(); line != lines.end(); ++line) 
        process_line(line_number++, line);
    
    if (file.is_open()) file.close(); 
}

template <typename TP>
void DataConverter::process_line(size_t line_number, TP token_pointer) 
{
    // Put your stuff here
    
    // Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> wspace_separator{" "}; 

    // Create a string from the line token and tokenize it
    std::string line(*token_pointer);
    tokenizer   elements{line, wspace_separator};
    
    for (auto& e : elements) {
        // Tokenizer creates a string, but because of the way we tokenized it
        // we know that it only has 1 element, so convert to char
        switch (e[0]) {
            case '0':
                break;
            case '1':
                break;
            case '-':
                break;
            default:
                std::cerr << "Error reading input data - exiting =(\n";
                exit(1);
        }
    }
}

}               // End namespcea haplo
