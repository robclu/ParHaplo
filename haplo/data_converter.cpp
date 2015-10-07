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
    
#define ZERO    0x00
#define ONE     0x01
#define TWO     0x02
#define THREE   0x03

DataConverter::DataConverter(const char* data_file_1, const char* data_file_2, const char* data_file_3):
_data(0), _rows(0), _columns(0), _base_a(0), _base_c(0), _base_t(0), _base_g(0), _chromosome(0)
{
    process_ground_truth(data_file_1);
    read_in_simulated_data(data_file_2);
    write_simulated_data_to_file(data_file_3);
}
    
void DataConverter::printMap() const
{
    int counter = 1;
    int num_chromosomes = 1;
    
    for(auto& chromosome: _ref_alt_chromosomes){
        std::cout << "Chromosome: " << num_chromosomes << std::endl;
        
        for (auto& elements: chromosome){
            
            for(auto& y: elements.second){
                std::cout << std::endl << counter << " : " <<  elements.first << std::endl;
                std::cout << ": " << std::endl;
                y.print();
                std::cout << std::endl;
                counter++;

            }
        }
        num_chromosomes++;
        
    }
    
    counter = 0;
    num_chromosomes = 1;
    
    for(auto& chromosome: _simulated_chromosomes){
        std::cout << "Chromosome: " << num_chromosomes << std::endl;
        for (auto& x: chromosome){
            std::cout << std::endl << counter << " : " <<  x.first << std::endl;
            std::cout << ": " << std::endl;
            for(auto& y: x.second){
                      y.print();
            }
            counter++;
            
        }
        num_chromosomes++;
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
    
    _base_a.resize(_columns);
    _base_c.resize(_columns);
    _base_t.resize(_columns);
    _base_g.resize(_columns);
    
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
    
void DataConverter::process_ground_truth(const char* data_file)
{
    io::mapped_file_source file(data_file);
    if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
    
    std::string data(file.data(), file.size());
    
    // Create a tokenizer to tokenize by newline character
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> nwline_separator{"\n"};
    
    // Tokenize the data into lines
    tokenizer lines{data, nwline_separator};
    
    size_t header_length = 5;
    
    for (const auto& line : lines) {
        if(header_length == 0) {
            determine_dataset_ref_sequence(line);
            _rows++;
        }
        else {
            header_length--;
        }

    }
    _chromosome = 0;
    if (file.is_open()) file.close();

    
}
    
void DataConverter::read_in_simulated_data(const char* data_file)
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
    
    
        for (const auto& line : lines) {
            convert_dataset_to_binary(line);
            _rows++;
        }
    
        _chromosome = 0;
    
        process_each_read();
        
        if (file.is_open()) file.close();
}

    
template <typename TP>    
void DataConverter::find_base_occurrance(const TP& line)
{
    // count occurances of the bases (a,c,t,g) in each column of a line
    for(size_t i = 0; i < _columns; ++i) {
        if(line[i] == 'a')
            _base_a.push_back(_base_a.at(i)++);
        else if(line[i] == 'c')
            _base_c.push_back(_base_c.at(i)++);
        else if(line[i] == 't')
            _base_t.push_back(_base_t.at(i)++);
        else if(line[i] == 'g')
            _base_g.push_back(_base_g.at(i)++);
    }
}

void DataConverter::determine_simulated_ref_sequence()
{
    size_t base_occurance[4];
    
    for (size_t i = 0 ; i < _columns; i++) {
        base_occurance[0] = {_base_a.at(i)};
        base_occurance[1] = {_base_c.at(i)};
        base_occurance[2] = {_base_t.at(i)};
        base_occurance[3] = {_base_g.at(i)};
        
        std::sort(base_occurance, base_occurance + 4);
        
        if(_base_a.at(i) == base_occurance[3])
            _ref_seq.push_back('a');
        else if(_base_c.at(i) == base_occurance[3])
            _ref_seq.push_back('c');
        else if(_base_t.at(i) == base_occurance[3])
            _ref_seq.push_back('t');
        else if(_base_g.at(i) == base_occurance[3])
            _ref_seq.push_back('g');
        
        if(_base_a.at(i) == base_occurance[2])
            _alt_seq.push_back('a');
        else if(_base_c.at(i) == base_occurance[2])
            _alt_seq.push_back('c');
        else if(_base_t.at(i) == base_occurance[2])
            _alt_seq.push_back('t');
        else if(_base_g.at(i) == base_occurance[2])
            _alt_seq.push_back('g');
    }

    
}
    
template <typename TP>
void DataConverter::determine_dataset_ref_sequence(const TP& token_pointer)
{
        //Create a tokenizer to tokenize by whitespace
        using tokenizer = boost::tokenizer<boost::char_separator<char>>;
        boost::char_separator<char> wspace_separator{" "};

        tokenizer   elements{token_pointer, wspace_separator};
        size_t column_counter = 0;
        size_t position = 0;
        char ref_base;
        char alt_base;
        size_t haplo_one = 2; // if no haplotype present, defaults to 2
        size_t haplo_two = 2;
    
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
        _ref_alt_chromosomes[_chromosome - 1][position].push_back(Base(ref_base, alt_base, true, haplo_one, haplo_two));

}
    
    
template <typename TP>
void DataConverter::convert_dataset_to_binary(const TP& token_pointer)
{
        //Create a tokenizer to tokenize by newline character and another by whitespace
        using tokenizer = boost::tokenizer<boost::char_separator<char>>;
        boost::char_separator<char> wspace_separator{" "};
        
        tokenizer   elements{token_pointer, wspace_separator};
        size_t column_counter = 0;
        size_t start_position = 0;
        size_t end_position = 0;
        std::string cigar_value;
        std::string sequence;
    
        size_t temp_chromosome = _chromosome;
    
    for (auto& e : elements) {

        if(column_counter == 2){
            
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
        else if(_chromosome >=1 && column_counter > 2){
            
            //position column
            if(column_counter == 3){
                start_position = std::stoul(e);
                
                if(temp_chromosome != _chromosome)
                    _start_of_chromosome_reads.push_back(start_position);
            }
            //cigar value column
            else if(column_counter == 5){
                cigar_value = e;
            }
            //read sequence column
            else if(column_counter == 9){
                sequence = e;
            }
            
        }
        
        column_counter++;
    }
    end_position = start_position + sequence.length();
    process_cigar_value(start_position, end_position, cigar_value, sequence);
    store_simulated_data(_chromosome, start_position, end_position, sequence);
    
}
    
void DataConverter::store_simulated_data(size_t chromosome, size_t start_position, size_t end_position, std::string sequence)
{
            if (_simulated_chromosomes.at(chromosome - 1).find(start_position) != _simulated_chromosomes.at(chromosome - 1).end()) {
                // Key is in the list, so add to the vector
                auto& existing_reads = _simulated_chromosomes.at(chromosome - 1)[start_position];
                existing_reads.push_back(Read(end_position, sequence));
            }
            else {
                std::vector<Read> first_read = {Read(end_position, sequence)};
                _simulated_chromosomes.at(chromosome - 1)[start_position] = first_read;
            }

        
}
    
void DataConverter::process_cigar_value(size_t& start_position, size_t& end_position, std::string& cigar_value, std::string& sequence)
{
    //7 operations for CIGAR notation
    //M - a match i.e. consistent reads
    //D - deletion i.e. missing base in sequence
    //I - insertion i.e. missing base in reference
    //S - soft clipping i.e. unaligned bases (can be removed)
    //H - hard clipping i.e. can be removed
    //P - padding i.e. represented by a gap
    //N - skipped bases on the reference i.e. represented by gaps
    //this vector stores the positions of characters in the cigar value
    std::vector<size_t> operations;
    std::vector<size_t> operation_positions;
    int counter = 0;
    size_t operation_value_start = 0;
    for(size_t i = 0; i < cigar_value.length(); ++i){
        if(isalpha(cigar_value[i])){
            counter++;
            operation_positions.push_back(i);
            if(counter == 1)
                operations.push_back(std::stoul(cigar_value.substr(0, i)));
            else{
                operation_value_start = operation_positions.rbegin()[1] + 1;
                operations.push_back(std::stoul(cigar_value.substr(operation_value_start, i - operation_value_start)));
            }

        }
        
    }
    
    
    std::vector<std::string> operated_portions;
    size_t current_position_in_read = 0;
    std::string temp;
    char value_of_operation = 0;
    size_t num_of_operations = 0;
    for(size_t i = 0; i < operation_positions.size(); ++i){
        
            value_of_operation = cigar_value.at(operation_positions.at(i));
            num_of_operations = operations.at(i);
        
            if(value_of_operation == 'M'){
                operated_portions.push_back(sequence.substr(current_position_in_read,num_of_operations));
                current_position_in_read = current_position_in_read + operations.at(i);
            }
            else if (value_of_operation== 'D' || value_of_operation == 'N' || value_of_operation == 'P' || value_of_operation == 'S' || value_of_operation == 'H' ){
                temp.clear();
                for(int index = 0; index < num_of_operations; ++index){
                    temp.insert(0,"2");
                }
                
                if(i == 0){
                    start_position = start_position + current_position_in_read + num_of_operations;
                    current_position_in_read = current_position_in_read + num_of_operations;
                }
                else if(i == operation_positions.size() - 1){
                    end_position-=num_of_operations;
                }
                else {
                    end_position+=num_of_operations;
                    operated_portions.push_back(temp);
                }
            }
            else if(value_of_operation == 'I'){
                operated_portions.push_back(sequence.substr(current_position_in_read,num_of_operations));

            }
       
    }
    sequence.clear();
    for(int i = 0; i < operated_portions.size(); ++i){
        sequence+=operated_portions.at(i);
    }
    
    operated_portions.clear();
    operations.clear();
    operation_positions.clear();
    
    
}
    
void DataConverter::process_each_read()
{
    std::string zero = "0";
    std::string one = "1";
    std::string two = "2";
    size_t no_of_chromosomes = 0;
    
    for(auto& chromosome: _simulated_chromosomes){
        
        for (auto& all_reads: chromosome){
            //read = each read in all_reads.second vector
            for(auto& read: all_reads.second){
                for(int i = 0; i < read._sequence.length(); ++i){
                    
                    if (_ref_alt_chromosomes.at(no_of_chromosomes).find(all_reads.first + i) == _ref_alt_chromosomes.at(no_of_chromosomes).end()) {
                        // Position is not in references bases, add new "assumed" base
                        _ref_alt_chromosomes[no_of_chromosomes][all_reads.first + i].push_back(Base(read._sequence.at(i), 'A', 0, 2, 2));
                    }
                    Base& base_at_position = _ref_alt_chromosomes[no_of_chromosomes][all_reads.first + i].at(0);
                    
                    if (read._sequence.at(i) == '2'){
                        read._binary_sequence+=two;
                    }
                    else if(read._sequence.at(i) == base_at_position._ref_base){
                        read._binary_sequence+=one;
                    }
                    else {
                        read._binary_sequence+=zero;
                    }
                }
                
            }
            
        }
        no_of_chromosomes++;
        
    }
    
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
    for (size_t i = 0; i < _columns; ++i) {
        if(line[i] == _ref_seq.at(i)) {
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
    
    
void DataConverter::write_simulated_data_to_file(const char* filename)
{
    std::ofstream myfile;
    std::vector<size_t> positions;
    myfile.open (filename);
    
    
    size_t start_position = 0;
    size_t num_of_elements_in_map = 0;
    for(int chromosome = 0; chromosome < 1; ++chromosome){
        start_position = _start_of_chromosome_reads.at(chromosome);
        num_of_elements_in_map = _simulated_chromosomes.at(chromosome).size();
        
        for (size_t i = start_position;; ++i) {
                if (_simulated_chromosomes.at(chromosome).find(i) != _simulated_chromosomes.at(chromosome).end()) {
                    for(int vec = 0; vec < _simulated_chromosomes.at(chromosome).at(i).size(); ++vec){
                        //std::cout << i << " " << _chr1_simulated_data.at(i).at(vec)._end_position << " " << _chr1_simulated_data.at(i).at(vec)._binary_sequence << std::endl;
                        myfile << i << " " << _simulated_chromosomes.at(chromosome).at(i).at(vec)._end_position << " " << _simulated_chromosomes.at(chromosome).at(i).at(vec)._binary_sequence << std::endl;
                    }
                    num_of_elements_in_map--;
                    
                }
                if(num_of_elements_in_map == 0)
                    break;
                
        }
        
        
    }
    
    myfile.close();

  
}
    
}// End namespcea haplo
