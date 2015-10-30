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

// All outputs to the command line are used for debugging
DataConverter::DataConverter(const char* data_file_1, const char* data_file_2, const char* data_file_3):
_data(0), _rows(0), _columns(0), _base_a(0), _base_c(0), _base_t(0), _base_g(0), _chromosome(0), _start_reading(0)
{
    convert_dataset_to_binary(data_file_1, data_file_2);
    //std::cout << "I converted to binary" << std::endl;
    write_dataset_to_file(data_file_3);
    //std::cout << "I outputted to file" << std::endl;
}
    
DataConverter::DataConverter(const char* data_file_1, const char* data_file_2):
_data(0), _rows(0), _columns(0), _base_a(0), _base_c(0), _base_t(0), _base_g(0), _chromosome(0)
{
    convert_simulated_data_to_binary(data_file_1);
    write_simulated_data_to_file(data_file_2);
    _total_num_elements = 0;
}
    
void DataConverter::print_dataset() const
{
    int counter = 1;
    int num_chromosomes = 1;
    /*
    // read through all chromosomes
    for(auto& chromosome: _ref_alt_chromosomes){
        std::cout << "Chromosome: " << num_chromosomes << std::endl;
        // all elements in the map per chromosome
        for (auto& elements: chromosome){
            // all elements in the vector of bases
            for(auto& base: elements.second){
                std::cout << counter << std::endl << "start position: " <<  elements.first << std::endl;
                base.print(); // print members of Base
                std::cout << std::endl;
                counter++;
            }
        }
        num_chromosomes++;
        
    }
    
    counter = 0;
    num_chromosomes = 1;
    
    // read through all chromosomes
    for(auto& chromosome: _simulated_chromosomes){
        std::cout << "Chromosome: " << num_chromosomes << std::endl;
        // all elements in the map per chromosome
        for (auto& elements: chromosome){
            std::cout << counter << std::endl << "start position: " <<  elements.first << std::endl;
            for(auto& read: elements.second){
                read.print(); // print members of read
                std::cout << std::endl;
            }
            counter++;
        }
        num_chromosomes++;
    }*/
}


void DataConverter::print_simulated() const
{
    //for (const auto& element : _data) std::cout << element;
    
    //std::cout << "Ref Sequence: " << std::endl;
    //for (const auto& element : _ref_seq) std::cout << element;
    
    //std::cout << std::endl;
    
    //std::cout << "Alt Sequence: " << std::endl;
    //for (const auto& element : _alt_seq) std::cout << element;
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
    
    // Resize vectors to incorporate number of columns in file
    _base_a.resize(_columns);
    _base_c.resize(_columns);
    _base_t.resize(_columns);
    _base_g.resize(_columns);
    
    for (const auto& line : lines) {
        find_base_occurrance(line);
        _rows++;
    }
    
    // Use processed data to formulate a ref sequence
    determine_simulated_ref_sequence();
    
    // Do actual ACTG > binary conversion
    for (const auto& line : lines) {
        process_each_line(line);
    }
    
    //size_t total_num_elements_value = 0;
    
    for(const auto& elements: _elements_per_line){
        _total_num_elements+=elements;
    }
    
    // If want to put total number of elements at the end of the file
    /*_total_num_elements = std::to_string(total_num_elements_value);
    
    for(int i = 0; i < total_num_elements.length(); ++i){
        _data.push_back(total_num_elements[i]);
    }*/

    
    if (file.is_open()) file.close();
}
    
size_t DataConverter::getTotalNumElements()
{
    return _total_num_elements;
}
    
void DataConverter::convert_dataset_to_binary(const char* data_file_1, const char* data_file_2)
{
    // Determine reference sequence by processing ground truth file
    process_ground_truth(data_file_1);
    //std::cout << "I processed ground truth" << std::endl;
    
    // Open file and convert to string (for tokenizer)
    io::mapped_file_source file(data_file_2);
    if (!file.is_open()) throw std::runtime_error("Could not open input file =(!\n");
        
    std::string data(file.data(), file.size());
        
    // Create a tokenizer to tokenize by newline character
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> nwline_separator{"\n"};
        
    // Tokenize the data into lines
    tokenizer lines{data, nwline_separator};
    
    // Process and store important data
    for (const auto& line : lines) {
        process_dataset(line);
        _rows++;
    }
    //std::cout << "I processed dataset" << std::endl;
    
    // Reset chromosome count in case new files need to be processed
    _chromosome = 0;
    
    // Do actual ACTG > binary conversion
    process_each_read();
    //std::cout << "I converted to binary" << std::endl;
    
    store_haplotype_answers();
        
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
    
    // Set header length to remove header when reading data
    size_t header_length = 5;
        
    for (const auto& line : lines) {
        if(header_length == 0) {
            // Use data to construct a map of reference bases
            determine_dataset_ref_sequence(line);
            _rows++;
        }
        else {
            header_length--;
        }
            
    }
    //print_dataset();
    // Reset chromosome count in case new files need to be processed
    _chromosome = 0;
    _start_reading = 0;
    
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
        
        // Find two most occurring bases
        std::sort(base_occurance, base_occurance + 4);
        
        // Assign most occurring base to reference sequence
        if(_base_a.at(i) == base_occurance[3])
            _ref_seq.push_back('a');
        else if(_base_c.at(i) == base_occurance[3])
            _ref_seq.push_back('c');
        else if(_base_t.at(i) == base_occurance[3])
            _ref_seq.push_back('t');
        else if(_base_g.at(i) == base_occurance[3])
            _ref_seq.push_back('g');
        
        // Assign second most occurring base to alternate sequence
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
void DataConverter::determine_dataset_ref_sequence(const TP& line)
{
    //Create a tokenizer to tokenize by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> wspace_separator{" "};
    
    
    // Tokenize the data into elements
    tokenizer   elements{line, wspace_separator};
    
    
    size_t column_counter = 0;  // Used to keep track of elements to store the correct information
    size_t position = 0;        // Store position of base
    char ref_base;
    char alt_base;
    size_t haplo_one = 2;       // if no haplotype present for base, defaults to 2
    size_t haplo_two = 2;
    
    for (auto& e : elements) {
        
        if(column_counter == 0){
            if(e.length() == 5){
                //store only chromosome 20, uncomment below block for all chromosomes
                if(std::stoul(e.substr(3,2)) == 20){
                    _start_reading = 1;
                }
                else{
                    _start_reading = 0;
                }
                
                
            }
            else{
                _start_reading = 0;
            }
            
            if(_start_reading == 1){
                _chromosome = std::stoul(e.substr(3,2));
            }
            else{
                break;
            }
        }

        
        // store all chromosomes
        /*if(column_counter == 0){
            if(_chromosome == 0){
                if(e[3]-48 == 1){
                    _start_reading = 1;
                }
            }
            if(_start_reading == 1 && _chromosome >=0 && _chromosome < 9){
                _chromosome = e[3] - 48;
            }
            else if(_chromosome >= 9 && _chromosome < 22){
                
                _chromosome = std::stoul(e.substr(3,2));
            }
            else{//if(_chromosome > 22)
                //std::cout << "Invalid chromosome" << std::endl;
                break;
            }
        }*/
        // Should not store chromosome if not chr1 - chr22
        else if(_chromosome >=1 && _chromosome <= 22 && column_counter > 0){
            
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
            
        } column_counter++;
    }
    // store base information in map
    if(_chromosome >= 1 && _chromosome <= 22)
        _ref_alt_chromosomes[_chromosome - 1][position].push_back(Base(ref_base, alt_base, true, haplo_one, haplo_two));

}
    
    
template <typename TP>
void DataConverter::process_dataset(const TP& token_pointer)
{
    //Create a tokenizer to tokenize by newline character and another by whitespace
    using tokenizer = boost::tokenizer<boost::char_separator<char>>;
    boost::char_separator<char> wspace_separator{"\t"};
    
    //std::cout << token_pointer << std::endl;
    
    tokenizer elements{token_pointer, wspace_separator};
    
    size_t column_counter = 0;              // Used to keep track of elements to store the correct information
    size_t start_position = 0;              // Store start position of read
    size_t end_position = 0;                // Store calculated end position of read
    std::string cigar_value;
    std::string sequence;
    size_t temp_chromosome = _chromosome;   // Store temporary indicator of chromosome to notify changes
    
    for (auto& e : elements) {

        if(column_counter == 2){
                if(e.length() == 5){
                    // Store only chromosome 20
                    if(std::stoul(e.substr(3,2)) == 20){
                        _start_reading = 1;
                    }
                    else{
                        _start_reading = 0;
                    }

                        
                }
                else{
                        _start_reading = 0;
                }

                
        // To store data for each chromosome separately
        /*if(_start_reading == 1 && _chromosome >=0 && _chromosome < 9){
                _chromosome = e[3] - 48;
            }
            else if(_chromosome >= 9 && _chromosome < 22){
                
                _chromosome = std::stoul(e.substr(3,2));
            }*/
                   if(_start_reading == 1){
                       _chromosome = std::stoul(e.substr(3,2));
                   }
            else{
                break;
            }
        }
        // Should not store chromosome if not chr1 - chr22
        else if(_chromosome >=1 && _chromosome <= 22 && column_counter > 2){
            
            // Position column
            if(column_counter == 3){
                start_position = std::stoul(e);
                
                // If chromosome indicator changes, store start position
                if(temp_chromosome != _chromosome){
                    _start_of_chromosome_reads.push_back(start_position);
                    std::cout << "chromosome: " << _chromosome << std::endl;
                }

            }
            // Cigar value column
            else if(column_counter == 5){
                cigar_value = e;
            }
            // Read sequence column
            else if(column_counter == 9){
                sequence = e;
            }
            
        }
        column_counter++;
    }
    // Process informaton and store data
    if(_chromosome >= 1 && _chromosome <= 22 ){
        end_position = start_position + sequence.length() - 1;
        process_cigar_value(start_position, end_position, cigar_value, sequence);
        store_read_data(_chromosome, start_position, end_position, sequence);
    }

    
}
    
void DataConverter::store_read_data(size_t chromosome, size_t start_position, size_t end_position, std::string sequence)
{
    // Check if a read exists for a particular position
    if (_simulated_chromosomes.at(chromosome - 1).find(start_position) != _simulated_chromosomes.at(chromosome - 1).end()) {
        // Position is in the map, so add read to the vector
        auto& existing_reads = _simulated_chromosomes.at(chromosome - 1)[start_position];
        existing_reads.push_back(Read(end_position, sequence));
    }
    else {
        // Position not in map, create and store new vector with read
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

    std::vector<size_t> operation_positions;            // Stores the positions of characters in the cigar value
    std::vector<size_t> operations;                     // Stores the number of operations that need to be carried out (indexed against operation_positions)
    int counter = 0;                                    // To check whether first operation or the remainder
    size_t operation_value_start = 0;
    
    // Read through characters in cigar string
    for(size_t i = 0; i < cigar_value.length(); ++i){
        if(isalpha(cigar_value[i])){
            counter++;
            // If character, store position in cigar string
            operation_positions.push_back(i);
            if(counter == 1)
                operations.push_back(std::stoul(cigar_value.substr(0, i)));
            else{
                // Find position of last operation and extract string between each operation
                operation_value_start = operation_positions.rbegin()[1] + 1;
                operations.push_back(std::stoul(cigar_value.substr(operation_value_start, i - operation_value_start)));
            }

        }
        
    }
    
    
    std::vector<std::string> operated_portions;        // Store parts of sequence that have been operated on
    size_t current_position_in_read = 0;               // Keeps track of position in sequence
    std::string temp;                                  // For adding gaps (represented by 2s) into the sequence
    char value_of_operation = 0;                       // Which operation to be applied (M, S, H... etc)
    size_t num_of_operations = 0;                      // How many operatiosn to be applied
    
    for(size_t i = 0; i < operation_positions.size(); ++i){
        
            value_of_operation = cigar_value.at(operation_positions.at(i));
            num_of_operations = operations.at(i);
        
            // Extract part of sequence based on num_of_operations
            if(value_of_operation == 'M'){
                
                operated_portions.push_back(sequence.substr(current_position_in_read,num_of_operations));
                current_position_in_read = current_position_in_read + operations.at(i);
                
            }
            // Add string of gaps (2s) based on num_of_operations
            else if (value_of_operation== 'D' || value_of_operation == 'N' || value_of_operation == 'P' || value_of_operation == 'S' || value_of_operation == 'H' ){
                
                // assume S only occurs at the beginning or end
                if(value_of_operation == 'S'){
                                // If occurs at the beginning, move start position of sequence (ignoring sequence elements)
                    if(i == 0){
                        current_position_in_read = current_position_in_read + num_of_operations;
                        end_position-=num_of_operations;
                    }
                    // If occurs at the end, move end position
                    else if(i == operation_positions.size() - 1){
                        end_position-=num_of_operations;
                        
                    }
                    
                }
                else if(value_of_operation == 'H'){
                    //do nothing
                }
                
                // Otherwise, store extra string to be added to the sequence
                else{
                    
                    temp.clear(); // Re-allocate space for string
                    for(int index = 0; index < num_of_operations; ++index){
                        temp.insert(0,"2");
                    }

                    end_position+=num_of_operations;
                    operated_portions.push_back(temp);
                    
                }
               
            }
            // Operation to be done only on the reference therefore extract part of sequence based on num_of_operations
            else if(value_of_operation == 'I'){
                operated_portions.push_back(sequence.substr(current_position_in_read,num_of_operations));
                current_position_in_read = current_position_in_read + operations.at(i);
                

            }
       
    }
    sequence.clear();
    // Add the operated portions of the sequence together to make the new sequence
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
    Base base_at_position;
    size_t start_position = 0;
    size_t num_of_elements_in_map = 0;
    
    // if number of elements in map needed
    /*for(int chromosome = 19; chromosome < 20; ++chromosome){
        num_of_elements_in_map = _simulated_chromosomes.at(chromosome).size();
        std::cout << "chromosome: " << chromosome + 1 << " " << num_of_elements_in_map << std::endl;
    }*/
    num_of_elements_in_map = 0;
    
    for(int chromosome = 19; chromosome < 20; ++chromosome){
        //std::cout << "Start converting chromosome:" << chromosome + 1 << std::endl;
        start_position = _start_of_chromosome_reads.at(chromosome-19);
        num_of_elements_in_map = _simulated_chromosomes.at(chromosome).size();
        
        // get all read sequences and convert to binary based on bases at that position
        for (size_t i = start_position;; ++i) {
            if(_simulated_chromosomes.at(chromosome).count(i) > 0){

                for(size_t reads = 0; reads < _simulated_chromosomes.at(chromosome).at(i).size(); ++reads){
                    Read& read = _simulated_chromosomes.at(chromosome).at(i).at(reads);
                    for(size_t len = 0; len < read._sequence.length(); ++len){
                        
                        if (_ref_alt_chromosomes.at(chromosome).count(i + len) > 0) {
                            base_at_position = _ref_alt_chromosomes.at(chromosome).at(i+len).at(0);
                            if (read._sequence.at(len) == '2'){
                                read._binary_sequence+=two;
                            }
                            else if(read._sequence.at(len) == base_at_position._ref_base){
                                read._binary_sequence+=one;
                            }
                            else{
                                read._binary_sequence+=zero;
                            }

                            
                        } // end if
                        else {
                            
                            if (read._sequence.at(len) == '2'){
                                read._binary_sequence+=two;
                            }
                            else{
                                read._binary_sequence+=one;
                            }
                            
                            // Position is not in references bases, add new "assumed" base
                            _ref_alt_chromosomes[chromosome][i + len].push_back(Base(read._sequence.at(len), 'A', 0, 2, 2));
                        } // end else
                        
                        
                    }//for 4
                    
                }// for 3
                
                    num_of_elements_in_map--;


           }
           if(num_of_elements_in_map == 1){
               std::cout << "Finished chromosome:" << chromosome + 1 << std::endl;
               break;
           }
            
        }//for 2
    }//for 1
}


template <typename TP>
void DataConverter::process_each_line(const TP& line)
{
    std::string start_position;
    size_t start_position_value;
    std::string end_position;
    size_t end_position_value;
    
    // Find start position of line
    for(size_t i = 0; i < line.length(); ++i){
        if(line[i] != '-'){
            start_position_value = i;
            start_position = std::to_string(start_position_value);
            break;
        }
    }
    
    // Fine end position of line
    for(size_t i = 0; i < line.length(); ++i){
        if(line[line.length() - 1 - i] != '-'){
            end_position_value = line.length() - i - 1;
            end_position = std::to_string(end_position_value);
            break;
        }
    }
    
    _elements_per_line.push_back(end_position_value - start_position_value + 1);
    
    for(int i = 0; i < start_position.length(); ++i){
        _data.push_back(start_position[i]);
    }
    
    _data.push_back(' ');
    
    for(int i = 0; i < end_position.length(); ++i){
        _data.push_back(end_position[i]);
    }
    
    _data.push_back(' ');

    for (size_t i = start_position_value; i < end_position_value + 1; ++i) {
        if(line[i] == _ref_seq.at(i)) {
            _data.push_back('1');
        }
        else if(line[i] == '-') {
            _data.push_back('-');
        }
        else {
            _data.push_back('0');
        }
    }
    // Add in the newline character
    _data.push_back('\n');
}
    
// uncomment if require matrix format output
/*template <typename TP>
void DataConverter::process_each_line(const TP& line)
{
    for (size_t i = 0; i < line.length(); ++i) {
        
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
}*/
    
std::vector<size_t> DataConverter::convert_data_to_binary(std::vector<char> input)
{
    std::vector<size_t> output;
    for(size_t i = 0; i < input.size(); ++i){
        
        if(input.at(i) == _ref_seq.at(i))
            output.push_back(ONE);
        // i.e. if ZERO
        else
            output.push_back(ZERO);
    }
    
    return output;
    
}

void DataConverter::write_simulated_data_to_file(const char* filename)
{
    // Create the parameters for the output file
    std::string filename_1 = filename;
    filename_1 = filename + std::string("_") + std::to_string(_total_num_elements) + std::string(".txt");
    //std::cout << filename_1 << std::endl;
    io::mapped_file_params file_params(filename_1);
    
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
    
void DataConverter::store_haplotype_answers()
{
    Base base_at_position;
    size_t start_position = 0;
    size_t num_of_elements_in_map = 0;
    std::vector<const char*> filenames = {"answers/chr_1.txt", "answers/chr_2.txt", "answers/chr_3.txt", "answers/chr_4.txt", "answers/chr_5.txt", "answers/chr_6.txt", "answers/chr_7.txt", "answers/chr_8.txt", "answers/chr_9.txt", "answers/chr_10.txt", "answers/chr_11.txt", "answers/chr_12.txt", "answers/chr_13.txt", "answers/chr_14.txt", "answers/chr_15.txt", "answers/chr_16.txt", "answers/chr_17.txt", "answers/chr_18.txt", "answers/chr_19.txt", "answers/chr_20.txt", "answers/chr_21.txt", "answers/chr_22.txt"};

    // store only chromosome 20, change chromosome to start at the beginning of the range required
    for(int chromosome = 19; chromosome < 20; ++chromosome){
        std::cout << "Start storing haplotype chromosome:" << chromosome + 1 << std::endl;
        start_position = _start_of_chromosome_reads.at(chromosome-19);
        num_of_elements_in_map = _simulated_chromosomes.at(chromosome).size();
        std::ofstream myfile;
        myfile.open(filenames.at(chromosome));
        
        for (size_t i = start_position;; ++i) {
            if(_simulated_chromosomes.at(chromosome).count(i) > 0){
                //std::cout << chromosome << ": " << i << std::endl;
                for(size_t reads = 0; reads < _simulated_chromosomes.at(chromosome).at(i).size(); ++reads){
                    Read& read = _simulated_chromosomes.at(chromosome).at(i).at(reads);
                    for(size_t len = 0; len < read._sequence.length(); ++len){
                        
                        if (_ref_alt_chromosomes.at(chromosome).count(i + len) > 0) {
                            base_at_position = _ref_alt_chromosomes.at(chromosome).at(i+len).at(0);
                            myfile << i + len << " " << base_at_position._real << " " << base_at_position._haplotype_one << " " << base_at_position._haplotype_two << std::endl;
                            
                            
                        } // end if
                        else {
                           
                        } // end else
                        
                        
                    }//for 4
                    
                }// for 3
                num_of_elements_in_map--;
                
            }//if
            if(num_of_elements_in_map == 1){
                //std::cout << "Finished chromosome:" << chromosome + 1 << std::endl;
                break;
            }
            
        }//for 2
        myfile.close();
    }//for 1
 
}
    
    
void DataConverter::write_dataset_to_file(const char* filename_test)
{
    size_t start_position = 0;              // To indicate lowest read in map
    size_t num_of_elements_in_map = 0;      // To keep track of how many remaining reads left to print
    size_t num_of_element_in_file = 0;
    std::vector<const char*> filenames = {"output_files/chr_1.txt", "output_files/chr_2.txt", "output_files/chr_3.txt", "output_files/chr_4.txt", "output_files/chr_5.txt", "output_files/chr_6.txt", "output_files/chr_7.txt", "output_files/chr_8.txt", "output_files/chr_9.txt", "output_files/chr_10.txt", "output_files/chr_11.txt", "output_files/chr_12.txt", "output_files/chr_13.txt", "output_files/chr_14.txt", "output_files/chr_15.txt", "output_files/chr_16.txt", "output_files/chr_17.txt", "output_files/chr_18.txt", "output_files/chr_19.txt", "output_files/chr_20.txt", "output_files/chr_21.txt", "output_files/chr_22.txt"};
    
    // store only chromosome 20, change chromosome to start at the beginning of the range required
    for(int chromosome = 19; chromosome < 20; ++chromosome){
        //std::cout << "Print chromosome: " << chromosome + 1 << std::endl;
        
        start_position = _start_of_chromosome_reads.at(chromosome - 19);
        num_of_elements_in_map = _simulated_chromosomes.at(chromosome).size();
        std::ofstream myfile;
        myfile.open (filenames.at(chromosome));
        
        for (size_t i = start_position;; ++i) {
                if (_simulated_chromosomes.at(chromosome).find(i) != _simulated_chromosomes.at(chromosome).end()) {
                    for(int vec = 0; vec < _simulated_chromosomes.at(chromosome).at(i).size(); ++vec){
                        //std::cout << i << " " << _chr1_simulated_data.at(i).at(vec)._end_position << " " << _chr1_simulated_data.at(i).at(vec)._binary_sequence << std::endl;
                        num_of_element_in_file+=(_simulated_chromosomes.at(chromosome).at(i).at(vec)._binary_sequence.length());
                        myfile << i - 1 << " " << _simulated_chromosomes.at(chromosome).at(i).at(vec)._end_position - 1 << " " << _simulated_chromosomes.at(chromosome).at(i).at(vec)._binary_sequence << std::endl;
                        
                    }
                    num_of_elements_in_map--;
                    
                }
            if(num_of_elements_in_map == 1){
                //std::cout << "End Chromosome: " << chromosome + 1 << std::endl;
                break;
            }
            
        }
        myfile << num_of_element_in_file << std::endl;
        num_of_element_in_file = 0;
        myfile.close();
    }
    
  
}
    
}// End namespcea haplo
