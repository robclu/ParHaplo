// ----------------------------------------------------------------------------------------------------------
/// @file   data_coverter.hpp
/// @brief  Header file for the data coverting class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_DATA_CONVERTER_HPP
#define PARAHAPLO_DATA_CONVERTER_HPP

#include "small_containers.hpp"

#include <array>
#include <string>
#include <vector>
#include <unordered_map>

#ifndef ONE
    #define     ZERO    0x00
    #define     ONE     0x01
    #define     TWO     0x02
#endif


namespace haplo {
    
    // ------------------------------------------------------------------------------------------------------
    /// @struct     Base
    /// @brief      Stores the relevant data per base for a particular position of the reference sequence in the ground truth file
    /// @param[in]  _ref_base          The reference base (occurrences represent by a 1)
    /// @param[in]  _alt_base          The alternate base (occurrences represent by a 0)
    /// @param[in]  _real              Boolean indicating whether a base exists (1) or is assumed (0) for the dataset
    /// @param[in]  _haplotype_one     Occurrence of base along one chromosome
    /// @param[in]  _haplotype_two     Occurrence of base along other chromosome
    // ------------------------------------------------------------------------------------------------------
    struct Base {
        Base()
        {}
        Base(char ref, char alt, bool r) : _ref_base(ref), _alt_base(alt), _real(r)
        {}
        Base(char ref, char alt, bool r, size_t hap_one, size_t hap_two) : _ref_base(ref), _alt_base(alt), _real(r), _haplotype_one(hap_one), _haplotype_two(hap_two)
        {}
        void print() const {
            std::cout << "ref: " << _ref_base << std::endl;
            std::cout << "alt: " << _alt_base << std::endl;
            std::cout << "real: " << _real << std::endl;
            std::cout << "hap 1: " << _haplotype_one << std::endl;
            std::cout << "hap 2: " << _haplotype_two << std::endl;
            
        }
        
        char     _ref_base;  // a, t, c, g
        char     _alt_base;  // a, t, c, g
        bool     _real;
        size_t   _haplotype_one;
        size_t   _haplotype_two;

    };
    
    // ------------------------------------------------------------------------------------------------------------------------
    /// @struct     Read
    /// @brief      Stores the relevant data per read for a particular dataset
    /// @param[in]  _end_position          The end position of the read
    /// @param[in]  _sequence              The sequenced bases for the read
    /// @param[in]  _binary_sequence       The converted binary sequence for the read (with respect to the reference sequence)
    // ------------------------------------------------------------------------------------------------------------------------
    struct Read {
        Read()
        {}
        Read(size_t end_position, std::string seq) : _end_position(end_position), _sequence(seq)
        {}
        void print() const {
            std::cout << "end: " << _end_position << std::endl;
            std::cout << "seq: " << _sequence << std::endl;
            std::cout << "bin: " << _binary_sequence << std::endl;
        }
        
        size_t          _end_position;
        std::string     _sequence; // a, t, c, g
        std::string     _binary_sequence; // 0, 1
        
    };

    using umap_read = std::unordered_map<size_t, std::vector<Read>>;
    using umap_base = std::unordered_map<size_t, std::vector<Base>>;
    using chromo_array_reads = std::array<umap_read, 22>;
    using chromo_array_bases = std::array<umap_base, 22>;

    
// ----------------------------------------------------------------------------------------------------------
/// @class      DataConverter
/// @class      Processes the reference sequence and dataset files and converts from ACTG to binary
// ----------------------------------------------------------------------------------------------------------
class DataConverter {
private:
    std::vector<char>       _data;                      //!< The converted data for smaller input files
    std::vector<size_t>     _elements_per_line;          //!< Number of elements in a single line of smaller input files
    size_t                  _rows;                      //!< Number of rows in the input file
    size_t                  _columns;                   //!< Number of columns in the input file
    size_t                  _chromosome;                //!< Chromosome number identifier
    std::vector<char>       _ref_seq;                   //!< Reference sequence for smaller input files
    std::vector<char>       _alt_seq;                   //!< Alternate sequence for smaller input files
    std::vector<size_t>     _base_a;                    //!< Counter for occurrences of base a, indexed by columns for smaller input files
    std::vector<size_t>     _base_c;                    //!< Counter for occurrences of base c, indexed by columns for smaller input files
    std::vector<size_t>     _base_t;                    //!< Counter for occurrences of base t, indexed by columns for smaller input files
    std::vector<size_t>     _base_g;                    //!< Counter for occurrences of base g, indexed by columns for smaller input files
   
    chromo_array_bases      _ref_alt_chromosomes;       //!< Array of maps to store all the bases for a specific ground truth file of references, indexed by chromosome number
    chromo_array_reads      _simulated_chromosomes;     //!< Array of maps to store all the reads for a specific dataset file, indexed by chromosome number

    std::vector<size_t>     _start_of_chromosome_reads; //!< Vector of starting positions for each chromosome, indexed by chromosome number


public:
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Class default constructor
    // ------------------------------------------------------------------------------------------------------
    DataConverter() {};

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts input data for bigger datasets (bam) and outputs it to a file
    /// @param[in]  data_file_1     Stores the reference sequence (ground truth file)
    /// @param[in]  data_file_2     Stores all the reads for a dataset (bam file)
    /// @param[in]  data_file_3     Stores the processed output data (.txt file)
    // ------------------------------------------------------------------------------------------------------
    DataConverter(const char* data_file_1, const char* data_file_2, const char* data_file_3);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts input data for smaller datasets and outputs it to a file
    /// @param[in]  data_file_1     Stores the smaller dataset (.txt file)
    /// @param[in]  data_file_3     Stores the processed output data (.txt file)
    // ------------------------------------------------------------------------------------------------------
    DataConverter(const char* data_file_1, const char* data_file_2);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Writes the converted data to a file for smaller input files
    /// @param[in]  data_file       Stores the processed output data
    // ------------------------------------------------------------------------------------------------------
    void write_simulated_data_to_file(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Writes the converted simulated data to a file for bigger input files
    /// @param[in]  data_file       Stores the processed output data
    // ------------------------------------------------------------------------------------------------------
    void write_dataset_to_file(const char* data_file);
    
    // DEBUGGING
    void print_simulated() const;
    
    //
    // DEBUGGING
    void print_dataset() const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Takes in BinaryArray elements and uses the reference sequence to convert back to ACTG
    /// @param[in]  input       Stores the binary elements
    /// @return     A vector of characters (ACTG)
    // ------------------------------------------------------------------------------------------------------
    template <size_t length>
    std::vector<char> convert_data_from_binary(BinaryArray<length, 2> input);
    
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Takes in ACTG elements and uses the reference sequence to convert to binary
    /// @param[in]  input       Stores the characters of the sequence
    /// @return     A vector of binary elements
    // ------------------------------------------------------------------------------------------------------
    std::vector<size_t> convert_data_to_binary(std::vector<char> input);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts a character to a byte for mapping
    /// @param[in]  input       Stores the character (ACTG)
    /// @return     The byte equivalent of the character
    // ------------------------------------------------------------------------------------------------------
    byte convert_char_to_byte(char input);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts a byte to a character for mapping
    /// @param[in]  input       Stores the byte value to be converted
    /// @return     The character equivalent of the byte
    // ------------------------------------------------------------------------------------------------------
    char convert_byte_to_char(byte input);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts the CIGAR value of a bam file in order to apply operations to the sequence of each read
    /// @param[in]  start_position       The start position of a read
    /// @param[in]  end_position         The end position of a read
    /// @param[in]  cigar_value          The cigar value of a read
    /// @param[in]  sequence             The sequence of a read
    // ------------------------------------------------------------------------------------------------------
    void process_cigar_value(size_t& start_position, size_t& end_position, std::string& cigar_value, std::string& sequence);

    
private:
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts the simulated input data into a binary equivalent
    /// @param[in]  data_file       Stores the simulated input data
    // ------------------------------------------------------------------------------------------------------
    void convert_simulated_data_to_binary(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes the simulated input and determines a reference sequence
    // ------------------------------------------------------------------------------------------------------
    void determine_simulated_ref_sequence();
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the occurrences of the 4 bases (ACTG) in each column of the simulated data
    /// @param[in]  line       Stores each line passed to be processed
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void find_base_occurrance(const TP& line);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Uses reference sequence to convert processed simulated data to binary
    /// @param[in]  line       Stores each line passed to be processed
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_each_line(const TP& line);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts the dataset input data into a binary equivalent using a given reference sequence
    /// @param[in]  data_file_1       Stores the reference sequence (ground truth)
    /// @param[in]  data_file_2       Stores the dataset (bam)
    // ------------------------------------------------------------------------------------------------------
    void convert_dataset_to_binary(const char* data_file_1, const char* data_file_2);

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes the ground truth file to be stored in a map of bases for each chromosome
    /// @param[in]  data_file       Stores the ground truth file
    // ------------------------------------------------------------------------------------------------------
    void process_ground_truth(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes the ground truth file to be stored in a map of bases for each chromosome
    /// @param[in]  data_file       Stores the ground truth file
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void determine_dataset_ref_sequence(const TP& token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Processes the dataset file and extracts the important informationt to be stored in a map of reads
    /// @param[in]  line       Stores each line passed to be processed
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_dataset(const TP& line);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts the CIGAR value of a bam file in order to apply operations to the sequence of each read
    /// @param[in]  start_position       The start position of a read
    /// @param[in]  end_position         The end position of a read
    /// @param[in]  cigar_value          The cigar value of a read
    // ------------------------------------------------------------------------------------------------------
    void store_read_data(size_t chromosome, size_t start_position, size_t end_position, std::string sequence);

    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Compares datasets reads to created reference bases in order to convert the sequence to  binary sequence
    // ------------------------------------------------------------------------------------------------------
    void process_each_read();
    
    
};

template <size_t length>
std::vector<char> DataConverter::convert_data_from_binary(BinaryArray<length, 2> input)
{
    std::vector<char> output;
    for(size_t i = 0; i < length; ++i){
        
        if(input.get(i) == ONE) {
            output.push_back(_ref_seq.at(i));
        }
        // i.e. if ZERO
        else {
            output.push_back(_alt_seq.at(i));
        }
    }
    
    return output;
}

}               // End namespace haplo
#endif          // PARAHAPLO_INPUT_CONVERTER_HPP
