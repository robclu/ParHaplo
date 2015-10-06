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

#ifndef ONE
    #define     ZERO    0x00
    #define     ONE     0x01
    #define     TWO     0x02
#endif

// !!!! Please do not modift Makefile -- ccompile with make converter_tests to compile only the converter
// tests

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      InputConverter
/// @class      Converts the input from ATCG to binary
// --------------------------------------------------- ------------------------------------------------------
class DataConverter {
private:
    std::vector<char>       _data;                  //!< The converted data
    size_t                  _rows;                  //!< Number of rows in the input file
    size_t                  _columns;               //!< Number of columns in the input file
    size_t                  _no_of_chromosomes;     //!< Chromosome number identifier
    std::string             _chromosomes;           //!< Chromosome letter identifier

    // Rename these all as _ref_seq .. _base_a, unless _a_base makes more
    // sense, but i don't understand what they mean
    std::vector<char>       _refSeq;
    std::vector<char>       _altSeq;
    std::vector<size_t>     _aBase;
    std::vector<size_t>     _cBase;
    std::vector<size_t>     _tBase;
    std::vector<size_t>     _gBase;
    
    //BinaryVector<2>         _chr1_ref_seq;
    std::vector<char>           _chr1_ref_seq;
    BinaryVector<2>         _chr2_ref_seq;
    BinaryVector<2>         _chr3_ref_seq;
    BinaryVector<2>         _chr4_ref_seq;
    BinaryVector<2>         _chr5_ref_seq;
    BinaryVector<2>         _chr6_ref_seq;
    BinaryVector<2>         _chr7_ref_seq;
    BinaryVector<2>       _chr8_ref_seq;
    /*std::vector<char>       _chr9_ref_seq;
    std::vector<char>       _chr10_ref_seq;
    std::vector<char>       _chr11_ref_seq;
    std::vector<char>       _chr12_ref_seq;
    std::vector<char>       _chr13_ref_seq;
    std::vector<char>       _chr14_ref_seq;
    std::vector<char>       _chr15_ref_seq;
    std::vector<char>       _chr16_ref_seq;
    std::vector<char>       _chr17_ref_seq;
    std::vector<char>       _chr18_ref_seq;
    std::vector<char>       _chr19_ref_seq;
    std::vector<char>       _chr20_ref_seq;
    std::vector<char>       _chr21_ref_seq;
    std::vector<char>       _chr22_ref_seq;*/
    
    std::vector<char>             _chr1_alt_seq;
    
    BinaryVector<2>         _haplotype_one;
    BinaryVector<2>         _haplotype_two;
    
public:
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    DataConverter() {};
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    DataConverter(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Writes the converted data to a file
    /// @param[in]  filename    The name of the file to write the data to
    // ------------------------------------------------------------------------------------------------------
    void write_data_to_file(const char* filename ); 
    
    // DEBUGGING
    void print() const;
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <size_t length>
    std::vector<char> convert_data_from_binary(BinaryArray<length, 2> input);
    
private:
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    void convert_simulated_data_to_binary(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    void convert_dataset_to_binary(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void find_base_occurrance(const TP& token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    void determine_simulated_ref_sequence();
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void determine_dataset_ref_sequence(const TP& token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_line(const TP& token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    byte convert_char_to_byte(char input);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    char convert_byte_to_char(byte input);
    
    
    
};

template <size_t length>
std::vector<char> DataConverter::convert_data_from_binary(BinaryArray<length, 2> input)
{
    std::vector<char> output;
    for(size_t i = 0; i < length; ++i){
        
        if(input.get(i) == ONE) {
            output.push_back(_refSeq.at(i));
        }
        else {
            output.push_back(_altSeq.at(i));
        }
    }
    
    return output;
}

}               // End namespcea haplo
#endif          // PARAHAPLO_INPUT_CONVERTER_HPP
