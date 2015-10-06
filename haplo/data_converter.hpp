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

// !!!! Please do not modift Makefile -- ccompile with make converter_tests to compile only the converter
// tests

namespace haplo {
    
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
        
        char     _ref_base;  // a t c g
        char     _alt_base;
        bool     _real;
        size_t   _haplotype_one;
        size_t   _haplotype_two;

    };

    
// ----------------------------------------------------------------------------------------------------------
/// @class      InputConverter
/// @class      Converts the input from ATCG to binary
// --------------------------------------------------- ------------------------------------------------------
class DataConverter {
private:
    std::vector<char>       _data;                  //!< The converted data
    size_t                  _rows;                  //!< Number of rows in the input file
    size_t                  _columns;               //!< Number of columns in the input file
    size_t                  _chromosome;     //!< Chromosome number identifier
    // Rename these all as _ref_seq .. _base_a, unless _a_base makes more
    // sense, but i don't understand what they mean
    std::vector<char>       _refSeq;
    std::vector<char>       _altSeq;
    std::vector<size_t>     _aBase;
    std::vector<size_t>     _cBase;
    std::vector<size_t>     _tBase;
    std::vector<size_t>     _gBase;
   
    // Makybe make this an array of unordered maps -- this looks pretty ugly
    //  
    //  using umap = std::unorderd_map;
    //  using chromo_array = std::array<umap, 22>
    //  
    //  Then you can declare it above as
    //  
    //      chromo_array    _ref_alt_chromosomes -- for example
    std::unordered_map<size_t, Base> _chr1_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr2_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr3_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr4_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr5_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr6_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr7_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr8_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr9_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr10_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr11_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr12_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr13_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr14_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr15_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr16_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr17_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr18_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr19_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr20_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr21_ref_and_alt_seq;
    std::unordered_map<size_t, Base> _chr22_ref_and_alt_seq;

public:
    
    // ----------------------------------------------------------------------------------------------------
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
    
    void printMap() const;
    
    template <size_t length>
    std::vector<char> convert_data_from_binary(BinaryArray<length, 2> input);
    
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    byte convert_char_to_byte(char input);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    char convert_byte_to_char(byte input);

    
private:
    void storeBaseData(size_t chromosome, size_t position, char ref_base, char alt_base, bool real, size_t haplo_one, size_t haplo_two);

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
