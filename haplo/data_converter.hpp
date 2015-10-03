// ----------------------------------------------------------------------------------------------------------
/// @file   data_coverter.hpp
/// @brief  Header file for the data coverting class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_DATA_CONVERTER_HPP
#define PARAHAPLO_DATA_CONVERTER_HPP

#include <vector>
#include <string>

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      InputConverter
/// @class      Converts the input from ATCG to binary
// --------------------------------------------------- ------------------------------------------------------
class DataConverter {
private:
    std::vector<char>       _data;      //!< The converted data
    size_t                  _rows;
    size_t                  _columns;
    std::vector<char>       _refSeq;
    std::vector<char>       _altSeq;
    std::vector<int>        _aBase;
    std::vector<int>        _cBase;
    std::vector<int>        _tBase;
    std::vector<int>        _gBase;
public:
    DataConverter(const char* data_file);
    
    // Public interface should have a save function to save to file

    // DEBUGGING
    void print() const;
    
private:
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    void convert_data(const char* data_file);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void find_base_occurrance(size_t line_number, TP token_pointer);
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    void determine_ref_sequence();
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    template <typename TP>
    void process_line(size_t line, TP token_pointer);
    
    
    // ------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------
    void write_data_to_file();
};

}               // End namespcea haplo
#endif          // PARAHAPLO_INPUT_CONVERTER_HPP
