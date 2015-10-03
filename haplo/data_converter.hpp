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
    void process_line(size_t line, TP token_pointer); 
};

}               // End namespcea haplo
#endif          // PARAHAPLO_INPUT_CONVERTER_HPP
