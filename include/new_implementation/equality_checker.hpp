// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo equality checker -- determines is rows or columns are equivalent
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_CPP_EQUALITY_CHECKER_HPP
#define PARAHAPLO_CPP_EQUALITY_CHECKER_HPP

namespace haplo {
namespace check {
    
static constexpr bool ROWS      = false;
static constexpr bool COLUMNS   = true;

}       // End namespace check

// ----------------------------------------------------------------------------------------------------------
/// @struct     EqualityChecker 
/// @brief      Functor to check if rows or columns are equivalent, where either rows or columns are checked
///             based on the context
/// @tparam     Type    If the check must be for rows or columns 
// ----------------------------------------------------------------------------------------------------------
template <bool Type>
struct EqualityChecker;

// Specialize for row checking 
template <> 
struct EqualityChecker<check::ROWS> {

// The type of the checker
static constexpr bool type = check::ROWS;

// ----------------------------------------------------------------------------------------------------------
/// @brief      Checks if two rows are equivalent (assumes that data is stored row-major)
/// @param[in]  data            The data where the rows are
/// @param[in]  row_one_idx     A pointer to the start of the first row to compare
/// @param[in]  row_two_idx     A pointer to the start of the second row to compare
/// @param[in]  row_length      The number of elements in the rows
/// @param[in]  stride          The number of elements in data between successive column elements
/// @tparam     DataType        The type of data container (must have a get and set opertator)
/// @return     If the two rows are quivalent
// ---------------------------------------------------------------------------------------------------------- 
template <typename DataType>
inline bool operator()(const  DataType& data    ,  
                       size_t row_one_idx       ,
                       size_t row_two_idx       , 
                       size_t row_length        ,
                       size_t stride            )
{
    bool are_equal = true;                          // Assume equal
    for (size_t col_idx = 0; col_idx < row_length && are_equal; ++col_idx) 
        are_equal = (data.get(row_one_idx * stride + col_idx) == data.get(row_two_idx * stride + col_idx))
                  ? true                                                            
                  : false;
    return are_equal;
}

};

// Specialize for column checking
template <>
struct EqualityChecker<check::COLUMNS> {

// The type of the checker
static constexpr bool type = check::COLUMNS;

// ----------------------------------------------------------------------------------------------------------
/// @brief      Checks if two columns are equivalent
/// @param[in]  data            The data where the columns are stored
/// @param[in]  col_one_idx     A pointer to the start of the first column to compare
/// @param[in]  col_two_idx     A pointer to the start of the second column to compare
/// @param[in]  col_length      The number of elements in the columns
/// @param[in]  stride          The number of elements in data between successive column elements
/// @tparam     DataType        The type of data container (must have a get funciton for a data element)
/// @return     If the two columns are quivalent
// ----------------------------------------------------------------------------------------------------------
template <typename DataType>
inline bool operator()(const  DataType& data    , 
                       size_t col_one_idx       ,
                       size_t col_two_idx       ,
                       size_t col_length        ,       
                       size_t stride            )
{
    bool are_equal = true;                          // Assume equal
    for (size_t row_idx = 0; row_idx < col_length && are_equal; ++row_idx) 
        are_equal = (data.get(col_one_idx + row_idx * stride) == data.get(col_two_idx + row_idx * stride)) 
                  ? true 
                  : false;
    return are_equal;
}

};

}           // End namespace haplo

#endif      // PARAHAPLO_CPP_EQUALITY_CHECKER_HPP
