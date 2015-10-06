// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for parahaplo search operations class
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_SEARCH_OPERATIONS_CPU_HPP
#define PARHAPLO_SEARCH_OPERATIONS_CPU_HPP

namespace haplo {
            
// ----------------------------------------------------------------------------------------------------------
/// @class 
// ----------------------------------------------------------------------------------------------------------

// Specialization for cpu selection operator
template <>
class NodeSelector<devices::cpu> {
public:
    
private:
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Constructor
    // ------------------------------------------------------------------------------------------------------
    NodeSelector() noexcept {}
       

    // ------------------------------------------------------------------------------------------------------
    /// @brief      Sorts some of the nodes so that the ones most likely to have the greatest impact on the
    ///             solution are evaluated first
    // ------------------------------------------------------------------------------------------------------
    void operator()();

};

template <>
void NodeSelector<cpu::devices>::operator()()
{
}

}               // End namespace haplo
#endif          // PARHAPLO_SEARCH_OPERATIONS_CPU_HPP

