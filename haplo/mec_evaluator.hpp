// ----------------------------------------------------------------------------------------------------------
/// @file   mec_evaluator.hpp
/// @brief  Header file for the evaluator class which determines the MEC score of a solution
// ----------------------------------------------------------------------------------------------------------

#ifndef PARAHAPLO_MEC_EVALUATOR_HPP
#define PARAHAPLO_MEC_EVALUATOR_HPP

#include "small_containers.hpp"
#include <array>

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class      Evaluator
/// @brief      Evaluates the MEC score for a haplotype
// ----------------------------------------------------------------------------------------------------------
class MecEvaluator {
public:
    // --------------------------------------------- ALIAS'S ------------------------------------------------
    using binary_vector     = BinaryVector<1>;
    using haplo_container   = std::array<binary_vector, 2>;
    // ------------------------------------------------------------------------------------------------------
public:
    // ------------------------------------------------------------------------------------------------------
    /// @brief  Constructor -- takes in the 
    // ------------------------------------------------------------------------------------------------------
    Evaluator() {};
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Gets the reconstruction rate
    /// @param[in]  ref_haplotype   The file for the reference haplotype 
    /// @param[in]  sol_haplotype   The file for the solution halotype
    // ------------------------------------------------------------------------------------------------------
    double operator()(const char* ref_haplotype, const char* sol_haplotype) const;
private:
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Converts an input file into a haplotype
    /// @param      input_file  The dfile to get the haplotypes from
    /// @param      haplo_arary The array to put the haplotypes into
    // ------------------------------------------------------------------------------------------------------
    void get_haplotypes(const char* inut_file, haplo_container& haplotypes) const;
    
    // ------------------------------------------------------------------------------------------------------
    /// @brief      Determines the reconstruction rate
    /// @param[in]  ref_haplos  The reference haplotypes
    /// @param[in]  sol_haplos  The solution haplotypes
    // ------------------------------------------------------------------------------------------------------
    double recon_rate(const haplo_container& ref_haplos, const haplo_container& sol_haplos) const;
};

}               // End namespace haplo
#endif          // PARAHAPLO_EVALUATOR_HPP
