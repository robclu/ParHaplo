// ----------------------------------------------------------------------------------------------------------
/// @file   evaluator_main.cpp
/// @brief  Main file for the parahaplo evaluator
// ----------------------------------------------------------------------------------------------------------

#include <iostream>

#include "../haplo/evaluator.hpp"

static constexpr char* reference     = "input_files/5268_reference.txt";
static constexpr char* solution      = "input_files/5268_solution.txt";

int main(int argc, char** argv)
{
    haplo::Evaluator evaluator;
 
    std::cout << "RR: " << evaluator(reference, solution) << "\n";
}

