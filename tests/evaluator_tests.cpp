// ----------------------------------------------------------------------------------------------------------
/// @file   evaluator_tests.cpp
/// @brief  Test suite for parahaplo evaluator class
// ----------------------------------------------------------------------------------------------------------

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
    #define BOOST_TEST_MODULE EvaluatorTests
#endif
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../haplo/evaluator.hpp"

static constexpr char* reference = "eval_files/eval_input_zero_ref.txt";
static constexpr char* solution =  "eval_files/eval_input_zero_sol.txt";

BOOST_AUTO_TEST_SUITE( EvaluatorSuite )
    
BOOST_AUTO_TEST_CASE( canCorrectlyEvaluateSolution )
{
    haplo::Evaluator evaluator;
    BOOST_CHECK( evaluator(reference, solution) == 1 );
}

BOOST_AUTO_TEST_SUITE_END()
