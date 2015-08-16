// ----------------------------------------------------------------------------------------------------------
/// @file   Header file for Block expressions
// ----------------------------------------------------------------------------------------------------------

#ifndef PARHAPLO_CPP_BLOCK_EXPRESSION_HPP
#define PARHAPLO_CPP_BLOCK_EXPRESSION_HPP

#include <vector>

namespace haplo {
    
// ----------------------------------------------------------------------------------------------------------
/// @class  BlockExpression
/// @brief  Base class from which all other Block classes derive. However, this uses static polymorphism 
///         rather than dynamic polymorphism so that the overhead of virtual function calls is removed, 
///         and additionally provides very clear expressions on Blocks when used in the interface.           \n
///                                                                                                          \n
///         The CRTP (Curiously Recurring Template Pattern) is used, see
///         https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern, additionally see:            \n
///                                                                                                          \n
///         Expression Templates:                                                                            \n
///         Wikipedia: https://en.wikipedia.org/wiki/Expression_templates                                    \n
///         Dr. Dobbs: http://www.drdobbs.com/c-expression-templates/184401627
/// @tparam T           The type of data used by the expressions 
/// @tapram Expression  The expression (Block, BlockSlice ...)
// ----------------------------------------------------------------------------------------------------------
template <typename T, typename Expression>
class BlockExpression {
public:
    // -------------------------------------- Typedefs ------------------------------------------------------
    using container_type    = std::vector<T>;
    using size_type         = typename container_type::size_type;
    using value_type        = typename container_type::value_type;           
    using reference         = typename container_type::reference;
    // ------------------------------------------------------------------------------------------------------
    
    // ------------------------------------------------------------------------------------------------------
    //! @brief     Returns the size of the Expression
    //! @return    The size of the Expression
    // ------------------------------------------------------------------------------------------------------
    size_type size() const { return static_cast<Expression const&>(*this).size(); }

    // ------------------------------------------------------------------------------------------------------
    //! @brief     Gets and element from the Expression data.
    //! @param[in] i   The element in the Expression which must be fetched.
    //! @return    The value of the element at position i of the Expression data.
    // ------------------------------------------------------------------------------------------------------
    value_type operator[](size_type i) const { return static_cast<Expression const&>(*this)[i]; }

    // ------------------------------------------------------------------------------------------------------
    //! @brief     Gets a reference to the Expression.
    //! @return    A reference to the Expression.
    // ------------------------------------------------------------------------------------------------------
    operator Expression&() { return static_cast<Expression&>(*this); }

    // ------------------------------------------------------------------------------------------------------
    //! @brief     Gets a constant reference to the Expression
    //! @return    A constant reference to the Expression
    // ------------------------------------------------------------------------------------------------------
    operator Expression const&() const   { return static_cast<const  Expression&>(*this); }
};

}       // End namespace haplo

#endif  // PARAHAPLO_CPP_BLOCK_EXPRESSION_HPP
