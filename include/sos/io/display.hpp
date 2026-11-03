/**
 * @file monomial.hpp
 * @brief Defines the differnt views of code classes.
 * @author Augusto Guerra
 * @date 2026-11-02
 * @copyright Copyright (c) 2026.
 */

#ifndef DISPLAY_HPP
#define DISPLAY_HPP

#include <iostream>

#include "../algebra/monomial.hpp"
#include "../algebra/mvpolynomial.hpp"

namespace sos::algebra {

    std::ostream& operator<<(std::ostream& os, const monomialT& mon);

    std::ostream& operator<<(std::ostream& os, const mvPolynomialT& poly);

}

#endif // DISPLAY_HPP
