/**
 * @file monomial.hpp
 * @brief Define the class of monomials, a monomial
 *        is represented as a vector of exponents
 *        for the variables.
 * @author Augusto Guerra
 * @date 2026-11-02
 * @copyright Copyright (c) 2026.
 */

#ifndef MONOMIAL_HPP
#define MONOMIAL_HPP

#include "types.hpp"

namespace sos::algebra {

    class monomialT {
            exponents exps;
            degree  deg;

        public:
            monomialT(exponents& exps);

            degree get_degree() const;
            const exponents& get_exponents() const;

            bool operator< (const monomialT& other) const;

            void display() const;
    };

}

#endif
