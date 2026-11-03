/**
 * @file mvpolynomial.hpp
 * @brief Define the class of multivariate
 *        polynomials, they are mappings from
 *        monomials to rational coefficients.
 * @author Augusto Guerra
 * @date 2026-11-02
 * @copyright Copyright (c) 2026.
 */

#ifndef MVPOLYNOMIAL_HPP
#define MVPOLYNOMIAL_HPP

#include <map>

#include <cstdint>

#include "types.hpp"
#include "monomial.hpp"


namespace sos::algebra {

    using monomialToCoeff = std::map<monomialT, coefficient>;

    class mvPolynomialT {
            monomialToCoeff terms;
            degree deg;

            uint64_t nvars;
            bool homogeneous;

            void updateData();

        public:
            mvPolynomialT();
            void addTerm(const monomialT& mon, const coefficient& coef);

            degree get_degree() const;
            uint64_t get_nvars() const;
            const monomialToCoeff& get_terms() const;
            bool is_homogeneous() const;

            void display() const;
    };

}

#endif // MVPOLYNOMIAL_HPP
