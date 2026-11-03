#include <cstdint>
#include <gmp.h>
#include <iostream>

#include "../../include/sos/algebra/mvpolynomial.hpp"

namespace sos::algebra {

    mvPolynomialT::mvPolynomialT() : deg(0), nvars(0), homogeneous(true) {
        std::cin >> nvars;

        coefficient c;
        while (std::cin >> c) {
            mpq_canonicalize(c.get_mpq_t());

            exponents e(nvars);
            for (uint64_t i = 0; i < nvars ; ++i) {
                std::cin >> e[i];
            }
            monomialT mon(e);
            addTerm(mon, c);
        }
        updateData();
    }

    void
    mvPolynomialT::addTerm(const monomialT& mon, const coefficient& c) {
        if (c == 0) { return; }

        auto it = terms.find(mon);
        if (it != terms.end()) {
            it->second += c;

            if (it->second == 0) {
                terms.erase(it);
            }
        }
        else {
            terms[mon] = c;
        }
    }

    void
    mvPolynomialT::updateData() {
        if (terms.empty()) {
            deg = 0;
            homogeneous = true;
            return;
        }

       auto it = terms.begin();

       deg = it->first.get_degree();
       homogeneous = true;

       ++it;
       for ( ; it != terms.end(); ++it) {
            degree d = it->first.get_degree();

            if (d != deg) {
                homogeneous = false;
            }
            if (d > deg) {
                deg = d;
            }
       }
    }

    degree mvPolynomialT::get_degree() const {
        return(deg);
    }

    uint64_t mvPolynomialT::get_nvars() const {
        return(nvars);
    }

    bool mvPolynomialT::is_homogeneous() const {
        return(homogeneous);
    }

    const monomialToCoeff& mvPolynomialT::get_terms() const {
        return(terms);
    }

} // end namespace sos::algebra
