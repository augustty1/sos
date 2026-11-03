#include "../../include/sos/algebra/monomial.hpp"

namespace sos::algebra {

    monomialT::monomialT(exponents& exps): exps(exps), deg(0) {
       for (exponent e : exps) {
            deg += e;
       } 
    }

    degree
    monomialT::get_degree() const {
        return(deg);
    }

    const exponents&
    monomialT::get_exponents() const {
        return(exps);
    }

    bool
    monomialT::operator< (const monomialT& other) const {
        return(exps < other.exps);
    }

} // end namespace sos::algebra
