#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include <gmpxx.h>

#include "../../include/sos/algebra/types.hpp"
#include "../../include/sos/algebra/monomial.hpp"
#include "../../include/sos/algebra/mvpolynomial.hpp"
#include "../../include/sos/io/display.hpp"

namespace sos::algebra {

    namespace {
        struct termBlock {
            std::string top;
            std::string mid;
            std::string bot;
        };

        std::string
        center(const std::string& s, size_t width) {
            size_t left = (width - s.length()) / 2;
            size_t right = width - s.length() - left;
            return(std::string(left, ' ') + s + std::string(right, ' '));
        }

        termBlock
        makeBlock(const std::string& sign, const coefficient& abs_c, const monomialT& mon) {
            bool constant = (mon.get_degree() == 0);

            std::ostringstream oss;
            if (!constant) { oss << mon; }
            std::string mon_str = oss.str();

            std::string num = abs_c.get_num().get_str();
            std::string den = abs_c.get_den().get_str();

            std::string coef, above, below;
            if (abs_c.get_den() == 1) {
                if (abs_c != 1 || constant) { coef = num; }
                above = std::string(coef.length(), ' ');
                below = above;
            }
            else {
                size_t width = std::max(num.length(), den.length()) + 2;
                coef  = std::string(width, '-');
                above = center(num, width);
                below = center(den, width);
            }

            std::string rest = mon_str;
            if (!coef.empty() && !rest.empty()) { rest = " " + rest; }

            std::string indent(sign.length(), ' ');
            std::string gap(rest.length(), ' ');

            return { indent + above + gap, sign + coef + rest, indent + below + gap };
        }

    } // end namespace

    std::ostream& operator<<(std::ostream& os, const monomialT& mon) {
        if (mon.get_degree() == 0) {
            os << "1";
            return os;
        }

        bool first = true;
        uint64_t index = 0;
        for (exponent e : mon.get_exponents()) {
            if (e != 0) {
                if (!first) { os << " "; }
                os << "x_" << index << "^" << e;
                first = false;
            }
            ++index;
        }
        return os;
    }

    std::ostream& operator<<(std::ostream& os, const mvPolynomialT& poly) {
        const auto& terms = poly.get_terms();

        if (terms.empty()) {
            os << "0";
            return os;
        }

        std::string top, mid, bot;
        bool first = true;
        bool has_fraction = false;

        for (const auto& t : terms) {
            const coefficient& c = t.second;
            bool neg = (c < 0);

            std::string sign;
            if (first) { sign = neg ? "- " : ""; }
            else       { sign = neg ? " - " : " + "; }

            coefficient abs_c = abs(c);
            termBlock b = makeBlock(sign, abs_c, t.first);

            top += b.top;
            mid += b.mid;
            bot += b.bot;

            if (abs_c.get_den() != 1) { has_fraction = true; }
            first = false;
        }

        if (has_fraction) {
            os << top << "\n" << mid << "\n" << bot;
        }
        else {
            os << mid;
        }
        return os;
    }

    void monomialT::display() const {
        std::cout << *this;
    }

    void mvPolynomialT::display() const {

        std::cout << "[Input multivariate polynomial]" << std::endl;
        std::cout << *this << std::endl;

        std::cout << '\n';

        std::cout << "[Total degree] " << this->get_degree() << std::endl;
        std::cout << "[Number of variables] " << this->get_nvars() << std::endl;
        std::cout << "[Homogeneous] " << (this->is_homogeneous() ? "True" : "False" ) << std::endl;

    }

} // end namespace sos::algebra
