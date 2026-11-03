#include <exception>
#include <iostream>

#include "../include/sos/algebra/mvpolynomial.hpp"

int main() {
    try {
        sos::algebra::mvPolynomialT poly;
        poly.display();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return(1);
    }
    return(0);
}
