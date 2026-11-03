/**
 * @file monomial.hpp
 * @brief Defines the elementary types to be used
 *        primaly within the algebra namespace.
 * @author Augusto Guerra
 * @date 2026-11-02
 * @copyright Copyright (c) 2026.
 */

#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>
#include <cstdint>
#include <gmpxx.h>

namespace sos::algebra {

    using coefficient = mpq_class;
    using exponent = std::uint32_t;
    using exponents = std::vector<exponent>;
    using degree = std::uint64_t;

}

#endif // TYPES_HPP
