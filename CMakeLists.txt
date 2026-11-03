cmake_minimum_required(VERSION 3.12)
project(sos_polynomial CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

set(LIBRARY_SOURCES
    src/algebra/monomial.cpp
    src/algebra/mvpolynomial.cpp
    src/io/display.cpp
)

add_executable(sos app/main.cpp ${LIBRARY_SOURCES})

target_link_libraries(sos PRIVATE gmpxx gmp)
