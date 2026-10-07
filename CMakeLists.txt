cmake_minimum_required(VERSION 3.15)
project(mvpoly_sos LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- SCIP -------------------------------------------------------------
# Mesmo padrão que o próprio SCIPpp usa: pacote CONFIG do SCIP.
# Se o CMake não achar, passe -DSCIP_DIR=/usr/lib/cmake/scip (ou onde estiver o scip-config.cmake)
find_package(SCIP CONFIG REQUIRED)
if(TARGET libscip AND NOT TARGET SCIP::SCIP)
  add_library(SCIP::SCIP ALIAS libscip)
endif()

# --- SCIPpp ------------------------------------------------------------
# O install do SCIPpp NAO instala a biblioteca estatica (so os headers), entao
# compilamos os 5 .cpp dele junto com o projeto. De onde vem o codigo:
#   1) -DSCIPPP_SOURCE_DIR=/caminho/do/clone/SCIPpp   (se voce ja clonou)
#   2) senao, o CMake baixa do GitHub automaticamente (precisa de internet)
set(SCIPPP_SOURCE_DIR "" CACHE PATH "Clone local do SCIPpp (opcional)")

if(NOT SCIPPP_SOURCE_DIR)
  include(FetchContent)
  # SOURCE_SUBDIR inexistente: baixa o codigo mas NAO roda o CMakeLists do SCIPpp
  FetchContent_Declare(scippp_src
    GIT_REPOSITORY https://github.com/scipopt/SCIPpp.git
    GIT_TAG        1.4.0
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  nao_existe)
  FetchContent_MakeAvailable(scippp_src)
  set(SCIPPP_SOURCE_DIR ${scippp_src_SOURCE_DIR})
endif()

if(NOT EXISTS ${SCIPPP_SOURCE_DIR}/include/scippp/model.hpp)
  message(FATAL_ERROR "SCIPPP_SOURCE_DIR invalido: ${SCIPPP_SOURCE_DIR}")
endif()

set(scippp_version "local")
configure_file(${SCIPPP_SOURCE_DIR}/source/version.cpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/scippp_version.cpp)
file(GLOB SCIPPP_SRC ${SCIPPP_SOURCE_DIR}/source/*.cpp)

add_library(ScipPP STATIC ${SCIPPP_SRC} ${CMAKE_CURRENT_BINARY_DIR}/scippp_version.cpp)
target_include_directories(ScipPP PUBLIC ${SCIPPP_SOURCE_DIR}/include)
target_link_libraries(ScipPP PUBLIC SCIP::SCIP)

# --- SCIP-SDP (libscipsdp; headers ficam em /usr/include/scip/) --------
find_library(SCIPSDP_LIB NAMES scipsdp REQUIRED)

message(STATUS "SCIPpp src   : ${SCIPPP_SOURCE_DIR}")
message(STATUS "SCIP-SDP lib : ${SCIPSDP_LIB}")

# --- GMP (mpq_class, usado em mvpolynomial) ---------------------------
find_library(GMPXX_LIB NAMES gmpxx REQUIRED)
find_library(GMP_LIB   NAMES gmp   REQUIRED)

# --- Executavel ---------------------------------------------------------
# Layout esperado:  src/main.cpp  src/mvpolynomial.cpp  include/mvpolynomial.hpp
# (os .cpp incluem "../include/mvpolynomial.hpp", entao nao precisa de include dir)
add_executable(sos
  src/main.cpp
  src/mvpolynomial.cpp
)
target_link_libraries(sos PRIVATE
  ScipPP
  ${SCIPSDP_LIB}   # antes do SCIP, por causa da ordem de resolucao de simbolos
  SCIP::SCIP
  ${GMPXX_LIB}     # gmpxx antes de gmp
  ${GMP_LIB}
)
