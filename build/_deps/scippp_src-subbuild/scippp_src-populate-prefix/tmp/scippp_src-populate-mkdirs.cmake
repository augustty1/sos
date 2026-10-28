# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-src")
  file(MAKE_DIRECTORY "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-src")
endif()
file(MAKE_DIRECTORY
  "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-build"
  "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix"
  "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix/tmp"
  "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix/src/scippp_src-populate-stamp"
  "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix/src"
  "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix/src/scippp_src-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix/src/scippp_src-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/augusto/desktop/research/project/sos/build/_deps/scippp_src-subbuild/scippp_src-populate-prefix/src/scippp_src-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
