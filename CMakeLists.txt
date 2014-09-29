cmake_minimum_required (VERSION 2.6)
project (deconflict CXX)

find_program (INTEL_COMPILER icc)
if (INTEL_COMPILER)
  include (CMakeForceCompiler)
  cmake_force_c_compiler (icc "Intel C Compiler")
  cmake_force_cxx_compiler (icpc "Intel C++ Compiler")
endif (INTEL_COMPILER)

# For the mac with openMP
# find_program (LLVM_GCC_4_2 llvm-gcc-4.2)
# if (LLVM_GCC_4_2)
#   include (CMakeForceCompiler)
#   cmake_force_c_compiler (llvm-gcc-4.2 "OpenMP compliant compiler in Mac OS")
#   cmake_force_c_compiler (llvm-g++-4.2 "OpenMP compliant compiler in Mac OS")
# endif (LLVM_GCC_4_2)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../ocaml-cmake/cmake")

find_package (OCaml REQUIRED)
include (UseOCaml)

include_directories (
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_OCaml_STD_LIBRARY_PATH}
  )

set_source_files_properties(src/ml_geodesy.c
  PROPERTIES LANGUAGE CXX)

add_library (mlgeodesy src/ml_geodesy.c)

add_ocaml_library (geodesy
  NATIVE
  SOURCES ocaml/geodesy
  C_LIBRARIES mlgeodesy)

install_ocaml_targets (geodesy
  DESTINATION lib/ocaml/site-lib/geodesy)

install_ocaml_interfaces (geodesy
  geodesy
  DESTINATION lib/ocaml/site-lib/geodesy)

install_ocaml_exports (geodesy
  DESTINATION lib/ocaml/site-lib/geodesy
  FILE geodesy-config-mli.cmake)

install (
  TARGETS mlgeodesy
  EXPORT geodesy-config
  DESTINATION lib/ocaml/site-lib/geodesy)

install (
  EXPORT geodesy-config
  DESTINATION lib/ocaml/site-lib/geodesy)

enable_testing()

# These functions are embedded into INTEL_COMPILER
# Those checks should be performed with other compilers though!
if (INTEL_COMPILER)
else (INTEL_COMPILER)
  set_source_files_properties(tests/test_sse2_math.c
    PROPERTIES LANGUAGE CXX)
  add_executable (tests_sse2_math tests/test_sse2_math.c)
  add_test (sse2_math tests_sse2_math)
endif (INTEL_COMPILER)

add_ocaml_executable (tests_geodesy NATIVE
  SOURCES tests/test_geodesy
  PACKAGES oUnit
  LIBRARIES geodesy
  C_LIBRARIES mlgeodesy)

add_test (geodesy tests_geodesy)

add_ocaml_executable (bench_geodesy NATIVE
  SOURCES tests/bench_geodesy
  LIBRARIES geodesy
  C_LIBRARIES mlgeodesy)



