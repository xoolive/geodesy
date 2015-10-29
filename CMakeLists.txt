cmake_minimum_required (VERSION 2.6)
project (geodesy CXX)

option (BUILD_TESTS "Build tests" OFF)

find_package (OCaml REQUIRED)
include (UseOCaml)

include_directories (
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_OCaml_STD_LIBRARY_PATH}
  )

add_library (mlgeodesy ocaml/ml_geodesy.cpp)

if (WIN32)
  set (OCAML_LINKER ${CMAKE_Flexlink_EXECUTABLE})
  execute_process(
    COMMAND         ${CMAKE_CXX_COMPILER} -print-file-name=libstdc++.a
    OUTPUT_VARIABLE LIBSTDCXX_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  execute_process(
    COMMAND         ${CMAKE_CXX_COMPILER} -print-file-name=libgcc_eh.a
    OUTPUT_VARIABLE LIBGCC_EH_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  set (OCAML_LINKER "${OCAML_LINKER} ${LIBSTDCXX_PATH} ${LIBGCC_EH_PATH}")
else (WIN32)
  set (OCAML_LINKER ${CMAKE_CXX_COMPILER})
endif (WIN32)

configure_file (
  ${CMAKE_SOURCE_DIR}/ocaml/META.in
  ${CMAKE_BINARY_DIR}/META)

add_ocaml_library (geodesy
  NATIVE
  SOURCES ocaml/geodesy
  C_LIBRARIES mlgeodesy)

add_ocaml_library (geodesy-byte
  BYTECODE
  SOURCES ocaml/geodesy
  C_LIBRARIES mlgeodesy)

if (BUILD_TESTS)

  enable_testing()

  find_ocaml_package(oUnit)

  add_ocaml_executable (tests_geodesy NATIVE
    SOURCES tests/test_geodesy
    PACKAGES oUnit
    LIBRARIES geodesy)

  add_test (geodesy tests_geodesy)

endif (BUILD_TESTS)

# If you don't want to install with ocamlfind:

install_ocaml_targets (geodesy
  DESTINATION lib/ocaml/site-lib/geodesy)

install_ocaml_targets (geodesy-byte
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
