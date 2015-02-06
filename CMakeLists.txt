cmake_minimum_required (VERSION 2.6)
project (geodesy CXX)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../ocaml-cmake/cmake")

option (BUILD_TESTS "Build tests" OFF)

find_package (OCaml REQUIRED)
include (UseOCaml)

include_directories (
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_OCaml_STD_LIBRARY_PATH}
  )

add_library (mlgeodesy ocaml/ml_geodesy.cpp)

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

if (BUILD_TESTS)
  enable_testing()

  find_ocaml_package(oUnit)

  # TODO -cclib "-Wl,-no_compact_unwind"
  add_ocaml_executable (tests_geodesy NATIVE
    SOURCES tests/test_geodesy
    PACKAGES oUnit
    LIBRARIES geodesy
    C_LIBRARIES mlgeodesy)

  add_test (geodesy tests_geodesy)

endif (BUILD_TESTS)
