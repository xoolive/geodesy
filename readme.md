A pragmatic and efficient geodesy toolkit
=========================================

The current project exposes interfaces in different languages to a pragmatic and efficient geodesy toolkit. Both spherical and WGS84 models are supported, though there are more functions implemented for the spherical model.

# C++

The source code is common to all languages and implemented with C++ templates for simple and double precision floating point values. When possible, SSE2 vectorisation is provided. (AVX expected) Interfaces are then provided for Python (through cython) and for OCaml.

The `src` folder contains the C++ implementation. Simply include resp. `spherical_geodesy.h` or `wgs84_geodesy.h` according to your needs. `*.cpp` files are just provided for aesthetical reasons, if you want to compile the library in advance.

# Python

The `geodesy` folder contains the Python interface and exposes both `geodesy.sphere` and `geodesy.wgs84` modules. You can install the module with the `setup.py` file provided. `numpy` is necessary to compile the library.

Recommended installation procedure:
- `virtualenv` (optional, but a good idea nonetheless):
```
 virtualenv geo_test
 source geo_test/bin/activate
```
- installation:
```
 pip install cython numpy
 pip install git+git://github.com/xolive/geodesy.git
```
You can then read the manual embedded in `geodesy.sphere` and `geodesy.wgs84`

# OCaml

The `ocaml` folder contains the OCaml interface. Refer to `geodesy.mli` for the interface.

Recommended installation procedure using `ocamlfind` and `opam` (from version 1.2), provided your `PATH` yields access to `cmake`:
```
 opam pin add geodesy git://github.com/xoolive/geodesy.git
```
You can then compile your OCaml source file with `ocamlfind`, e.g.:
```
 ocamlfind ocamlc -linkpkg -package oUnit,geodesy tests/test_geodesy.ml
```

The `CMakeLists.txt` also offers a decent installation procedure provided you want to embed the library in bigger projects. You will need to clone the [ocaml-cmake](https://github.com/ocaml-cmake/ocaml-cmake) project first and set the `CMAKE_MODULE_PATH` accordingly.

You can also use the `make install` target for a project installation that you may later import in a different project. (`find_package (geodesy REQUIRED)` in your new `CMakeLists.txt`)


