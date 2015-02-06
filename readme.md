A pragmatic and efficient geodesy toolkit (C++, Python, OCaml)
==============================================================

The current project exposes interfaces in different languages to a pragmatic and efficient geodesy toolkit. Both spherical and WGS84 models are supported, though there are more functions implemented for the spherical model.

The source code is common to all languages and implemented with C++ templates for simple and double precision floating point values. When possible, SSE2 vectorisation is provided. (AVX expected) Interfaces are then provided for Python (through cython) and for OCaml.

- The `src` folder contains the C++ implementation. Simply include resp. `spherical_geodesy.h` or `wgs84_geodesy.h` according to your needs. `*.cpp` files are just provided for aesthetical reasons, if you want to compile the library in advance.

- The `geodesy` folder contains the Python interface and exposes both `geodesy.sphere` and `geodesy.wgs84` modules. You can install the module with the `setup.py` file provided. `numpy` is necessary to compile the library.

Recommended installation procedure :
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

- The `ocaml` folder contains the OCaml interface. Refer to `geodesy.mli` for the interface.

The `CMakeLists.txt` offers a decent installation procedure provided you want to embed the library in bigger projects. You will need to clone the [ocaml-cmake](https://github.com/ocaml-cmake/ocaml-cmake) project first (and possibly edit the `CMAKE_MODULE_PATH` variable in the `CMakeLists.txt` accordingly).

An automatic installation procedure via opam/ocamlfind is considered. Any contribution is welcome!


