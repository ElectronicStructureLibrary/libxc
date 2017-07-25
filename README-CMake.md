# libxc [![Build Status](https://travis-ci.org/loriab/libxc.svg?branch=master)](https://travis-ci.org/loriab/libxc)

Miguel A.L. Marques's Libxc (http://www.tddft.org/programs/Libxc) wrapped in CMake for Psi4 (https://github.com/psi4/psi4)

### History

This is the Libxc project (http://www.tddft.org/programs/Libxc) by
Prof. Miguel A.L. Marques of Martin-Luther-Universität Halle-Wittenberg.

Libxc is written in C. It has source and manual are available at the above
website and, as distributed, builds with `make`.

### This Repository

Libxc has been in the *ab initio* quantum chemistry package Psi4
(http://psicode.org/, https://github.com/psi4/psi4) since Novenber 2016. In Psi4,
it builds with `cmake`. This repository contains an unpacked tarball of Libxc 3.0.0
source, which is otherwise untouched, and the files below

* [cmake/](cmake) directory
* [CMakeLists.txt](CMakeLists.txt) top-level
* [testsuite/CMakeLists.txt](testsuite/CMakeLists.txt) tests
* [config.h.cmake.in](config.h.cmake.in) a dummy config file
* this README-CMake.md

Those files will be preserved in this repository (may move to psi4/libxc), but the
unpacked tarball may be changed to a download command of upstream Libxc. The stress
here is that the CMake build system is the only value added by this repository.

#### Caveats

* tested on Linux and Mac, static and shared lib, namespaced and non-namespaced headers, but really only to the extent that it works for Psi4
* all the fancy libtool options and Fortran interface _not_ tested
* test suite executed after build via `ctest`. But it has always totally passed or totally failed, which doesn't inspire confidence
* The generated `libxc_docs.txt` is large, and the generation step sometimes balks on it, leading to `xc_funcs.h` not found errors. Just execute again.

#### Version

This codebase was copied from upstream (above website) at 3.0.0.

#### Building

```bash
cmake -H. -Bobjdir
cd objdir && make
make install
```

The build is also responsive to

* static/shared toggle `BUILD_SHARED_LIBS`
* install location `CMAKE_INSTALL_PREFIX`
* namespacing of headers `NAMESPACE_INSTALL_INCLUDEDIR`
* of course, `CMAKE_C_COMPILER`, `BUILD_TESTING`, and `CMAKE_C_FLAGS`

See [CMakeLists.txt](CMakeLists.txt) for options details. All these build options should be passed as `cmake -DOPTION`.

#### Detecting

This project installs with `LibxcConfig.cmake`, `LibxcConfigVersion.cmake`, and `LibxcTargets.cmake` files suitable for use with CMake [`find_package()`](https://cmake.org/cmake/help/v3.2/command/find_package.html) in `CONFIG` mode.

* `find_package(Libxc)` - find any xc libraries and headers
* `find_package(Libxc 3.0.0 EXACT CONFIG REQUIRED COMPONENTS static)` - find Libxc exactly version 3.0.0 built with static libraries or die trying

See [cmake/LibxcConfig.cmake.in](cmake/LibxcConfig.cmake.in) for details of how to detect the Config file and what CMake variables and targets are exported to your project.

#### Using

After `find_package(Libxc ...)`,

* test if package found with `if(${Libxc_FOUND})` or `if(TARGET Libxc::xc)`
* link to library (establishes dependency), including header and definitions configuration with `target_link_libraries(mytarget Libxc::xc)`
* include header files using `target_include_directories(mytarget PRIVATE $<TARGET_PROPERTY:Libxc::xc,INTERFACE_INCLUDE_DIRECTORIES>)`
* compile target applying `-DUSING_Libxc` definition using `target_compile_definitions(mytarget PRIVATE $<TARGET_PROPERTY:Libxc::xc,INTERFACE_COMPILE_DEFINITIONS>)`
