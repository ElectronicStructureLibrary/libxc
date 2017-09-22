# Libxc

Libxc is a library of exchange-correlation functionals for
density-functional theory. The aim is to provide a portable, well
tested and reliable set of exchange and correlation functionals that
can be used by a variety of programs.

For more information, please check the manual at
http://www.tddft.org/programs/Libxc

## INSTALLATION

### Autotools

The recommended way to install the library is by using GNU Autotools.

To install the library, just use the standard procedure:
```
./configure --prefix=PATH/TO/LIBXC
make
make check
make install
```

If you're not using a stable release tarball, you'll first need to
generate ```configure``` with ```àutoreconf -i```.


### CMake

Support for CMake has also been recently contributed by Lori Burns.

The CMake file has the following caveats
* tested on Linux and Mac, static and shared lib, namespaced and non-namespaced headers, but really only to the extent that it works for Psi4
* all the fancy libtool options and Fortran interface _not_ tested
* test suite executed after build via `ctest`. But it has always totally passed or totally failed, which doesn't inspire confidence
* The generated `libxc_docs.txt` is large, and the generation step sometimes balks on it, leading to `xc_funcs.h` not found errors. Just execute again.

#### Building with CMake

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

#### Detecting with CMake

CMake builds install with `LibxcConfig.cmake`, `LibxcConfigVersion.cmake`, and `LibxcTargets.cmake` files suitable for use with CMake [`find_package()`](https://cmake.org/cmake/help/v3.2/command/find_package.html) in `CONFIG` mode.

* `find_package(Libxc)` - find any xc libraries and headers
* `find_package(Libxc 3.0.0 EXACT CONFIG REQUIRED COMPONENTS static)` - find Libxc exactly version 3.0.0 built with static libraries or die trying

See [cmake/LibxcConfig.cmake.in](cmake/LibxcConfig.cmake.in) for details of how to detect the Config file and what CMake variables and targets are exported to your project.

#### Use with CMake

After `find_package(Libxc ...)`,

* test if package found with `if(${Libxc_FOUND})` or `if(TARGET Libxc::xc)`
* link to library (establishes dependency), including header and definitions configuration with `target_link_libraries(mytarget Libxc::xc)`
* include header files using `target_include_directories(mytarget PRIVATE $<TARGET_PROPERTY:Libxc::xc,INTERFACE_INCLUDE_DIRECTORIES>)`
* compile target applying `-DUSING_Libxc` definition using `target_compile_definitions(mytarget PRIVATE $<TARGET_PROPERTY:Libxc::xc,INTERFACE_COMPILE_DEFINITIONS>)`
