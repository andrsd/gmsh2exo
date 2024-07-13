gmsh2exo
========

[![qa](https://github.com/andrsd/gmsh2exo/actions/workflows/qa.yml/badge.svg)](https://github.com/andrsd/gmsh2exo/actions/workflows/qa.yml)
[![build](https://github.com/andrsd/gmsh2exo/actions/workflows/build.yml/badge.svg)](https://github.com/andrsd/gmsh2exo/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/andrsd/gmsh2exo/graph/badge.svg?token=JHXYRN6N2X)](https://codecov.io/gh/andrsd/gmsh2exo)
[![License](http://img.shields.io/:license-mit-blue.svg)](https://andrsd.mit-license.org/)


Convert gmsh mesh files into exodusII files.

## Features

- Conversion:
  - 3D: physical volumes are converted to blocks, physical surfaces are converted into side sets
  - 2D: physical surfaces are converted into blocks, physical curves are converted into side sets
- Support for first-order elements: LINE2, TRI3, QUAD4, TET4 and HEX8
- Simple command line interface

## Requirements

- C++17 compiler
- [fmt 9.1](https://github.com/fmtlib/fmt)
- [exodusIIcpp](https://github.com/andrsd/exodusIIcpp)
- [gmshparsercpp](https://github.com/andrsd/gmshparsercpp)

NOTE: `exodusIIcpp` has additional requirements like netCDF4, HDF5 and exodusii from SEACAS.

## How to Build 

1. Make sure you have [SEACAS](https://github.com/sandialabs/seacas/) and [fmt](https://github.com/fmtlib/fmt) installed.
   `exodusIIcpp` and `gmshparser` will be build as part of the build process.
2. Clone the `gmsh2exo` repository
3. Go to the location where you cloned `gmsh2exo`
4. Create a build directory and go into it:
   ```shell
   mkdir build
   cd build
   ```
5. Configure the build
   ```shell
   cmake ..
   ```
6. Build `gmsh2exo`
   ```shell
   make
   ```
7. If you plan on installing the package, supply `-DCMAKE_INSTALL_PREFIX=/install/prefix` during the configure step and do
   ```shell
   make install
   ```

### Notes

- `SEACAS` can be installed [as a whole package](https://github.com/sandialabs/seacas/?tab=readme-ov-file#build-instructions)
  or just with [exdousII](https://github.com/sandialabs/seacas/?tab=readme-ov-file#exodus).

  Debian-based linux distributions offer exodusII as a pre-built package: `libexodusii-dev`

- `exodusIIcpp` and `gmshparsercpp` will build as dynamic libraries by default.
  If you want to build them as static ones, supply `-DEXODUSIICPP_LIBRARY_TYPE=STATIC` and/or `-DGMSHPARSERCPP_LIBRARY_TYPE=STATIC` during the configre step.

## How to Use

```shell
$ gmsh2exo <msh file> <exo file>
```

## Support

If you are having issues and/or looking for help, please let us know in the [discussions](https://github.com/andrsd/gmsh2exo/discussions) located here on github.

## License

The project is licensed under the MIT license.
