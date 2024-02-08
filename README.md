gmsh2exo
========

[![qa](https://github.com/andrsd/gmsh2exo/actions/workflows/qa.yml/badge.svg)](https://github.com/andrsd/gmsh2exo/actions/workflows/qa.yml)
[![build](https://github.com/andrsd/gmsh2exo/actions/workflows/build.yml/badge.svg)](https://github.com/andrsd/gmsh2exo/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/andrsd/gmsh2exo/graph/badge.svg?token=JHXYRN6N2X)](https://codecov.io/gh/andrsd/gmsh2exo)
[![License](http://img.shields.io/:license-mit-blue.svg)](https://andrsd.mit-license.org/)


Convert gmsh mesh files into exodusII files.

## Features

- Converts individual "blocks" respecting physical entities assignment
- Generates side sets from DIM-1 physical entities
- Simple command line interface
- Support for first-order elements: LINE2, TRI3, QUAD4, TET4 and HEX8

## Requirements

- C++17 compiler
- [fmt 9.1](https://github.com/fmtlib/fmt)
- [exodusIIcpp](https://github.com/andrsd/exodusIIcpp)
- [gmshparsercpp](https://github.com/andrsd/gmshparsercpp)

NOTE: `exodusIIcpp` has additional requirements like netCDF4, HDF5 and exodusii from SEACAS.

## How to use it

```shell
$ gmsh2exo <msh file> <exo file>
```
