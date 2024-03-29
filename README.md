# OpenDSP

OpenDSP is an open source, cross platform, C++ Implementation of common digital
signal processing algorithms.

## Documentation

The available functions are documented with Doxygen style comments in the C++
header files. For now the easiest way to view the documentation is by simply
opening the header files and reading the comments, but I plan on finding a way
to actually build and formally publish the documentation in the near future.

## Disclaimer

This library is still a work in progress, but since it can already calculate
fast Fourier transforms of arbitrary sizes relatively quickly I am choosing to
make it public now because others may find it useful.

# Build Instructions

OpenDSP can be built as a standalone CMake project.

    cmake -H. -Bbuild
    cmake --build build

## Requirements

* CMake v3.14.0 (or later)
* C++11 compliant compiler (or later)
* Conan package manager v1.18.1 (or later)

## Dependencies

This project requires the Conan package manager to automatically manage all of
its dependencies. To install Conan see:
https://docs.conan.io/en/latest/installation.html
