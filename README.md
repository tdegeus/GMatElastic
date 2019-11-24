
# LinearElastic

[![Travis](https://travis-ci.com/tdegeus/GMatElastic.svg?branch=master)](https://travis-ci.com/tdegeus/GMatElastic)
[![Build status](https://ci.appveyor.com/api/projects/status/c4c1l93dmejo76ym?svg=true)](https://ci.appveyor.com/project/tdegeus/gmatelastic)

Linear elastic material model. An overview of the theory can be found in `docs/theory` in particular in this [PDF](docs/readme.pdf).

>   **Disclaimer**
>   
>   This library is free to use under the [MIT license](https://github.com/tdegeus/GMatElastic/blob/master/LICENSE). Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word-of-mouth advertisement, etc. Bug reports or feature requests can be filed on [GitHub](https://github.com/tdegeus/GMatElastic). As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.
>   
>   Download: [.zip file](https://github.com/tdegeus/GMatElastic/zipball/master) | [.tar.gz file](https://github.com/tdegeus/GMatElastic/tarball/master).
>   
>   (c - [MIT](https://github.com/tdegeus/GMatElastic/blob/master/LICENSE)) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | [github.com/tdegeus/GMatElastic](https://github.com/tdegeus/GMatElastic)

# Contents

<!-- MarkdownTOC -->

- [Implementation](#implementation)
- [Installation](#installation)
    - [C++ headers](#c-headers)
        - [Using conda](#using-conda)
    - [From source](#from-source)
    - [Python module](#python-module)
        - [Using conda](#using-conda-1)
    - [From source](#from-source-1)
- [Compiling](#compiling)
    - [By hand](#by-hand)
    - [Using pkg-config](#using-pkg-config)
    - [Using `CMakeLists.txt`](#using-cmakeliststxt)

<!-- /MarkdownTOC -->

# Implementation

The headers are meant to be self-explanatory, please check them out:

* [Cartesian3d.h](include/LinearElastic/Cartesian3d.h)

Only a tiny example is presented here, that is meant to understand the code's structure:

```cpp
#include <LinearElastic/Cartesian3d.h>

int main()
{
    // a single material point
    // - create class
    LinearElastic::Cartesian3d::Elastic elastic(K, G);
    // - compute stress [allocate result]
    Sig = elastic.Stress(Eps);
    ...
    // - compute stress [no allocation]
    elastic.stress(Eps, Sig); 
    ...

    // a "matrix" of material points
    // - create class
    LinearElastic::Cartesian3d::Elastic matrix(nelem, nip);
    // - set material
    matrix.setElastic(I, K, G);
    // - compute stress [allocate result]
    Sig = matrix.Stress(Eps);
    ...
    // - compute stress [no allocation]
    matrix.stress(Eps, Sig); 
    ...
}
```

# Installation

## C++ headers

### Using conda

```bash
conda install -c conda-forge gmatelastic
```

## From source

```bash
# Download GMatElastic
git checkout https://github.com/tdegeus/GMatElastic.git
cd GMatElastic

# Install headers, CMake and pkg-config support
cmake .
make install
```

## Python module

### Using conda

> Warning: this has the disadvantage of xsimd optimisation being switched off

```bash
conda install -c conda-forge python-gmatelastic
```

## From source

> To get the prerequisites you can use conda
> 
> ```bash
> conda install -c conda-forge pyxtensor
> conda install -c conda-forge xsimd
> ```

```bash
# Download GMatElastic
git checkout https://github.com/tdegeus/GMatElastic.git
cd GMatElastic

# Compile and install the Python module
python setup.py build
python setup.py install
```

# Compiling

## By hand

Presuming that the compiler is `c++`, compile using:

```
c++ -I/path/to/GMatElastic/include ...
```

## Using pkg-config

Presuming that the compiler is `c++`, compile using:

```
c++ `pkg-config --cflags GMatElastic` ...
```

## Using `CMakeLists.txt`

Using *GMatElastic* the `CMakeLists.txt` can be as follows

```cmake
cmake_minimum_required(VERSION 3.1)

project(example)

find_package(xtensor REQUIRED)
find_package(GMatElastic REQUIRED)

add_executable(example example.cpp)

target_link_libraries(example
    PRIVATE
    xtensor
    GMatElastic)
```

Compilation can then proceed using 

```bash
cmake .
make
```
