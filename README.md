# GMatLinearElastic

[![Travis](https://travis-ci.com/tdegeus/GMatElastic.svg?branch=master)](https://travis-ci.com/tdegeus/GMatElastic)

Linear elastic material model.
An overview of the theory can be found in `docs/readme.tex` 
conveniently compiled to this [PDF](docs/readme.pdf).

# Contents

<!-- MarkdownTOC levels="1,2,3" -->

- [Disclaimer](#disclaimer)
- [Implementation](#implementation)
    - [C++ and Python](#c-and-python)
    - [Cartesian3d](#cartesian3d)
        - [Overview](#overview)
        - [Function names](#function-names)
        - [Storage](#storage)
        - [Example](#example)
    - [Debugging](#debugging)
- [Installation](#installation)
    - [C++ headers](#c-headers)
        - [Using conda](#using-conda)
        - [From source](#from-source)
    - [Python module](#python-module)
        - [Using conda](#using-conda-1)
        - [From source](#from-source-1)
- [Compiling](#compiling)
    - [Using CMake](#using-cmake)
        - [Example](#example-1)
        - [Targets](#targets)
        - [Optimisation](#optimisation)
    - [By hand](#by-hand)
    - [Using pkg-config](#using-pkg-config)
- [References / Credits](#references--credits)

<!-- /MarkdownTOC -->

# Disclaimer

This library is free to use under the
[MIT license](https://github.com/tdegeus/GMatElastic/blob/master/LICENSE).
Any additions are very much appreciated, in terms of suggested functionality, code,
documentation, testimonials, word-of-mouth advertisement, etc.
Bug reports or feature requests can be filed on
[GitHub](https://github.com/tdegeus/GMatElastic).
As always, the code comes with no guarantee.
None of the developers can be held responsible for possible mistakes.

Download: 
[.zip file](https://github.com/tdegeus/GMatElastic/zipball/master) |
[.tar.gz file](https://github.com/tdegeus/GMatElastic/tarball/master).

(c - [MIT](https://github.com/tdegeus/GMatElastic/blob/master/LICENSE))
T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me |
[github.com/tdegeus/GMatElastic](https://github.com/tdegeus/GMatElastic)

# Implementation

## C++ and Python

The code is a C++ header-only library (see [installation notes](#c-headers)), 
but a Python module is also provided (see [installation notes](#python-module)).
The interfaces are identical except:

+   All *xtensor* objects (`xt::xtensor<...>`) are *NumPy* arrays in Python. 
    Overloading based on rank is also available in Python.
+   The Python module cannot change output objects in-place: 
    only functions whose name starts with a capital letter are included, see below.
+   All `::` in C++ are `.` in Python.

## Cartesian3d

[Cartesian3d.h](include/GMatLinearElastic/Cartesian3d.h)

### Overview

At the material point level to model is implemented in the class:

+   `Elastic`: linear elastic material model.

There is a `Matrix` class that allows you to combine 
have a single API for a matrix of material points. 

### Function names

+   Functions whose name starts with a capital letter (e.g. `Stress`) 
    return their result (allocating it internally).
+   Functions whose name starts with a small letter (e.g. `stress`) 
    write to the, fully allocated, last input argument(s) 
    (avoiding re-allocation, but making the user responsible to do it properly).

### Storage

+   Scalar
    ```cpp
    double
    ```

+   2nd-order tensor
    ```cpp
    xt::xtensor_fixed<double, xt::xshape<3, 3>> = 
    GMatElastic::Cartesian3d::Tensor2
    ```

+   List *(i)* of second order tensors *(x,y)* : *A(i,x,y)*
    ```cpp
    xt::xtensor<double,3>
    ```
    Note that the shape is `[I, 3, 3]`.

+   Matrix *(i,j)* of second order tensors *(x,y)* : *A(i,j,x,y)*
    ```cpp
    xt::xtensor<double,4>
    ```
    Note that the shape is `[I, J, 3, 3]`.

### Example

Only a partial examples are presented here, meant to understand the code's structure.

#### Individual material point

```cpp
#include <GMatLinearElastic/Cartesian3d.h>

namespace GMat = GMatLinearElastic::Cartesian3d;

int main()
{
    // a single material point
    GMat::Elastic elastic(K, G);
    ...
    
    // set strain tensor (follows e.g. from FEM discretisation)
    GMat::Tensor2 Eps;
    ...
    
    // compute stress (including allocation of the result)
    GMat::Tensor2 Sig = elastic.Stress(Eps);
    // OR compute stress without (re)allocating the results
    // in this case "Sig" has to be of the correct type and shape
    elastic.stress(Eps, Sig); 
    ...

    return 0;
}
```

#### Matrix of material points

```cpp
#include <GMatLinearElastic/Cartesian3d.h>

namespace GMat = GMatLinearElastic::Cartesian3d;

int main()
{
    // a matrix, of shape [nelem, nip], of material points
    GMat::Elastic matrix(nelem, nip);

    // set materials:
    // points where I(x,y) == 1 are assigned, points where I(x,y) == 0 are skipped
    // all points can only be assigned once
    matrix.setElastic(I, K, G);
    ...

    // set strain tensor (follows e.g. from FEM discretisation)
    xt::xtensor<double,4> eps = xt::empty<double>({nelem, nip, 3ul, 3ul});
    ... 

    // compute stress (allocate result)
    xt::xtensor<double,4> sig = matrix.Stress(eps);
    // OR compute stress without (re)allocating the results
    // in this case "sig" has to be of the correct type and shape
    matrix.stress(eps, sig); 
    ...

    return 0;
}
```

## Debugging

To enable assertions define `GMATELASTIC_ENABLE_ASSERT` 
**before** including *GMatElastic* for the first time. 

Using *CMake* this can be done using the `GMatElastic::assert` target 
(see [below](#using-cmake)).

>   To also enable assertions of *xtensor* also define `XTENSOR_ENABLE_ASSERT` 
>   **before** including *xtensor* (and *GMatElastic*) for the first time. 
>   
>   Using *CMake* all assertions are enabled using the `GMatElastic::debug` target 
>   (see [below](#using-cmake)).

# Installation

## C++ headers

### Using conda

```bash
conda install -c conda-forge gmatelastic
```

### From source

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

```bash
conda install -c conda-forge python-gmatelastic
```

Note that *xsimd* and hardware optimisations are **not enabled**. 
To enable them you have to compile on your system, as is discussed next.

### From source

>   You need *xtensor*, *pyxtensor* and optionally *xsimd* as prerequisites. 
>   Additionally, Python needs to know how to find them. 
>   The easiest is to use *conda* to get the prerequisites:
> 
>   ```bash
>   conda install -c conda-forge pyxtensor
>   conda install -c conda-forge xsimd
>   ```
>   
>   If you then compile and install with the same environment 
>   you should be good to go. 
>   Otherwise, a bit of manual labour might be needed to
>   treat the dependencies.

```bash
# Download GMatElastic
git checkout https://github.com/tdegeus/GMatElastic.git
cd GMatElastic

# Compile and install the Python module
python setup.py build
python setup.py install
# OR you can use one command (but with less readable output)
python -m pip install .
```

# Compiling

## Using CMake

### Example

Using *GMatElastic* your `CMakeLists.txt` can be as follows

```cmake
cmake_minimum_required(VERSION 3.1)
project(example)
find_package(GMatElastic REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE GMatElastic)
```

### Targets

The following targets are available:

*   `GMatElastic`
    Includes *GMatElastic* and the *xtensor* dependency.

*   `GMatElastic::assert`
    Enables assertions by defining `GMATELASTIC_ENABLE_ASSERT`.

*   `GMatElastic::debug`
    Enables all assertions by defining 
    `GMATELASTIC_ENABLE_ASSERT` and `XTENSOR_ENABLE_ASSERT`.

*   `GMatElastic::compiler_warings`
    Enables compiler warnings (generic).

### Optimisation

It is advised to think about compiler optimization and enabling *xsimd*. 
Using *CMake* this can be done using the `xtensor::optimize` and `xtensor::use_xsimd` targets. 
The above example then becomes:

```cmake
cmake_minimum_required(VERSION 3.1)
project(example)
find_package(GMatElastic REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE 
    GMatElastic 
    xtensor::optimize 
    xtensor::use_xsimd)
```

See the [documentation of xtensor](https://xtensor.readthedocs.io/en/latest/) concerning optimization.

## By hand

Presuming that the compiler is `c++`, compile using:

```
c++ -I/path/to/GMatElastic/include ...
```

Note that you have to take care of the *xtensor* dependency, the C++ version, optimization, 
enabling *xsimd*, ...

## Using pkg-config

Presuming that the compiler is `c++`, compile using:

```
c++ `pkg-config --cflags GMatElastic` ...
```

Note that you have to take care of the *xtensor* dependency, the C++ version, optimization, 
enabling *xsimd*, ...

# References / Credits

+   [xtensor](https://github.com/QuantStack/xtensor) is used under the hood.
