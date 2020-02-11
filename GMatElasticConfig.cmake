# GMatElastic cmake module
#
# This module sets the target:
#
#     GMatElastic
#
# In addition, it sets the following variables:
#
#     GMatElastic_FOUND - true if the library is found
#     GMatElastic_VERSION - the library's version
#     GMatElastic_INCLUDE_DIRS - directory containing the library's headers
#
# The following support targets are defined to simplify things:
#
#     GMatElastic::compiler_warnings - enable compiler warnings
#     GMatElastic::assert - enable library assertions
#     GMatElastic::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatElastic"

if(NOT TARGET GMatElastic)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatElasticTargets.cmake")
endif()

# Define "GMatElastic_INCLUDE_DIRS"

get_target_property(
    GMatElastic_INCLUDE_DIRS
    GMatElastic
    INTERFACE_INCLUDE_DIRECTORIES)

# Find dependencies

find_dependency(xtensor)

# Define support target "GMatElastic::compiler_warnings"

if(NOT TARGET GMatElastic::compiler_warnings)
    add_library(GMatElastic::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET GMatElastic::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET GMatElastic::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "GMatElastic::assert"

if(NOT TARGET GMatElastic::assert)
    add_library(GMatElastic::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastic::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTIC_ENABLE_ASSERT)
endif()

# Define support target "GMatElastic::debug"

if(NOT TARGET GMatElastic::debug)
    add_library(GMatElastic::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastic::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT GMATELASTIC_ENABLE_ASSERT)
endif()
