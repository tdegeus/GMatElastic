
cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatElastic)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatElastic/config.h"
             _GMAT_version_defines
     REGEX "#define GMATELASTIC_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${_GMAT_version_defines})
    if(ver MATCHES "#define GMATELASTIC_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(_GMAT_VERSION_${CMAKE_MATCH_1}
            "${CMAKE_MATCH_2}"
            CACHE INTERNAL "")
    endif()
endforeach()

set(_GMAT_VERSION ${_GMAT_VERSION_MAJOR}.${_GMAT_VERSION_MINOR}.${_GMAT_VERSION_PATCH})

message(STATUS "Building GMatElastic v${_GMAT_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)

add_library(GMatElastic INTERFACE)

target_include_directories(GMatElastic INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatElastic INTERFACE xtensor)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)

install(TARGETS GMatElastic EXPORT GMatElastic-targets)

install(
    EXPORT GMatElastic-targets
    FILE GMatElasticTargets.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastic")

set(_GMAT ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GMatElasticConfigVersion.cmake"
    VERSION ${_GMAT_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMAT})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GMatElasticConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/GMatElasticConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastic")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GMatElastic.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/GMatElastic.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GMatElastic.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("GMatElasticConfig.cmake")

if(BUILD_TESTS)
    add_subdirectory(test)
endif()
