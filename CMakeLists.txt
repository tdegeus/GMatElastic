
cmake_minimum_required(VERSION 3.0)

# Options
# =======

project(GMatElastic)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Internals
# =========

set(GMATELASTIC_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Version
# =======

file(
    STRINGS "${GMATELASTIC_INCLUDE_DIR}/GMatElastic/config.h" GMatElastic_version_defines
    REGEX "#define GMATELASTIC_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${GMatElastic_version_defines})
    if(ver MATCHES "#define GMATELASTIC_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(GMATELASTIC_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()

set(GMATELASTIC_VERSION
    ${GMATELASTIC_VERSION_MAJOR}.${GMATELASTIC_VERSION_MINOR}.${GMATELASTIC_VERSION_PATCH})

message(STATUS "Building GMatElastic v${GMATELASTIC_VERSION}")

# Build
# =====

set(GMATELASTIC_HEADERS
    ${GMATELASTIC_INCLUDE_DIR}/GMatElastic/Cartesian3d.h
    ${GMATELASTIC_INCLUDE_DIR}/GMatElastic/Cartesian3d.hpp
    ${GMATELASTIC_INCLUDE_DIR}/GMatElastic/Cartesian3d_Elastic.hpp
    ${GMATELASTIC_INCLUDE_DIR}/GMatElastic/Cartesian3d_Matrix.hpp
    ${GMATELASTIC_INCLUDE_DIR}/GMatElastic/config.h)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
    ${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${GMATELASTIC_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets)

export(
    EXPORT ${PROJECT_NAME}-targets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

install(
    FILES ${GMATELASTIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

set(GMATELASTIC_CMAKECONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for ${PROJECT_NAME}Config.cmake")

configure_package_config_file(
    "GMatElasticConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION ${GMATELASTIC_CMAKECONFIG_INSTALL_DIR})

set(_GMATELASTIC_CMAKE_SIZEOF_VOID_P ${CMAKE_SIZEOF_VOID_P})

unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${GMATELASTIC_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMATELASTIC_CMAKE_SIZEOF_VOID_P})

install(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION ${GMATELASTIC_CMAKECONFIG_INSTALL_DIR})

install(
    EXPORT ${PROJECT_NAME}-targets
    FILE "${PROJECT_NAME}Targets.cmake"
    DESTINATION ${GMATELASTIC_CMAKECONFIG_INSTALL_DIR})

configure_file(
    "GMatElastic.pc.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
     @ONLY)

install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

if(BUILD_TESTS)
    add_subdirectory(test)
endif()
