cmake_minimum_required(VERSION 3.16)

project(gmsh2exo
    VERSION 0.5
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(FetchContent)

option(GMSH2EXO_BUILD_TESTS "Build tests" NO)

find_package(fmt REQUIRED)

FetchContent_Declare(
    exodusIIcpp
    GIT_REPOSITORY https://github.com/andrsd/exodusIIcpp.git
    GIT_TAG        v1.5
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_EXODUSIICPP
    FETCHCONTENT_UPDATES_DISCONNECTED_EXODUSIICPP
    EXODUSIICPP_BUILD_TESTS
)

FetchContent_Declare(
    gmshparsercpp
    GIT_REPOSITORY https://github.com/andrsd/gmshparsercpp.git
    GIT_TAG        v1.2
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_GMSHPARSERCPP
    FETCHCONTENT_UPDATES_DISCONNECTED_GMSHPARSERCPP
    GMSHPARSERCPP_BUILD_TESTS
)

FetchContent_Declare(
    cxxopts
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG        v3.1.1
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_GMSHPARSERCPP
    FETCHCONTENT_UPDATES_DISCONNECTED_GMSHPARSERCPP
    CXXOPTS_BUILD_EXAMPLES
    CXXOPTS_BUILD_TESTS
    CXXOPTS_ENABLE_INSTALL
    CXXOPTS_ENABLE_WARNINGS
    CXXOPTS_USE_UNICODE_HELP
)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        release-1.11.0
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
mark_as_advanced(FORCE
    BUILD_GMOCK
    INSTALL_GTEST
    FETCHCONTENT_SOURCE_DIR_GOOGLETEST
    FETCHCONTENT_UPDATES_DISCONNECTED_GOOGLETEST
)

FetchContent_MakeAvailable(exodusIIcpp gmshparsercpp cxxopts googletest)
mark_as_advanced(FORCE
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_QUIET
    FETCHCONTENT_UPDATES_DISCONNECTED
)

#
configure_file(include/Gmsh2ExoConfig.h.in Gmsh2ExoConfig.h)

add_subdirectory(src)
