cmake_minimum_required(VERSION 3.16)

project(gmsh2exo
    VERSION 0.7
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(FetchContent)
include(${CMAKE_SOURCE_DIR}/cmake/CodeCoverage.cmake)

option(GMSH2EXO_BUILD_TESTS "Build tests" NO)

find_package(fmt 9.1 REQUIRED)

FetchContent_Declare(
    exodusIIcpp
    GIT_REPOSITORY https://github.com/andrsd/exodusIIcpp.git
    GIT_TAG        v1.9
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_EXODUSIICPP
    FETCHCONTENT_UPDATES_DISCONNECTED_EXODUSIICPP
    EXODUSIICPP_BUILD_TESTS
)

FetchContent_Declare(
    gmshparsercpp
    GIT_REPOSITORY https://github.com/andrsd/gmshparsercpp.git
    GIT_TAG        v1.3
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_GMSHPARSERCPP
    FETCHCONTENT_UPDATES_DISCONNECTED_GMSHPARSERCPP
    GMSHPARSERCPP_BUILD_TESTS
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

if(EXODUSIICPP_LIBRARY_TYPE STREQUAL "STATIC")
    FetchContent_GetProperties(exodusIIcpp)
    if(NOT exodusIIcpp_POPULATED)
        FetchContent_Populate(exodusIIcpp)
        add_subdirectory(${exodusiicpp_SOURCE_DIR} ${exodusiicpp_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
else()
    FetchContent_MakeAvailable(exodusIIcpp)
endif()

if(GMSHPARSERCPP_LIBRARY_TYPE STREQUAL "STATIC")
    FetchContent_GetProperties(gmshparsercpp)
    if(NOT gmshparsercpp_POPULATED)
        FetchContent_Populate(gmshparsercpp)
        add_subdirectory(${gmshparsercpp_SOURCE_DIR} ${gmshparsercpp_BINARY_DIR} EXCLUDE_FROM_ALL)
    endif()
else()
    FetchContent_MakeAvailable(gmshparsercpp)
endif()

if (GMSH2EXO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

mark_as_advanced(FORCE
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_QUIET
    FETCHCONTENT_UPDATES_DISCONNECTED
)

#
file(GLOB SRCS CONFIGURE_DEPENDS src/*.cpp)
add_executable(${PROJECT_NAME} ${SRCS})
target_code_coverage(${PROJECT_NAME})

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/contrib/cxxopts
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        fmt::fmt
        exodusIIcpp
        gmshparsercpp
)

install(
    TARGETS gmsh2exo
    RUNTIME DESTINATION bin
)

#
configure_file(include/Gmsh2ExoConfig.h.in Gmsh2ExoConfig.h)
