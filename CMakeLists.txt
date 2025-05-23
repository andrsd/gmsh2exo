cmake_minimum_required(VERSION 3.16)

project(gmsh2exo
    VERSION 0.9
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
    GIT_TAG        b4172a0999e7393d629e8d52872512e3e4368a89
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_EXODUSIICPP
    FETCHCONTENT_UPDATES_DISCONNECTED_EXODUSIICPP
    EXODUSIICPP_BUILD_TESTS
)

FetchContent_Declare(
    gmshparsercpp
    GIT_REPOSITORY https://github.com/andrsd/gmshparsercpp.git
    GIT_TAG        7ed3ae3d8c251e79cee6de21ba73e83f4ac45759
)
mark_as_advanced(FORCE
    FETCHCONTENT_SOURCE_DIR_GMSHPARSERCPP
    FETCHCONTENT_UPDATES_DISCONNECTED_GMSHPARSERCPP
    GMSHPARSERCPP_BUILD_TESTS
)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
mark_as_advanced(FORCE
    BUILD_GMOCK
    INSTALL_GTEST
    FETCHCONTENT_SOURCE_DIR_GOOGLETEST
    FETCHCONTENT_UPDATES_DISCONNECTED_GOOGLETEST
)

FetchContent_MakeAvailable(exodusIIcpp)
FetchContent_MakeAvailable(gmshparsercpp)

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

target_compile_definitions(${PROJECT_NAME}
    PRIVATE
        -DGMSH2EXO_VERSION="${CMAKE_PROJECT_VERSION}"
)

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
