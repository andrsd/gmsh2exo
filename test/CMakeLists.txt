project(gmsh2exo-test)

include(TesterExodusII)

set(GMSH2EXO_ASSETS_DIR ${PROJECT_SOURCE_DIR}/assets CACHE PATH "" FORCE)
mark_as_advanced(FORCE GMSH2EXO_ASSETS_DIR)
set(GMSH2EXO_GOLD_DIR ${PROJECT_SOURCE_DIR}/gold CACHE PATH "" FORCE)
mark_as_advanced(FORCE GMSH2EXO_GOLD_DIR)

add_test_exodiff(
    NAME   2blk
    BIN    gmsh2exo
    INPUT  ${GMSH2EXO_ASSETS_DIR}/2blk.msh
    OUTPUT ${PROJECT_BINARY_DIR}/2blk.exo
    GOLD   ${GMSH2EXO_GOLD_DIR}/2blk.exo
)
set_tests_properties(
    2blk
    PROPERTIES
        ENVIRONMENT LLVM_PROFILE_FILE=2blk-test.profraw
)

add_test_exodiff(
    NAME   square
    BIN    gmsh2exo
    INPUT  ${GMSH2EXO_ASSETS_DIR}/square.msh
    OUTPUT ${PROJECT_BINARY_DIR}/square.exo
    GOLD   ${GMSH2EXO_GOLD_DIR}/square.exo
)
set_tests_properties(
    square
    PROPERTIES
        ENVIRONMENT LLVM_PROFILE_FILE=square-test.profraw
)

add_test_exodiff(
    NAME   cube-2x2x2
    BIN    gmsh2exo
    INPUT  ${GMSH2EXO_ASSETS_DIR}/cube-2x2x2.msh
    OUTPUT ${PROJECT_BINARY_DIR}/cube-2x2x2.exo
    GOLD   ${GMSH2EXO_GOLD_DIR}/cube-2x2x2.exo
)
set_tests_properties(
    cube-2x2x2
    PROPERTIES
        ENVIRONMENT LLVM_PROFILE_FILE=cube-2x2x2-test.profraw
)
