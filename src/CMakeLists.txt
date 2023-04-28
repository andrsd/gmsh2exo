add_executable(${PROJECT_NAME}
    main.cpp
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${PROJECT_BINARY_DIR}
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        exodusIIcpp
        gmshparsercpp
        cxxopts
)

install(
    TARGETS gmsh2exo
    RUNTIME DESTINATION bin
)
