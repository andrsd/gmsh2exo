add_executable(${PROJECT_NAME}
    main.cpp
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        exodusIIcpp
        gmshparsercpp
)

install(
    TARGETS gmsh2exo
    RUNTIME DESTINATION bin
)
