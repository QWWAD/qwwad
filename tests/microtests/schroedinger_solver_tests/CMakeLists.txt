if( VERBOSE )
    message( "    /schroedinger_solver_tests" )
endif()

include_directories( ${PROJECT_SOURCE_DIR}/src ${GTEST_INCLUDE_DIR} )

set( TARGET_SRC
        abstract_schroedinger_solver_tests.cpp
    )

add_executable( schroedinger_solver_tests ${TARGET_SRC} )
target_link_libraries( schroedinger_solver_tests
                       libqwwad
                       ${GTEST_LIBS}
                     )

add_test_to_make_check( schroedinger_solver_tests )
