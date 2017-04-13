if( VERBOSE )
    message( "  /integration-tests" )
endif()

include_directories( ${PROJECT_SOURCE_DIR}/src ${GTEST_INCLUDE_DIR} )

set( TARGET_SRC
        qwwad-schroedinger-infinite-well-tests.cpp
   )

add_executable( integration_tests ${TARGET_SRC} )
target_link_libraries( integration_tests
                       libqwwad
                       ${GTEST_LIBS}
                     )

add_test_to_make_check( integration_tests )
