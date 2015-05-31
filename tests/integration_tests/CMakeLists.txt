if( VERBOSE )
    message( "  /integration-tests" )
endif()

include_directories( ${GSL_INCLUDE_DIR} ${QWWAD_SOURCE_DIR} ${GTEST_INCLUDE_DIR} )
set( TARGET_SRC
        qwwad-schroedinger-infinite-well-tests.cpp
   )

add_executable( integration_tests ${TARGET_SRC} )
target_link_libraries( integration_tests
                       qwwad
                       ${GTEST_LIBS}
                       ${GSL_LIBRARIES} 
                     )

add_test_to_make_check( integration_tests )
