if( VERBOSE )
    message( "  /integration-tests" )
endif()

#include_directories( ${PROJECT_SOURCE_DIR}/src ${GTEST_INCLUDE_DIR} )

add_qwwad_test(qwwad-schroedinger-infinite-well-tests)
