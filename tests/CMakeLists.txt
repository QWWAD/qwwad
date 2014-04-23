if( VERBOSE )
    message( "tests" )
endif()

project( tests )


set( TARGET_SRC
        some_tests.cpp
    )

add_executable( tests ${TARGET_SRC} )
include_directories( ${src_SOURCE_DIR} ${GTEST_INCLUDE_DIR} )
target_link_libraries( tests gtest )

install( TARGETS tests
         RUNTIME DESTINATION test
       )

