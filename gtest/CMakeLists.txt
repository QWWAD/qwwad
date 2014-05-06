if( VERBOSE )
    message( "-gtest" )
endif()

project( gtest )

# The GTest include directory will need to be added to the search path for all tests.
# Therefore it needs to be exported at parent scope so it can be accessed from other
# CMakeLists files.
set( GTEST_INCLUDE_DIR ${gtest_SOURCE_DIR}/include PARENT_SCOPE )
# The relevant thread options needed to be added when using gtest...
# Add these to a standard variable so they are not forgotten!
set( GTEST_LIBS gtest ${CMAKE_THREAD_LIBS_INIT} PARENT_SCOPE )

# Since only the default main function (from gtest_main.cc) is used within the tests,
# gtest and gtest_main are compiled into one library to link against. If in the future
# a different main function is required these will need to be split into two.
include_directories( ${gtest_SOURCE_DIR} ${gtest_SOURCE_DIR}/include )
add_library( gtest src/gtest_main.cc src/gtest-all.cc )
