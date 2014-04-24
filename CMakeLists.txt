cmake_minimum_required( VERSION 2.8.11 )
list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeScripts" )

project( QWWAD )

option( VERBOSE "Show information about CMake build configuration." )

find_package( LAPACK REQUIRED )
find_package( Boost 1.35 REQUIRED )
find_package( GSL REQUIRED )

if( VERBOSE )
    message( "Adding folders to build:" )
endif()

add_subdirectory( src )
add_subdirectory( gtest )
add_subdirectory( tests )
