cmake_minimum_required( VERSION 2.8.11 )

include(CMakeScripts/ConfigEnv.cmake)

list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeScripts" )

project( "Quantum Wells, Wires and Dots" )
set( QWWAD_VERSION 1.5+devel )
set( QWWAD_URL "https://launchpad.net/qwwad" )
set( QWWAD_BUGREPORT "https://bugs.launchpad.net/qwwad" )

option( VERBOSE "Show information about CMake build configuration." )

# Enable C++11 builds
# TODO: With CMake > 3.1, can use target_compile_features() to do this automatically
set(CMAKE_CXX_FLAGS "-std=c++11")

find_package( LAPACK REQUIRED )
find_package( Armadillo 4.000 REQUIRED )
find_package( Boost 1.42 COMPONENTS program_options REQUIRED )
find_package( GSL REQUIRED )
include( CMakeScripts/make_check_macros.cmake )

if( VERBOSE )
    message( "Generating config.h" )
endif()


configure_file( ${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_BINARY_DIR}/config.h )
include_directories( ${CMAKE_BINARY_DIR} )

add_definitions( -DHAVE_CONFIG_H )

if( VERBOSE )
    message( "Adding folders to build:" )
endif()

add_subdirectory( src )
add_subdirectory( examples )

setup_make_check()
add_subdirectory( gtest EXCLUDE_FROM_ALL )
add_subdirectory( tests EXCLUDE_FROM_ALL )
