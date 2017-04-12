cmake_minimum_required( VERSION 3.0 )

include(CMakeScripts/ConfigEnv.cmake)

list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeScripts" )

# set a default build configuration
set(CMAKE_BUILD_TYPE_INIT "Release")

# Allow the version to be set in the project command
cmake_policy(SET CMP0048 NEW)

# Define project metadata
# Use even minor-version numbers (e.g., 1.6.x) for development versions
# Use odd  minor-version numbers (e.g., 1.7.x) for "stable" releases
project(qwwad
	VERSION 1.6
       )

message("-- Building ${PROJECT_NAME} version: ${qwwad_VERSION}")
set( QWWAD_URL "https://launchpad.net/qwwad" )
set( QWWAD_BUGREPORT "https://bugs.launchpad.net/qwwad" )

# Define user-configurable build options
option( VERBOSE "Show information about CMake build configuration." )

# Enable C++11 builds
# TODO: With CMake > 3.1, can use target_compile_features() to do this automatically
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Define minimum versions of dependencies here
set(ARMADILLO_REQUIRED_VERSION 4.000)
set(BOOST_REQUIRED_VERSION     1.42)
set(LIBXMLPP_REQUIRED_VERSION  2.30)

# Check that required libraries are available
find_package( Armadillo ${ARMADILLO_REQUIRED_VERSION} REQUIRED )
find_package( Boost ${BOOST_REQUIRED_VERSION}
	      COMPONENTS program_options
	      REQUIRED )
find_package( GSL REQUIRED )
find_package( LAPACK REQUIRED )

include( CMakeScripts/make_check_macros.cmake )

# Define useful directory variables
set(INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INC_DIR "${CMAKE_INSTALL_PREFIX}/include" CACHE PATH "Installation directory for headers")

# Create and install pkg-config file to help 3rd party programs discover the QWWAD
# library
if(VERBOSE)
    message("Generating pkg-config file")
endif()

configure_file( ${CMAKE_SOURCE_DIR}/qwwad.pc.cmake ${CMAKE_BINARY_DIR}/qwwad.pc @ONLY )
install(FILES ${CMAKE_BINARY_DIR}/qwwad.pc DESTINATION share/pkgconfig)

# Create a shared header with platform-specific config data.  This is only used during
# building, and isn't intended to be installed
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
add_subdirectory( doc )
add_subdirectory( examples )

setup_make_check()
add_subdirectory( gtest EXCLUDE_FROM_ALL )
add_subdirectory( tests EXCLUDE_FROM_ALL )

# Configure tarball packaging
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Quantum wells, wires and dots")
set(CPACK_PACKAGE_VENDOR              "Alexander Valavanis and Paul Harrison")
set(CPACK_GENERATOR                   "TGZ")
set(CPACK_SOURCE_GENERATOR            "TGZ")
set(CPACK_PACKAGE_DESCRIPTION_FILE    "${CMAKE_SOURCE_DIR}/README")
set(CPACK_RESOURCE_FILE_LICENSE       "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_VERSION             "${qwwad_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR       "${qwwad_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR       "${qwwad_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH       "${qwwad_VERSION_PATCH}")
set(CPACK_PACKAGE_NAME                "${PROJECT_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME    "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")

include(CPack)
