cmake_minimum_required( VERSION 3.0 )

include(CMakeScripts/ConfigEnv.cmake)
include(FindPkgConfig)

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

pkg_check_modules( LIBXMLPP REQUIRED "libxml++-2.6 >= ${LIBXMLPP_REQUIRED_VERSION}" )
include_directories(SYSTEM ${LIBXMLPP_INCLUDE_DIRS})

include( CMakeScripts/make_check_macros.cmake )

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

# A function to simplify building and installing QWWAD programs
macro(add_qwwad_program program_name description)
	list(APPEND qwwad_programs             ${program_name})
	list(APPEND qwwad_program_descriptions ${description})
endmacro(add_qwwad_program)

# Define a list of all QWWAD programs, and their descriptions here
add_qwwad_program(qwwad_charge_density           "charge density in a heterostructure")
add_qwwad_program(qwwad_cs_single_spiral         "atomic positions in single-spiral of zinc blende crystal")
add_qwwad_program(qwwad_cs_zinc_blende           "atomic positions in a zinc blende crystal")
add_qwwad_program(qwwad_density_of_states        "density of states in 1D, 2D and 3D systems")
add_qwwad_program(qwwad_diffuse                  "solve diffusion equation for a nominal heterostructure")
add_qwwad_program(qwwad_ef_band_edge             "band-edge potential for a heterostructure")
add_qwwad_program(qwwad_ef_cylindrical_wire      "eigenstates for a cylindrical quantum wire")
add_qwwad_program(qwwad_ef_cylindrical_wire_wf   "eigenstates for a cylindrical quantum wire (wave functions)")
add_qwwad_program(qwwad_ef_dispersion            "dispersion relation for a set of energy subbands")
add_qwwad_program(qwwad_ef_donor_generic         "eigenstates for a donor in a heterostructure (generic search)")
add_qwwad_program(qwwad_ef_donor_specific        "eigenstates for a donor in a heterostructure (specific form-factor)")
add_qwwad_program(qwwad_ef_exciton               "eigenstates for an exciton in a heterostructure")
add_qwwad_program(qwwad_ef_generic               "eigenstates in a generic heterostructure")
add_qwwad_program(qwwad_ef_infinite_well         "eigenstates in an infinite quantum well")
add_qwwad_program(qwwad_ef_infinite_wire         "eigenstates in an infinite quantum wire")
add_qwwad_program(qwwad_ef_parabolic_well        "eigenstates in a parabolic well")
add_qwwad_program(qwwad_ef_plot                  "translate wavefunction data into plottable form")
add_qwwad_program(qwwad_ef_plot_3d               "generate 3D wavefunction plotting script for MATLAB")
add_qwwad_program(qwwad_ef_poeschl_teller        "eigenstates in a Poeschl-Teller well")
add_qwwad_program(qwwad_ef_spherical_dot         "eigenstates in a spherical quantum dot")
add_qwwad_program(qwwad_ef_spherical_dot_wf      "eigenstates in a spherical quantum dot (wavefunctions)")
add_qwwad_program(qwwad_ef_square_well           "eigenstates in a finite square quantum well")
add_qwwad_program(qwwad_ef_superlattice          "eigenstates of a Kronig-Penney superlattice")
add_qwwad_program(qwwad_ef_zeeman                "Zeeman-splitting contribution to potential profile")
add_qwwad_program(qwwad_fermi_distribution       "Fermi-Dirac distributions for a set of subbands")
add_qwwad_program(qwwad_material_property        "look up property for a given material")
add_qwwad_program(qwwad_mesh                     "generate 1D mesh for numerical simulations")
add_qwwad_program(qwwad_poisson                  "space-charge potential from Poission equation")
add_qwwad_program(qwwad_population_init          "initial estimate of subband populations")
add_qwwad_program(qwwad_pp_charge_density        "charge-density from pseudopotential calculations")
add_qwwad_program(qwwad_pp_dispersion            "dispersion relation from pseudopotential calculations")
add_qwwad_program(qwwad_pp_form_factor           "form-factor for pseudopotential calculations")
add_qwwad_program(qwwad_pp_large_basis           "large-basis pseudopotential calculation")
add_qwwad_program(qwwad_pp_large_basis_so        "large-basis pseudopotential calculation with spin-orbit splitting")
add_qwwad_program(qwwad_pp_lattice_vector_table  "sort reciprocal lattice vectors in acending magnitude")
add_qwwad_program(qwwad_pp_superlattice          "pseudopotential calculation of states in superlattice")
add_qwwad_program(qwwad_reciprocal_fcc           "reciprocal lattice vectors for FCC crystal")
add_qwwad_program(qwwad_reciprocal_cube          "reciprocal lattice vectors for simple cubic crystal")
add_qwwad_program(qwwad_reciprocal_single_spiral "reciprocal lattice vectors for single spiral of FCC crystal")
add_qwwad_program(qwwad_specific_heat_capacity   "specific heat capacity")
add_qwwad_program(qwwad_spin_flip_raman          "spin-flip Raman spectrum")
add_qwwad_program(qwwad_sr_acoustic_phonon       "acoustic phonon scattering rate")
add_qwwad_program(qwwad_sr_alloy_disorder        "alloy disorder scattering rate")
add_qwwad_program(qwwad_sr_carrier_carrier       "carrier-carrier scattering rate")
add_qwwad_program(qwwad_sr_interface_roughness   "interface roughness scattering rate")
add_qwwad_program(qwwad_sr_impurity              "impurity scattering rate")
add_qwwad_program(qwwad_sr_lo_phonon             "LO-phonon scattering rate")
add_qwwad_program(qwwad_sr_radiative             "radiative scattering rate")
add_qwwad_program(qwwad_superlattice_k           "wave-vectors for superlattice pseudopotential model")
add_qwwad_program(qwwad_thermal_1d               "temperature profile using a 1D numerical simulation")
add_qwwad_program(qwwad_thermal_rc               "temperature profile using a 1D R-C model")
add_qwwad_program(qwwad_tx_double_barrier        "transmission through a double barrier")
add_qwwad_program(qwwad_tx_double_barrier_iv     "current-voltage relation for a double barrier")
add_qwwad_program(qwwad_tx_single_barrier        "transmission through a single barrier")
add_qwwad_program(qwwad_uncertainty              "uncertainty product for position and momentum")
add_qwwad_program(qwwad_critical_thickness       "critical thickness for a strained film")

add_subdirectory( src )
add_subdirectory( doc )
add_subdirectory( examples )

setup_make_check()
include_directories(src src/googletest/googletest/include ${CMAKE_CURRENT_BINARY_DIR})
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
