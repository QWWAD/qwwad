if( VERBOSE )
    message( "-src" )
endif()

include_directories( ${ARMADILLO_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} )

set( wf_options wf_options.cpp wf_options.h )

add_subdirectory( qwwad )

add_executable( qwwad_ef_donor_generic
                qwwad_ef_donor_generic.cpp )

target_link_libraries( qwwad_ef_donor_generic qwwad )

add_executable( qwwad_ef_dispersion
                qwwad_ef_dispersion.cpp )

target_link_libraries( qwwad_ef_dispersion qwwad )

add_executable( qwwad_density_of_states
                qwwad_density_of_states.cpp )

target_link_libraries( qwwad_density_of_states qwwad )

install( PROGRAMS
           ${CMAKE_BINARY_DIR}/src/qwwad_ef_dispersion
           ${CMAKE_BINARY_DIR}/src/qwwad_ef_donor_generic
           ${CMAKE_BINARY_DIR}/src/qwwad_density_of_states
         DESTINATION
           ${CMAKE_INSTALL_PREFIX}/bin )
