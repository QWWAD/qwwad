if( VERBOSE )
    message( "-src" )
endif()

include_directories( ${ARMADILLO_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} )

set( wf_options wf_options.cpp wf_options.h )

add_subdirectory( qwwad )

# Compile all the QWWAD programs
add_executable( qwwad_ef_donor_generic  qwwad_ef_donor_generic.cpp )
add_executable( qwwad_ef_dispersion     qwwad_ef_dispersion.cpp )
add_executable( qwwad_ef_infinite_well  qwwad_ef_infinite_well.cpp )
add_executable( qwwad_density_of_states qwwad_density_of_states.cpp )

# Link QWWAD programs to the QWWAD library
target_link_libraries( qwwad_ef_donor_generic  qwwad )
target_link_libraries( qwwad_ef_dispersion     qwwad )
target_link_libraries( qwwad_ef_infinite_well  qwwad )
target_link_libraries( qwwad_density_of_states qwwad )

install( TARGETS
           qwwad_ef_dispersion
           qwwad_ef_donor_generic
	   qwwad_ef_infinite_well
           qwwad_density_of_states
         RUNTIME DESTINATION bin )
