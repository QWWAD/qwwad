if( VERBOSE )
    message( "-src" )
endif()

include_directories( ${ARMADILLO_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} )

set( wf_options wf_options.cpp wf_options.h )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

add_subdirectory( qwwad )

# A convenience function for adding an executable program
function(add_qwwad_program progname)
	set(_CPPSRC ${progname}.cpp)
	add_executable(${progname} ${_CPPSRC})
	target_link_libraries(${progname} qwwad)
	install(TARGETS ${progname} RUNTIME DESTINATION bin )
endfunction()

add_qwwad_program(qwwad_density_of_states)
add_qwwad_program(qwwad_ef_dispersion)
add_qwwad_program(qwwad_ef_donor_generic)
add_qwwad_program(qwwad_ef_infinite_well)
