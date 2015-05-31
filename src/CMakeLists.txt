if( VERBOSE )
    message( "-src" )
endif()

include_directories( ${ARMADILLO_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} )

set( wf_options wf_options.cpp wf_options.h )

add_subdirectory( qwwad )

add_executable( qwwad_find_donor_state qwwad_find_donor_state.cpp )
target_link_libraries( qwwad_find_donor_state qwwad )

if( NOT APPLE )
add_executable( dispersion_relation dispersion_relation.cpp ${wf_options} )
target_link_libraries( dispersion_relation qwwad ${Boost_LIBRARIES}
	${GSL_LIBRARIES} )

add_executable( dos dos.cpp ${qwwad_options} )
target_link_libraries( dos qwwad ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
endif()
