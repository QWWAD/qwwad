if( VERBOSE )
    message( "  /qwwad" )
endif()

set( qwwad_src
	data-checker.cpp
        donor-energy-minimiser.cpp
        dos-functions.cpp
        double-barrier.cpp
        eigenstate.cpp
        fermi.cpp
        file-io.cpp
        file-io-deprecated.cpp
        linear-algebra.cpp
        maths-helpers.cpp
        options.cpp
        subband.cpp
        mesh.cpp
        schroedinger-solver.cpp
        schroedinger-solver-donor.cpp
        schroedinger-solver-donor-2D.cpp
        schroedinger-solver-donor-3D.cpp
        schroedinger-solver-donor-variable.cpp
        schroedinger-solver-finite-well.cpp
        schroedinger-solver-full.cpp
        schroedinger-solver-infinite-well.cpp
        schroedinger-solver-kronig-penney.cpp
        schroedinger-solver-taylor.cpp
        schroedinger-solver-tridiagonal.cpp
	wf_options.cpp
    )

set( qwwad_h
	data-checker.h
        donor-energy-minimiser.h
        dos-functions.h
        double-barrier.h
        eigenstate.h
        fermi.h
        file-io.h
        file-io-deprecated.h
        linear-algebra.h
        maths-helpers.h
        options.h
        subband.h
        constants.h
        mesh.h
        schroedinger-solver.h
        schroedinger-solver-donor.h
        schroedinger-solver-donor-2D.h
        schroedinger-solver-donor-3D.h
        schroedinger-solver-donor-variable.h
        schroedinger-solver-finite-well.h
        schroedinger-solver-full.h
        schroedinger-solver-infinite-well.h
	schroedinger-solver-kronig-penney.h
        schroedinger-solver-taylor.h
        schroedinger-solver-tridiagonal.h
	wf_options.h
    )

add_library( qwwad SHARED ${qwwad_src} ${qwwad_h} )
set_target_properties( qwwad
	               PROPERTIES
		       SOVERSION 1
		       VERSION   ${QWWAD_VERSION} )
target_link_libraries( qwwad ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${ARMADILLO_LIBRARIES} )

# Install the shared QWWAD library
install(TARGETS qwwad
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib)

# Install the QWWAD API headers in the "qwwad" subfolder
install(FILES ${qwwad_h}
	DESTINATION ${INSTALL_INC_DIR}/qwwad)
