if( VERBOSE )
    message( "-src" )
endif()

project( QWWAD )

include_directories( ${ARMADILLO_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} )

set( d0_helpers d0-helpers.c d0-helpers.h )
set( dos_funcs dos-functions.cpp dos-functions.h )
set( qwwad_options qwwad-options.cpp qwwad-options.h )
set( wf_options wf_options.cpp wf_options.h ${qwwad_options} )
set( qwwad_fileio qwwad-fileio.cpp qwwad-fileio.h )

set( qwwad_src
        dos-functions.cpp
        double-barrier.cpp
        qclsim-fermi.cpp
        qclsim-fileio.cpp
        qclsim-linalg.cpp
        qclsim-maths.cpp
        qclsim-subband.cpp
        qwwad-heterostructure.cpp
        qwwad-schroedinger.cpp
        qwwad-schroedinger-finite-well.cpp
        qwwad-schroedinger-full.cpp
        qwwad-schroedinger-infinite-well.cpp
        qwwad-schroedinger-taylor.cpp
        qwwad-schroedinger-tridiagonal.cpp
    )

set( qwwad_h
        dos-functions.h
        double-barrier.h
        qclsim-fermi.h
        qclsim-fileio.h
        qclsim-linalg.h
        qclsim-maths.h
        qclsim-subband.h
        qwwad-constants.h
        qwwad-heterostructure.h
        qwwad-schroedinger.h
        qwwad-schroedinger-finite-well.h
        qwwad-schroedinger-full.h
        qwwad-schroedinger-infinite-well.h
        qwwad-schroedinger-taylor.h
        qwwad-schroedinger-tridiagonal.h
    )

add_library( qwwad ${qwwad_src} )

add_executable( d0 d0.c ${d0_helpers} )
target_link_libraries( d0 ${GSL_LIBRARIES} )
add_executable( d02D d02D.c ${d0_helpers} )
target_link_libraries( d02D ${GSL_LIBRARIES} )
add_executable( d03D d03D.c ${d0_helpers} )
target_link_libraries( d03D ${GSL_LIBRARIES} )

if( NOT APPLE )
add_executable( dispersion_relation dispersion_relation.cpp ${wf_options} )
target_link_libraries( dispersion_relation qwwad ${Boost_LIBRARIES}
	${GSL_LIBRARIES} )

add_executable( dos dos.cpp ${qwwad_options} ${dos_funcs} ${qwwad_fileio} )
target_link_libraries( dos qwwad ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
endif()
