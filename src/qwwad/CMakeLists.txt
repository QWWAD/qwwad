if( VERBOSE )
    message( "  /qwwad" )
endif()

set( qwwad_src
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
        heterostructure.cpp
        schroedinger-solver.cpp
        schroedinger-solver-donor.cpp
        schroedinger-solver-donor-2D.cpp
        schroedinger-solver-donor-3D.cpp
        schroedinger-solver-donor-variable.cpp
        schroedinger-solver-finite-well.cpp
        schroedinger-solver-full.cpp
        schroedinger-solver-infinite-well.cpp
        schroedinger-solver-taylor.cpp
        schroedinger-solver-tridiagonal.cpp
    )

set( qwwad_h
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
        qwwad-constants.h
        heterostructure.h
        schroedinger-solver.h
        schroedinger-solver-donor.h
        schroedinger-solver-donor-2D.h
        schroedinger-solver-donor-3D.h
        schroedinger-solver-donor-variable.h
        qwwad-schroedinger-finite-well.h
        qwwad-schroedinger-full.h
        qwwad-schroedinger-infinite-well.h
        qwwad-schroedinger-taylor.h
        qwwad-schroedinger-tridiagonal.h
    )

add_library( qwwad ${qwwad_src} )
target_link_libraries( qwwad ${Boost_LIBRARIES} ${GSL_LIBRARIES} )
