if( VERBOSE )
    message( "src" )
endif()

set( d0_helpers d0-helpers.c d0-helpers.h )
set( dos_funcs dos-functions.cpp dos-functions.h )
set( qwwad_options qwwad-options.cpp qwwad-options.h )
set( wf_options wf_options.cpp wf_options.h ${qwwad_options} )
set( qwwad_fileio qwwad-fileio.cpp qwwad-fileio.h )

#set( qclsim_src
#        qclsim-fermi.cpp
#        qclsim-fileio.cpp
#        qclsim-linalg.cpp
#        qclsim-maths.cpp
#        qclsim-schroedinger.cpp
#        qclsim-subband.cpp
#    )
#
#set( qclsim_h
#        qclsim-constants.h
#        qclsim-fermi.h
#        qclsim-fileio.h
#        qclsim-linalg.h
#        qclsim-maths.h
#        qclsim-schroedinger.h
#        qclsim-subband.h
#    )
#
#add_library( qclsim ${qclsim_src} )

#add_executable( d0 d0.c ${d0_helpers} )
#add_executable( d02D d02D.c ${d0_helpers} )
#add_executable( d03D d03D.c ${d0_helpers} )

#add_executable( dispersion_relation dispersion_relation.cpp ${wf_options} )
#target_link_libraries( dispersion_relation qclsim )
#
#add_executable( dos dos.cpp ${qwwad_options} ${dos_funcs} ${qwwad_fileio} )
#target_link_libraries( dos qclsim )
