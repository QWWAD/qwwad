if(VERBOSE)
	message("-doc")
endif()

# Check that the help2man program is available
find_program(HELP2MAN help2man DOC "Location of the help2man program")

# Only build manpages if the required tools are available
if( NOT HELP2MAN )
	message(WARNING "Unable to find help2man program.  UNIX manpages will not be installed")
else()
	message("help2man program found at: ${HELP2MAN}")

	# Define a macro to make it easier to generate manpages
	# - progname is the name of the program
	# - summary is a 1-line summary of the program behaviour
	function(add_qwwad_manpage progname summary)
		set(_MANPAGENAME ${CMAKE_CURRENT_BINARY_DIR}/${progname}.1)
		set(_MANPAGESRC  ${CMAKE_CURRENT_SOURCE_DIR}/${progname}.h2m)
		set(_PROGPATH ${CMAKE_BINARY_DIR}/src/${progname})

		add_custom_command(OUTPUT ${_MANPAGENAME}
			COMMAND ${HELP2MAN} --output=${_MANPAGENAME} --name=${summary} --include=${_MANPAGESRC} ${_PROGPATH}
			DEPENDS ${_MANPAGESRC} ${_PROGPATH}
			COMMENT "Building manpage for ${progname}"
			VERBATIM
			)
		add_custom_target(manpage_${progname} ALL DEPENDS ${_MANPAGENAME})
		install(FILES ${_MANPAGENAME} DESTINATION share/man/man1)
	endfunction(add_qwwad_manpage)

	# List all manpages to be installed here, along with descriptions
	add_qwwad_manpage(qwwad_density_of_states  "compute density of states in 1D, 2D and 3D systems")
	add_qwwad_manpage(qwwad_ef_dispersion      "compute the dispersion relation for a set of energy subbands")
	add_qwwad_manpage(qwwad_ef_donor_generic   "find eigenstates for a donor in a heterostructure (generic search)")
	add_qwwad_manpage(qwwad_ef_infinite_well   "compute eigenstates in an infinite quantum well")
	add_qwwad_manpage(qwwad_ef_plot            "translate wavefunction data into plottable form")
	add_qwwad_manpage(qwwad_ef_square_well     "compute eigenstates in a finite square quantum well")
	add_qwwad_manpage(qwwad_ef_superlattice    "compute eigenstates of a Kronig-Penney superlattice")
	add_qwwad_manpage(qwwad_fermi_distribution "find Fermi-Dirac distributions for a set of subbands")
endif()
