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
	function(add_manpage progname summary)
		message("Generating manpage for ${progname}")
		message("Using summary: ${summary}")

		set(_MANPAGENAME ${CMAKE_CURRENT_BINARY_DIR}/${progname}.1)
		set(_MANPAGESRC  ${CMAKE_CURRENT_SOURCE_DIR}/${progname}.h2m)
		set(_PROGPATH ${CMAKE_BINARY_DIR}/src/${progname})
		message("Manpage name: ${_MANPAGENAME}")
		message("Program path: ${_PROGPATH}")

		add_custom_command(OUTPUT ${_MANPAGENAME}
			COMMAND ${HELP2MAN} --output=${_MANPAGENAME} --name=${summary} --include=${_MANPAGESRC} ${_PROGPATH}
			DEPENDS ${_MANPAGESRC} ${_PROGPATH}
			COMMENT "Building manpage for ${progname}"
			VERBATIM
			)
		add_custom_target(manpage_${progname} ALL DEPENDS ${_MANPAGENAME})
		install(FILES ${_MANPAGENAME} DESTINATION share/man/man1)
	endfunction(add_manpage)

	# List all manpages to be installed here, along with descriptions
	add_manpage("qwwad_ef_infinite_well" "compute eigenstates in an infinite quantum well")
	add_manpage("qwwad_ef_donor_generic" "find state of electron orbiting a donor ion in a heterostructure using a generic search")
endif()
