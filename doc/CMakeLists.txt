if(VERBOSE)
	message("-doc")
endif()

if(BUILD_MANPAGES)
	# Check that the help2man program is available
	find_program(HELP2MAN help2man DOC "Location of the help2man program")

	# Only build manpages if the required tools are available
	if( NOT HELP2MAN )
		message(WARNING "Unable to find help2man program.  UNIX manpages will not be installed")
	elseif(WIN32)
		message(STATUS "Configuring Windows build.  UNIX manpages will not be installed")
	else()
		message("help2man program found at: ${HELP2MAN}")

		# Store number of programs in the _nprog variable
		list(LENGTH qwwad_programs _nprog)
		math(EXPR _iprog_max "${_nprog} -1")

		foreach(_iprog RANGE ${_iprog_max})
			list(GET qwwad_programs ${_iprog} _progname)
			list(GET qwwad_program_descriptions ${_iprog} _summary)

			set(_MANPAGENAME ${CMAKE_CURRENT_BINARY_DIR}/${_progname}.1)
			set(_MANPAGESRC  ${CMAKE_CURRENT_SOURCE_DIR}/${_progname}.h2m)
			set(_PROGPATH ${CMAKE_BINARY_DIR}/src/${_progname})

			# Only build manpage if the relevant input file exists
			if(EXISTS ${_MANPAGESRC})
				add_custom_command(OUTPUT ${_MANPAGENAME}
					COMMAND ${HELP2MAN} --output=${_MANPAGENAME} --name=${_summary} --include=${_MANPAGESRC} ${_PROGPATH}
					DEPENDS ${_MANPAGESRC} ${_PROGPATH}
					COMMENT "Building manpage for ${_progname}"
					VERBATIM
					)
				add_custom_target(manpage_${_progname} ALL DEPENDS ${_MANPAGENAME})
				install(FILES ${_MANPAGENAME} DESTINATION share/man/man1)
			else()
				message(WARNING "Cannot find manpage input file for ${_progname}")
			endif()
		endforeach(_iprog)
	endif()
endif(BUILD_MANPAGES)
