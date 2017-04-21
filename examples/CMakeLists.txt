if(VERBOSE)
	message("-examples")
endif()

# Define a list of all the examples subfolders
set( qwwad_example_folders
	schroedinger-solutions
	numerical-solutions
	diffusion
	impurities
	excitons
	simple-wires-and-dots
	carrier-scattering
	#empirical-pseudopotential
	#pseudopotential-het
	#waveguide
	)

foreach(_folder IN LISTS qwwad_example_folders)
	# Define installation path
	set(_SS_INSTALL "${QWWAD_SHARE_INSTALL}/examples/${_folder}")

	message("Installing scripts to ${_SS_INSTALL}")

	# Install index file
	install(FILES ${_folder}/index DESTINATION ${_SS_INSTALL})

	# Find and install all scripts
	file(GLOB _SS_SCRIPTS "${_folder}/*.sh")
	install(PROGRAMS ${_SS_SCRIPTS} DESTINATION ${_SS_INSTALL})
endforeach()


