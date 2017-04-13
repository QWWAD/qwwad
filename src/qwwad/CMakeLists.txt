if( VERBOSE )
    message( "  /qwwad" )
endif()

# Convenience macro to add library modules
macro(add_libqwwad_module modname)
	list(APPEND qwwad_src ${modname}.cpp)
	list(APPEND qwwad_h   ${modname}.h)
endmacro()

add_libqwwad_module(data-checker)
add_libqwwad_module(debye)
add_libqwwad_module(donor-energy-minimiser)
add_libqwwad_module(donor-energy-minimiser-fast)
add_libqwwad_module(donor-energy-minimiser-linear)
add_libqwwad_module(dos-functions)
add_libqwwad_module(double-barrier)
add_libqwwad_module(eigenstate)
add_libqwwad_module(fermi)
add_libqwwad_module(file-io)
add_libqwwad_module(file-io-deprecated)
add_libqwwad_module(intersubband-transition)
add_libqwwad_module(linear-algebra)
add_libqwwad_module(material)
add_libqwwad_module(material-library)
add_libqwwad_module(material-property)
add_libqwwad_module(material-property-constant)
add_libqwwad_module(material-property-interp)
add_libqwwad_module(material-property-numeric)
add_libqwwad_module(material-property-poly)
add_libqwwad_module(material-property-string)
add_libqwwad_module(maths-helpers)
add_libqwwad_module(mesh)
add_libqwwad_module(options)
add_libqwwad_module(poisson-solver)
add_libqwwad_module(ppff)
add_libqwwad_module(pplb-functions)
add_libqwwad_module(ppsop)
add_libqwwad_module(subband)
add_libqwwad_module(scattering-calculator-LO)
add_libqwwad_module(schroedinger-solver)
add_libqwwad_module(schroedinger-solver-donor)
add_libqwwad_module(schroedinger-solver-donor-2D)
add_libqwwad_module(schroedinger-solver-donor-3D)
add_libqwwad_module(schroedinger-solver-donor-variable)
add_libqwwad_module(schroedinger-solver-finite-well)
add_libqwwad_module(schroedinger-solver-full)
add_libqwwad_module(schroedinger-solver-infinite-well)
add_libqwwad_module(schroedinger-solver-kronig-penney)
add_libqwwad_module(schroedinger-solver-poeschl-teller)
add_libqwwad_module(schroedinger-solver-shooting)
add_libqwwad_module(schroedinger-solver-taylor)
add_libqwwad_module(schroedinger-solver-tridiagonal)
add_libqwwad_module(wf_options)

add_library( libqwwad SHARED ${qwwad_src} ${qwwad_h} )
set_target_properties( libqwwad
	               PROPERTIES
		       OUTPUT_NAME ${PROJECT_NAME}
		       SOVERSION   1
		       VERSION     ${qwwad_VERSION}
		     )

target_link_libraries( libqwwad
	${GSL_LIBRARIES}
	${Boost_LIBRARIES}
	${LAPACK_LIBRARIES}
	${ARMADILLO_LIBRARIES}
	${LIBXMLPP_LIBRARIES} )

# Install the shared QWWAD library
install(TARGETS libqwwad
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib)

# Install the QWWAD API headers in the "qwwad" subfolder
install(FILES ${qwwad_h}
	DESTINATION ${INSTALL_INC_DIR}/qwwad)
