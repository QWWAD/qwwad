if( VERBOSE )
    message( "-src" )
endif()

include_directories( ${ARMADILLO_INCLUDE_DIRS} ${GSL_INCLUDE_DIR} )

set( wf_options wf_options.cpp wf_options.h )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

add_subdirectory( qwwad )

# Loop through all the programs in the qwwad_program list
# and specify their build, link & install "recipes"
foreach(_prog IN LISTS qwwad_programs)
	set(_CPPSRC ${_prog}.cpp)
	add_executable(${_prog} ${_CPPSRC})
	target_link_libraries(${_prog} libqwwad)
	install(TARGETS ${_prog} RUNTIME DESTINATION bin )
endforeach()

add_library(gtest SHARED googletest/googletest/src/gtest-all.cc)
add_library(gtest_main STATIC googletest/googletest/src/gtest_main.cc)
target_include_directories(gtest PRIVATE googletest/googletest)
