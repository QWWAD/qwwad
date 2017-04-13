# CMake macros for setting up the 'make check' behavour for Makefiles.
# This uses the CTest program (part of CMake) as a driver for all the tests.
# Doing so adds a 'RUN_TESTS' project for IDEs (e.g. Visual Studio) or a 'test'
# project in Makefiles (run with 'make test') which will run CTest.
# (Note that these projects require the tests to be build before they are built
#  as the test projects themselves are not dependencies.)
# The test may also be run by running the CTest program manually from within the
# build directory.
# The 'check' project is added manually using the add_custom_target command which
# is essentially just a wrapper to run the CTest program.
# Unlike the RUN_TESTS/test projects however, the test projects themselves are
# added as dependencies to the check project such that running 'make check'
# before the test executables are run will cause the relevent projects to be built.

macro( add_qwwad_test test_name )
	add_executable(${test_name} ${test_name}.cpp)
	target_link_libraries(${test_name} libqwwad gtest gtest_main)
	add_test(NAME ${test_name} COMMAND ${test_name})
	add_dependencies( check ${test_name} )
endmacro()
