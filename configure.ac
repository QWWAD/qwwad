#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([qwwad], [1.0], [a.valavanis@leeds.ac.uk])
AC_CONFIG_SRCDIR([src/csss.c])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

# The number in the options list here refers to the minimum version
# of Automake that can be accepted
AM_INIT_AUTOMAKE([-W -Wall 1.11])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor modf sqrt])

AC_CONFIG_FILES([Makefile
		 config/Makefile
		 examples/Makefile
		 m4/Makefile
		 src/Makefile])
AC_OUTPUT
