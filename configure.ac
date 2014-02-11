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

dnl Quieten the Make output
AM_SILENT_RULES([yes])

# Checks for programs.
AC_PROG_CC

dnl Enable some compiler flags by default
CFLAGS="-Wall $CFLAGS"

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

dnl Check for a CBLAS implementation
AC_SEARCH_LIBS([cblas_dgemm], [cblas gslcblas],,
	       AC_MSG_ERROR([Cannot find a CBLAS implementation]))

dnl Check for the GSL library
AC_CHECK_LIB([gsl], [gsl_blas_dgemm],,
	     AC_MSG_ERROR([Cannot locate the GSL library.  Please
			  install it and try again.]))

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h strings.h complex.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor modf sqrt])

AC_CONFIG_FILES([Makefile
		 config/Makefile
		 examples/Makefile
		 m4/Makefile
		 src/Makefile])
AC_OUTPUT
