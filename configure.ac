#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

AC_INIT([Quantum Wells, Wires and Dots], [1.0+devel], [https://bugs.launchpad.net/qwwad], [qwwad], [https://launchpad.net/qwwad])

dnl If we have an old version of Autoconf, then we need to set the PACKAGE_URL variable
dnl manually.
dnl TODO: Get rid of this and bump the version to 2.64 in AC_PREREQ when all reasonable
dnl       OSes have Autoconf >= 2.64 available
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]),[2.64]), [-1],
  [AC_MSG_WARN([Using an old version of Autoconf])
   AC_DEFINE([PACKAGE_BUGREPORT], ["https://bugs.launchpad.net/qwwad"], [The address where bugs for the package should be sent])
   AC_DEFINE([PACKAGE_URL], ["https://launchpad.net/qwwad"], [The URL for the package])
  ]
)

AC_CONFIG_SRCDIR([src/csss.c])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])

# The number in the options list here refers to the minimum version
# of Automake that can be accepted
AM_INIT_AUTOMAKE([-W -Wall 1.11 subdir-objects])

dnl Quieten the Make output
AM_SILENT_RULES([yes])

AC_LANG([C++])

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AM_PROG_CC_C_O
AC_PROG_CXX

dnl Libtool needs access to the "ar" archiver so check that it exists.
dnl TODO: This macro is new in Automake 1.12 so we can apply this
dnl       check by default once we update our Automake dependency.
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Set up libtool
AC_PROG_LIBTOOL
LT_PREREQ([2.2])
LT_INIT

dnl Enable some compiler flags by default
CFLAGS="-Wall $CFLAGS"
CXXFLAGS="-Wall $CXXFLAGS"

dnl Set up very strict build for development purposes
CPPFLAGS="-W -Wall -Werror"

# Checks for libraries.
AC_CHECK_LIB([m], [pow])

# Check for Boost Program Options

# This hack automatically checks for boost in the BOOST_HOME directory if
# it has been defined.  Care is probably needed to ensure that this doesn't
# cause unwanted side-effects.
if [[ "x$BOOST_HOME" != "x" ]]; then
	ac_boost_path=${BOOST_HOME}
fi

AX_BOOST_BASE([1.35],,
	      AC_MSG_ERROR([Boost C++ extensions unavailable.  Please install and try again]))

AX_BOOST_PROGRAM_OPTIONS

# Check for linear algebra libraries
AX_BLAS(,AC_MSG_ERROR([BLAS library unavailable.  Please install it and try again]))
AX_LAPACK(,AC_MSG_ERROR([LAPACK library unavailable.  Please install it and try again]))
LIBS="$BLAS_LIBS $LAPACK_LIBS $LIBS"

dnl Check for GNU Scientific Library
PKG_CHECK_MODULES(GSL, gsl)
CPPFLAGS="$GSL_CFLAGS $CPPFLAGS"
LIBS="$GSL_LIBS $LIBS"

# Check if we have the new "official" C header for LAPACK.  Also, note that
# there is a bug in LAPACK that causes incorrect linkage with C++11. See
# http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=4221
AC_CHECK_HEADER([lapacke.h],
		AC_CHECK_LIB([lapacke], [LAPACKE_dlamch], with_lapacke=yes, with_lapacke=no),
		with_lapacke=no)

if test "x$with_lapacke" = "xyes"; then
	AC_MSG_RESULT([Using official C bindings for LAPACK])
	AC_DEFINE(HAVE_LAPACKE,1,[Build with lapacke.h])
	LIBS="$LIBS -llapacke"
else
	AC_MSG_WARN([Official C bindings for LAPACK unavailable.  Using our own bindings instead!])
fi

dnl Check for an armadillo implementation
AC_CHECK_HEADER([armadillo])
AC_MSG_CHECKING([whether we can link to the Armadillo library])
LIBS="$LIBS -larmadillo"
AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <armadillo>]],
		                [[arma::mat A = arma::randu<arma::mat>(4,5)]])],
               [have_armadillo=yes], [have_armadillo=no])
AC_MSG_RESULT([$have_armadillo])

if test "x$have_armadillo" = "xno"; then
	AC_MSG_ERROR([Armadillo library is not usable])
fi

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h strings.h complex.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Configure pthreads.
have_pthreads=no
AX_PTHREAD([have_pthreads=yes], [have_pthreads=no])
AM_CONDITIONAL([HAVE_PTHREADS],[test "x$have_pthreads" = "xyes"])
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_LIBS)

# Checks for library functions.
AC_CHECK_FUNCS([floor modf sqrt])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
