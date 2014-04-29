#if HAVE_DOXYGEN
# OPTIONAL_DOCS = doc
#endif

ACLOCAL_AMFLAGS=-I m4
#SUBDIRS=$(OPTIONAL_DOCS)

EXTRA_DIST=autogen.sh                     \
	   config/install-sh              \
	   config/depcomp                 \
	   config/missing                 \
	   m4/ax_blas.m4                  \
	   m4/ax_boost_base.m4            \
	   m4/ax_boost_program_options.m4 \
	   m4/ax_lapack.m4

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
LDADD = $(BOOST_LDFLAGS)             \
	$(BOOST_PROGRAM_OPTIONS_LIB)


noinst_HEADERS = src/d0-helpers.h    \
		 src/dos-functions.h \
		 src/dox.h           \
		 src/ef-helpers.h    \
		 src/maths.h         \
		 src/ppff.h          \
		 src/struct.h 

bin_PROGRAMS =  src/csss    \
		src/cszb    \
		src/d0      \
		src/d02D    \
		src/d03D    \
		src/dispersion_relation \
		src/dos     \
		src/ebe     \
		src/efcwire \
		src/efcwwf  \
	       	src/effv    \
		src/efiw    \
		src/efiwire \
		src/efkpsl  \
	       	src/efmfv   \
		src/efpqw   \
		src/efsdot  \
		src/efsdwf  \
	       	src/efshoot \
		src/efsqw   \
	       	src/efsx    \
	       	src/efwf    \
		src/efxm    \
		src/efxv    \
		src/find_heterostructure \
		src/gde     \
		src/hup     \
		src/ivdb    \
		src/ovl     \
		src/ppcd    \
		src/ppee    \
		src/pplb    \
		src/pplbso  \
		src/ppsg    \
		src/ppsl    \
		src/ppvfq   \
		src/pth     \
		src/rlv-fcc \
		src/rlv-sc  \
		src/rlv-ss  \
	       	src/sbp     \
		src/scp     \
		src/sfr     \
		src/slk     \
		src/sradp   \
		src/srcc    \
		src/srelo   \
	       	src/srmpr   \
	       	src/srrad   \
	       	src/tdb     \
	       	src/tsb

d0_helpers    = src/d0-helpers.c src/d0-helpers.h
qwwad_options = src/qwwad-options.cpp src/qwwad-options.h
wf_options    = src/wf_options.cpp src/wf_options.h $(qwwad_options)
qwwad_fileio  = src/qwwad-fileio.cpp src/qwwad-fileio.h

# QCLsim infrastructure (released from Leeds code)
noinst_LIBRARIES = src/libqclsim.a
src_libqclsim_a_SOURCES =  \
			   src/dos-functions.cpp       \
			   src/dos-functions.h         \
			   src/double-barrier.cpp      \
			   src/double-barrier.h        \
			   src/qclsim-constants.h      \
			   src/qclsim-fermi.cpp        \
			   src/qclsim-fermi.h          \
			   src/qclsim-fileio.cpp       \
			   src/qclsim-fileio.h         \
			   src/qclsim-linalg.cpp       \
			   src/qclsim-linalg.h         \
			   src/qclsim-maths.cpp        \
			   src/qclsim-maths.h          \
			   src/qclsim-schroedinger.cpp \
			   src/qclsim-schroedinger.h   \
			   src/qclsim-schroedinger-finite-well.cpp \
			   src/qclsim-schroedinger-finite-well.h   \
			   src/qclsim-subband.cpp      \
			   src/qclsim-subband.h        \
			   src/qwwad-heterostructure.cpp \
			   src/qwwad-heterostructure.h

src_d0_SOURCES   = src/d0.c $(d0_helpers)

src_d02D_SOURCES = src/d02D.c $(d0_helpers)
src_d03D_SOURCES = src/d03D.c $(d0_helpers)

src_dispersion_relation_SOURCES = src/dispersion_relation.cpp $(wf_options)
src_dispersion_relation_LDADD   = src/libqclsim.a $(LDADD)

src_dos_SOURCES  = src/dos.cpp $(qwwad_options) $(qwwad_fileio)
src_dos_LDADD    = src/libqclsim.a $(LDADD)

src_efiw_SOURCES = src/efiw.cpp $(qwwad_options)
src_efiw_LDADD   = src/libqclsim.a $(LDADD)

src_efkpsl_SOURCES = src/efkpsl.cpp $(qwwad_options)
src_efkpsl_LDADD   = src/libqclsim.a $(LDADD)

src_efsqw_SOURCES = src/efsqw.cpp $(qwwad_options)
src_efsqw_LDADD   = src/libqclsim.a $(LDADD)

src_find_heterostructure_SOURCES = src/find_heterostructure.cpp $(qwwad_options)
src_find_heterostructure_LDADD   = src/libqclsim.a $(LDADD)

src_ivdb_SOURCES = src/ivdb.cpp $(qwwad_options)
src_ivdb_LDADD   = src/libqclsim.a $(LDADD)

src_ppvfq_SOURCES = src/ppvfq.c src/ppff.h
src_ppvfq_LDADD   = src/libqclsim.a

src_pplbso_SOURCES = src/pplbso.c src/ppsop.h

src_sbp_SOURCES  = src/sbp.cpp $(qwwad_options) $(qwwad_fileio)
src_sbp_LDADD    = src/libqclsim.a $(LDADD)

src_tdb_SOURCES  = src/tdb.cpp $(qwwad_options) $(qwwad_fileio)
src_tdb_LDADD    = src/libqclsim.a $(LDADD)

src_tsb_SOURCES  = src/tsb.cpp $(qwwad_options) $(qwwad_fileio)
src_tsb_LDADD    = src/libqclsim.a $(LDADD)

# Build rule for installing manual pages
dist_man1_MANS = man/dos.1

man/dos.1: src/dos.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/dos$(EXEEXT)
	help2man --output=man/dos.1 --name="compute density of states in 1D, 2D and 3D systems" src/dos$(EXEEXT)

# Install examples
exdatadir = $(datadir)/qwwad/examples
 
ex2datadir = $(exdatadir)/solutions-to-schroedingers-equation
ex2srcdir = examples/solutions-to-schroedingers-equation

EXTRA_DIST += examples

ex2data_DATA = $(ex2srcdir)/infinite-well-energy-vs-width.sh      \
               $(ex2srcdir)/infinite-well-energy-vs-width.agr     \
               $(ex2srcdir)/infinite-well-energy-vs-width.eps     \
               $(ex2srcdir)/infinite-well-wavefunctions.sh        \
               $(ex2srcdir)/infinite-well-wavefunctions.agr       \
               $(ex2srcdir)/infinite-well-wavefunctions.eps       \
               $(ex2srcdir)/double-barrier-transmission-width.sh  \
               $(ex2srcdir)/double-barrier-transmission-width.agr \
               $(ex2srcdir)/double-barrier-transmission-width.eps

### COMPILE THE GTEST LIBRARY ITSELF ###

EXTRA_DIST += \
  gtest/include/gtest/gtest-param-test.h.pump \
  gtest/include/gtest/internal/gtest-param-util-generated.h.pump \
  gtest/include/gtest/internal/gtest-tuple.h.pump \
  gtest/include/gtest/internal/gtest-type-util.h.pump \
  gtest/src/gtest-death-test.cc \
  gtest/src/gtest-filepath.cc \
  gtest/src/gtest-internal-inl.h \
  gtest/src/gtest-port.cc \
  gtest/src/gtest-printers.cc \
  gtest/src/gtest-test-part.cc \
  gtest/src/gtest-typed-test.cc \
  gtest/src/gtest.cc \
  gtest/include/gtest/internal/gtest-param-util-generated.h \
  gtest/include/gtest/internal/gtest-type-util.h \
  gtest/include/gtest/internal/gtest-param-util.h \
  gtest/include/gtest/internal/gtest-port.h \
  gtest/include/gtest/internal/gtest-linked_ptr.h \
  gtest/include/gtest/internal/gtest-internal.h \
  gtest/include/gtest/internal/gtest-death-test-internal.h \
  gtest/include/gtest/internal/gtest-filepath.h \
  gtest/include/gtest/internal/gtest-string.h \
  gtest/include/gtest/internal/gtest-tuple.h \
  gtest/include/gtest/gtest.h \
  gtest/include/gtest/gtest-typed-test.h \
  gtest/include/gtest/gtest-message.h \
  gtest/include/gtest/gtest-test-part.h \
  gtest/include/gtest/gtest_pred_impl.h \
  gtest/include/gtest/gtest-printers.h \
  gtest/include/gtest/gtest-death-test.h \
  gtest/include/gtest/gtest-param-test.h \
  gtest/include/gtest/gtest-spi.h \
  gtest/include/gtest/gtest_prod.h 

# Modifies compiler and linker flags for pthreads compatibility.
if HAVE_PTHREADS
  GTEST_CXXFLAGS = @PTHREAD_CFLAGS@ -DGTEST_HAS_PTHREAD=1
  GTEST_LIBS = @PTHREAD_LIBS@
else
  GTEST_CXXFLAGS = -DGTEST_HAS_PTHREAD=0
endif

# We define the global AM_CPPFLAGS as everything we compile includes from these
# directories.
AM_CPPFLAGS += -isystem $(srcdir)/gtest/include -I$(srcdir)/gtest

# Build rules for libraries.
check_LTLIBRARIES = gtest/lib/libgtest.la gtest/lib/libgtest_main.la
#check_CPPFLAGS = -I$(srcdir)/gtest/include

gtest_lib_libgtest_la_SOURCES = gtest/src/gtest-all.cc
gtest_lib_libgtest_la_CXXFLAGS = $(GTEST_CXXFLAGS) -Wno-missing-field-initializers
gtest_lib_libgtest_la_LIBADD = $(GTEST_LIBS)

gtest_lib_libgtest_main_la_SOURCES = gtest/src/gtest_main.cc
gtest_lib_libgtest_main_la_CXXFLAGS = $(GTEST_CXXFLAGS)
gtest_lib_libgtest_main_la_LIBADD = gtest/lib/libgtest.la $(GTEST_LIBS)

### TESTS ###
TESTS = tests/tests
check_PROGRAMS = tests/tests
tests_tests_SOURCES = tests/some_tests.cpp
tests_tests_CXXFLAGS = $(GTEST_CXXFLAGS)
tests_tests_LDADD = gtest/lib/libgtest_main.la \
                    gtest/lib/libgtest.la \
	            $(GTEST_LIBS)

