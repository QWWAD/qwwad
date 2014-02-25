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

noinst_HEADERS = src/const.h         \
		 src/d0-helpers.h    \
		 src/dos-functions.h \
		 src/dox.h           \
		 src/ef-helpers.h    \
		 src/maths.h         \
		 src/ppff.h          \
		 src/qclsim-constants.h \
		 src/qclsim-schroedinger.h \
		 src/struct.h 

bin_PROGRAMS =  src/csss    \
		src/cszb    \
		src/d0      \
		src/d02D    \
		src/d03D    \
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
		src/gde     \
		src/hup     \
		src/ivdb    \
		src/ovl     \
		src/ppcd    \
		src/ppee    \
		src/ppsg    \
		src/ppvfq   \
		src/pth     \
		src/rlv-fcc \
		src/rlv-sc  \
		src/rlv-ss  \
	       	src/sbp     \
		src/scp     \
		src/sfr     \
		src/slk     \
		src/srcc    \
		src/srelo   \
	       	src/srmpr   \
	       	src/srrad   \
	       	src/tdb     \
	       	src/tsb

d0_helpers    = src/d0-helpers.c src/d0-helpers.h
dos_funcs     = src/dos-functions.cpp src/dos-functions.h
qwwad_options = src/qwwad-options.cpp src/qwwad-options.h

# QCLsim infrastructure (released from Leeds code)
noinst_LIBRARIES = src/libqclsim.a
src_libqclsim_a_SOURCES =  \
			   src/qclsim-fileio.cpp       \
			   src/qclsim-fileio.h         \
			   src/qclsim-linalg.cpp       \
			   src/qclsim-linalg.h         \
			   src/qclsim-maths.cpp        \
			   src/qclsim-maths.h          \
			   src/qclsim-schroedinger.cpp \
			   src/qclsim-schroedinger.h

src_d0_SOURCES   = src/d0.c $(d0_helpers)
src_d02D_SOURCES = src/d02D.c $(d0_helpers)
src_d03D_SOURCES = src/d03D.c $(d0_helpers)

src_dos_SOURCES  = src/dos.cpp $(qwwad_options) $(dos_funcs)
src_dos_LDADD    = src/libqclsim.a $(LDADD)

src_efiw_SOURCES = src/efiw.cpp $(qwwad_options)
src_efiw_LDADD   = src/libqclsim.a $(LDADD)

src_sbp_SOURCES  = src/sbp.cpp $(qwwad_options)
src_sbp_LDADD    = src/libqclsim.a $(LDADD)
