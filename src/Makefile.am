# Makefile.am - Automake input file for QWWAD project (source code)

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
LDADD = $(BOOST_LDFLAGS)             \
	$(BOOST_PROGRAM_OPTIONS_LIB)

noinst_HEADERS = const.h      \
		 d0-helpers.h \
		 dox.h        \
		 ef-helpers.h \
		 maths.h      \
		 ppff.h       \
		 qclsim-constants.h \
		 qclsim-schroedinger.h \
		 struct.h 

bin_PROGRAMS =  csss    \
		cszb    \
		d0      \
		d02D    \
		d03D    \
		dos     \
		ebe     \
		efcwire \
		efcwwf  \
	       	effv    \
		efiw    \
		efiwire \
		efkpsl  \
	       	efmfv   \
		efpqw   \
		efsdot  \
		efsdwf  \
	       	efshoot \
		efsqw   \
	       	efsx    \
	       	efwf    \
		efxm    \
		efxv    \
		gde     \
		hup     \
		ivdb    \
		ovl     \
		ppcd    \
		ppee    \
		ppsg    \
		ppvfq   \
		pth     \
		rlv-fcc \
		rlv-sc  \
		rlv-ss  \
	       	sbp     \
		scp     \
		sfr     \
		slk     \
		srcc    \
		srelo   \
	       	srmpr   \
	       	srrad   \
	       	tdb     \
	       	tsb

d0_helpers    = d0-helpers.c d0-helpers.h
creaddata     = creaddata.cpp creaddata.h
qwwad_options = qwwad-options.cpp qwwad-options.h

# QCLsim infrastructure (released from Leeds code)
qclsim_maths        = qclsim-maths.cpp qclsim-maths.h
qclsim_linalg       = qclsim-linalg.cpp qclsim-linalg.h $(qclsim_maths)
qclsim_schroedinger = qclsim-schroedinger.cpp qclsim-schroedinger.h

d0_SOURCES   = d0.c $(d0_helpers)
d02D_SOURCES = d02D.c $(d0_helpers)
d03D_SOURCES = d03D.c $(d0_helpers)
dos_SOURCES  = dos.cpp
efiw_SOURCES = efiw.cpp qclsim-schroedinger.cpp $(qwwad_options) $(qclsim_linalg) $(qclsim_schoedinger)
