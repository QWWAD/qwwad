#if HAVE_DOXYGEN
# OPTIONAL_DOCS = doc
#endif

ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS  = $(GSL_CFLAGS)
AM_LDFLAGS   = $(GSL_LIBS) \
               $(BOOST_LDFLAGS) \
               $(BOOST_PROGRAM_OPTIONS_LIB)

EXTRA_DIST=autogen.sh                     \
	   config/install-sh              \
	   config/depcomp                 \
	   config/missing                 \
	   m4/ax_blas.m4                  \
	   m4/ax_boost_base.m4            \
	   m4/ax_boost_program_options.m4 \
	   m4/ax_lapack.m4

noinst_HEADERS = src/dos-functions.h \
		 src/ef-helpers.h    \
		 src/maths.h         \
		 src/ppff.h          \
		 src/struct.h 

bin_PROGRAMS =  chargedensity \
                csss    \
		cszb    \
		densityinput \
		dispersion_relation \
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
		efsqw   \
		efss    \
	       	efsx    \
		efxv    \
		find_heterostructure \
		find_poisson_potential \
		gde     \
		hup     \
		i0      \
		ivdb    \
		ppcd    \
		ppee    \
		ppsg    \
		pth     \
		qwwad_find_donor_state \
		rlv-fcc \
		rlv-sc  \
		rlv-ss  \
	       	sbp     \
		sfr     \
		slk     \
		sradp   \
		srcc    \
		srelo   \
	       	srmpr   \
	       	srrad   \
	       	tdb     \
	       	tsb     \
		wfplot

# Build most functionality into library infrastructure
lib_LTLIBRARIES = libqwwad.la
libqwwad_la_SOURCES =  \
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
			   src/qclsim_poisson_solver.cpp \
			   src/qclsim_poisson_solver.h   \
			   src/qclsim-subband.cpp      \
			   src/qclsim-subband.h        \
			   src/qwwad-donor-energy-minimiser.cpp \
			   src/qwwad-donor-energy-minimiser.h   \
			   src/qwwad-heterostructure.cpp \
			   src/qwwad-heterostructure.h   \
			   src/qwwad-fileio.cpp          \
			   src/qwwad-fileio.h            \
			   src/qwwad-options.cpp         \
			   src/qwwad-options.h           \
			   src/qwwad-schroedinger.cpp    \
			   src/qwwad-schroedinger.h      \
			   src/qwwad-schroedinger-donor.h      \
			   src/qwwad-schroedinger-donor.cpp    \
			   src/qwwad-schroedinger-donor-2D.h      \
			   src/qwwad-schroedinger-donor-2D.cpp    \
			   src/qwwad-schroedinger-donor-3D.h      \
			   src/qwwad-schroedinger-donor-3D.cpp    \
			   src/qwwad-schroedinger-donor-variable.h   \
			   src/qwwad-schroedinger-donor-variable.cpp \
			   src/qwwad-schroedinger-finite-well.cpp \
			   src/qwwad-schroedinger-finite-well.h   \
			   src/qwwad-schroedinger-full.h          \
			   src/qwwad-schroedinger-full.cpp        \
			   src/qwwad-schroedinger-infinite-well.cpp \
			   src/qwwad-schroedinger-infinite-well.h   \
			   src/qwwad-schroedinger-poeschl-teller.cpp \
			   src/qwwad-schroedinger-poeschl-teller.h   \
			   src/qwwad-schroedinger-shooting.cpp      \
			   src/qwwad-schroedinger-shooting.h        \
			   src/qwwad-schroedinger-taylor.cpp        \
			   src/qwwad-schroedinger-taylor.h          \
			   src/qwwad-schroedinger-tridiagonal.cpp   \
			   src/qwwad-schroedinger-tridiagonal.h     \
			   src/wf_options.cpp                       \
			   src/wf_options.h

libqwwad_la_CFLAGS = $(GSL_CFLAGS) \
		     $(BOOST_CPPFLAGS)

libqwwad_la_LDFLAGS = $(BOOST_LDFLAGS)

libqwwad_la_LIBADD = $(GSL_LIBS) \
                     $(BOOST_PROGRAM_OPTIONS_LIB) \
		     $(LAPACK_LIBS) \
		     $(BLAS_LIBS)   \
		     $(FLIBS)

chargedensity_SOURCES = src/chargedensity.cpp
chargedensity_LDADD   = libqwwad.la

csss_SOURCES = src/csss.c

cszb_SOURCES = src/cszb.c

densityinput_SOURCES = src/densityinput.cpp
densityinput_LDADD   = libqwwad.la

dispersion_relation_SOURCES = src/dispersion_relation.cpp
dispersion_relation_LDADD   = libqwwad.la

dos_SOURCES  = src/dos.cpp
dos_LDADD    = libqwwad.la

ebe_SOURCES = src/ebe.c

efcwire_SOURCES = src/efcwire.c

efcwwf_SOURCES = src/efcwwf.c

effv_SOURCES = src/effv.c

efiw_SOURCES = src/efiw.cpp
efiw_LDADD   = libqwwad.la

efiwire_SOURCES = src/efiwire.c

efkpsl_SOURCES = src/efkpsl.cpp
efkpsl_LDADD   = libqwwad.la

efmfv_SOURCES = src/efmfv.cpp
efmfv_LDADD   = libqwwad.la

efpqw_SOURCES = src/efpqw.cpp
efpqw_LDADD   = libqwwad.la

efsdot_SOURCES = src/efsdot.c

efsdwf_SOURCES = src/efsdwf.c

efss_SOURCES = src/efss.cpp
efss_LDADD   = libqwwad.la

efsqw_SOURCES = src/efsqw.cpp
efsqw_LDADD   = libqwwad.la

efsx_SOURCES = src/efsx.c

efxv_SOURCES = src/efxv.cpp
efxv_LDADD   = libqwwad.la

find_heterostructure_SOURCES = src/find_heterostructure.cpp
find_heterostructure_LDADD   = libqwwad.la

find_poisson_potential_SOURCES = src/find_poisson_potential.cpp
find_poisson_potential_LDADD   = libqwwad.la

gde_SOURCES = src/gde.cpp
gde_LDADD   = libqwwad.la

hup_SOURCES = src/hup.cpp
hup_LDADD   = libqwwad.la

i0_SOURCES  = src/i0.cpp
i0_LDADD    = libqwwad.la

ivdb_SOURCES = src/ivdb.cpp
ivdb_LDADD   = libqwwad.la

ppcd_SOURCES = src/ppcd.c
ppcd_CFLAGS  = $(GSL_CFLAGS)
ppcd_LDADD   = $(GSL_LIBS)

ppee_SOURCES = src/ppee.c

#pplb_SOURCES  = src/pplb.c
#pplb_CPPFLAGS = $(PTHREAD_CFLAGS)
#pplb_LDADD    = $(PTHREAD_LIBS)

ppsg_SOURCES = src/ppsg.c

# ppsl_SOURCES = src/ppsl.c
# ppsl_CFLAGS = $(GSL_CFLAGS)
# ppsl_LDADD  = $(GSL_LIBS)

# ppvfq_SOURCES = src/ppvfq.c src/ppff.h
# ppvfq_LDADD   = libqwwad.la

# pplbso_SOURCES = src/pplbso.c src/ppsop.h

pth_SOURCES = src/pth.cpp
pth_LDADD   = libqwwad.la

qwwad_find_donor_state_SOURCES = src/qwwad_find_donor_state.cpp
qwwad_find_donor_state_LDADD   = libqwwad.la

rlv_fcc_SOURCES = src/rlv-fcc.c

rlv_sc_SOURCES = src/rlv-sc.c

rlv_ss_SOURCES = src/rlv-ss.c

sbp_SOURCES  = src/sbp.cpp
sbp_LDADD    = libqwwad.la

sfr_SOURCES  = src/sfr.c

slk_SOURCES  = src/slk.c

sradp_SOURCES = src/sradp.c
sradp_CFLAGS  = $(GSL_CFLAGS)
sradp_LDADD   = $(GSL_LIBS)

srcc_SOURCES  = src/srcc.c

srelo_SOURCES  = src/srelo.c

srmpr_SOURCES  = src/srmpr.c

srrad_SOURCES  = src/srrad.c

tdb_SOURCES  = src/tdb.cpp
tdb_LDADD    = libqwwad.la

tsb_SOURCES  = src/tsb.cpp
tsb_LDADD    = libqwwad.la

wfplot_SOURCES = src/wfplot.cpp
wfplot_LDADD   = libqwwad.la

# Build rule for installing manual pages
dist_man1_MANS = man/dispersion_relation.1 \
		 man/dos.1                 \
		 man/efiw.1                \
		 man/sbp.1                 \
		 man/efsqw.1               \
		 man/wfplot.1              \
		 man/efkpsl.1              \
		 man/tsb.1                 \
		 man/tdb.1                 \
		 man/ivdb.1                \
		 man/find_heterostructure.1 \
		 man/efxv.1                 \
		 man/efss.1                 \
		 man/efpqw.1                \
		 man/pth.1                  \
		 man/qwwad_find_donor_state.1 \
		 man/find_poisson_potential.1 \
		 man/hup.1                    \
		 man/densityinput.1           \
		 man/chargedensity.1          \
		 man/gde.1

EXTRA_DIST += doc/find_heterostructure.h2m   \
	      doc/efxv.h2m                   \
	      doc/efss.h2m                   \
	      doc/efpqw.h2m                  \
	      doc/pth.h2m                    \
	      doc/find_poisson_potential.h2m \
              doc/efiw.h2m                   \
	      doc/dispersion_relation.h2m    \
	      doc/dos.h2m                    \
	      doc/sbp.h2m                    \
	      doc/efsqw.h2m                  \
	      doc/wfplot.h2m                 \
	      doc/efkpsl.h2m                 \
	      doc/tsb.h2m                    \
	      doc/tdb.h2m                    \
	      doc/ivdb.h2m                   \
	      doc/hup.h2m                    \
	      doc/densityinput.h2m           \
	      doc/chargedensity.h2m          \
	      doc/gde.h2m                    \
	      doc/qwwad_find_donor_state.h2m 

man/dos.1: src/dos.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/dos.h2m | dos$(EXEEXT) 
	$(MKDIR_P) man
	help2man --output=$@ --name="compute density of states in 1D, 2D and 3D systems" --include $(top_srcdir)/doc/dos.h2m ./dos$(EXEEXT)

man/efiw.1: src/efiw.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/efiw.h2m | efiw$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates in an infinite quantum well" --include $(top_srcdir)/doc/efiw.h2m ./efiw$(EXEEXT)

man/dispersion_relation.1: src/dispersion_relation.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/dispersion_relation.h2m | dispersion_relation$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute the dispersion relation for an energy subband" --include $(top_srcdir)/doc/dispersion_relation.h2m ./dispersion_relation$(EXEEXT)

man/sbp.1: src/sbp.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/sbp.h2m | sbp$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find Fermi-Dirac distribution functions for a set of subbands" --include $(top_srcdir)/doc/sbp.h2m ./sbp$(EXEEXT)

man/efsqw.1: src/efsqw.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/efsqw.h2m | efsqw$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates of a finite square quantum well" --include $(top_srcdir)/doc/efsqw.h2m ./efsqw$(EXEEXT)

man/wfplot.1: src/wfplot.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/wfplot.h2m | wfplot$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="translate wavefunction data into a prettier plottable form." --include $(top_srcdir)/doc/wfplot.h2m ./wfplot$(EXEEXT)

man/efkpsl.1: src/efkpsl.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/efkpsl.h2m | efkpsl$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates of a Kronig-Penney superlattice." --include $(top_srcdir)/doc/efkpsl.h2m ./efkpsl$(EXEEXT)

man/tsb.1: src/tsb.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/tsb.h2m | tsb$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute transmission coefficient through single tunnelling barrier." --include $(top_srcdir)/doc/tsb.h2m ./tsb$(EXEEXT) 

man/tdb.1: src/tdb.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/tdb.h2m | tdb$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute transmission coefficient through double tunnelling barrier." --include $(top_srcdir)/doc/tdb.h2m ./tdb$(EXEEXT)

man/ivdb.1: src/ivdb.cpp $(top_srcdir)/configure.ac $(top_srcdir)/doc/ivdb.h2m | ivdb$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute current--voltage characteristics for a double tunnelling barrier." --include $(top_srcdir)/doc/ivdb.h2m ./ivdb$(EXEEXT)

man/find_heterostructure.1: src/find_heterostructure.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/find_heterostructure.h2m | find_heterostructure$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate spatial mesh and output alloy & doping profiles." --include $(top_srcdir)/doc/find_heterostructure.h2m ./find_heterostructure$(EXEEXT)

man/efxv.1: src/efxv.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/efxv.h2m | efxv$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="tabulate material parameters for heterostructure." --include $(top_srcdir)/doc/efxv.h2m ./efxv$(EXEEXT)

man/efss.1: src/efss.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/efss.h2m | efss$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="solve the Schroedinger equation numerically." --include $(top_srcdir)/doc/efss.h2m ./efss$(EXEEXT)

man/efpqw.1: src/efpqw.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/efpqw.h2m | efpqw$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate a parabolic alloy profile." --include $(top_srcdir)/doc/efpqw.h2m ./efpqw$(EXEEXT)

man/pth.1: src/pth.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/pth.h2m | pth$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates of a Poeschl-Teller potential hole." --include $(top_srcdir)/doc/pth.h2m ./pth$(EXEEXT) 

man/find_poisson_potential.1: src/find_poisson_potential.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/find_poisson_potential.h2m | find_poisson_potential$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find the Poisson potential induced by a charge profile." --include $(top_srcdir)/doc/find_poisson_potential.h2m ./find_poisson_potential$(EXEEXT) 

man/hup.1: src/hup.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/hup.h2m | hup$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find the uncertainty product for position and momentum." --include $(top_srcdir)/doc/hup.h2m ./hup$(EXEEXT) 

man/densityinput.1: src/densityinput.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/densityinput.h2m | densityinput$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate an appoximate carrier energy distribution." --include $(top_srcdir)/doc/densityinput.h2m ./densityinput$(EXEEXT) 

man/chargedensity.1: src/chargedensity.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/chargedensity.h2m | chargedensity$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate a charge density profile for a 2D system." --include $(top_srcdir)/doc/chargedensity.h2m ./chargedensity$(EXEEXT) 

man/gde.1: src/gde.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/gde.h2m | gde$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate a 1D diffusion profile." --include $(top_srcdir)/doc/gde.h2m ./gde$(EXEEXT) 

man/qwwad_find_donor_state.1: src/qwwad_find_donor_state.cpp $(top_srcdir)/configure.ac \
    $(top_srcdir)/doc/qwwad_find_donor_state.h2m | qwwad_find_donor_state$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find ground state of electron orbiting a donor ion in a heterostructure." --include $(top_srcdir)/doc/qwwad_find_donor_state.h2m ./qwwad_find_donor_state$(EXEEXT) 

# Install examples
exdatadir = $(datadir)/qwwad/examples

exdata_DATA = examples/README \
	      examples/template.sh
 
ex2datadir = $(exdatadir)/schroedinger-solutions
ex2srcdir = examples/schroedinger-solutions

EXTRA_DIST += examples

ex2data_DATA = $(ex2srcdir)/index                                 \
               $(ex2srcdir)/infinite-well-energy-vs-width.sh      \
               $(ex2srcdir)/infinite-well-energy-vs-width.agr     \
               $(ex2srcdir)/infinite-well-energy-vs-width.eps     \
               $(ex2srcdir)/infinite-well-wavefunctions.sh        \
               $(ex2srcdir)/infinite-well-wavefunctions.agr       \
               $(ex2srcdir)/infinite-well-wavefunctions.eps       \
	       $(ex2srcdir)/non-parabolic-dispersion.sh           \
	       $(ex2srcdir)/non-parabolic-dispersion.agr          \
	       $(ex2srcdir)/non-parabolic-dispersion.eps          \
	       $(ex2srcdir)/density-of-states.sh                  \
	       $(ex2srcdir)/density-of-states-bulk.agr            \
	       $(ex2srcdir)/density-of-states-bulk.eps            \
	       $(ex2srcdir)/density-of-states-2D.agr              \
	       $(ex2srcdir)/density-of-states-2D.eps              \
	       $(ex2srcdir)/fermi-distribution-T-sweep.sh         \
	       $(ex2srcdir)/fermi-distribution-T-sweep.agr        \
	       $(ex2srcdir)/fermi-distribution-T-sweep.eps        \
	       $(ex2srcdir)/fermi-energy-T-sweep.agr              \
	       $(ex2srcdir)/fermi-energy-T-sweep.eps              \
	       $(ex2srcdir)/fermi-distribution-concentration-sweep.sh  \
	       $(ex2srcdir)/fermi-distribution-concentration-sweep.agr \
	       $(ex2srcdir)/fermi-distribution-concentration-sweep.eps \
	       $(ex2srcdir)/fermi-energy-concentration-sweep.agr \
	       $(ex2srcdir)/fermi-energy-concentration-sweep.eps \
	       $(ex2srcdir)/thermalised-distribution.sh          \
	       $(ex2srcdir)/thermalised-distribution.agr         \
	       $(ex2srcdir)/thermalised-distribution.eps         \
	       $(ex2srcdir)/finite-well-matching.sh              \
	       $(ex2srcdir)/finite-well-matching.agr             \
	       $(ex2srcdir)/finite-well-matching.eps             \
	       $(ex2srcdir)/finite-well-width.sh                 \
	       $(ex2srcdir)/finite-well-E-L.agr                  \
	       $(ex2srcdir)/finite-well-E-L.eps                  \
	       $(ex2srcdir)/finite-well-wf.sh                    \
	       $(ex2srcdir)/finite-well-wf.agr                   \
	       $(ex2srcdir)/finite-well-wf.eps                   \
	       $(ex2srcdir)/effective-mass-mismatch.sh           \
	       $(ex2srcdir)/effective-mass-mismatch-E-L.agr      \
	       $(ex2srcdir)/effective-mass-mismatch-E-L.eps      \
	       $(ex2srcdir)/effective-mass-mismatch-DE-L.agr     \
	       $(ex2srcdir)/effective-mass-mismatch-DE-L.eps     \
	       $(ex2srcdir)/infinite-barrier-height-limit.sh     \
	       $(ex2srcdir)/infinite-barrier-height-limit.agr    \
	       $(ex2srcdir)/infinite-barrier-height-limit.eps    \
	       $(ex2srcdir)/infinite-barrier-mass-limit.sh       \
	       $(ex2srcdir)/infinite-barrier-mass-limit.agr      \
	       $(ex2srcdir)/infinite-barrier-mass-limit.eps      \
	       $(ex2srcdir)/large-barrier-mass-wf.sh             \
	       $(ex2srcdir)/large-barrier-mass-wf.agr            \
	       $(ex2srcdir)/large-barrier-mass-wf.eps            \
	       $(ex2srcdir)/infinite-SL-40-40-E-k.sh             \
	       $(ex2srcdir)/infinite-SL-40-40-E-k.agr            \
	       $(ex2srcdir)/infinite-SL-40-40-E-k.eps            \
	       $(ex2srcdir)/infinite-SL-60-60-E-k.sh             \
	       $(ex2srcdir)/infinite-SL-miniband-width.sh        \
	       $(ex2srcdir)/infinite-SL-miniband-width.agr       \
	       $(ex2srcdir)/infinite-SL-miniband-width.eps       \
	       $(ex2srcdir)/infinite-SL-mass.sh                  \
	       $(ex2srcdir)/infinite-SL-mass.agr                 \
	       $(ex2srcdir)/infinite-SL-mass.eps                 \
	       $(ex2srcdir)/single-barrier-tx-width.sh           \
	       $(ex2srcdir)/single-barrier-tx-width.agr          \
	       $(ex2srcdir)/single-barrier-tx-width.eps          \
	       $(ex2srcdir)/single-barrier-tx-height.sh          \
	       $(ex2srcdir)/single-barrier-tx-height.agr         \
	       $(ex2srcdir)/single-barrier-tx-height.eps         \
               $(ex2srcdir)/double-barrier-tx-width.sh           \
               $(ex2srcdir)/double-barrier-tx-width.agr          \
               $(ex2srcdir)/double-barrier-tx-width.eps          \
               $(ex2srcdir)/double-barrier-tx-height.sh          \
               $(ex2srcdir)/double-barrier-tx-height.agr         \
               $(ex2srcdir)/double-barrier-tx-height.eps         \
               $(ex2srcdir)/double-barrier-IV-temperature.sh     \
               $(ex2srcdir)/double-barrier-IV-temperature.agr    \
               $(ex2srcdir)/double-barrier-IV-temperature.eps    \
               $(ex2srcdir)/double-barrier-IV-width.sh           \
               $(ex2srcdir)/double-barrier-IV-width.agr          \
               $(ex2srcdir)/double-barrier-IV-width.eps          \
               $(ex2srcdir)/landau-quantisation.sh               \
               $(ex2srcdir)/landau-quantisation.agr              \
               $(ex2srcdir)/landau-quantisation.eps


exnsdatadir = $(exdatadir)/numerical-solutions
exnssrcdir = examples/numerical-solutions

exnsdata_DATA = $(exnssrcdir)/index \
		$(exnssrcdir)/bisect-example.sh \
		$(exnssrcdir)/bisect-example.agr \
		$(exnssrcdir)/bisect-example.eps \
		$(exnssrcdir)/shooting-method-wavefunctions.sh    \
		$(exnssrcdir)/shooting-method-wavefunction-even.agr    \
		$(exnssrcdir)/shooting-method-wavefunction-even.eps    \
		$(exnssrcdir)/shooting-method-wavefunction-odd.agr    \
		$(exnssrcdir)/shooting-method-wavefunction-odd.eps    \
		$(exnssrcdir)/shooting-convergence.sh             \
		$(exnssrcdir)/shooting-convergence.eps             \
		$(exnssrcdir)/shooting-convergence.agr            \
		$(exnssrcdir)/shooting-convergence-LB.sh           \
		$(exnssrcdir)/shooting-convergence-LB-E.agr           \
		$(exnssrcdir)/shooting-convergence-LB-E.eps           \
		$(exnssrcdir)/shooting-convergence-LB-wf.agr           \
		$(exnssrcdir)/shooting-convergence-LB-wf.eps           \
		$(exnssrcdir)/poeschl-teller-hole-vary-lambda.sh  \
		$(exnssrcdir)/poeschl-teller-hole-E-lambda.dat    \
		$(exnssrcdir)/poeschl-teller-hole-V.agr           \
		$(exnssrcdir)/poeschl-teller-hole-V.eps           \
		$(exnssrcdir)/poeschl-teller-hole-convergence.sh  \
		$(exnssrcdir)/poeschl-teller-hole-convergence.dat \
		$(exnssrcdir)/variable-mass-0.2-alloy.sh          \
		$(exnssrcdir)/variable-mass-0.2-alloy.dat         \
		$(exnssrcdir)/variable-mass-0.75-alloy-well-width.sh  \
		$(exnssrcdir)/variable-mass-0.75-alloy-well-width.dat \
		$(exnssrcdir)/variable-mass-0.75-alloy-dz.sh  \
		$(exnssrcdir)/variable-mass-0.75-alloy-dz.dat \
		$(exnssrcdir)/double-quantum-well.sh          \
		$(exnssrcdir)/double-quantum-well-E-vs-LB.agr \
		$(exnssrcdir)/double-quantum-well-E-vs-LB.eps \
		$(exnssrcdir)/double-quantum-well-wf.agr      \
		$(exnssrcdir)/double-quantum-well-wf.eps      \
		$(exnssrcdir)/multiple-qw-40-40.sh            \
		$(exnssrcdir)/multiple-qw-40-40-E-N.agr       \
		$(exnssrcdir)/multiple-qw-40-40-E-N.eps       \
		$(exnssrcdir)/multiple-qw-40-40-wf.agr        \
		$(exnssrcdir)/multiple-qw-40-40-wf.eps        \
		$(exnssrcdir)/multiple-qw-100-100.sh          \
		$(exnssrcdir)/multiple-qw-100-100-wf.agr      \
		$(exnssrcdir)/multiple-qw-100-100-wf.eps      \
		$(exnssrcdir)/quantum-confined-stark-effect.sh \
		$(exnssrcdir)/quantum-confined-stark-effect-E-F.agr \
		$(exnssrcdir)/quantum-confined-stark-effect-E-F.eps \
		$(exnssrcdir)/field-induced-anti-crossings.sh \
		$(exnssrcdir)/field-induced-anti-crossings-E-F.agr \
		$(exnssrcdir)/field-induced-anti-crossings-E-F.eps \
		$(exnssrcdir)/symmetry-and-selection-rules.sh \
		$(exnssrcdir)/symmetry-and-selection-rules.agr \
		$(exnssrcdir)/symmetry-and-selection-rules.eps \
		$(exnssrcdir)/field-effect-variable-epsilon.sh  \
		$(exnssrcdir)/field-effect-variable-epsilon.eps \
		$(exnssrcdir)/field-effect-variable-epsilon.agr \
		$(exnssrcdir)/uncertainty-principle-finite-well.sh  \
		$(exnssrcdir)/uncertainty-principle-finite-well.eps \
		$(exnssrcdir)/uncertainty-principle-finite-well.agr \
		$(exnssrcdir)/uncertainty-principle-mass-limit.sh  \
		$(exnssrcdir)/uncertainty-principle-mass-limit.eps \
		$(exnssrcdir)/uncertainty-principle-mass-limit.agr \
		$(exnssrcdir)/nonparabolic-dE-vs-width.sh \
		$(exnssrcdir)/nonparabolic-dE-vs-width.agr \
		$(exnssrcdir)/nonparabolic-dE-vs-width.eps \
		$(exnssrcdir)/nonparabolic-dE-vs-width-inalgaas.sh \
		$(exnssrcdir)/nonparabolic-dE-vs-width-inalgaas.agr \
		$(exnssrcdir)/nonparabolic-dE-vs-width-inalgaas.eps \
		$(exnssrcdir)/poisson-schroedinger.sh \
		$(exnssrcdir)/poisson-schroedinger-E-I.agr \
		$(exnssrcdir)/poisson-schroedinger-E-I.eps \
		$(exnssrcdir)/poisson-schroedinger-field.agr \
		$(exnssrcdir)/poisson-schroedinger-field.eps \
		$(exnssrcdir)/poisson-schroedinger-sigma.agr \
		$(exnssrcdir)/poisson-schroedinger-sigma.eps \
		$(exnssrcdir)/poisson-schroedinger-Vp.agr \
		$(exnssrcdir)/poisson-schroedinger-Vp.eps \
		$(exnssrcdir)/poisson-schroedinger-Vt.agr \
		$(exnssrcdir)/poisson-schroedinger-Vt.eps \
		$(exnssrcdir)/poisson-schroedinger-comp.sh \
		$(exnssrcdir)/poisson-schroedinger-comp-E-I.agr \
		$(exnssrcdir)/poisson-schroedinger-comp-E-I.eps \
		$(exnssrcdir)/modulation-doping.sh \
		$(exnssrcdir)/modulation-doping-V.agr \
		$(exnssrcdir)/modulation-doping-V.eps \
		$(exnssrcdir)/high-electron-mobility-transistor.sh \
		$(exnssrcdir)/high-electron-mobility-transistor-v-I.agr \
		$(exnssrcdir)/high-electron-mobility-transistor-v-I.eps \
		$(exnssrcdir)/high-electron-mobility-transistor-wf-I.agr \
		$(exnssrcdir)/high-electron-mobility-transistor-wf-I.eps \
		$(exnssrcdir)/parabolic-well-solutions.sh \
		$(exnssrcdir)/parabolic-well-solutions.agr \
		$(exnssrcdir)/parabolic-well-solutions.eps \
		$(exnssrcdir)/parabolic-well-energies.agr \
		$(exnssrcdir)/parabolic-well-energies.eps \
		$(exnssrcdir)/parabolic-well-E-10.sh \
		$(exnssrcdir)/parabolic-well-E-10.dat

exdiffdatadir = $(exdatadir)/diffusion
exdiffsrcdir = examples/diffusion

exdiffdata_DATA = $(exdiffsrcdir)/index \
		  $(exdiffsrcdir)/constant-step.sh \
		  $(exdiffsrcdir)/constant-step.eps \
		  $(exdiffsrcdir)/constant-step.agr \
		  $(exdiffsrcdir)/error-function.sh \
		  $(exdiffsrcdir)/error-function.eps \
		  $(exdiffsrcdir)/error-function.agr \
		  $(exdiffsrcdir)/constant-Dt.sh     \
		  $(exdiffsrcdir)/constant-Dt.eps    \
		  $(exdiffsrcdir)/constant-Dt.agr    \
		  $(exdiffsrcdir)/constant-D.sh      \
		  $(exdiffsrcdir)/constant-D.eps     \
		  $(exdiffsrcdir)/constant-D.agr     \
		  $(exdiffsrcdir)/constant-D-E-t.eps \
		  $(exdiffsrcdir)/constant-D-E-t.agr \
		  $(exdiffsrcdir)/concentration-dependent.sh  \
		  $(exdiffsrcdir)/concentration-dependent.agr \
		  $(exdiffsrcdir)/concentration-dependent.eps \
		  $(exdiffsrcdir)/depth-dependent.sh          \
		  $(exdiffsrcdir)/depth-dependent.agr         \
		  $(exdiffsrcdir)/depth-dependent.eps         \
		  $(exdiffsrcdir)/time-dependent.sh           \
		  $(exdiffsrcdir)/time-dependent.agr          \
		  $(exdiffsrcdir)/time-dependent.eps          \
		  $(exdiffsrcdir)/delta-doping.sh             \
		  $(exdiffsrcdir)/delta-doping-d.agr          \
		  $(exdiffsrcdir)/delta-doping-d.eps          \
		  $(exdiffsrcdir)/delta-doping-vp.agr         \
		  $(exdiffsrcdir)/delta-doping-vp.eps         \
		  $(exdiffsrcdir)/delta-doping-E.agr          \
		  $(exdiffsrcdir)/delta-doping-E.eps

eximpdatadir = $(exdatadir)/impurities
eximpsrcdir = examples/impurities

eximpdata_DATA = $(eximpsrcdir)/index                \
		 $(eximpsrcdir)/E-variational-2D.sh  \
		 $(eximpsrcdir)/E-variational-2D.agr \
		 $(eximpsrcdir)/E-variational-2D.eps \
		 $(eximpsrcdir)/E-donor-2D.sh        \
		 $(eximpsrcdir)/E-donor-2D-total.eps \
		 $(eximpsrcdir)/E-donor-2D-total.agr \
		 $(eximpsrcdir)/E-donor-2D-binding.eps \
		 $(eximpsrcdir)/E-donor-2D-binding.agr \
		 $(eximpsrcdir)/E-donor-2D-lambda.eps  \
		 $(eximpsrcdir)/E-donor-2D-lambda.agr  \
		 $(eximpsrcdir)/E-binding-2D.sh        \
		 $(eximpsrcdir)/E-binding-2D.eps       \
		 $(eximpsrcdir)/E-binding-2D.agr       \
		 $(eximpsrcdir)/E-binding-2D-ED0.eps   \
		 $(eximpsrcdir)/E-binding-2D-ED0.agr   \
		 $(eximpsrcdir)/E-donor-3D-vs-2D.sh    \
		 $(eximpsrcdir)/E-donor-3D-vs-2D.agr   \
		 $(eximpsrcdir)/E-donor-3D-vs-2D.eps   \
		 $(eximpsrcdir)/E-donor-3D-wf.agr      \
		 $(eximpsrcdir)/E-donor-3D-wf.eps      \
		 $(eximpsrcdir)/E-donor-3D-wf-vs-2D.agr \
		 $(eximpsrcdir)/E-donor-3D-wf-vs-2D.eps \
		 $(eximpsrcdir)/E-binding-3D.sh         \
		 $(eximpsrcdir)/E-binding-3D.agr        \
		 $(eximpsrcdir)/E-binding-3D.eps        \
		 $(eximpsrcdir)/E-binding-3D-zoomed.agr \
		 $(eximpsrcdir)/E-binding-3D-zoomed.eps \
		 $(eximpsrcdir)/E-binding-3D-lambda.agr \
		 $(eximpsrcdir)/E-binding-3D-lambda.eps \
		 $(eximpsrcdir)/E-zeta.sh               \
		 $(eximpsrcdir)/E-zeta.agr              \
		 $(eximpsrcdir)/E-zeta.eps              \
		 $(eximpsrcdir)/E-donor-variable.sh     \
		 $(eximpsrcdir)/E-donor-variable.eps    \
		 $(eximpsrcdir)/E-donor-variable.agr    \
		 $(eximpsrcdir)/E-donor-variable-zeta-lambda.eps \
		 $(eximpsrcdir)/E-donor-variable-zeta-lambda.agr \
		 $(eximpsrcdir)/E-donor-comparison.sh            \
		 $(eximpsrcdir)/E-donor-comparison.dat           \
		 $(eximpsrcdir)/spin-flip-raman.sh 
#		 $(eximpsrcdir)/E-binding-alternative.sh

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
  GTEST_CXXFLAGS = @PTHREAD_CFLAGS@ -DGTEST_HAS_PTHREAD=1 -I $(srcdir)/src/
  GTEST_LIBS = @PTHREAD_LIBS@
else
  GTEST_CXXFLAGS = -DGTEST_HAS_PTHREAD=0 -I $(srcdir)/src
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
TESTS = tests/run-all-tests

check_PROGRAMS = tests/run-all-tests

tests_run_all_tests_SOURCES = tests/integration_tests/qwwad-schroedinger-infinite-well-tests.cpp \
			      tests/integration_tests/qwwad-schroedinger-tridiagonal-tests.cpp \
			      tests/microtests/schroedinger_solver_tests/testable_abstract_schroedinger_solver.h \
			      tests/microtests/schroedinger_solver_tests/abstract_schroedinger_solver_tests.cpp

tests_run_all_tests_CXXFLAGS = $(GTEST_CXXFLAGS)
tests_run_all_tests_LDADD = gtest/lib/libgtest_main.la \
			    gtest/lib/libgtest.la \
			    $(GTEST_LIBS) \
			    libqwwad.la
