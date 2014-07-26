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

bin_PROGRAMS =  src/chargedensity \
                src/csss    \
		src/cszb    \
		src/d0      \
		src/d02D    \
		src/d03D    \
		src/densityinput \
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
		src/efsqw   \
		src/efss    \
	       	src/efsx    \
		src/efxv    \
		src/find_heterostructure \
		src/find_poisson_potential \
		src/gde     \
		src/hup     \
		src/ivdb    \
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
		src/scps    \
		src/sfr     \
		src/slk     \
		src/sradp   \
		src/srcc    \
		src/srelo   \
	       	src/srmpr   \
	       	src/srrad   \
	       	src/tdb     \
	       	src/tsb     \
		src/wfplot

d0_helpers    = src/d0-helpers.c src/d0-helpers.h
qwwad_options = src/qwwad-options.cpp src/qwwad-options.h
wf_options    = src/wf_options.cpp src/wf_options.h $(qwwad_options)
qwwad_fileio  = src/qwwad-fileio.cpp src/qwwad-fileio.h

# QCLsim infrastructure (released from Leeds code)
noinst_LIBRARIES = src/libqwwad.a
src_libqwwad_a_SOURCES =  \
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
			   src/qwwad-heterostructure.cpp \
			   src/qwwad-heterostructure.h   \
			   src/qwwad-schroedinger.cpp    \
			   src/qwwad-schroedinger.h      \
			   src/qwwad-schroedinger-finite-well.cpp \
			   src/qwwad-schroedinger-finite-well.h   \
			   src/qwwad-schroedinger-full.h          \
			   src/qwwad-schroedinger-full.cpp        \
			   src/qwwad-schroedinger-infinite-well.cpp \
			   src/qwwad-schroedinger-infinite-well.h   \
			   src/qwwad-schroedinger-shooting.cpp      \
			   src/qwwad-schroedinger-shooting.h        \
			   src/qwwad-schroedinger-taylor.cpp        \
			   src/qwwad-schroedinger-taylor.h          \
			   src/qwwad-schroedinger-tridiagonal.cpp   \
			   src/qwwad-schroedinger-tridiagonal.h

src_chargedensity_SOURCES = src/chargedensity.cpp $(wf_options)
src_chargedensity_LDADD   = src/libqwwad.a $(LDADD)

src_densityinput_SOURCES = src/densityinput.cpp $(wf_options)
src_densityinput_LDADD   = src/libqwwad.a $(LDADD)

src_d0_SOURCES   = src/d0.c $(d0_helpers)

src_d02D_SOURCES = src/d02D.c $(d0_helpers)
src_d03D_SOURCES = src/d03D.c $(d0_helpers)

src_dispersion_relation_SOURCES = src/dispersion_relation.cpp $(wf_options)
src_dispersion_relation_LDADD   = src/libqwwad.a $(LDADD)

src_dos_SOURCES  = src/dos.cpp $(qwwad_options) $(qwwad_fileio)
src_dos_LDADD    = src/libqwwad.a $(LDADD)

src_efiw_SOURCES = src/efiw.cpp $(qwwad_options)
src_efiw_LDADD   = src/libqwwad.a $(LDADD)

src_efkpsl_SOURCES = src/efkpsl.cpp $(qwwad_options)
src_efkpsl_LDADD   = src/libqwwad.a $(LDADD)

src_efpqw_SOURCES = src/efpqw.cpp $(qwwad_options)
src_efpqw_LDADD   = src/libqwwad.a $(LDADD)

src_efss_SOURCES = src/efss.cpp $(qwwad_options)
src_efss_LDADD   = src/libqwwad.a $(LDADD)

src_efsqw_SOURCES = src/efsqw.cpp $(qwwad_options)
src_efsqw_LDADD   = src/libqwwad.a $(LDADD)

src_efxv_SOURCES = src/efxv.cpp $(qwwad_options)
src_efxv_LDADD   = src/libqwwad.a $(LDADD)

src_find_heterostructure_SOURCES = src/find_heterostructure.cpp $(qwwad_options)
src_find_heterostructure_LDADD   = src/libqwwad.a $(LDADD)

src_find_poisson_potential_SOURCES = src/find_poisson_potential.cpp $(qwwad_options)
src_find_poisson_potential_LDADD   = src/libqwwad.a $(LDADD)

src_hup_SOURCES = src/hup.cpp $(wf_options)
src_hup_LDADD   = src/libqwwad.a $(LDADD)

src_ivdb_SOURCES = src/ivdb.cpp $(qwwad_options)
src_ivdb_LDADD   = src/libqwwad.a $(LDADD)

src_ppvfq_SOURCES = src/ppvfq.c src/ppff.h
src_ppvfq_LDADD   = src/libqwwad.a

src_pplbso_SOURCES = src/pplbso.c src/ppsop.h

src_pth_SOURCES = src/pth.cpp $(qwwad_options)
src_pth_LDADD   = src/libqwwad.a $(LDADD)

src_sbp_SOURCES  = src/sbp.cpp $(qwwad_options) $(qwwad_fileio)
src_sbp_LDADD    = src/libqwwad.a $(LDADD)

src_tdb_SOURCES  = src/tdb.cpp $(qwwad_options) $(qwwad_fileio)
src_tdb_LDADD    = src/libqwwad.a $(LDADD)

src_tsb_SOURCES  = src/tsb.cpp $(qwwad_options) $(qwwad_fileio)
src_tsb_LDADD    = src/libqwwad.a $(LDADD)

src_wfplot_SOURCES = src/wfplot.cpp $(wf_options)
src_wfplot_LDADD   = src/libqwwad.a $(LDADD)

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
		 man/ivdb.1

man/dos.1: src/dos.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/dos$(EXEEXT)
	help2man --output=$@ --name="compute density of states in 1D, 2D and 3D systems" src/dos$(EXEEXT)

man/efiw.1: src/efiw.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/efiw$(EXEEXT)
	help2man --output=$@ --name="compute eigenstates in an infinite quantum well" src/efiw$(EXEEXT)

man/dispersion_relation.1: src/dispersion_relation.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/dispersion_relation$(EXEEXT)
	help2man --output=$@ --name="compute the dispersion relation for an energy subband" src/dispersion_relation$(EXEEXT)

man/sbp.1: src/sbp.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/sbp$(EXEEXT)
	help2man --output=$@ --name="find Fermi-Dirac distribution functions for a set of subbands" src/sbp$(EXEEXT)

man/efsqw.1: src/efsqw.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/efsqw$(EXEEXT)
	help2man --output=$@ --name="compute eigenstates of a finite square quantum well" src/efsqw$(EXEEXT)

man/wfplot.1: src/wfplot.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/wfplot$(EXEEXT)
	help2man --output=$@ --name="translate wavefunction data into a prettier plottable form." src/wfplot$(EXEEXT)

man/efkpsl.1: src/efkpsl.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/efkpsl$(EXEEXT)
	help2man --output=$@ --name="compute eigenstates of a Kronig-Penney superlattice." src/efkpsl$(EXEEXT)

man/tsb.1: src/tsb.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/tsb$(EXEEXT)
	help2man --output=$@ --name="compute transmission coefficient through single tunnelling barrier." src/tsb$(EXEEXT)

man/tdb.1: src/tdb.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/tdb$(EXEEXT)
	help2man --output=$@ --name="compute transmission coefficient through double tunnelling barrier." src/tdb$(EXEEXT)

man/ivdb.1: src/ivdb.cpp $(top_srcdir)/configure.ac
	$(MKDIR_P) man
	$(MAKE) $(AM_MAKEFLAGS) src/ivdb$(EXEEXT)
	help2man --output=$@ --name="compute current--voltage characteristics for a double tunnelling barrier." src/ivdb$(EXEEXT)

# Install examples
exdatadir = $(datadir)/qwwad/examples

exdata_DATA = examples/README \
	      examples/template.sh
 
ex2datadir = $(exdatadir)/solutions-to-schroedingers-equation
ex2srcdir = examples/solutions-to-schroedingers-equation

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
		$(exnssrcdir)/poisson-schroedinger-V.agr \
		$(exnssrcdir)/poisson-schroedinger-V.eps \
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
		$(exnssrcdir)/high-electron-mobility-transistor-wf-I.eps


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
			      tests/microtests/schroedinger_solver_tests/abstract_schroedinger_solver_tests.cpp

tests_run_all_tests_CXXFLAGS = $(GTEST_CXXFLAGS)
tests_run_all_tests_LDADD = gtest/lib/libgtest_main.la \
			    gtest/lib/libgtest.la \
			    $(GTEST_LIBS) \
			    src/libqwwad.a

