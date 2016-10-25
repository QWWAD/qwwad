#if HAVE_DOXYGEN
# OPTIONAL_DOCS = doc
#endif

ACLOCAL_AMFLAGS=-I m4
AM_CPPFLAGS  = $(GSL_CFLAGS) -DQWWAD_PKGDATADIR='"$(pkgdatadir)"'
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

dist_pkgdata_DATA = src/material-library.xml

# Install pkg-config support file
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = qwwad.pc

# Deprecated header files that:
# 1. will never be part of libqwwad
# 2. are used by multiple programs
noinst_HEADERS = src/ef-helpers.h    \
		 src/maths.h         \
		 src/struct.h 

# The full set of QWWAD programs, to be installed into the user's bin folder
bin_PROGRAMS = \
		qwwad_charge_density           \
		qwwad_cs_single_spiral         \
		qwwad_cs_zinc_blende           \
		qwwad_density_of_states        \
		qwwad_diffuse                  \
		qwwad_ef_band_edge             \
		qwwad_ef_cylindrical_wire      \
		qwwad_ef_cylindrical_wire_wf   \
		qwwad_ef_dispersion            \
		qwwad_ef_donor_generic         \
		qwwad_ef_donor_specific        \
		qwwad_ef_exciton               \
		qwwad_ef_generic               \
		qwwad_ef_infinite_well         \
		qwwad_ef_infinite_wire         \
		qwwad_ef_parabolic_well        \
		qwwad_ef_plot                  \
		qwwad_ef_poeschl_teller        \
		qwwad_ef_spherical_dot         \
		qwwad_ef_spherical_dot_wf      \
		qwwad_ef_square_well           \
		qwwad_ef_superlattice          \
	       	qwwad_ef_zeeman                \
	       	qwwad_fermi_distribution       \
		qwwad_material_property        \
		qwwad_mesh                     \
		qwwad_poisson                  \
		qwwad_population_init          \
		qwwad_pp_charge_density        \
		qwwad_pp_dispersion            \
		qwwad_pp_form_factor           \
		qwwad_pp_large_basis           \
		qwwad_pp_large_basis_so        \
		qwwad_pp_lattice_vector_table  \
		qwwad_pp_superlattice          \
		qwwad_reciprocal_fcc           \
		qwwad_reciprocal_cube          \
		qwwad_reciprocal_single_spiral \
		qwwad_specific_heat_capacity   \
		qwwad_spin_flip_raman          \
		qwwad_sr_acoustic_phonon       \
		qwwad_sr_alloy_disorder        \
		qwwad_sr_carrier_carrier       \
		qwwad_sr_interface_roughness   \
		qwwad_sr_impurity              \
		qwwad_sr_lo_phonon             \
		qwwad_sr_radiative             \
		qwwad_superlattice_k           \
		qwwad_thermal_1d               \
		qwwad_thermal_rc               \
	       	qwwad_tx_double_barrier        \
		qwwad_tx_double_barrier_iv     \
	       	qwwad_tx_single_barrier        \
		qwwad_uncertainty              \
		qwwad_critical_thickness       \
	       	effv                           \
	       	efsx

# Executable scripts to be installed into the user's
# bin folder
dist_bin_SCRIPTS = src/xyz2pdb

# Build most functionality into library infrastructure
lib_LTLIBRARIES = libqwwad.la

# Version info for the shared library.
# Note that the numbers must be set as follows:
# (1) "current"  - the latest MAJOR ABI version. Increase this if it is incompatible with the last version
# (2) "revision" - the latest MINOR ABI version. Increase this if minor,
#                  compatible changes have been introduced.
#                  Reset to zero if the current ABI version increases.
# (3) "age"      - the number of historical major ABI versions that are backward-compatible with the current
#                  version. This should stay as zero unless the current version is guaranteed to
#                  support older interfaces
libqwwadversionflags = -version-info 1:0:0

# Shared library header files that will be installed into the user's system
includeqwwaddir=$(includedir)/qwwad
includeqwwad_HEADERS = \
			   src/qwwad/constants.h                          \
			   src/qwwad/data-checker.h                       \
			   src/qwwad/donor-energy-minimiser-fast.h        \
			   src/qwwad/donor-energy-minimiser-linear.h      \
			   src/qwwad/donor-energy-minimiser.h             \
			   src/qwwad/dos-functions.h                      \
			   src/qwwad/double-barrier.h                     \
			   src/qwwad/eigenstate.h                         \
			   src/qwwad/fermi.h                              \
			   src/qwwad/file-io.h                            \
			   src/qwwad/lapack-declarations.h                \
			   src/qwwad/linear-algebra.h                     \
			   src/qwwad/maths-helpers.h                      \
			   src/qwwad/mesh.h                               \
			   src/qwwad/options.h                            \
			   src/qwwad/schroedinger-solver.h                \
			   src/qwwad/schroedinger-solver-donor.h          \
			   src/qwwad/schroedinger-solver-donor-2D.h       \
			   src/qwwad/schroedinger-solver-donor-3D.h       \
			   src/qwwad/schroedinger-solver-donor-variable.h \
			   src/qwwad/schroedinger-solver-finite-well.h    \
			   src/qwwad/schroedinger-solver-full.h           \
			   src/qwwad/schroedinger-solver-infinite-well.h  \
			   src/qwwad/schroedinger-solver-kronig-penney.h  \
			   src/qwwad/schroedinger-solver-poeschl-teller.h \
			   src/qwwad/schroedinger-solver-shooting.h       \
			   src/qwwad/schroedinger-solver-taylor.h         \
			   src/qwwad/schroedinger-solver-tridiagonal.h    \
			   src/qwwad/subband.h

# Source files for libqwwad that will be compiled into the binary library
libqwwad_la_SOURCES =  \
		           src/qwwad/data-checker.cpp                       \
			   src/qwwad/dos-functions.cpp                      \
			   src/qwwad/double-barrier.cpp                     \
			   src/qwwad/fermi.cpp                              \
			   src/qwwad/donor-energy-minimiser.cpp             \
			   src/qwwad/donor-energy-minimiser-linear.cpp      \
			   src/qwwad/donor-energy-minimiser-fast.cpp        \
			   src/qwwad/mesh.cpp                               \
			   src/qwwad/file-io.cpp                            \
			   src/qwwad/file-io-deprecated.cpp                 \
			   src/qwwad/file-io-deprecated.h                   \
			   src/qwwad/options.cpp                            \
			   src/qwwad/intersubband-transition.h              \
			   src/qwwad/intersubband-transition.cpp            \
			   src/qwwad/linear-algebra.cpp                     \
			   src/qwwad/maths-helpers.cpp                      \
			   src/qwwad/poisson-solver.cpp                     \
			   src/qwwad/poisson-solver.h                       \
			   src/qwwad/scattering-calculator-LO.h             \
			   src/qwwad/scattering-calculator-LO.cpp           \
			   src/qwwad/schroedinger-solver.cpp                \
			   src/qwwad/schroedinger-solver-donor.cpp          \
			   src/qwwad/schroedinger-solver-donor-2D.cpp       \
			   src/qwwad/schroedinger-solver-donor-3D.cpp       \
			   src/qwwad/schroedinger-solver-donor-variable.cpp \
			   src/qwwad/schroedinger-solver-finite-well.cpp    \
			   src/qwwad/schroedinger-solver-full.cpp           \
			   src/qwwad/schroedinger-solver-infinite-well.cpp  \
			   src/qwwad/schroedinger-solver-kronig-penney.cpp  \
			   src/qwwad/schroedinger-solver-poeschl-teller.cpp \
			   src/qwwad/schroedinger-solver-shooting.cpp       \
			   src/qwwad/schroedinger-solver-taylor.cpp         \
			   src/qwwad/schroedinger-solver-tridiagonal.cpp    \
			   src/qwwad/eigenstate.cpp                         \
			   src/qwwad/subband.cpp                            \
			   src/qwwad/wf_options.cpp                         \
			   src/qwwad/wf_options.h                           \
			   src/qwwad-material.cpp                           \
			   src/qwwad-material.h                             \
			   src/qwwad-material-property.cpp                  \
			   src/qwwad-material-property.h                    \
			   src/qwwad-material-property-constant.cpp         \
			   src/qwwad-material-property-constant.h           \
			   src/qwwad-material-property-numeric.cpp          \
			   src/qwwad-material-property-numeric.h            \
			   src/qwwad-material-property-interp.cpp           \
			   src/qwwad-material-property-interp.h             \
			   src/qwwad-material-property-poly.cpp             \
			   src/qwwad-material-property-poly.h               \
			   src/qwwad-material-property-string.cpp           \
			   src/qwwad-material-property-string.h             \
			   src/material_library.cpp                         \
			   src/material_library.h                           \
			   src/qwwad-debye.cpp                              \
			   src/qwwad-debye.h

libqwwad_la_CFLAGS = $(GSL_CFLAGS) \
		     $(BOOST_CPPFLAGS)

libqwwad_la_LDFLAGS = $(BOOST_LDFLAGS) \
		      $(libqwwadversionflags)

libqwwad_la_LIBADD = $(GSL_LIBS) \
                     $(BOOST_PROGRAM_OPTIONS_LIB) \
		     $(LAPACK_LIBS) \
		     $(BLAS_LIBS)   \
		     $(FLIBS)

LDADD = libqwwad.la

# Deprecated programs
# TODO: Remove these when they are no longer needed
effv_SOURCES = src/effv.c
efsx_SOURCES = src/efsx.c

qwwad_sr_alloy_disorder_SOURCES = src/qwwad_sr_alloy_disorder.cpp
qwwad_charge_density_SOURCES    = src/qwwad_charge_density.cpp
qwwad_critical_thickness_SOURCES = src/qwwad_critical_thickness.cpp
qwwad_cs_single_spiral_SOURCES = src/qwwad_cs_single_spiral.cpp
qwwad_cs_zinc_blende_SOURCES = src/qwwad_cs_zinc_blende.cpp
qwwad_population_init_SOURCES = src/qwwad_population_init.cpp
qwwad_ef_dispersion_SOURCES = src/qwwad_ef_dispersion.cpp
qwwad_density_of_states_SOURCES  = src/qwwad_density_of_states.cpp
qwwad_ef_exciton_SOURCES = src/qwwad_ef_exciton.cpp
qwwad_ef_cylindrical_wire_SOURCES = src/qwwad_ef_cylindrical_wire.cpp
qwwad_ef_cylindrical_wire_wf_SOURCES = src/qwwad_ef_cylindrical_wire_wf.cpp
qwwad_ef_infinite_well_SOURCES = src/qwwad_ef_infinite_well.cpp
qwwad_ef_infinite_wire_SOURCES = src/qwwad_ef_infinite_wire.cpp
qwwad_ef_superlattice_SOURCES = src/qwwad_ef_superlattice.cpp
qwwad_ef_zeeman_SOURCES = src/qwwad_ef_zeeman.cpp
qwwad_ef_parabolic_well_SOURCES = src/qwwad_ef_parabolic_well.cpp
qwwad_ef_spherical_dot_SOURCES = src/qwwad_ef_spherical_dot.cpp
qwwad_ef_spherical_dot_wf_SOURCES = src/qwwad_ef_spherical_dot_wf.cpp
qwwad_ef_generic_SOURCES = src/qwwad_ef_generic.cpp
qwwad_ef_square_well_SOURCES = src/qwwad_ef_square_well.cpp
qwwad_ef_band_edge_SOURCES = src/qwwad_ef_band_edge.cpp
qwwad_mesh_SOURCES = src/qwwad_mesh.cpp
qwwad_poisson_SOURCES = src/qwwad_poisson.cpp
qwwad_diffuse_SOURCES = src/qwwad_diffuse.cpp
qwwad_uncertainty_SOURCES = src/qwwad_uncertainty.cpp
qwwad_sr_interface_roughness_SOURCES = src/qwwad_sr_interface_roughness.cpp
qwwad_ef_donor_generic_SOURCES  = src/qwwad_ef_donor_generic.cpp
qwwad_sr_impurity_SOURCES  = src/qwwad_sr_impurity.cpp
qwwad_tx_double_barrier_iv_SOURCES = src/qwwad_tx_double_barrier_iv.cpp
qwwad_material_property_SOURCES = src/qwwad_material_property.cpp
qwwad_pp_charge_density_SOURCES = src/qwwad_pp_charge_density.cpp src/ppff.h src/ppff.cpp
qwwad_pp_dispersion_SOURCES = src/qwwad_pp_dispersion.cpp
qwwad_pp_large_basis_SOURCES  = src/qwwad_pp_large_basis.cpp src/pplb-functions.cpp src/pplb-functions.h src/ppff.h src/ppff.cpp
qwwad_pp_lattice_vector_table_SOURCES = src/qwwad_pp_lattice_vector_table.cpp
qwwad_pp_superlattice_SOURCES = src/qwwad_pp_superlattice.cpp src/ppff.h src/ppff.cpp
qwwad_pp_form_factor_SOURCES = src/qwwad_pp_form_factor.cpp src/ppff.h src/ppff.cpp
qwwad_pp_large_basis_so_SOURCES = src/qwwad_pp_large_basis_so.cpp src/pplb-functions.h src/pplb-functions.cpp src/ppff.h src/ppff.cpp src/ppsop.h src/ppsop.cpp
qwwad_ef_poeschl_teller_SOURCES = src/qwwad_ef_poeschl_teller.cpp
qwwad_ef_donor_specific_SOURCES = src/qwwad_ef_donor_specific.cpp
qwwad_reciprocal_fcc_SOURCES = src/qwwad_reciprocal_fcc.cpp
qwwad_reciprocal_cube_SOURCES = src/qwwad_reciprocal_cube.cpp
qwwad_reciprocal_single_spiral_SOURCES = src/qwwad_reciprocal_single_spiral.cpp
qwwad_fermi_distribution_SOURCES  = src/qwwad_fermi_distribution.cpp
qwwad_spin_flip_raman_SOURCES  = src/qwwad_spin_flip_raman.cpp
qwwad_superlattice_k_SOURCES  = src/qwwad_superlattice_k.cpp
qwwad_sr_acoustic_phonon_SOURCES = src/qwwad_sr_acoustic_phonon.cpp
qwwad_sr_carrier_carrier_SOURCES  = src/qwwad_sr_carrier_carrier.cpp
qwwad_sr_lo_phonon_SOURCES  = src/qwwad_sr_lo_phonon.cpp
qwwad_sr_radiative_SOURCES  = src/qwwad_sr_radiative.cpp
qwwad_specific_heat_capacity_SOURCES = src/qwwad_specific_heat_capacity.cpp
qwwad_tx_double_barrier_SOURCES  = src/qwwad_tx_double_barrier.cpp
qwwad_thermal_1d_SOURCES = src/qwwad_thermal_1d.cpp
qwwad_thermal_rc_SOURCES = src/qwwad_thermal_rc.cpp
qwwad_tx_single_barrier_SOURCES  = src/qwwad_tx_single_barrier.cpp
qwwad_ef_plot_SOURCES = src/qwwad_ef_plot.cpp

# Build rule for installing manual pages
dist_man1_MANS = \
		 man/qwwad_density_of_states.1  \
		 man/qwwad_ef_dispersion.1      \
		 man/qwwad_ef_infinite_well.1   \
		 man/qwwad_ef_square_well.1     \
		 man/qwwad_fermi_distribution.1 \
		 man/qwwad_ef_plot.1              \
		 man/qwwad_ef_superlattice.1              \
		 man/qwwad_tx_single_barrier.1                 \
		 man/qwwad_tx_double_barrier.1                 \
		 man/qwwad_tx_double_barrier_iv.1                \
		 man/qwwad_mesh.1 \
		 man/qwwad_ef_band_edge.1                 \
		 man/qwwad_ef_generic.1                 \
		 man/qwwad_ef_parabolic_well.1                \
		 man/qwwad_ef_poeschl_teller.1                  \
		 man/qwwad_ef_donor_specific.1 \
		 man/qwwad_ef_donor_generic.1  \
		 man/qwwad_ef_exciton.1        \
		 man/qwwad_poisson.1 \
		 man/qwwad_uncertainty.1                    \
		 man/qwwad_population_init.1           \
		 man/qwwad_charge_density.1          \
		 man/qwwad_diffuse.1                 \
		 man/qwwad_ef_zeeman.1               \
		 man/qwwad_spin_flip_raman.1

CLEANFILES = $(dist_man1_MANS)

EXTRA_DIST += \
	      doc/qwwad_density_of_states.h2m  \
	      doc/qwwad_ef_dispersion.h2m      \
              doc/qwwad_ef_infinite_well.h2m   \
	      doc/qwwad_ef_square_well.h2m     \
	      doc/qwwad_fermi_distribution.h2m \
	      doc/qwwad_mesh.h2m   \
	      doc/qwwad_ef_band_edge.h2m                   \
	      doc/qwwad_ef_generic.h2m                   \
	      doc/qwwad_ef_parabolic_well.h2m                  \
	      doc/qwwad_ef_poeschl_teller.h2m                    \
	      doc/qwwad_poisson.h2m \
	      doc/qwwad_ef_plot.h2m                 \
	      doc/qwwad_ef_superlattice.h2m                 \
	      doc/qwwad_tx_single_barrier.h2m                    \
	      doc/qwwad_tx_double_barrier.h2m                    \
	      doc/qwwad_tx_double_barrier_iv.h2m                   \
	      doc/qwwad_uncertainty.h2m                    \
	      doc/qwwad_population_init.h2m           \
	      doc/qwwad_charge_density.h2m          \
	      doc/qwwad_diffuse.h2m                    \
	      doc/qwwad_ef_donor_specific.h2m      \
	      doc/qwwad_ef_donor_generic.h2m       \
	      doc/qwwad_ef_exciton.h2m             \
	      doc/qwwad_ef_zeeman.h2m              \
	      doc/qwwad_spin_flip_raman.h2m

man/qwwad_density_of_states.1: src/qwwad_density_of_states.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_density_of_states.h2m | qwwad_density_of_states$(EXEEXT) 
	$(MKDIR_P) man
	help2man --output=$@ --name="compute density of states in 1D, 2D and 3D systems" --include $(top_srcdir)/doc/qwwad_density_of_states.h2m ./qwwad_density_of_states$(EXEEXT)

man/qwwad_ef_infinite_well.1: src/qwwad_ef_infinite_well.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_ef_infinite_well.h2m | qwwad_ef_infinite_well$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates in an infinite quantum well" --include $(top_srcdir)/doc/qwwad_ef_infinite_well.h2m ./qwwad_ef_infinite_well$(EXEEXT)

man/qwwad_ef_dispersion.1: src/qwwad_ef_dispersion.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_ef_dispersion.h2m | qwwad_ef_dispersion$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute the dispersion relation for a set of energy subbands" --include $(top_srcdir)/doc/qwwad_ef_dispersion.h2m ./qwwad_ef_dispersion$(EXEEXT)

man/qwwad_fermi_distribution.1: src/qwwad_fermi_distribution.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_fermi_distribution.h2m | qwwad_fermi_distribution$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find Fermi-Dirac distribution functions for a set of subbands" --include $(top_srcdir)/doc/qwwad_fermi_distribution.h2m ./qwwad_fermi_distribution$(EXEEXT)

man/qwwad_ef_square_well.1: src/qwwad_ef_square_well.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_ef_square_well.h2m | qwwad_ef_square_well$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates of a finite square quantum well" --include $(top_srcdir)/doc/qwwad_ef_square_well.h2m ./qwwad_ef_square_well$(EXEEXT)

man/qwwad_ef_plot.1: src/qwwad_ef_plot.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_ef_plot.h2m | qwwad_ef_plot$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="translate wavefunction data into a prettier plottable form." --include $(top_srcdir)/doc/qwwad_ef_plot.h2m ./qwwad_ef_plot$(EXEEXT)

man/qwwad_ef_superlattice.1: src/qwwad_ef_superlattice.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_ef_superlattice.h2m | qwwad_ef_superlattice$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates of a Kronig-Penney superlattice." --include $(top_srcdir)/doc/qwwad_ef_superlattice.h2m ./qwwad_ef_superlattice$(EXEEXT)

man/qwwad_tx_single_barrier.1: src/qwwad_tx_single_barrier.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_tx_single_barrier.h2m | qwwad_tx_single_barrier$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute transmission coefficient through single tunnelling barrier." --include $(top_srcdir)/doc/qwwad_tx_single_barrier.h2m ./qwwad_tx_single_barrier$(EXEEXT) 

man/qwwad_tx_double_barrier.1: src/qwwad_tx_double_barrier.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_tx_double_barrier.h2m | qwwad_tx_double_barrier$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute transmission coefficient through double tunnelling barrier." --include $(top_srcdir)/doc/qwwad_tx_double_barrier.h2m ./qwwad_tx_double_barrier$(EXEEXT)

man/qwwad_tx_double_barrier_iv.1: src/qwwad_tx_double_barrier_iv.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp $(top_srcdir)/doc/qwwad_tx_double_barrier_iv.h2m | qwwad_tx_double_barrier_iv$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute current--voltage characteristics for a double tunnelling barrier." --include $(top_srcdir)/doc/qwwad_tx_double_barrier_iv.h2m ./qwwad_tx_double_barrier_iv$(EXEEXT)

man/qwwad_mesh.1: src/qwwad_mesh.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_mesh.h2m | qwwad_mesh$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate spatial mesh and output alloy & doping profiles." --include $(top_srcdir)/doc/qwwad_mesh.h2m ./qwwad_mesh$(EXEEXT)

man/qwwad_ef_band_edge.1: src/qwwad_ef_band_edge.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_band_edge.h2m | qwwad_ef_band_edge$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="tabulate material parameters for heterostructure." --include $(top_srcdir)/doc/qwwad_ef_band_edge.h2m ./qwwad_ef_band_edge$(EXEEXT)

man/qwwad_ef_generic.1: src/qwwad_ef_generic.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_generic.h2m | qwwad_ef_generic$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="solve the Schroedinger equation numerically." --include $(top_srcdir)/doc/qwwad_ef_generic.h2m ./qwwad_ef_generic$(EXEEXT)

man/qwwad_ef_parabolic_well.1: src/qwwad_ef_parabolic_well.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_parabolic_well.h2m | qwwad_ef_parabolic_well$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate a parabolic alloy profile." --include $(top_srcdir)/doc/qwwad_ef_parabolic_well.h2m ./qwwad_ef_parabolic_well$(EXEEXT)

man/qwwad_ef_poeschl_teller.1: src/qwwad_ef_poeschl_teller.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_poeschl_teller.h2m | qwwad_ef_poeschl_teller$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="compute eigenstates of a Poeschl-Teller potential hole." --include $(top_srcdir)/doc/qwwad_ef_poeschl_teller.h2m ./qwwad_ef_poeschl_teller$(EXEEXT) 

man/qwwad_poisson.1: src/qwwad_poisson.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_poisson.h2m | qwwad_poisson$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find the Poisson potential induced by a charge profile." --include $(top_srcdir)/doc/qwwad_poisson.h2m ./qwwad_poisson$(EXEEXT) 

man/qwwad_uncertainty.1: src/qwwad_uncertainty.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_uncertainty.h2m | qwwad_uncertainty$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find the uncertainty product for position and momentum." --include $(top_srcdir)/doc/qwwad_uncertainty.h2m ./qwwad_uncertainty$(EXEEXT) 

man/qwwad_population_init.1: src/qwwad_population_init.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_population_init.h2m | qwwad_population_init$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate an appoximate carrier energy distribution." --include $(top_srcdir)/doc/qwwad_population_init.h2m ./qwwad_population_init$(EXEEXT) 

man/qwwad_charge_density.1: src/qwwad_charge_density.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_charge_density.h2m | qwwad_charge_density$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate a charge density profile for a 2D system." --include $(top_srcdir)/doc/qwwad_charge_density.h2m ./qwwad_charge_density$(EXEEXT) 

man/qwwad_diffuse.1: src/qwwad_diffuse.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_diffuse.h2m | qwwad_diffuse$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="generate a 1D diffusion profile." --include $(top_srcdir)/doc/qwwad_diffuse.h2m ./qwwad_diffuse$(EXEEXT) 

man/qwwad_ef_donor_specific.1: src/qwwad_ef_donor_specific.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_donor_specific.h2m | qwwad_ef_donor_specific$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find ground state of electron orbiting a donor ion in a heterostructure." --include $(top_srcdir)/doc/qwwad_ef_donor_specific.h2m ./qwwad_ef_donor_specific$(EXEEXT) 

man/qwwad_ef_zeeman.1: src/qwwad_ef_zeeman.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_zeeman.h2m | qwwad_ef_zeeman$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find the Zeeman potential for a 2D system." --include $(top_srcdir)/doc/qwwad_ef_zeeman.h2m ./qwwad_ef_zeeman$(EXEEXT) 

man/qwwad_spin_flip_raman.1: src/qwwad_spin_flip_raman.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_spin_flip_raman.h2m | qwwad_spin_flip_raman$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find the spin-flip Raman spectrum." --include $(top_srcdir)/doc/qwwad_spin_flip_raman.h2m ./qwwad_spin_flip_raman$(EXEEXT) 

man/qwwad_ef_donor_generic.1: src/qwwad_ef_donor_generic.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_donor_generic.h2m | qwwad_ef_donor_generic$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find state of electron orbiting a donor ion in a heterostructure using a generic search." --include $(top_srcdir)/doc/qwwad_ef_donor_generic.h2m ./qwwad_ef_donor_generic$(EXEEXT) 

man/qwwad_ef_exciton.1: src/qwwad_ef_exciton.cpp $(top_srcdir)/configure.ac src/qwwad/options.cpp \
    $(top_srcdir)/doc/qwwad_ef_exciton.h2m | qwwad_ef_exciton$(EXEEXT)
	$(MKDIR_P) man
	help2man --output=$@ --name="find ground state of an exciton within a heterostructure." --include $(top_srcdir)/doc/qwwad_ef_exciton.h2m ./qwwad_ef_exciton$(EXEEXT) 

# Install examples
exdatadir = $(datadir)/qwwad/examples

exdata_DATA = examples/README \
	      examples/template.sh
 
ex2datadir = $(exdatadir)/schroedinger-solutions
ex2srcdir = examples/schroedinger-solutions

EXTRA_DIST += examples

ex2data_DATA = $(ex2srcdir)/index                                 \
               $(ex2srcdir)/infinite-well-energy-vs-width.sh      \
               $(ex2srcdir)/infinite-well-wavefunctions.sh        \
	       $(ex2srcdir)/non-parabolic-dispersion.sh           \
	       $(ex2srcdir)/density-of-states.sh                  \
	       $(ex2srcdir)/fermi-distribution-T-sweep.sh         \
	       $(ex2srcdir)/fermi-distribution-concentration-sweep.sh  \
	       $(ex2srcdir)/thermalised-distribution.sh          \
	       $(ex2srcdir)/finite-well-matching.sh              \
	       $(ex2srcdir)/finite-well-width.sh                 \
	       $(ex2srcdir)/finite-well-wf.sh                    \
	       $(ex2srcdir)/effective-mass-mismatch.sh           \
	       $(ex2srcdir)/infinite-barrier-height-limit.sh     \
	       $(ex2srcdir)/infinite-barrier-mass-limit.sh       \
	       $(ex2srcdir)/large-barrier-mass-wf.sh             \
	       $(ex2srcdir)/infinite-SL-40-40-E-k.sh             \
	       $(ex2srcdir)/infinite-SL-60-60-E-k.sh             \
	       $(ex2srcdir)/infinite-SL-miniband-width.sh        \
	       $(ex2srcdir)/infinite-SL-mass.sh                  \
	       $(ex2srcdir)/single-barrier-tx-width.sh           \
	       $(ex2srcdir)/single-barrier-tx-height.sh          \
               $(ex2srcdir)/double-barrier-tx-width.sh           \
               $(ex2srcdir)/double-barrier-tx-height.sh          \
               $(ex2srcdir)/double-barrier-IV-temperature.sh     \
               $(ex2srcdir)/double-barrier-IV-width.sh           \
               $(ex2srcdir)/landau-quantisation.sh

exnsdatadir = $(exdatadir)/numerical-solutions
exnssrcdir = examples/numerical-solutions

exnsdata_DATA = $(exnssrcdir)/index                                  \
		$(exnssrcdir)/bisect-example.sh                      \
		$(exnssrcdir)/shooting-method-wavefunctions.sh       \
		$(exnssrcdir)/shooting-convergence.sh                \
		$(exnssrcdir)/shooting-convergence-LB.sh             \
		$(exnssrcdir)/poeschl-teller-hole-vary-lambda.sh     \
		$(exnssrcdir)/poeschl-teller-hole-convergence.sh     \
		$(exnssrcdir)/variable-mass-0.2-alloy.sh             \
		$(exnssrcdir)/variable-mass-0.75-alloy-well-width.sh \
		$(exnssrcdir)/variable-mass-0.75-alloy-dz.sh         \
		$(exnssrcdir)/double-quantum-well.sh                 \
		$(exnssrcdir)/multiple-qw-40-40.sh                   \
		$(exnssrcdir)/multiple-qw-100-100.sh                 \
		$(exnssrcdir)/quantum-confined-stark-effect.sh       \
		$(exnssrcdir)/field-induced-anti-crossings.sh        \
		$(exnssrcdir)/symmetry-and-selection-rules.sh        \
		$(exnssrcdir)/field-effect-variable-epsilon.sh       \
		$(exnssrcdir)/uncertainty-principle-finite-well.sh   \
		$(exnssrcdir)/uncertainty-principle-mass-limit.sh    \
		$(exnssrcdir)/nonparabolic-dE-vs-width.sh            \
		$(exnssrcdir)/nonparabolic-dE-vs-width-inalgaas.sh   \
		$(exnssrcdir)/poisson-schroedinger.sh                \
		$(exnssrcdir)/poisson-schroedinger-comp.sh           \
		$(exnssrcdir)/modulation-doping.sh                   \
		$(exnssrcdir)/high-electron-mobility-transistor.sh   \
		$(exnssrcdir)/parabolic-well-solutions.sh            \
		$(exnssrcdir)/parabolic-well-E-10.sh

exdiffdatadir = $(exdatadir)/diffusion
exdiffsrcdir = examples/diffusion

exdiffdata_DATA = $(exdiffsrcdir)/index \
		  $(exdiffsrcdir)/constant-step.sh \
		  $(exdiffsrcdir)/error-function.sh \
		  $(exdiffsrcdir)/constant-Dt.sh     \
		  $(exdiffsrcdir)/constant-D.sh      \
		  $(exdiffsrcdir)/concentration-dependent.sh  \
		  $(exdiffsrcdir)/depth-dependent.sh          \
		  $(exdiffsrcdir)/time-dependent.sh           \
		  $(exdiffsrcdir)/delta-doping.sh

eximpdatadir = $(exdatadir)/impurities
eximpsrcdir = examples/impurities

eximpdata_DATA = $(eximpsrcdir)/index                    \
		 $(eximpsrcdir)/E-variational-2D.sh      \
		 $(eximpsrcdir)/E-donor-2D.sh            \
		 $(eximpsrcdir)/E-binding-2D.sh          \
		 $(eximpsrcdir)/E-donor-3D-vs-2D.sh      \
		 $(eximpsrcdir)/E-binding-3D.sh          \
		 $(eximpsrcdir)/E-zeta.sh                \
		 $(eximpsrcdir)/E-donor-variable.sh      \
		 $(eximpsrcdir)/E-donor-comparison.sh    \
		 $(eximpsrcdir)/spin-flip-raman.sh       \
		 $(eximpsrcdir)/E-binding-alternative.sh \
		 $(eximpsrcdir)/E-new-vs-old.sh          \
		 $(eximpsrcdir)/E-excited-states.sh

exexcdatadir = $(exdatadir)/excitons
exexcsrcdir = examples/excitons

exexcdata_DATA = $(exexcsrcdir)/index              \
		 $(exexcsrcdir)/infinite-well.sh   \
		 $(exexcsrcdir)/lambda-beta-map.sh \
		 $(exexcsrcdir)/p-vs-a.sh          \
		 $(exexcsrcdir)/finite-well.sh     \
		 $(exexcsrcdir)/multiple-qw-wf.sh  \
		 $(exexcsrcdir)/multiple-qw.sh

exqwddatadir = $(exdatadir)/simple-wires-and-dots
exqwdsrcdir  = examples/simple-wires-and-dots

exqwddata_DATA = $(exqwdsrcdir)/index \
		 $(exqwdsrcdir)/infinite-rectangular-wire.sh \
		 $(exqwdsrcdir)/infinite-rectangular-wire-wf.gnuplot

exsrdatadir = $(exdatadir)/carrier-scattering
exsrsrcdir  = examples/carrier-scattering

exsrdata_DATA = $(exsrsrcdir)/index           \
		$(exsrsrcdir)/cc-screening.sh \
		$(exsrsrcdir)/cc-avg.sh     \
		$(exsrsrcdir)/cc-avg-N.sh   \
		$(exsrsrcdir)/cc-avg-T.sh   \
		$(exsrsrcdir)/cc-intra-N.sh \
		$(exsrsrcdir)/cc-intra-T.sh \
		$(exsrsrcdir)/cc-auger.sh   \
		$(exsrsrcdir)/cc-asymmetric.sh \
		$(exsrsrcdir)/imp-screening.sh \
		$(exsrsrcdir)/imp-profiles.sh \
		$(exsrsrcdir)/imp-profiles-N.sh \
		$(exsrsrcdir)/imp-N.sh \
		$(exsrsrcdir)/imp-E.sh \
		$(exsrsrcdir)/ado-x-well.sh \
		$(exsrsrcdir)/ado-x-barrier.sh    \
		$(exsrsrcdir)/ado-x-double-well.sh \
		$(exsrsrcdir)/ado-LB-double-well.sh \
		$(exsrsrcdir)/ifr-x-barrier.sh \
		$(exsrsrcdir)/ifr-Ek.sh \
		$(exsrsrcdir)/ifr-Delta.sh \
		$(exsrsrcdir)/ifr-Lambda.sh \
		$(exsrsrcdir)/ifr-LW.sh \
		$(exsrsrcdir)/LO-phonon.sh \
		$(exsrsrcdir)/LO-phonon-avg.sh \
		$(exsrsrcdir)/LO-phonon-LW-screen.sh \
		$(exsrsrcdir)/LO-phonon-Te-screen.sh \
		$(exsrsrcdir)/ac-E-T.sh \
		$(exsrsrcdir)/ac-E-em-ab.sh \
		$(exsrsrcdir)/ac-Ek.sh      \
		$(exsrsrcdir)/ac-T.sh 

exeppdatadir = $(exdatadir)/empirical-pseudopotential
exeppsrcdir  = examples/empirical-pseudopotential

#exeppdata_DATA = $(exeppsrcdir)/

expphetdatadir = $(exdatadir)/pseudopotential-het
expphetsrcdir  = examples/pseudopotential-het

expphetdata_DATA = $(expphetsrcdir)/dot-array-transport.sh

exwgdatadir = $(exdatadir)/waveguide
exwgsrcdir  = examples/waveguide

exwgdata_DATA = $(exwgsrcdir)/specific-heat.sh

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
