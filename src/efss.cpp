/**
 * \file   efss.cpp
 * \brief  Use Schroedinger library to solve 2D SE for QW system
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \author Jonathan Cooper <el06jdc@leeds.ac.uk>
 *
 * \details  This program solves the Schrodinger equation for any 
 *            two-dimensional nonrelativistic quantum-well system (QWS).
 *
 *            It takes inputs in the form of the conduction band profile,
 *            effective mass profile, and non-parabolicity (alpha) profile
 *            all in SI units, and all taken from the files Ve.dat, Me.dat,
 *            and alpha.dat respectively.
 *
 *            The output consists of the energy solutions being outputted to the
 *            command line (for inital inspection) and the following output files
 *            being generated:
 *            - Energy solutions: Ee.r (in meV)
 *            - WF's in numbered files i.e:
 *              - wf_e1.r
 *              - wf_e2.r...etc.
 *
 *            The effective mass can be treated as either a constant value,
 *            a spatially-varying value, or a nonparabolic effective masses.
 */

#include <iostream>
#include <cstdlib>
#include <error.h>
#include "qwwad-options.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-linalg.h"
#include "qwwad-schroedinger-full.h"
#include "qwwad-schroedinger-shooting.h"
#include "qwwad-schroedinger-taylor.h"
#include "qwwad-schroedinger-tridiagonal.h"

using namespace Leeds;
using namespace Leeds::constants;

/** 
 * \brief The type of solver to use
 */
enum SolverType {
    MATRIX_CONSTANT_MASS,  ///< Matrix method using constant effective mass across entire structure
    MATRIX_VARIABLE_MASS,  ///< Matrix method using Spatially-varying effective mass through structure (parabolic)

    /**
     * \brief   Energy and spatially-dependent effective mass (nonparabolic dispersion)
     *
     * \details This uses the full matrix method, after Cooper et al., J. Appl. Phys. (2010)
     *          It gives a direct and exact solution but is extremely slow!
     */
    MATRIX_FULL_NONPARABOLIC,

    /**
     * \brief	Approximate method for solving nonparabolic Schroedinger equation
     *
     * \details	This approximate method uses a Taylor expansion of the nonparabolic effective
     * 		mass to simplify the eigenvalue problem obtained from the Schroedinger equation.
     * 		(Details of which can be found in Alharbi, Opt. Quant. Electron., 40, 551-559 (2008).
     *
     * 		This approximation is only valid for states which are energetically close the
     * 		conduction band edge. For states whose energy above the conduction band edge is
     * 		comparable with the band gap the effective mass is overestimated causing the energy of
     * 		the state to be lower than expected. Furthermore, the effective mass is increasingly
     * 		overestimated with increasing energy to a point where, at an energy approximately
     * 		one band gap above the conduction band edge, many states bunch together and become
     * 		degenerate.
     *
     * 		Users that are simulating deep quantum wells, where the conduction band offset is
     * 		comparable to that of the band gap, are advised to use the 'nonparabolic' option with
     * 		'fwf' in order to stop this situation occurring.
     */
    MATRIX_TAYLOR_NONPARABOLIC,

    SHOOTING_CONSTANT_MASS, ///< Shooting method using constant mass
    SHOOTING_VARIABLE_MASS, ///< Shooting-method using spatially-varying mass
    SHOOTING_NONPARABOLIC   ///< Shooting-method using nonparabolic dispersion
};

/** 
 * \brief Store for command line inputs
 * \todo  Use qclsim options class instead
 */
class FwfOptions : public Options {
    private:
        SolverType type; ///< The type of Schroedinger solver to use

    public:
        void print() const {};
        SolverType get_type() const {return type;}
        double get_mass() const {return vm["mass"].as<double>() * me;}
        std::string get_mass_filename() const {return vm["mass-file"].as<std::string>();}
        std::string get_alpha_filename() const {return vm["alpha-file"].as<std::string>();}
        std::string get_potential_filename() const {return vm["v-file"].as<std::string>();}
        size_t      get_nst_max() const {return vm["nst-max"].as<size_t>();}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the minimum energy spacing between states [J]
        double get_dE() const {return vm["dE"].as<double>()*1e-3*e;}

        /// \returns true if we want to just calculate a trial wavefunction at a given energy
        bool try_energy() const {return (vm.count("try-energy") == 1);}

        /// \returns the energy at which to calculate a trial wavefunction [J]
        double get_trial_energy() const {return vm["try-energy"].as<double>() * 1e-3*e;}

        FwfOptions(int argc, char* argv[]) :
            type(MATRIX_VARIABLE_MASS)
        {
            // No default can be set here... we want the confining potential to be used
            // by default rather than a manually-specified number!
            add_numeric_option("E-cutoff",          "Cut-off energy for solutions [meV]");

            program_specific_options->add_options()
                ("solver", po::value<std::string>()->default_value("matrix-variable-mass"),
                 "Set the way in which the Schroedinger equation is solved.  You can "
                 "use one of the following methods:\n"
                 "\n"
                 "    matrix-variable-mass       \tMatrix solver, with spatially varying, energy-independent effective mass\n\n"
                 "    matrix-constant-mass       \tMatrix solver, with constant effective mass\n\n"
                 "    matrix-full-nonparabolic   \tMatrix solver, using slow, but very accounting for nonparabolic dispersion\n\n"
                 "    matrix-taylor-nonparabolic \tMatrix solver, using fast accounting for nonparabolic dispersion, but breaks down as E(state) - E(band edge) approaches the bandgap\n\n"
                 "    shooting-variable-mass     \tShooting solver, with spatially varying, energy-independent effective mass\n\n"
                 "    shooting-constant-mass     \tShooting solver, using constant effective mass\n\n"
                 "    shooting-nonparabolic      \tShooting solver, using nonparabolic dispersion\n\n")

                ("mass",
                 po::value<double>()->default_value(0.067),
                 "The constant effective mass to use across the entire structure. "
                 "This option only has an effect when used with the matrix-constant-mass "
                 "or shooting-constant-mass solvers.")

                ("dE,d", po::value<double>()->default_value(1e-3),
                 "Minimum separation (in energy) between states [meV]. "
                 "This is only used with the shooting-method solvers.")

                ("mass-file", 
                 po::value<std::string>()->default_value("m.r"),
                 "Set filename from which effective mass profile is read. "
                 "This is only needed if you are not using constant effective "
                 "mass.")

                ("particle,p", po::value<char>()->default_value('e'),
                 "Particle to be used: 'e', 'h' or 'l'")

                ("alpha-file",
                 po::value<std::string>()->default_value("alpha.r"),
                 "Set filename from which nonparabolicity parameter profile "
                 "is read.")

                ("v-file",
                 po::value<std::string>()->default_value("v.r"),
                 "Set filename from which the confining potential is read.")

                ("nst-max",
                 po::value<size_t>()->default_value(0),
                 "Set maximum number of subbands to find.  The default (0) means "
                 "that all states will be found up to maximum confining potential, "
                 "or the cut-off energy (if specified).")

                ("try-energy", po::value<double>(),
                 "Calculate a trial wavefunction at a given energy [meV] and write to "
                 "file.")
                ;

            std::string doc = "Solve the 1D Schroedinger equation with "
                "the effective mass/envelope function approximations.";

            add_prog_specific_options_and_parse(argc, argv, doc);

            // Parse the calculation type
            std::string solver_arg(vm["solver"].as<std::string>());

            if     (!strcmp(solver_arg.c_str(), "matrix-constant-mass"))
                type = MATRIX_CONSTANT_MASS;
            else if(!strcmp(solver_arg.c_str(), "matrix-variable-mass"))
                type = MATRIX_VARIABLE_MASS;
            else if(!strcmp(solver_arg.c_str(), "matrix-full-nonparabolic"))
                type = MATRIX_FULL_NONPARABOLIC;
            else if(!strcmp(solver_arg.c_str(), "matrix-taylor-nonparabolic"))
                type = MATRIX_TAYLOR_NONPARABOLIC;
            else if(!strcmp(solver_arg.c_str(), "shooting-constant-mass"))
                type = SHOOTING_CONSTANT_MASS;
            else if(!strcmp(solver_arg.c_str(), "shooting-variable-mass"))
                type = SHOOTING_VARIABLE_MASS;
            else if(!strcmp(solver_arg.c_str(), "shooting-nonparabolic"))
                type = SHOOTING_NONPARABOLIC;
            else
            {
                std::ostringstream oss;
                oss << "Cannot parse solver type: " << solver_arg;
                throw std::runtime_error(oss.str());
            }

            if(get_verbose()) print();	
        }
};

/**
 * Function to format and output solutions
 *
 * \param[in] solutions The set of states to output
 * \param[in] z         Spatial positions [m]
 * \param[in] opt       User options
 *
 * \details   Outputs energy solutions to the command line
 *            before dumping them to the file 'Ee.dat'.
 *            Also dumps WFs to separate (numbered) files:
 *            IE: wf_e1.dat
 *                wf_e2.dat... etc.
 */
static void output(const std::vector<State>    &solutions, 
                   const std::valarray<double> &z,
                   const FwfOptions            &opt)
{
    // Check solutions were found
    if(solutions.empty())
        printf("No solutions found!!!\n");
    else
    {
        if(opt.get_verbose())
        {
            printf("Energy solutions:\n");

            // Output all solutions to screen
            for(unsigned int ist=0; ist < solutions.size(); ist++)
                printf("%u\t%9.3f meV\n", ist, solutions[ist].get_E()/(1e-3*e));
        }

        const char   p       = opt.get_particle();
        char energy_filename[9];
        sprintf(energy_filename,"E%c.r",p);
        char wf_prefix[9];
        sprintf(wf_prefix,"wf_%c",p);

        State::write_to_file(energy_filename,
                             wf_prefix,
                             ".r",
                             solutions,
                             z,
                             true);
    }
}

int main(int argc, char *argv[]){
    const FwfOptions opt(argc, argv);

    // Read data from file
    std::valarray<double> z; // Spatial locations [m]
    std::valarray<double> V; // Potential profile [J]
    read_table_xy(opt.get_potential_filename().c_str(), z, V);

    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    std::valarray<double> z_tmp;
    std::valarray<double> m(nz); // Band-edge effective mass [kg]
    std::valarray<double> alpha(nz); // Nonparabolicity parameter [1/J]

    // Decide whether to read mass and nonparabolicity data depending on
    // the solver we're using
    switch(opt.get_type()) {
        case MATRIX_TAYLOR_NONPARABOLIC:
        case MATRIX_FULL_NONPARABOLIC:
        case SHOOTING_NONPARABOLIC:
            read_table_xy(opt.get_alpha_filename().c_str(), z_tmp, alpha);
        case MATRIX_VARIABLE_MASS:
        case SHOOTING_VARIABLE_MASS:
            read_table_xy(opt.get_mass_filename().c_str(), z_tmp, m);
            break;
        case MATRIX_CONSTANT_MASS:
        case SHOOTING_CONSTANT_MASS:
            m += opt.get_mass();
            break;
        default:
            error(EXIT_FAILURE, 0, "Unrecognised effective mass type");
    }

    // By default, we set the number of states automatically
    // within the range of the potential profile
    const size_t nst_max = opt.get_nst_max();

    // If we have a flat potential, we need to take care!
    if(gsl_fcmp(V.max(), V.min(), 1e-6*e) == 0 && nst_max == 0)
    {
        if(opt.get_verbose())
            std::cout << "Flat potential in " << opt.get_potential_filename()
                      << ".  Will assume this is an infinite quantum well " << std::endl;
    }

    // Print out some information about the calculation if requested
    if(opt.get_verbose())
    {
        if (nst_max == 0)
            std::cout << "Searching for solutions between " << V.min()/(e*1e-3)
                      << " meV and " << V.max()/(e*1e-3) << std::endl;
        else
            std::cout << "Searching for " << nst_max
                      << " solutions above the band-edge." << std::endl;

        std::cout << nz << " points in spatial profile with spatial step of "
                    << dz*1e9 << "nm." << std::endl;
    }

    SchroedingerSolver *se = NULL; // Solver for Schroedinger equation

    switch(opt.get_type())
    {
        case MATRIX_CONSTANT_MASS:
        case MATRIX_VARIABLE_MASS:
            se = new SchroedingerSolverTridiag(m,
                                               V,
                                               z,
                                               nst_max);
            break;
        case MATRIX_FULL_NONPARABOLIC:
            se = new SchroedingerSolverFull(m,
                                            alpha,
                                            V,
                                            z,
                                            nst_max);
            break;
        case MATRIX_TAYLOR_NONPARABOLIC:
            se = new SchroedingerSolverTaylor(m,
                                              alpha,
                                              V,
                                              z,
                                              nst_max);
            break;
        case SHOOTING_CONSTANT_MASS:
        case SHOOTING_VARIABLE_MASS:
        case SHOOTING_NONPARABOLIC:
            se = new SchroedingerSolverShooting(m,
                                                alpha,
                                                V,
                                                z,
                                                opt.get_dE(),
                                                nst_max);
    }

    // Set cut-off energy if desired
    if(opt.vm.count("E-cutoff") > 0)
        se->set_E_cutoff(opt.get_numeric_option("E-cutoff") * e/1000);

    if (opt.try_energy() && (opt.get_type() == SHOOTING_CONSTANT_MASS || opt.get_type() == SHOOTING_VARIABLE_MASS))
    {
        const double E_trial = opt.get_trial_energy();
        const std::valarray<double> psi = dynamic_cast<SchroedingerSolverShooting *>(se)->trial_wavefunction(E_trial);

        char wf_filename[9];
        sprintf(wf_filename, "wf_%cE.r", opt.get_particle());
        write_table_xy(wf_filename, z, psi);
    }
    else
    {
        std::vector<State> solutions = se->get_solutions(true);
        output(solutions, z, opt);
    }

    delete se;

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
