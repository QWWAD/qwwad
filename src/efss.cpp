/**
 * \file   efss.cpp
 * \brief  Use Schroedinger library to solve 2D SE for QW system
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \author Jonathan Cooper <el06jdc@leeds.ac.uk>
 *
 * \details  This program solves the nonrelativistic Schrodinger equation
 *           for any one-dimensional nonrelativistic potential profile.
 *
 *           It takes inputs in the form of the conduction band profile,
 *           effective mass profile, and non-parabolicity (alpha) profile
 *           all in SI units.
 */

#include <iostream>
#include <cstdlib>
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
    MATRIX_PARABOLIC,  ///< Matrix method (parabolic bands)

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

    SHOOTING_PARABOLIC,   ///< Shooting method (parabolic dispersion)
    SHOOTING_NONPARABOLIC ///< Shooting-method using nonparabolic dispersion
};

/** 
 * \brief Store for command line inputs
 * \todo  Use qclsim options class instead
 */
class FwfOptions : public Options {
    private:
        SolverType type; ///< The type of Schroedinger solver to use

    public:
        SolverType get_type() const {return type;}

        FwfOptions(int argc, char* argv[]) :
            type(MATRIX_PARABOLIC)
        {
            // No default can be set here... we want the confining potential to be used
            // by default rather than a manually-specified number!
            add_numeric_option("E-cutoff",              "Cut-off energy for solutions [meV]");
            add_numeric_option("mass",                  "The constant effective mass to use across the entire structure. "
                                                        "If unspecified, the mass profile will be read from file.");
            add_numeric_option("dE,d",       1e-3,      "Minimum separation (in energy) between states [meV]. "
                                                        "This is only used with the shooting-method solvers.");
            add_string_option ("mass-file",  "m.r",     "Filename from which effective mass profile is read. "
                                                        "This is only needed if you are not using constant effective "
                                                        "mass.");
            add_char_option   ("particle,p", 'e',       "Particle to be used: 'e', 'h' or 'l'");
            add_string_option ("alpha-file", "alpha.r", "Filename from which nonparabolicity parameter profile is read.");
            add_string_option ("v-file",     "v.r",     "Filename from which confining potential is read.");
            add_size_option   ("nst-max",     0,        "Maximum number of subbands to find.  The default (0) means "
                                                        "that all states will be found up to maximum confining potential, "
                                                        "or the cut-off energy (if specified).");
            add_numeric_option("try-energy",            "Calculate a trial wavefunction at a given energy [meV] and "
                                                        "write to file. "
                                                        "This only works with Shooting solvers.");
            add_string_option ("solver",     "matrix",  "Set the way in which the Schroedinger "
                                                        "equation is solved. See the manual for "
                                                        "a detailed list of the options");
            std::string doc = "Solve the 1D Schroedinger equation with "
                "the effective mass/envelope function approximations.";

            add_prog_specific_options_and_parse(argc, argv, doc);

            // Parse the calculation type
            std::string solver_arg(vm["solver"].as<std::string>());

            if     (!strcmp(solver_arg.c_str(), "matrix"))
                type = MATRIX_PARABOLIC;
            else if(!strcmp(solver_arg.c_str(), "matrix-full-nonparabolic"))
                type = MATRIX_FULL_NONPARABOLIC;
            else if(!strcmp(solver_arg.c_str(), "matrix-taylor-nonparabolic"))
                type = MATRIX_TAYLOR_NONPARABOLIC;
            else if(!strcmp(solver_arg.c_str(), "shooting"))
                type = SHOOTING_PARABOLIC;
            else if(!strcmp(solver_arg.c_str(), "shooting-nonparabolic"))
                type = SHOOTING_NONPARABOLIC;
            else
            {
                std::ostringstream oss;
                oss << "Cannot parse solver type: " << solver_arg;
                throw std::runtime_error(oss.str());
            }
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
        std::cerr << "No solutions found!" << std::endl;
    else
    {
        if(opt.get_verbose())
        {
            std::cout << "Energy solutions:" << std::endl;

            // Output all solutions to screen
            for(unsigned int ist=0; ist < solutions.size(); ist++)
                std::cout << ist << "\t" << std::fixed << solutions[ist].get_E() * 1000/e << " meV" << std::endl;
        }

        const char p = opt.get_char_option("particle");
        std::ostringstream energy_filename;
        energy_filename << "E" << p << ".r";
        std::ostringstream wf_prefix;
        wf_prefix << "wf_" << p;

        State::write_to_file(energy_filename.str(),
                             wf_prefix.str(),
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
    read_table_xy(opt.get_string_option("v-file").c_str(), z, V);

    const size_t nz = z.size();
    const double dz = z[1] - z[0];

    std::valarray<double> z_tmp;
    std::valarray<double> m(nz); // Band-edge effective mass [kg]
    std::valarray<double> alpha(nz); // Nonparabolicity parameter [1/J]

    // Read nonparabolicity data from file if needed
    if(opt.get_type() == MATRIX_TAYLOR_NONPARABOLIC ||
       opt.get_type() == MATRIX_FULL_NONPARABOLIC   ||
       opt.get_type() == SHOOTING_NONPARABOLIC)
        read_table_xy(opt.get_string_option("alpha-file").c_str(), z_tmp, alpha);

    // Set a constant effective mass if specified.
    // Read spatially-varying profile from file if not.
    if(opt.vm.count("mass"))
        m += opt.get_numeric_option("mass") * me;
    else
        read_table_xy(opt.get_string_option("mass-file").c_str(), z_tmp, m);

    // By default, we set the number of states automatically
    // within the range of the potential profile
    const size_t nst_max = opt.get_size_option("nst-max");

    // If we have a flat potential, the user needs to specify either
    // a cut-off energy or a number of states, otherwise we can't
    // proceed
    if(gsl_fcmp(V.max(), V.min(), 1e-6*e) == 0 && nst_max == 0 && opt.vm.count("E-cutoff") == 0)
    {
        std::cerr << "Flat potential detected in " << opt.get_string_option("v-file")
                  << ".  You must either specify a cut-off energy using --E-cutoff "
                  << "or a number of states using --nst-max" << std::endl;
        exit(EXIT_FAILURE);
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
        case MATRIX_PARABOLIC:
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
        case SHOOTING_PARABOLIC:
        case SHOOTING_NONPARABOLIC:
            se = new SchroedingerSolverShooting(m,
                                                alpha,
                                                V,
                                                z,
                                                opt.get_numeric_option("dE") * e/1000,
                                                nst_max);
    }

    // Set cut-off energy if desired
    if(opt.vm.count("E-cutoff") > 0)
        se->set_E_cutoff(opt.get_numeric_option("E-cutoff") * e/1000);

    // Output a single trial wavefunction
    if (opt.vm.count("try-energy") != 0 && (opt.get_type() == SHOOTING_PARABOLIC || opt.get_type() == SHOOTING_NONPARABOLIC))
    {
        const double E_trial = opt.get_numeric_option("try-energy") * e/1000;
        const std::valarray<double> psi = dynamic_cast<SchroedingerSolverShooting *>(se)->trial_wavefunction(E_trial);

        std::ostringstream wf_filename;
        wf_filename << "wf_" << opt.get_char_option("particle") << "E.r";
        write_table_xy(wf_filename.str().c_str(), z, psi);
    }
    else // Output all wavefunctions
    {
        std::vector<State> solutions = se->get_solutions(true);
        output(solutions, z, opt);
    }

    delete se;

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
