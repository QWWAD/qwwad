/**
 * \file   d02D.cpp
 * \brief  Calculates state of electron attached to donor in a user-specified potential
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details Implements a variational technique to calculate the
 *          uncorrelated one particle energies of an electron attatched to a 
 *          single donor at any position, in any user supplied potential.  
 *          The potential is read from the file v.r
 *
 *          This version is a single variational parameter calculation---the 
 *          two-dimensional (2D) approximation trial wavefunction:
 *
 *          \f[
 *            \Psi(z)=\chi(z) \exp(-r^{\prime\prime}/\lambda)
 *          \f]
 *
 *          where \f$r^{\prime\prime} = \sqrt(x^2+y^2)\f$
 *
 *   Input files:
 *  r_d.r		donor (or acceptor positions)
 *  v.r		one-dimensional potential
 *
 *  Output files:
 *  e.r		total energies for each r_d
 *  l.r		Bohr radii (lambda) for each r_d
 *  wfn.r		wave functions, both Psi and chi, n=0,1,2..
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"
#include "qwwad-donor-energy-minimiser.h"
#include "qwwad-schroedinger-donor-2D.h"
#include "qwwad-schroedinger-donor-3D.h"

using namespace Leeds;
using namespace constants;

int main(int argc,char *argv[])
{
    Options opt;
    std::string doc("Find state of electron attached to a donor in a 2D system");

    opt.add_numeric_option("dE,d",              1, "Energy step for Shooting solver [meV]");
    opt.add_numeric_option("epsilon,e",     13.18, "Bulk relative permittivity");
    opt.add_numeric_option("mass,m",        0.067, "Bulk effective mass (relative to free electron)");
    opt.add_numeric_option("lambdastart,s",    50, "Initial value for Bohr radius search [Angstrom]");
    opt.add_numeric_option("lambdastep,t",      1, "Step size for Bohr radius search [Angstrom]");
    opt.add_numeric_option("lambdastop,u",     -1, "Final value for Bohr radius search [Angstrom]");
    opt.add_string_option ("lambdasearch", "fast", "Method to use for locating Bohr radius (\"fast\" or \"linear\")");
    opt.add_string_option ("symmetry",       "2D", "Symmetry of hydrogenic wave function (\"2D\" or \"3D\")");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    /* computational default values */
    const double delta_E = opt.get_numeric_option("dE") * 1e-3*e;    // Energy increment [J]
    const double epsilon = opt.get_numeric_option("epsilon") * eps0; // Permittivity [F/m]
    const double mstar   = opt.get_numeric_option("mass") * me;      // Effective mass [kg]

    const double lambda_start = opt.get_numeric_option("lambdastart") * 1e-10; // Initial Bohr radius [m]
    const double lambda_step  = opt.get_numeric_option("lambdastep")  * 1e-10; // Bohr radius increment [m]
    const double lambda_stop  = opt.get_numeric_option("lambdastop")  * 1e-10; // Final Bohr radius [m]

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> V; // Confining potential [J]
    read_table_xy("v.r", z, V);

    // Read list of donor (or acceptor) positions
    std::valarray<double> r_d; // [m]
    read_table_x("r_d.r", r_d);

    // Solutions for each donor position
    std::valarray<double> E0(r_d.size());       // Binding energy [J]
    std::valarray<double> lambda_0(r_d.size()); // Bohr radius [m]

    // Perform variational calculation for each donor/acceptor position
    for(unsigned int i_d = 0; i_d < r_d.size(); ++i_d)
    {
        // Create an initial estimate of the Schroedinger solution using a guess at lambda
        SchroedingerSolverDonor *se = 0;

        if(opt.get_string_option("symmetry") == "2D")
            se = new SchroedingerSolverDonor2D(mstar, V, z, epsilon, r_d[i_d], lambda_0[i_d], delta_E, 1);
        else if(opt.get_string_option("symmetry") == "3D")
            se = new SchroedingerSolverDonor3D(mstar, V, z, epsilon, r_d[i_d], lambda_0[i_d], delta_E, 1);
        else
        {
            std::cerr << "Unrecognised symmetry type: " << opt.get_string_option("symmetry") << std::endl;
            exit(EXIT_FAILURE);
        }

        // Now, use a minimisation technique to correct the Bohr radius and find the minimum energy
        // solution
        DonorEnergyMinimiser minimiser(se, lambda_start, lambda_step, lambda_stop);

        if(opt.get_string_option("lambdasearch") == "linear")
            minimiser.minimise(MINIMISE_LINEAR);
        else if (opt.get_string_option("lambdasearch") == "fast")
            minimiser.minimise(MINIMISE_FAST);
        else
        {
            std::cerr << "Unrecognised search type: " << opt.get_string_option("lambdasearch") << std::endl;
            exit(EXIT_FAILURE);
        }

        // Read out the solutions now that we've minimised the energy
        std::vector<State> solutions = se->get_solutions();
        E0[i_d]                      = solutions[0].get_E();
        lambda_0[i_d]                = se->get_lambda();

        // Get the complete wavefunction
        std::valarray<double> psi(solutions[0].psi_array());

        // Get the wavefunction (without the hydrogenic factor)
        std::vector<State> solutions_chi = se->get_solutions_chi();
        std::valarray<double> chi(solutions_chi[0].psi_array());

        /* generate output filename (and open file for writing) using 
           the basis wf%i.r where the integer %i is the donor index i_d  */
        char   filename[9];     /* character string for wavefunction filename  */
        sprintf(filename,"wf%i.r",i_d);

        write_table_xyz(filename, z, psi, chi);
        delete se;
    }/* end loop over r_d */

    /* Output neutral dopant binding energies (E) and 
       Bohr radii (lambda) in meV and Angstrom respectively */
    const std::valarray<double> E_out = E0*1000.0/e;
    const std::valarray<double> r_d_out = r_d/1e-10;
    const std::valarray<double> lambda_out = lambda_0*1.0e10;
    write_table_xy("e.r", r_d_out, E_out);
    write_table_xy("l.r", r_d_out, lambda_out);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :