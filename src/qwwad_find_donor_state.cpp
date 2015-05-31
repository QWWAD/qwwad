/**
 * \file   qwwad-find-donor-state.cpp
 * \brief  Calculates state of electron attached to donor in a heterostructure potential
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
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"
#include "qwwad/donor-energy-minimiser.h"
#include "qwwad/schroedinger-solver-donor-2D.h"
#include "qwwad/schroedinger-solver-donor-3D.h"
#include "qwwad/schroedinger-solver-donor-variable.h"

using namespace QWWAD;
using namespace constants;

int main(int argc,char *argv[])
{
    Options opt;
    std::string doc("Find state of electron attached to a donor in a 2D system");

    opt.add_option<double>     ("dE,d",              1, "Energy step for Shooting solver [meV]");
    opt.add_option<double>     ("epsilon,e",     13.18, "Bulk relative permittivity");
    opt.add_option<double>     ("mass,m",        0.067, "Bulk effective mass (relative to free electron)");
    opt.add_option<double>     ("lambdastart,s",    50, "Initial value for Bohr radius search [Angstrom]");
    opt.add_option<double>     ("lambdastep,t",      1, "Step size for Bohr radius search [Angstrom]");
    opt.add_option<double>     ("lambdastop,u",     -1, "Final value for Bohr radius search [Angstrom]");
    opt.add_option<double>     ("zetastart,w",   0.001, "Initial value for symmetry parameter search");
    opt.add_option<double>     ("zetastep,x",     0.01, "Step size for symmetry parameter search");
    opt.add_option<double>     ("zetastop,y",       -1, "Final value for symmetry parameter search");
    opt.add_option<std::string>("searchmethod", "fast", "Method to use for locating parameters (\"fast\" or \"linear\")");
    opt.add_option<std::string>("symmetry",       "2D", "Symmetry of hydrogenic wave function (\"2D\", \"3D\" or \"variable\")");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    const auto delta_E = opt.get_option<double>("dE") * 1e-3*e;    // Energy increment [J]
    const auto epsilon = opt.get_option<double>("epsilon") * eps0; // Permittivity [F/m]
    const auto mstar   = opt.get_option<double>("mass") * me;      // Effective mass [kg]

    const auto lambda_start = opt.get_option<double>("lambdastart") * 1e-10; // Initial Bohr radius [m]
    const auto lambda_step  = opt.get_option<double>("lambdastep")  * 1e-10; // Bohr radius increment [m]
    const auto lambda_stop  = opt.get_option<double>("lambdastop")  * 1e-10; // Final Bohr radius [m]
    const auto zeta_start   = opt.get_option<double>("zetastart"); // Initial symmetry parameter
    const auto zeta_step    = opt.get_option<double>("zetastep");  // Symmetry parameter increment
    const auto zeta_stop    = opt.get_option<double>("zetastop");  // Final symmetry parameter

    std::valarray<double> z; // Spatial location [m]
    std::valarray<double> V; // Confining potential [J]
    read_table("v.r", z, V);

    // Read list of donor (or acceptor) positions
    std::valarray<double> r_d; // [m]
    read_table("r_d.r", r_d);

    // Solutions for each donor position
    std::valarray<double> E0(r_d.size());       // Binding energy [J]
    std::valarray<double> lambda_0(r_d.size()); // Bohr radius [m]
    std::valarray<double> zeta_0(r_d.size()); // symmetry parameter

    // Perform variational calculation for each donor/acceptor position
    for(unsigned int i_d = 0; i_d < r_d.size(); ++i_d)
    {
        // Create an initial estimate of the Schroedinger solution using a guess at lambda
        SchroedingerSolverDonor *se = 0;

        if(opt.get_option<std::string>("symmetry") == "2D")
            se = new SchroedingerSolverDonor2D(mstar, V, z, epsilon, r_d[i_d], lambda_0[i_d], delta_E);
        else if(opt.get_option<std::string>("symmetry") == "3D")
            se = new SchroedingerSolverDonor3D(mstar, V, z, epsilon, r_d[i_d], lambda_0[i_d], delta_E);
        else if(opt.get_option<std::string>("symmetry") == "variable")
        {
            se = new SchroedingerSolverDonorVariable(mstar, V, z, epsilon, r_d[i_d], lambda_0[i_d], zeta_0[i_d], delta_E);
        }
        else
        {
            std::cerr << "Unrecognised symmetry type: " << opt.get_option<std::string>("symmetry") << std::endl;
            exit(EXIT_FAILURE);
        }

        // Now, use a minimisation technique to correct the Bohr radius and find the minimum energy
        // solution
        DonorEnergyMinimiser minimiser(se, lambda_start, lambda_step, lambda_stop);
        minimiser.set_zeta_params(zeta_start, zeta_step, zeta_stop);

        const auto search_method = opt.get_option<std::string>("searchmethod");

        if(search_method == "linear")
            minimiser.minimise(MINIMISE_LINEAR);
        else if (search_method == "fast")
            minimiser.minimise(MINIMISE_FAST);
        else
        {
            std::cerr << "Unrecognised search type: " << search_method << std::endl;
            exit(EXIT_FAILURE);
        }

        // Read out the solutions now that we've minimised the energy
        const auto solutions = se->get_solutions();
        E0[i_d]              = solutions[0].get_energy();
        lambda_0[i_d]        = se->get_lambda();

        if(opt.get_option<std::string>("symmetry") == "variable")
            zeta_0[i_d] = dynamic_cast<SchroedingerSolverDonorVariable *>(se)->get_zeta();

        // Get the complete wavefunction
        const auto psi = solutions[0].get_wavefunction_samples();

        // Get the wavefunction (without the hydrogenic factor)
        const auto solutions_chi = se->get_solutions_chi();
        const auto chi = solutions_chi[0].get_wavefunction_samples();

        /* generate output filename (and open file for writing) using 
           the basis wf%i.r where the integer %i is the donor index i_d  */
        char   filename[9];     /* character string for wavefunction filename  */
        sprintf(filename,"wf%i.r",i_d);

        write_table(filename, z, psi, chi);
        delete se;
    }/* end loop over r_d */

    /* Output neutral dopant binding energies (E) and 
       Bohr radii (lambda) in meV and Angstrom respectively */
    const std::valarray<double> E_out = E0*1000.0/e;
    const std::valarray<double> r_d_out = r_d/1e-10;
    const std::valarray<double> lambda_out = lambda_0*1.0e10;
    write_table("e.r", r_d_out, E_out);
    write_table("l.r", r_d_out, lambda_out);

    if(opt.get_option<std::string>("symmetry") == "variable")
        write_table("zeta.r", r_d_out, zeta_0);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
