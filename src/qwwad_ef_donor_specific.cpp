/**
 * \file   qwwad_ef_donor_specific.cpp
 * \brief  Calculates state of electron attached to donor in a heterostructure potential
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details Implements a variational technique to calculate the
 *          uncorrelated one particle energies of an electron attatched to a 
 *          single donor at any position, in any user supplied potential.  
 */

#include <cstdlib>
#include <iostream>
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/options.h"
#include "qwwad/donor-energy-minimiser-linear.h"
#include "qwwad/donor-energy-minimiser-fast.h"
#include "qwwad/schroedinger-solver-donor-2D.h"
#include "qwwad/schroedinger-solver-donor-3D.h"
#include "qwwad/schroedinger-solver-donor-variable.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Configure command-line options
 */
auto configure_options(int argc, char** argv) -> Options
{
    Options opt;
    std::string doc("Find state of electron attached to a donor in a 2D system");

    opt.add_option<double>     ("dE,d",                   1, "Energy step for Shooting solver [meV]");
    opt.add_option<double>     ("dcpermittivity,e",   13.18, "Bulk relative permittivity");
    opt.add_option<double>     ("mass,m",             0.067, "Bulk effective mass (relative to free electron)");
    opt.add_option<double>     ("lambdastart,s",         50, "Initial value for Bohr radius search [Angstrom]");
    opt.add_option<double>     ("lambdastep,t",           1, "Step size for Bohr radius search [Angstrom]");
    opt.add_option<double>     ("lambdastop,u",          -1, "Final value for Bohr radius search [Angstrom]");
    opt.add_option<double>     ("donorposition,r",           "Location of donor ion [Angstrom]");
    opt.add_option<double>     ("zetastart,w",        0.001, "Initial value for symmetry parameter search");
    opt.add_option<double>     ("zetastep,x",          0.01, "Step size for symmetry parameter search");
    opt.add_option<double>     ("zetastop,y",            -1, "Final value for symmetry parameter search");
    opt.add_option<std::string>("searchmethod",      "fast", R"(Method to use for locating parameters ("fast" or "linear"))");
    opt.add_option<std::string>("symmetry",            "2D", R"(Symmetry of hydrogenic wave function ("2D", "3D" or "variable"))");
    opt.add_option<std::string>("totalpotentialfile", "v.r", "Filename from which the total potential is read.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto delta_E = opt.get_option<double>("dE") * 1e-3*e;    // Energy increment [J]
    const auto epsilon = opt.get_option<double>("dcpermittivity") * eps0; // Permittivity [F/m]
    const auto mstar   = opt.get_option<double>("mass") * me;      // Effective mass [kg]

    // Read symmetry and minimiser types
    const auto symmetry_string = opt.get_option<std::string>("symmetry");
    const auto search_method   = opt.get_option<std::string>("searchmethod");

    const auto lambda_start = opt.get_option<double>("lambdastart") * 1e-10; // Initial Bohr radius [m]
    const auto lambda_step  = opt.get_option<double>("lambdastep")  * 1e-10; // Bohr radius increment [m]
    const auto lambda_stop  = opt.get_option<double>("lambdastop")  * 1e-10; // Final Bohr radius [m]
    const auto zeta_start   = opt.get_option<double>("zetastart"); // Initial symmetry parameter
    const auto zeta_step    = opt.get_option<double>("zetastep");  // Symmetry parameter increment
    const auto zeta_stop    = opt.get_option<double>("zetastop");  // Final symmetry parameter

    arma::vec z; // Spatial location [m]
    arma::vec V; // Confining potential [J]
    const auto totalpotentialfile = opt.get_option<std::string>("totalpotentialfile");
    read_table(totalpotentialfile, z, V);

    // Get donor location [m].  If unspecified, assume it's in the middle
    auto r_d = 0.0;

    if (opt.get_argument_known("donorposition") > 0)
        r_d = opt.get_option<double>("donorposition") * 1e-10;
    else
        r_d = (z[z.size()-1] + z[0])/2.0;

    // Initial estimates of the orbital geometry
    auto lambda_0 = lambda_start; // Bohr radius [m]
    auto zeta_0   = zeta_start;   // symmetry parameter

    // Create an initial estimate of the Schroedinger solution
    SchroedingerSolverDonor *se = nullptr;

    if(symmetry_string == "2D")
        se = new SchroedingerSolverDonor2D(mstar, V, z, epsilon, r_d, lambda_0, delta_E);
    else if(symmetry_string == "3D")
        se = new SchroedingerSolverDonor3D(mstar, V, z, epsilon, r_d, lambda_0, delta_E);
    else if(symmetry_string == "variable")
        se = new SchroedingerSolverDonorVariable(mstar, V, z, epsilon, r_d, lambda_0, zeta_0, delta_E);
    else
    {
        std::cerr << "Unrecognised symmetry type: " << opt.get_option<std::string>("symmetry") << std::endl;
        exit(EXIT_FAILURE);
    }

    // Now, use a minimiser to correct the orbital and find the minimum energy solution
    DonorEnergyMinimiser *minimiser = nullptr;

    if(search_method == "linear")
        minimiser = new DonorEnergyMinimiserLinear(se, lambda_start, lambda_step, lambda_stop);
    else if (search_method == "fast")
        minimiser = new DonorEnergyMinimiserFast(se, lambda_start, lambda_step, lambda_stop);
    else
    {
        std::cerr << "Unrecognised search type: " << search_method << std::endl;
        exit(EXIT_FAILURE);
    }

    minimiser->set_zeta_params(zeta_start, zeta_step, zeta_stop);
    minimiser->minimise();

    // Read out the solutions now that we've minimised the energy
    // Note that we only calculate the ground state at the moment
    const auto solutions = se->get_solutions();
    const auto E0        = solutions[0].get_energy();
    lambda_0             = se->get_lambda();

    if(opt.get_option<std::string>("symmetry") == "variable")
        zeta_0 = dynamic_cast<SchroedingerSolverDonorVariable *>(se)->get_zeta();

    // Get the complete wavefunction
    const auto psi = solutions[0].get_wavefunction_samples();

    // Get the wavefunction (without the hydrogenic factor)
    const auto solutions_chi = se->get_solutions_chi();
    const auto chi = solutions_chi[0].get_wavefunction_samples();

    // Save the ground-state wavefunction
    write_table("wf_e1.r", z, psi);

    // Save the ground-state wavefunction (no hydrogenic factor)
    write_table("wf_chi_e1.r", z, chi);

    // Output the search log
    write_table("searchlog.r",
                minimiser->get_lambda_history(),
                minimiser->get_zeta_history(),
                minimiser->get_E_history());

    delete minimiser;
    delete se;

    // Output neutral dopant binding energies (E) and 
    // Bohr radii (lambda) in meV and Angstrom respectively
    std::vector<double> indices(1,1);
    std::vector<double> E_out(1,E0*1000.0/e);
    std::vector<double> lambda_out(1,lambda_0*1.0e10);
    write_table("Ee.r", indices, E_out);
    write_table("l.r",  indices, lambda_out);

    if(opt.get_option<std::string>("symmetry") == "variable")
    {
        std::vector<double> zeta_out(1,zeta_0);
        write_table("zeta.r", indices, zeta_out);
    }

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
