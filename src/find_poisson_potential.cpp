/**
 * \file    solve_poisson.cpp
 * \brief   Solves Poisson equation to calculate space-charge induced potential
 * \author  Jonathan Cooper
 * \date    2013-02-06
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <cstdlib>
#include <valarray>

#include "qwwad-options.h"
#include "qclsim_poisson_solver.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"

using namespace Leeds;
using namespace Leeds::constants;

/**
 * \brief Get user options
 *
 * \param[in] argc Number of command-line arguments
 * \param[in] argv Array of command-line arguments
 *
 * \return The user options
 */
Options get_options(int argc, char* argv[])
{
    Options opt;

    const std::string doc("Find the Poisson potential induced by a given charge profile");

    opt.add_switch        ("uncharged",                 "True if there is no charge in the structure");
    opt.add_switch        ("centred",                   "True if the potential should be pivoted "
                                                        "around the centre of the structure");
    opt.add_switch        ("mixed",                     "Use mixed boundary conditions.  By default, "
                                                        "the space-charge effect is assumed to give "
                                                        "zero-field boundary conditions.  By supplying "
                                                        "this option, nonzero boundary fields can exist.");
    opt.add_string_option ("Vbasefile",                 "File containing baseline potential to be added to Poisson potential");
    opt.add_string_option ("potential-file", "v_p.r",   "Filename to which the Poisson potential is written.");
    opt.add_string_option ("Vtotalfile",     "v.r",     "Filename to which the total potential is written.");
    opt.add_string_option ("charge-file",    "rho.r",   "Set filename from which to read charge density profile.");
    opt.add_numeric_option("field,E",                   "Set external electric field [kV/cm]. Only specify if "
                                                        "the voltage drop needs to be fixed. Otherwise will be "
                                                        "equal to inbuilt potential from zero-field Poisson solution.");
    opt.add_numeric_option("offset",              0,    "Set potential at spatial point closest to origin [meV].");
    opt.add_switch        ("ptype",                     "Dopants are to be treated as acceptors, and wavefunctions "
                                                        "treated as hole states");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc, char* argv[])
{
    Options opt = get_options(argc, argv);

    std::valarray<double> z;
    std::valarray<double> _eps; // Low-frequency permittivity [F/m]
    read_table("eps-dc.r", z, _eps);

    const size_t nz = z.size();

    std::valarray<double> z2(nz);
    std::valarray<double> rho(nz); // Charge-profile

    // Read space-charge profile, if desired
    if(!opt.get_switch("uncharged"))
        read_table(opt.get_string_option("charge-file").c_str(), z2, rho);

    // Convert charge density into S.I. units
    rho *= e;

    // If we're using a p-type system, invert the charge profile so we have
    // a positive energy scale
    if(opt.get_switch("ptype"))
        rho *= -1;

    // Calculate Poisson potential due to charge within structure using cyclic boundary conditions
    const double dz = z[1] - z[0];
    std::valarray<double> phi(nz);
    if(opt.get_switch("mixed"))
    {
        Poisson poisson(_eps, dz, MIXED);
        phi = poisson.solve(rho);

        // Only fix the voltage across the structure if an applied field is specified.
        // (Otherwise just return the cyclic solution!)
        if(opt.vm.count("field"))
        {
            // Now solve the Laplace equation to find the contribution due to applied bias
            // Find voltage drop per period and take off the voltage drop from the charge discontinuity
            // within the structure. This will ensure that the voltage drop is equal to that specified
            // rather than being the sum of applied bias and voltage due to charge which is an unknown
            // quantity.
            const double field  = opt.get_numeric_option("field") * 1000.0 * 100.0; // V/m
            const double V_drop = field * e * (z[nz-1] - z[0]) - phi[nz-1];

            if(opt.get_verbose())
                std::cout << "Voltage drop per period: " << V_drop << "V\n";

            // Instantiate Poisson class to solve the Laplace equation
            Poisson laplace(_eps, dz, DIRICHLET);
            phi += laplace.solve_laplace(V_drop);

            if(opt.get_switch("centred"))
                    phi -= V_drop/2.0;
        }
    }
    else
    {
        if(opt.vm.count("field"))
        {
            Poisson poisson(_eps, dz, DIRICHLET);
            const double field  = opt.get_numeric_option("field") * 1000.0 * 100.0; // V/m
            const double V_drop = field * e * (z[nz-1] - z[0])*(nz+2)/nz - phi[nz-1];
            if(opt.get_verbose())
                std::cout << "Voltage drop per period: " << V_drop << "V\n";
            phi = poisson.solve(rho, V_drop);

            if(opt.get_switch("centred"))
                    phi -= V_drop/2.0;
        }
        else
        {
            Poisson poisson(_eps, dz, ZERO_FIELD);
            phi = poisson.solve(rho);
        }

        phi -= opt.get_numeric_option("offset") * e/1000; // Minus offset since potential has not yet been inverted
    }

    // Invert potential as we output in electron potential instead of absolute potential.
    phi *= -1;

    // Get field profile [V/m]
    std::valarray<double> F(z.size());
    for(unsigned int iz = 1; iz < nz-1; ++iz)
        F[iz] = (phi[iz+1] - phi[iz-1])/(2*dz*e);

    write_table_xy("field.r", z, F);
    write_table_xy(opt.get_string_option("potential-file").c_str(), z, phi);

    // Calculate total potential and add on the baseline
    // potential if desired
    std::valarray<double> Vtotal = phi;
    std::valarray<double> Vbase(phi.size());

    if (opt.vm.count("Vbasefile"))
    {
        std::valarray<double> zbase(phi.size());
        read_table(opt.get_string_option("Vbasefile").c_str(), zbase, Vbase);

        // TODO: Add more robust checking of z, zbase identicality here
        if(zbase.size() != z.size())
        {
            std::ostringstream oss;
            oss << "Baseline and Poisson potential profiles have different lengths (" << zbase.size() << ") and (" << z.size() << " respectively";
            throw std::runtime_error(oss.str());
        }

        Vtotal = phi + Vbase;
    }

    write_table_xy(opt.get_string_option("Vtotalfile").c_str(), z, Vtotal);
    
    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
