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
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <valarray>

#include "qwwad-options.h"
#include "qclsim_poisson_solver.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"

using namespace Leeds;
using namespace Leeds::constants;

/**
 * \brief User-settings for Poisson solver
 */
class PoissonOptions : public Options 
{
    private:
        double E; ///< External field [J/m]
        double offset; ///< Value of potential at spatial point closest to origin [V]

    public:
        PoissonOptions(int argc, char* argv[]);

        double get_field() const {return E;}
        double get_offset() const {return offset;}
        bool get_mixed() const {return vm["mixed"].as<bool>();}
        std::string get_charge_density_filename() const {return vm["charge-file"].as<std::string>();}
        bool field_applied() const {return vm.count("field");}
};

/**
 * \brief Constructor: Define and parse all user options
 *
 * \param[in] argc Number of command-line arguments
 * \param[in] argv Array of command-line arguments
 */
PoissonOptions::PoissonOptions(int argc, char* argv[]) :
    E(0.0),
    offset(0.0)
{
    add_switch       ("uncharged",               "True if there is no charge in the structure");
    add_switch       ("centred",                 "True if the potential should be pivoted around the centre of the structure");
    add_string_option("Vbasefile",               "File containing baseline potential to be added to Poisson potential");
    add_string_option("potential-file", "v_p.r", "Filename to which the Poisson potential is written.");
    add_string_option("Vtotalfile",     "v.r",   "Filename to which the total potential is written.");

    program_specific_options->add_options()
        ("field,E", po::value<double>(),
         "Set external electric field [kV/cm]. Only specify if the voltage drop needs to be fixed. Otherwise will be equal to inbuilt potential from zero-field Poisson solution.")

        ("offset", po::value<double>(),
         "Set value of potential at spatial point closest to origin. Will be zero by default.")

        ("mixed", po::bool_switch(),
         "Use mixed boundary conditions.  By default, the space-charge effect is assumed to give zero-field boundary conditions.  By supplying this option, nonzero boundary fields can exist.")

        ("charge-file",
         po::value<std::string>()->default_value("sigma.r"),
         "Set filename from which to read charge density profile.")
        ;

    std::string doc = "Find the space-charge induced potential for a "
        "one-dimensional charge profile [C/m^3].  The charge density "
        "is read from an input file, and the Poisson potential [V] is "
        "written to the output.";

    add_prog_specific_options_and_parse(argc, argv, doc);

    // Rescale external field to J/m if field specified
    if(vm.count("field"))
        E = vm["field"].as<double>() * 1000.0 * 100.0;

    if(vm.count("offset"))
        offset = vm["offset"].as<double>();
}

/**
 * \brief     main function for program
 *
 * \param[in] argc The number of command-line arguments
 * \param[in] argv Array of command-line arguments
 * 
 * \returns   Exit code of 0 signifies successful completion
 */
int main(int argc, char* argv[])
{
    PoissonOptions opt(argc, argv);

    std::valarray<double> z;
    std::valarray<double> _eps; // Low-frequency permittivity [F/m]
    read_table_xy("eps-dc.r", z, _eps);

    const size_t nz = z.size();

    std::valarray<double> z2(nz);
    std::valarray<double> rho(nz); // Charge-profile

    // Read space-charge profile, if desired
    if(!opt.get_switch("uncharged"))
        read_table_xy(opt.get_charge_density_filename().c_str(), z2, rho);

    // Convert charge density into S.I. units
    rho *= e;

    // Calculate Poisson potential due to charge within structure using cyclic boundary conditions
    const double dz = z[1] - z[0];
    std::valarray<double> phi(nz);
    if(opt.get_mixed() == true)
    {
        Poisson poisson(_eps, dz, MIXED);
        phi = poisson.solve(rho);

        // Only fix the voltage across the structure if an applied field is specified.
        // (Otherwise just return the cyclic solution!)
        if(opt.field_applied())
        {
            // Now solve the Laplace equation to find the contribution due to applied bias
            // Find voltage drop per period and take off the voltage drop from the charge discontinuity
            // within the structure. This will ensure that the voltage drop is equal to that specified
            // rather than being the sum of applied bias and voltage due to charge which is an unknown
            // quantity.
            const double V_drop = opt.get_field() * e * (z[nz-1] - z[0]) - phi[nz-1];

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
        Poisson poisson(_eps, dz, DIRICHLET);
        if(opt.field_applied())
        {
            const double V_drop = opt.get_field() * e * (z[nz-1] - z[0])*(nz+2)/nz - phi[nz-1];
            if(opt.get_verbose())
                std::cout << "Voltage drop per period: " << V_drop << "V\n";
            phi = poisson.solve(rho, V_drop);

            if(opt.get_switch("centred"))
                    phi -= V_drop/2.0;
        }
        else
            phi = poisson.solve(rho);

        if(opt.vm.count("offset"))
            phi -= opt.get_offset(); // Minus offset since potential has not yet been inverted
    }

    // Invert potential as we output in electron potential instead of absolute potential.
    phi *= -1;

    // Get field profile [V/m]
    std::valarray<double> F(z.size());
    for(unsigned int iz = 1; iz < nz-1; ++iz)
        F[iz] = (phi[iz+1] - phi[iz-1])/(dz*e);

    write_table_xy("field.r", z, F);
    write_table_xy(opt.get_string_option("potential-file").c_str(), z, phi);

    // Calculate total potential and add on the baseline
    // potential if desired
    std::valarray<double> Vtotal = phi;
    std::valarray<double> Vbase(phi.size());

    if (opt.vm.count("Vbasefile"))
    {
        std::valarray<double> zbase(phi.size());
        read_table_xy(opt.get_string_option("Vbasefile").c_str(), zbase, Vbase);

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
