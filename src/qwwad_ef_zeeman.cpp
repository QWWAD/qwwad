/**
 * \file   qwwad_ef_zeeman.cpp
 * \brief  Adds Zeeman splitting term to an input potential
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavania@leeds.ac.uk>
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_math.h>
#include "struct.h"
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Computes the effective temperature T0(x)
 *
 * \details  (J.A.Gaj, C.Bodin-Deshayes, P.Peyla, J.Cibert, G.Feuillet,
 *           Y.Merle d'Aubigne, R.Romestain and A.Wasiela, Proc. 21st
 *           Int. Conf. Phys. Semiconductors 1936 (1992))
 */
static double Teff(const double x)
{
    const double F = 40.7;
    const double G = 4.6;

    return F*x/(1+G*x);
}

/**
 * \brief Compute the effective spin <Sz(x)>
 *
 * \details Data taken from 
 *
 *  J.A.Gaj, C.Bodin-Deshayes, P.Peyla, J.Cibert, G.Feuillet,
 *  Y.Merle d'Aubigne, R.Romestain and A.Wasiela, Proc. 21st
 *  Int. Conf. Phys. Semiconductors 1936 (1992)               
 *
 *  Note Gaj gives	Delta E(saturation)=6A+6B
 * 
 *  But Delta E(saturation)=x(N0alpha+N0beta)<Sz>
 *
 *  Hence <Sz> follows from Gaj's expression
 */
static double Seff(const double N0alpha,
                   const double N0beta,
                   const double x)
{
    const double A = 2488.0*e*1e-3;
    const double B = -57880.0*e*1e-3;
    const double C = 152.7;
    const double D = -20760.0*e*1e-3;
    const double E = 8.083;
    const double s = (A + B*x*x/(1.0+C*x*x) + D*x/(1.0+E*x))/(N0alpha+N0beta);

    return s;
}

/**
 * \brief Configure command-line options for the program
 */
Options configure_options(int argc, char** argv)
{
    Options opt;

    std::string doc("Find the Zeeman splitting contribution to a potential profile.");

    opt.add_option<double>     ("magneticfield,B",            0.0, "Magnetic field along growth axis [T].");
    opt.add_option<char>       ("particle,p",                 'e', "ID of particle to be used: 'e', 'h' or 'l', for "
                                                                   "electrons, heavy holes or light holes respectively.");
    opt.add_option<bool>       ("spinup",                          "Spin direction is 'up'. Down-spin assumed if this flag is not used.");
    opt.add_option<double>     ("Tl",                         1.8, "Temperature of crystal lattice [K]");
    opt.add_option<std::string>("material,M",            "cdmnte", "Material ID: Currently, only \"cdmnte\" for Cd(1-x)Mn(x)Te is supported");
    opt.add_option<std::string>("bandedgepotentialfile",  "v_b.r", "File containing baseline potential to be added to Zeeman potential");
    opt.add_option<std::string>("totalpotentialfile",       "v.r", "Filename to which the total potential is written.");
    opt.add_option<std::string>("zeemanpotentialfile",    "v_z.r", "Filename to which the Zeeman potential is written.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
};


int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto MF       = opt.get_option<double>("magneticfield"); // [T]
    const auto p        = opt.get_option<char>  ("particle");
    const auto spinup   = opt.get_option<bool>  ("spinup");        // True = 'up spin'; False = 'down spin'
    const auto Tl       = opt.get_option<double>("Tl");            // Lattice temperature [K]
    const auto material = opt.get_option<std::string>("material");

    if(material != "cdmnte")
    {
        std::cerr << "Only CdMnTe is supported at present" << std::endl;
        exit(EXIT_FAILURE);
    }

    double N0alpha = 0.220*e; // Magnetic parameter
    double N0beta  = 0.880*e; // Magnetic parameter
    double J       = 2.5;     // magnetic ion spin

    // Open baseline potential and alloy files
    arma::vec z;  // Spatial locations [m]
    arma::vec Vb; // Band-edge potential [J]

    const auto bandedgefilename = opt.get_option<std::string>("bandedgepotentialfile");
    read_table(bandedgefilename, z, Vb);

    const auto nz = z.size(); // Number of spatial points

    // TODO: Check that alloy file uses same coordinates
    arma::vec z_tmp;
    arma::vec x;     // Alloy fraction at each point
    read_table("x.r", z_tmp, x);

    arma::vec V_zeeman(nz);

    // Add field in the form +/- 3A, +/- 3B or +/- B depending on spin
    for(unsigned int iz=0; iz<nz; iz++)
    {
        const double T0 = Teff(x[iz]); // Effective temperature
        const double Sz = Seff(N0alpha,N0beta,x[iz]);

        const double y = 2*J*mu_b*MF/(kB*(Tl+T0));

        const double A=x[iz]*N0alpha*Sz*sf_brillouin(J,y)/6;
        const double B=x[iz]*N0beta *Sz*sf_brillouin(J,y)/6;

        // calculate change in potential due to magnetic field
        switch(p)
        {
            case 'e': 
                if(spinup) V_zeeman[iz] =  3*A;
                else       V_zeeman[iz] = -3*A;
                break;
            case 'h':
                if(spinup) V_zeeman[iz] =  3*B;
                else       V_zeeman[iz] = -3*B;
                break;
            case 'l': 
                if(spinup) V_zeeman[iz] = B;
                else       V_zeeman[iz] = B;
                break;
        }
    }

    const arma::vec V_total = Vb + V_zeeman;

    // Write data to file
    const auto totalpotentialfile  = opt.get_option<std::string>("totalpotentialfile");
    const auto zeemanpotentialfile = opt.get_option<std::string>("zeemanpotentialfile");
    write_table(totalpotentialfile,  z, V_total);
    write_table(zeemanpotentialfile, z, V_zeeman);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
