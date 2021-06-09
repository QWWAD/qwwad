/**
 * \file   ado.cpp
 * \brief  Interface roughness scattering rate solver
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/subband.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string doc("Find the alloy disorder scattering rate.");

    opt.add_option<bool>  ("noblocking,S",         "Disable final-state blocking.");
    opt.add_option<double>("delta",             3, "Interface roughness height [angstrom]");
    opt.add_option<double>("lambda",           50, "Interface roughness correlation length [angstrom]");
    opt.add_option<double>("mass,m",        0.067, "Band-edge effective mass (relative to free electron)");
    opt.add_option<char>  ("particle,p",      'e', "ID of particle to be used: 'e', 'h' or 'l', for "
                                                   "electrons, heavy holes or light holes respectively.");
    opt.add_option<double>("temperature,T",   300, "Temperature of carrier distribution.");
    opt.add_option<double>("Ecutoff",              "Cut-off energy for carrier distribution [meV]. If not specified, then 5kT above band-edge.");
    opt.add_option<size_t>("nki",             101, "Number of initial wave-vector samples.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto Delta   =  opt.get_option<double>("delta")*1e-10;  // Roughness height [m]
    const auto Lambda  =  opt.get_option<double>("lambda")*1e-10; // Roughness correlation length [m]
    const auto m       =  opt.get_option<double>("mass")*me;      // Band-edge effective mass [kg]
    const auto p       =  opt.get_option<char>  ("particle");     // Particle ID
    const auto T       =  opt.get_option<double>("temperature");  // Temperature [K]
    const auto b_flag  = !opt.get_option<bool>("noblocking");     // Include final-state blocking by default
    const auto nki     =  opt.get_option<size_t>("nki");          // number of ki calculations

    std::ostringstream E_filename; // Energy filename string
    E_filename << "E" << p << ".r";
    std::ostringstream wf_prefix;  // Wavefunction filename prefix
    wf_prefix << "wf_" << p;

    // Read data for all subbands from file
    std::vector<Subband> subbands = Subband::read_from_file(E_filename.str(),
                                                            wf_prefix.str(),
                                                            ".r",
                                                            m);

    // Read and set carrier distributions within each subband
    arma::vec  Ef;      // Fermi energies [J]
    arma::uvec indices; // Subband indices (garbage)
    read_table("Ef.r", indices, Ef);
    Ef *= e/1000.0; // Rescale to J

    for(unsigned int isb = 0; isb < subbands.size(); ++isb) {
        subbands[isb].set_distribution_from_Ef_Te(Ef[isb], T);
    }

    // Read potential profile
    arma::vec z;
    arma::vec V;
    read_table("v.r", z, V);

    // Read interface locations
    arma::uvec iz_I;
    read_table("interfaces.r", iz_I);

    // Read list of wanted transitions
    arma::uvec i_indices;
    arma::uvec f_indices;

    read_table("rrp.r", i_indices, f_indices);

    FILE *Favg=fopen("ifr-avg.dat","w"); // open file for output of weighted means

    // Loop over all desired transitions
    for(unsigned int itx = 0; itx < i_indices.size(); ++itx)
    {
        // State indices for this transition (NB., these are indexed from 1)
        unsigned int i = i_indices[itx];
        unsigned int f = f_indices[itx];

        // Convenience labels for each subband (NB., these are indexed from 0)
        const Subband isb = subbands[i-1];
        const Subband fsb = subbands[f-1];

        // Subband minima
        const double Ei = isb.get_E_min();
        const double Ef = fsb.get_E_min();

        // Find minimum initial wave-vector that allows scattering
        const double Efi = Ef - Ei;
        double kimin = 0.0;

        if(Efi > 0) {
            kimin = sqrt(2*m*Efi)/hBar;
        }

        double kimax = 0;
        double Ecutoff = 0.0; // Maximum kinetic energy in initial subband

        // Use user-specified value if given
        if(opt.get_argument_known("Ecutoff"))
        {
            Ecutoff = opt.get_option<double>("Ecutoff")*e/1000;

            if(Ecutoff+Ei < Ef)
            {
                std::cerr << "No scattering permitted from state " << i << "->" << f << " within the specified cut-off energy." << std::endl;
                std::cerr << "Extending range automatically" << std::endl;
                Ecutoff += Ef;
            }
        }
        // Otherwise, use a fixed, 5kT range
        else
        {
            kimax   = isb.get_k_max(T);
            Ecutoff = hBar*hBar*kimax*kimax/(2*m);

            if(Ecutoff+Ei < Ef) {
                Ecutoff += Ef;
            }
        }

        kimax = isb.get_k_at_Ek(Ecutoff);

        const double dki=(kimax-kimin)/((float)nki - 1); // step length for loop over ki

        arma::vec Wbar_integrand_ki(nki); // initialise integral for average scattering rate
        arma::vec Wif(nki);               // Scattering rate for a given initial wave vector
        arma::vec Ei_t(nki);              // Total energy of initial state (for output file) [meV]

        // Find interface-roughness scattering matrix element
        const double dz = z[1] - z[0];
        const size_t nz = z.size();
        arma::vec dV_dz(nz); // Derivative of potential profile [J/m]

        for (unsigned int iz = 1; iz < nz-1; ++iz) {
            dV_dz[iz] = (V[iz+1] - V[iz-1])/dz;
        }

        // Assume periodic boundary conditions
        dV_dz[0]    = (V[1] - V[nz-1])/dz;
        dV_dz[nz-1] = (V[0] - V[nz-2])/dz;

        const auto psi_i  = isb.psi_array();
        const auto psi_f  = fsb.psi_array();
        const arma::cx_vec psi_if = psi_i%psi_f;
        double F_if_sq = 0.0;

        // Get contributions from each interface
        // Note that we don't include the last one, since this is the edge of the system, where
        // psi = 0
        for (unsigned int I=0; I < iz_I.size()-1; ++I)
        {
            unsigned int iz_L = 0; // Lower bound of interface
            unsigned int iz_U = 0; // Upper bound of interface

            if(I != 0) {
                iz_L = (iz_I[I] + iz_I[I-1])/2;
            } else {
                iz_L = iz_I[0]/2;
            }

            iz_U = (iz_I[I] + iz_I[I+1])/2;

            arma::cx_vec F_integrand_dz(iz_U-iz_L);

            for (unsigned int iz = iz_L; iz < iz_U; ++iz) {
                F_integrand_dz[iz-iz_L] = psi_if[iz]*dV_dz[iz];
            }

            const auto F_if = integral(F_integrand_dz, dz);
            const auto abs_F_if = abs(F_if);
            F_if_sq = abs_F_if * abs_F_if;
        }

        // calculate scattering rate for all ki
        for(unsigned int iki=0;iki<nki;iki++)
        {
            const double ki=kimin+dki*iki; // carrier momentum
            const double ki_sqr = ki*ki;

            // Find energy-conserving final wave-vector
            // This should be positive if the kimin value is correct
            const double kf_sqr = ki_sqr + 2*m*(Ei - Ef)/(hBar*hBar);
            assert(kf_sqr >= 0.0);
            const double kf = sqrt(kf_sqr);
            assert(!std::isnan(kf));
            const double beta = exp(-(ki_sqr + kf_sqr)*Lambda*Lambda/4) * gsl_sf_bessel_I0(ki*kf*Lambda*Lambda/2);

            Wif[iki] = pi*m*Delta*Delta*Lambda*Lambda/(hBar*hBar*hBar) * beta * F_if_sq;

            // Include final-state blocking factor
            if (b_flag) {
                Wif[iki] *= (1 - fsb.get_occupation_at_k(kf));
            }

            Ei_t[iki] = isb.get_E_total_at_k(ki) * 1000/e;

            /* calculate Fermi-Dirac weighted mean of scattering rates over the 
               initial carrier states, note that the integral step length 
               dE=2*sqr(hBar)*ki*dki/(2m)					*/
            Wbar_integrand_ki[iki] = Wif[iki]*ki*isb.get_occupation_at_k(ki);
        } /* end ki	*/

        /* output scattering rate versus carrier energy=subband minima+in-plane
           kinetic energy						*/
        std::ostringstream filename; // output filename
        filename << "ifr" << i << f << ",r";
        write_table(filename.str(), Ei_t, Wif);

        const double Wbar = integral(Wbar_integrand_ki, dki)/(pi*isb.get_total_population());

        fprintf(Favg,"%i %i %20.17le\n", i,f,Wbar);
} /* end while over states */

fclose(Favg);	/* close weighted mean output file	*/

return EXIT_SUCCESS;
} /* end main */
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
