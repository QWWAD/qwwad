/*==================================================================
           sradp Scattering Rate Acoustic Deformation Potential
  ==================================================================*/

/* This program calculates the carrier-acoustic deformation potential
   (acoustic phonon) scattering rate for parallel parabolic
   intra- and intersubband events.  The required rates are provided
   by the user in the file `rrp.r'.  The other necessary inputs are listed
   below.


	Input files:		rrp.r	contains required rates
				wf_xy.r	x=particle y=state
				Ex.r	x=particle, energies

	Output files:		AC[a,e]if.r	absorption and emission rates


    Paul Harrison, March 2000 
    									*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <complex>
#include "qwwad/options.h"
#include "qwwad/file-io.h"
#include "qwwad/subband.h"
#include "qwwad/constants.h"
#include "qwwad/maths-helpers.h"

using namespace QWWAD;
using namespace constants;

static void ff_table(double         dKz,
                     const Subband &isb,
                     const Subband &fsb,
                     unsigned int   nKz,
                     arma::vec     &Kz,
                     arma::vec     &Gifsqr);

/* This function outputs the formfactors into files	*/
static void ff_output(const arma::vec &Kz,
                      const arma::vec &Gifsqr,
                      unsigned int        i,
                      unsigned int        f)
{
    std::ostringstream filename;	/* output filename				*/
    filename << "G" << i << f << ".r";	
    write_table(filename.str(), Kz, Gifsqr);
}

auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string doc("Find the acoustic-phonon deformation potential scattering rate.");

    opt.add_option<bool>  ("outputff,a",              "Output form-factors to file.");
    opt.add_option<bool>  ("noblocking,b",            "Disable final-state blocking.");
    opt.add_option<double>("latticeconst,A",    5.65, "Lattice constant in growth direction [angstrom]");
    opt.add_option<double>("Ephonon,E",         2.0,  "Energy of acoustic phonon [meV]");
    opt.add_option<double>("vs,v"      ,     5117.0,  "Speed of sound [m/s]");
    opt.add_option<double>("mass,m",          0.067,  "Band-edge effective mass (relative to free electron)");
    opt.add_option<double>("density,d",       5317.5, "Mass density [kg/m^3]");
    opt.add_option<double>("Da,D",               7.0, "Acoustic deformation potential [eV]");
    opt.add_option<char>  ("particle,p",        'e',  "ID of particle to be used: 'e', 'h' or 'l', for "
                                                      "electrons, heavy holes or light holes respectively.");
    opt.add_option<double>("Te",                300,  "Carrier temperature [K].");
    opt.add_option<double>("Tl",                300,  "Lattice temperature [K].");
    opt.add_option<double>("Ecutoff",                 "Cut-off energy for carrier distribution [meV]. If not specified, then 5kT above band-edge.");
    opt.add_option<size_t>("nki",               301,  "Number of initial wave-vector samples.");
    opt.add_option<size_t>("nkz",               301,  "Number of phonon wave-vector samples.");
    opt.add_option<size_t>("ntheta",            101,  "Number of strips in theta angle integration");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto ff_flag =  opt.get_option<bool>("outputff");               // True if formfactors are wanted
    const auto b_flag  = !opt.get_option<bool>("noblocking");             // Final state blocking on by default
    const auto A0      =  opt.get_option<double>("latticeconst") * 1e-10; // Lattice constant [m]
    const auto rho     =  opt.get_option<double>("density");              // Mass density [kg/m^3]
    const auto Ephonon =  opt.get_option<double>("Ephonon")*e/1000;       // Acoustic phonon energy [J]
    const auto Da      =  opt.get_option<double>("Da")*e;                 // Acoustic deformation potential [J]
    const auto m       =  opt.get_option<double>("mass")*me;              // Band-edge effective mass [kg]
    const auto p       =  opt.get_option<char>("particle");               // Particle ID
    const auto Te      =  opt.get_option<double>("Te");                   // Carrier temperature [K]
    const auto Tl      =  opt.get_option<double>("Tl");                   // Lattice temperature [K]
    const auto Vs      =  opt.get_option<double>("vs");                   // Speed of sound [m/s]
    const auto nki     =  opt.get_option<size_t>("nki");                  // number of ki calculations
    const auto nKz     =  opt.get_option<size_t>("nkz");                  // number of Kz calculations
    const auto ntheta  =  opt.get_option<size_t>("ntheta");               // number of samples over angle

    FILE	*FACa;		/* pointer to absorption output file		*/
    FILE	*FACe;		/* pointer to emission   output file		*/

    // calculate step length in phonon wave-vector
    const double dKz=2/(A0*nKz); // Taken range of phonon integration as 2/A0

    const double dtheta=pi/static_cast<double>(ntheta-1); // theta integration from 0 to pi

    // Can save a bit of time by calculating cosines in advance
    arma::vec cos_theta(ntheta);

    for(unsigned int itheta = 0; itheta < ntheta; ++itheta) {
        cos_theta[itheta]     = cos(itheta*dtheta);
    }

    // calculate often used constants
    const double N0=1/(exp(Ephonon/(kB*Tl))-1); // Bose-Einstein factor

    // Find pre-factors for scattering rates
    const double Upsilon_a = Da*Da*m*N0/(rho*Vs*4*pi*pi*hBar*hBar);
    const double Upsilon_e = Da*Da*m*(N0+1)/(rho*Vs*4*pi*pi*hBar*hBar);

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
        subbands[isb].set_distribution_from_Ef_Te(Ef[isb], Te);
    }

    // Read list of wanted transitions
    arma::uvec i_indices;
    arma::uvec f_indices;

    read_table("rrp.r", i_indices, f_indices);
    const size_t ntx = i_indices.size();
    arma::vec Wabar(ntx);
    arma::vec Webar(ntx);

    // Loop over all desired transitions
    for(unsigned int itx = 0; itx < ntx; ++itx)
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

        arma::vec Kz(nKz);
        arma::vec Gifsqr(nKz);
        ff_table(dKz,isb,fsb,nKz,Kz,Gifsqr);		/* generates formfactor table	*/
        arma::vec Kz_sqr(nKz);

        for(unsigned int iKz = 0; iKz < nKz; ++iKz) {
            Kz_sqr[iKz] = Kz[iKz]*Kz[iKz];
        }

        // Output formfactors if desired
        if(ff_flag) {
            ff_output(Kz, Gifsqr, i, f);
        }

        /* Generate filename for particular mechanism and open file	*/
        std::ostringstream ab_filename; // absorption rate filename
        ab_filename << "ACa" << i << f << ".r";
        FACa=fopen(ab_filename.str().c_str(),"w");

        std::ostringstream em_filename; // emission rate filename
        em_filename << "ACe" << i << f << ".r";
        FACe=fopen(em_filename.str().c_str(),"w");

        // As a zero energy phonon is assumed, no need to 
        // consider emission and absorption processes as in e-LO scattering
        double kimax = 0;
        double Ecutoff = 0.0; // Maximum kinetic energy in initial subband
        const double DeltaE = Ef - Ei;

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
            kimax   = isb.get_k_max(Te);
            Ecutoff = hBar*hBar*kimax*kimax/(2*m);

            if(Ecutoff+Ei < Ef) {
                Ecutoff += Ef;
            }
        }

        kimax = isb.get_k_at_Ek(Ecutoff);

        const double dki=kimax/((float)nki);
        arma::vec Waif(nki); // Absorption scattering rate at this wave-vector [1/s]
        arma::vec Weif(nki); // Emission scattering rate at this wave-vector [1/s]
        arma::vec Wabar_integrand_ki(nki); // Average scattering rate [1/s]
        arma::vec Webar_integrand_ki(nki); // Average scattering rate [1/s]
        const double tmp = 2*m*DeltaE/(hBar*hBar);

        for(unsigned int iki=0;iki<nki;iki++)       /* calculate e-AC rate for all ki	*/
        {
            const double ki=dki*(float)iki+dki/100;	/* second term avoids ki=0 pole	*/

            arma::vec Wif_integrand_dtheta(ntheta);

            /* Integral around angle theta	*/
            for(unsigned int itheta=0;itheta<ntheta;itheta++)
            {
                arma::vec Wif_integrand_dKz(nKz);

                /* Integral over phonon wavevector Kz	*/
                for(unsigned int iKz=0;iKz<nKz;iKz++)
                {
                    const double ki_cos_theta = ki*cos_theta[itheta];
                    const double arg = ki_cos_theta * ki_cos_theta - tmp;	// sqrt argument
                    if(arg>=0)
                    {
                        const double sqrt_arg = sqrt(arg);

                        // solutions for phonon wavevector Kz
                        const double alpha1 =  sqrt_arg - ki_cos_theta;
                        const double alpha2 = -sqrt_arg - ki_cos_theta;

                        const double alpha1_sqr = alpha1*alpha1;
                        const double alpha2_sqr = alpha2*alpha2;

                        /* alpha1 and alpha2 represent solutions for the in-plane polar
                           coordinate Kxy of the carrier momentum---they must be positive, hence
                           use Heaviside unit step function to ignore other contributions	*/
                        Wif_integrand_dKz[iKz] = Gifsqr[iKz]*
                            (alpha1*Theta(alpha1)*sqrt(alpha1_sqr + Kz_sqr[iKz])+
                             alpha2*Theta(alpha2)*sqrt(alpha2_sqr + Kz_sqr[iKz]))/
                            (alpha1-alpha2);
                    }
                } /* end integral over Kz	*/

                Wif_integrand_dtheta[itheta] = integral(Wif_integrand_dKz, dKz);
            } /* end integral over theta	*/

            const double Wif = integral(Wif_integrand_dtheta, dtheta);
            Waif[iki]=2*Upsilon_a*Wif;
            Weif[iki]=2*Upsilon_e*Wif;

            /* Now check for energy conservation!, would be faster with a nasty `if'
               statement just after the beginning of the ki loop!                 */
            const double Eki = isb.get_Ek_at_k(ki);
            const double Ef_em = Eki - DeltaE - Ephonon;
            const double Ef_ab = Eki - DeltaE + Ephonon;
            Weif[iki] *= Theta(Ef_em);
            Waif[iki] *= Theta(Ef_ab);

            // Include final-state blocking factor
            if (b_flag)
            {
                // Final wave-vector
                if(Ef_em >= 0)
                {
                    const double kf_em = sqrt(Ef_em*2*m)/hBar;
                    Weif[iki] *= (1.0 - fsb.get_occupation_at_k(kf_em));
                }

                if(Ef_ab >= 0)
                {
                    const double kf_ab = sqrt(Ef_ab*2*m)/hBar;
                    Waif[iki] *= (1.0 - fsb.get_occupation_at_k(kf_ab));
                }
            }

            Wabar_integrand_ki[iki] = Waif[iki]*ki*isb.get_occupation_at_k(ki);
            Webar_integrand_ki[iki] = Weif[iki]*ki*isb.get_occupation_at_k(ki);

            /* output scattering rate versus carrier energy=subband minima+in-plane
               kinetic energy						*/
            fprintf(FACa,"%20.17le %20.17le\n",(Ei + gsl_pow_2(hBar*ki)/(2*m))/
                    (1e-3*e),Waif[iki]);

            fprintf(FACe,"%20.17le %20.17le\n",(Ei + gsl_pow_2(hBar*ki)/(2*m))/
                    (1e-3*e),Weif[iki]);
        }
        Wabar[itx] = integral(Wabar_integrand_ki, dki)/(pi*isb.get_total_population());
        Webar[itx] = integral(Webar_integrand_ki, dki)/(pi*isb.get_total_population());

        fclose(FACa);	/* close output file for this mechanism	*/
        fclose(FACe);	/* close output file for this mechanism	*/
    } /* end while over states */

    write_table("ACa-if.r", i_indices, f_indices, Wabar);
    write_table("ACe-if.r", i_indices, f_indices, Webar);
    return EXIT_SUCCESS;
} /* end main */

/**
 * \brief calculates the overlap integral squared between the two states
 */
static auto Gsqr(double   Kz,
                 const Subband &isb,
                 const Subband &fsb) -> double
{
 const auto z = isb.z_array();
 const double dz = z[1] - z[0];
 const double nz = z.size();
 const auto psi_i = isb.psi_array();
 const auto psi_f = fsb.psi_array();

 std::complex<double> I(0,1); // Imaginary unit

 // Find form-factor integral
 arma::cx_vec G_integrand_dz(nz);

 for(unsigned int iz=0; iz<nz; ++iz) {
     G_integrand_dz[iz] = exp(Kz*z[iz]*I) * psi_i[iz] * psi_f[iz];
 }

 auto G = integral(G_integrand_dz, dz);

 return norm(G);
}

/**
 * \brief Computes the formfactor at a range of phonon wave-vectors
 */
static void ff_table(const double   dKz,
                     const Subband &isb,
                     const Subband &fsb,
                     unsigned int   nKz,
                     arma::vec     &Kz,
                     arma::vec     &Gifsqr)
{
    for(unsigned int iKz=0;iKz<nKz;iKz++)
    {
        Kz[iKz]     = iKz*dKz;                 // Magnitude of phonon wave vector
        Gifsqr[iKz] = Gsqr(Kz[iKz], isb, fsb); // Squared form-factor
    }
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
