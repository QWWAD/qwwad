/**
 * \file   imp.cpp
 * \brief  Impurity scattering rate solver
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "qwwad/constants.h"
#include "qwwad/subband.h"
#include "qwwad/file-io.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

static void output_ff(double                      W, // Arbitrary well width to generate q
                      const std::vector<Subband> &subbands,
                      unsigned int                i,
                      unsigned int                f,
                      const arma::vec            &d);

auto FF_table(double           epsilon,
              const Subband   &isb,
              const Subband   &fsb,
              const arma::vec &d,
              size_t           nq,
              bool             S_flag,
              double           E_cutoff) -> gsl_spline *;

auto configure_options(int argc, char** argv) -> Options
{
    Options opt;

    std::string doc("Find the impurity scattering rate.");

    opt.add_option<bool>  ("outputff,a",           "Output form-factors to file.");
    opt.add_option<bool>  ("noscreening,S",        "Disable screening of the Coulomb interaction.");
    opt.add_option<bool>  ("noblocking,b",         "Disable final-state blocking.");
    opt.add_option<double>("epsilon,e",     13.18, "Low-frequency dielectric constant");
    opt.add_option<double>("mass,m",        0.067, "Band-edge effective mass (relative to free electron)");
    opt.add_option<char>  ("particle,p",      'e', "ID of particle to be used: 'e', 'h' or 'l', for "
                                                   "electrons, heavy holes or light holes respectively.");
    opt.add_option<double>("temperature,T",   300, "Temperature of carrier distribution.");
    opt.add_option<double>("width,w",         250, "Width of quantum well [angstrom]. (Solely for output).");
    opt.add_option<double>("Ecutoff",              "Cut-off energy for carrier distribution [meV]. If not specified, then 5kT above band-edge.");
    opt.add_option<size_t>("nki",             101, "Number of initial wave-vector samples.");
    opt.add_option<size_t>("nq",              101, "Number of strips in scattering vector integration");
    opt.add_option<size_t>("ntheta",          101, "Number of strips in theta angle integration");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

auto main(int argc,char *argv[]) -> int
{
    const auto opt = configure_options(argc, argv);

    const auto epsilon =  opt.get_option<double>("epsilon")*eps0; // Low frequency dielectric constant [F/m]
    const auto ff_flag =  opt.get_option<bool>  ("outputff");     // True if formfactors are wanted
    const auto m       =  opt.get_option<double>("mass")*me;      // Band-edge effective mass [kg]
    const auto p       =  opt.get_option<char>  ("particle");     // Particle ID
    const auto T       =  opt.get_option<double>("temperature");  // Temperature [K]
    const auto W       =  opt.get_option<double>("width")*1e-10;  // a well width, same as Smet [angstrom]
    const auto S_flag  = !opt.get_option<bool>  ("noscreening");  // Include screening by default
    const auto b_flag  = !opt.get_option<bool>  ("noblocking");   // Include final-state blocking by default
    const auto nki     =  opt.get_option<size_t>("nki");          // number of ki calculations
    const auto ntheta  =  opt.get_option<size_t>("ntheta");       // number of strips in theta integration
    const auto nq      =  opt.get_option<size_t>("nq");           // number of q_perp values for lookup table

    /* calculate step lengths	*/
    const double dtheta=2*pi/((float)ntheta - 1); // step length for theta integration

    // Can save a bit of time by calculating cosines in advance
    arma::vec cos_theta(ntheta);

    for(unsigned int itheta = 0; itheta < ntheta; ++itheta) {
        cos_theta[itheta] = cos(itheta*dtheta);
    }

    std::ostringstream E_filename; // Energy filename string
    E_filename << "E" << p << ".r";
    std::ostringstream wf_prefix;  // Wavefunction filename prefix
    wf_prefix << "wf_" << p;

    // Read data for all subbands from file
    auto subbands = Subband::read_from_file(E_filename.str(),
                                            wf_prefix.str(),
                                            ".r",
                                            m);

    // Read and set carrier distributions within each subband
    arma::vec  Ef;      // Fermi energies [J]
    arma::uvec indices; // Subband indices (garbage)
    read_table("Ef.r", indices, Ef);
    Ef *= e/1000.0; // Rescale to J

    // Read doping profile
    arma::vec z_d; // Spatial location
    arma::vec d;   // Volume doping [m^{-3}]
    read_table("d.r", z_d, d);

    for(unsigned int isb = 0; isb < subbands.size(); ++isb) {
        subbands[isb].set_distribution_from_Ef_Te(Ef[isb], T);
    }

    // Read list of wanted transitions
    arma::uvec i_indices;
    arma::uvec f_indices;

    read_table("rrp.r", i_indices, f_indices);

    auto *Favg=fopen("imp-avg.dat","w");	/* open file for output of weighted means */

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

        // Output form-factors if desired
        if(ff_flag) {
            output_ff(W,subbands,i,f,d);
        }

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

        gsl_spline *FF = FF_table(epsilon, isb, fsb, d,nq,S_flag,Ecutoff); // Form factor table

        gsl_interp_accel *acc = gsl_interp_accel_alloc(); // Creates accelerator for interpolation of FF

        /* calculate maximum value of ki & kj and hence kj step length	*/
        const double dki=(kimax-kimin)/((float)nki - 1); // step length for loop over ki

        arma::vec Wbar_integrand_ki(nki); // initialise integral for average scattering rate
        arma::vec Wif(nki);               // Scattering rate for a given initial wave vector
        arma::vec Ei_t(nki);              // Total energy of initial state (for output file) [meV]

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
            const double two_kif = 2*ki*kf;
            const double ki_sqr_plus_kf_sqr = ki_sqr + kf_sqr;

            // Now perform innermost integral (over theta)
            arma::vec Wif_integrand_theta(ntheta);

            for(unsigned int itheta=0;itheta<ntheta;itheta++)
            {
                // Calculate scattering vector
                // Note that most of the terms here are computed before the loop, so we
                // only need to look up the cos(theta)
                const double q_sqr = ki_sqr_plus_kf_sqr + two_kif * cos_theta[itheta];
                const double q = sqrt(q_sqr);
                assert(!std::isnan(q));

                // Find the form-factor at this wave-vector by looking it up in the
                // spline we created earlier
                Wif_integrand_theta[itheta] = gsl_spline_eval(FF, q, acc);
            } /* end theta */

            Wif[iki] = integral(Wif_integrand_theta, dtheta);

            // Multiply by pre-factor
            Wif[iki] *= m*e*e*e*e / (4*pi*hBar*hBar*hBar*epsilon*epsilon);

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
        std::ostringstream filename;	/* character string for output filename		*/
        filename << "imp" << i << f << ".r";
        write_table(filename.str(), Ei_t, Wif);

        const double Wbar = integral(Wbar_integrand_ki, dki)/(pi*isb.get_total_population());

        fprintf(Favg,"%i %i %20.17le\n", i,f,Wbar);

        gsl_spline_free(FF);
        gsl_interp_accel_free(acc);
} /* end while over states */

fclose(Favg);	/* close weighted mean output file	*/

return EXIT_SUCCESS;
} /* end main */

/** Tabulate the matrix element defined as 
 *    C_if⁺(q,z') = ∫_{z'}^∞ dz ψ_i(z) ψ_f(z)/exp(qz)]
 *  for a given wavevector, with respect to position
 */
auto find_Cif_p(const arma::cx_vec &psi_if, 
                const arma::vec    &exp_qz,
                const arma::vec    &z) -> arma::cx_vec
{
    const size_t nz = z.size();
    arma::cx_vec Cif_p(nz);
    const double dz=z[1]-z[0];

    Cif_p[nz-1] = psi_if[nz-1] / exp_qz[nz-1] * dz;

    for(int iz = nz-2; iz >=0; iz--) {
        Cif_p[iz] = Cif_p[iz+1] + psi_if[iz] / exp_qz[iz] * dz;
    }

    return Cif_p;
}

/** 
 * \brief Tabulate scattering matrix element component.
 *
 * \details defined as:
 *    C_if⁻(q,z') = ∫_{-∞}^{z'} dz ψ_i(z) ψ_f(z) exp(qz)
 *  for a given wavevector, with respect to position
 *
 * Note that the upper limit has to be the point just BEFORE each z'
 * value so that we don't double count
 */
auto find_Cif_m(const arma::cx_vec &psi_if, 
                const arma::vec    &exp_qz,
                const arma::vec    &z) -> arma::cx_vec
{
    const size_t nz = z.size();
    arma::cx_vec Cif_m(nz);
    const double dz = z[1]-z[0];

    // Seed the first value as zero
    Cif_m[0] = 0;

    // Now, perform a block integration by summing on top of the previous
    // value in the array
    for(unsigned int iz = 1; iz < nz; iz++) {
        Cif_m[iz] = Cif_m[iz-1] + psi_if[iz-1] * exp_qz[iz-1] * dz;
    }

    return Cif_m;
}

/** 
 * \brief Create an array of exp(qz) with respect to position
 *
 * \param q[in]       Scattering vector [1/m]
 * \param exp_qz[out] Output array (should be initialised before calling)
 * \param z[in]       Spatial positions [m]
 *
 * \todo  This is also useful for e-e scattering. Push into library
 */
auto find_exp_qz(const double q, const arma::vec &z) -> arma::vec
{
    //const double Lp = z.max() - z.min();

    // Use the midpoint of the z array as the origin, so as to minimise the
    // magnitude of the exponential terms
    return exp(q * (z - z[0]));
}

/** Find the matrix element Iif at a given dopant location z'.
 *
 * The matrix element is defined as
 *  I_if(q,z') = ∫dz ψ_i(z) ψ_f(z) exp(-q|z-z'|),
 * where z is the carrier location.  The numerical solution
 * can however be speeded up by replacing the modulus function
 * with the sum of two integrals.  We can say that
 *
 *  I_if(q,z') = C_if⁻(q,z')/exp(qz') + C_if⁺(q,z') exp(qz')',
 *
 * Therefore, we have separated the z' dependence from the
 * z dependence of the matrix element.
 */
auto Iif(unsigned int iz0,
         const arma::cx_vec &Cif_p,
         const arma::cx_vec &Cif_m, 
         const arma::vec    &exp_qz) -> std::complex<double>
{
    return Cif_m[iz0]/exp_qz[iz0] + Cif_p[iz0]*exp_qz[iz0];
}

/* This function calculates the overlap integral
 */
auto J(const double     q_perp,
         const Subband   &isb,
         const Subband   &fsb,
         const arma::vec &d) -> double
{
 const auto z = isb.z_array();
 const auto nz = z.size();
 const auto dz = z[1] - z[0];

 // Convenience labels for wave-functions in each subband
 const auto psi_i = isb.psi_array();
 const auto psi_f = fsb.psi_array();

 // Products of wavefunctions can be computed in advance
 const auto psi_if = psi_i % psi_f;

 const auto expTerm   = find_exp_qz(q_perp, z);
 const auto Cif_plus  = find_Cif_p(psi_if, expTerm, z);
 const auto Cif_minus = find_Cif_m(psi_if, expTerm, z);

 arma::vec Jif_integrand(nz);

 // Integral of i(=0) and f(=2) over z
 for(unsigned int iz=0;iz<nz;iz++) {
     const auto _Iif = Iif(iz, Cif_plus, Cif_minus, expTerm);
     const auto abs_Iif = abs(_Iif);
     Jif_integrand[iz] = abs_Iif * abs_Iif * d[iz];
 }

 const double Jif = integral(Jif_integrand, dz);

 return Jif;
}

/**
 *  \brief Compute the form factor Jif/q^2
 */
auto FF_table(const double     epsilon,
                      const Subband   &isb,
                      const Subband   &fsb,
                      const arma::vec &d,
                      const size_t     nq,
                      const bool       S_flag,
                      const double     E_cutoff) -> gsl_spline *
{
    const double kimax = isb.get_k_at_Ek(E_cutoff*1.1); // Max value of ki [1/m]
    const double Ei = isb.get_E_min();
    const double Ef = fsb.get_E_min();
    const double m  = isb.get_effective_mass();
    const double kfmax = sqrt(kimax*kimax + 2*m*(Ei - Ef)/(hBar*hBar));

    // maximum in-plane wave vector
    const double q_max = sqrt(kimax*kimax + kfmax*kfmax + 2*kimax*kfmax);

    const double dq=q_max/((float)(nq-1));	// interval in q_perp

    arma::vec q(nq);
    arma::vec FF(nq);

    for(unsigned int iq=0;iq<nq;iq++)
    {
        q[iq] = iq*dq;

        // Scattering matrix element
        const double Jif = J(q[iq], isb, fsb, d);

        // Thomas--Fermi screening wave-vector
        double q_TF = 0.0;

        // Allow screening to be turned off
        if(S_flag) {
            q_TF = m*e*e/(2*pi*epsilon*hBar*hBar);
        }

        // Screening permittivity * wave vector
        // Note that the pole at q_perp=0 is avoided as long as screening is included
        FF[iq] = Jif / (q[iq]*q[iq] + q_TF*q_TF + 2*q[iq]*q_TF);
    }

    // Fix singularity by "clipping" the top off it:
    if(!S_flag) {
        FF[0] = FF[1];
    }

    // Pack the table of FF vs q into a cubic spline
    gsl_spline *q_FF = gsl_spline_alloc(gsl_interp_cspline, nq);
    gsl_spline_init(q_FF, &(q[0]), &(FF[0]), nq);

    return q_FF;
}

/* This function outputs the formfactors into files	*/
static void output_ff(const double        W, // Arbitrary well width to generate q
                      const std::vector<Subband> &subbands,
                      const unsigned int  i,
                      const unsigned int  f,
                      const arma::vec    &d)
{
    std::ostringstream filename;	/* output filename				*/

 /* First generate filename and then open file	*/
 filename << "J" << i << f << ".r";
 auto *FA=fopen(filename.str().c_str(),"w"); // output file for form factors versus q_perp
 if(FA == nullptr) {
     std::cerr << "Error: Cannot open input file '" << filename.str() << "'." << std::endl;
     exit(EXIT_FAILURE);
 }

 // Convenience labels for each subband (NB., these are indexed from 0)
 const Subband isb = subbands[i-1];
 const Subband fsb = subbands[f-1];

 for(unsigned int iq=0;iq<100;iq++)
 {
  const double q_perp=6*iq/(100*W); // In-plane scattering vector
  const double Jif = J(q_perp, isb, fsb, d);
  fprintf(FA,"%le %le\n",q_perp*W,gsl_pow_2(Jif));
 }

 fclose(FA);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
