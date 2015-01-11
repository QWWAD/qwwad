/**
 * \file   srcc.cpp
 * \brief  Carrier-carrier scattering rate solver
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
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
#include "qclsim-constants.h"
#include "qclsim-subband.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

static void output_ff(const double       W,
                      const std::vector<Subband> &subbands,
                      const unsigned int i,
                      const unsigned int j,
                      const unsigned int f,
                      const unsigned int g);

gsl_spline * FF_table(const double                 Deltak0sqr,
                      const double                 epsilon,
                      const Subband               &isb,
                      const Subband               &jsb,
                      const Subband               &fsb,
                      const Subband               &gsb,
                      const double                 T,
                      const size_t                 nq,
                      const bool                   S_flag,
                      const double                 E_cutoff = -1);

double PI(const Subband &isb,
          const size_t   q_perp,
          const double   T);

Options configure_options(int argc, char* argv[])
{
    Options opt;

    std::string doc("Find the carrier-carrier scattering rate.");

    opt.add_switch        ("outputff,a",           "Output form-factors to file.");
    opt.add_switch        ("noscreening,S",        "Disable screening of the Coulomb interaction.");
    opt.add_numeric_option("epsilon,e",     13.18, "Low-frequency dielectric constant");
    opt.add_numeric_option("mass,m",        0.067, "Band-edge effective mass (relative to free electron)");
    opt.add_char_option   ("particle,p",      'e', "ID of particle to be used: 'e', 'h' or 'l', for "
                                                   "electrons, heavy holes or light holes respectively.");
    opt.add_numeric_option("temperature,T",   300, "Temperature of carrier distribution.");
    opt.add_numeric_option("width,w",         250, "Width of quantum well [angstrom]. (Solely for output).");
    opt.add_numeric_option("Ecutoff",              "Cut-off energy for carrier distribution [meV]. If not specified, then 5kT above band-edge.");
    opt.add_size_option   ("nki",             101, "Number of initial wave-vector samples for first carrier");
    opt.add_size_option   ("nkj",             101, "Number of initial wave-vector samples for second carrier");
    opt.add_size_option   ("nq",              101, "Number of strips in scattering vector integration");
    opt.add_size_option   ("ntheta",          101, "Number of strips in alpha angle integration");
    opt.add_size_option   ("nalpha",          101, "Number of strips in theta angle integration");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const Options opt = configure_options(argc, argv);

    const double epsilon = opt.get_numeric_option("epsilon")*eps0; // Low frequency dielectric constant [F/m]
    const bool   ff_flag = opt.get_switch("outputff");             // True if formfactors are wanted
    const double m       = opt.get_numeric_option("mass")*me;      // Band-edge effective mass [kg]
    const char   p       = opt.get_char_option("particle");  	   // Particle ID
    const double T       = opt.get_numeric_option("temperature");  // Temperature [K]
    const double W       = opt.get_numeric_option("width")*1e-10;  // a well width, same as Smet [angstrom]
    const bool   S_flag  = !opt.get_switch("noscreening");	   // Include screening by default
    const size_t nki     = opt.get_size_option("nki");             // number of ki calculations
    const size_t nkj     = opt.get_size_option("nkj");             // number of strips in |kj| integration
    const size_t nalpha  = opt.get_size_option("nalpha");          // number of strips in alpha integration
    const size_t ntheta  = opt.get_size_option("ntheta");          // number of strips in theta integration
    const size_t nq      = opt.get_size_option("nq");              // number of q_perp values for lookup table

    /* calculate step lengths	*/
    const double dalpha=2*pi/((float)nalpha - 1); // step length for alpha integration
    const double dtheta=2*pi/((float)ntheta - 1); // step length for theta integration

    // Can save a bit of time by calculating cosines in advance
    std::valarray<double> cos_theta(ntheta);

    for(unsigned int itheta = 0; itheta < ntheta; ++itheta)
        cos_theta[itheta] = cos(itheta*dtheta);

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
    std::valarray<double>       Ef;      // Fermi energies [J]
    std::valarray<double>       N;       // Subband populations [m^{-2}]
    std::valarray<unsigned int> indices; // Subband indices (garbage)
    read_table_xy("Ef.r", indices, Ef);
    Ef *= e/1000.0; // Rescale to J
    read_table_xy("N.r", indices, N);	// read populations
    N *= 1e+10*1e+4; // convert from units of 10^10cm^-2->m^-2

    for(unsigned int isb = 0; isb < subbands.size(); ++isb)
        subbands[isb].set_distribution(Ef[isb], N[isb]);

    // Read list of wanted transitions
    std::valarray<unsigned int> i_indices;
    std::valarray<unsigned int> j_indices;
    std::valarray<unsigned int> f_indices;
    std::valarray<unsigned int> g_indices;

    read_table_xyzu("rr.r", i_indices, j_indices, f_indices, g_indices);

    FILE *FccABCD=fopen("ccABCD.r","w");	/* open file for output of weighted means */

    // Loop over all desired transitions
    for(unsigned int itx = 0; itx < i_indices.size(); ++itx)
    {
        // State indices for this transition (NB., these are indexed from 1)
        unsigned int i = i_indices[itx];
        unsigned int j = j_indices[itx];
        unsigned int f = f_indices[itx];
        unsigned int g = g_indices[itx];

        // Convenience labels for each subband (NB., these are indexed from 0)
        const Subband isb = subbands[i-1];
        const Subband jsb = subbands[j-1];
        const Subband fsb = subbands[f-1];
        const Subband gsb = subbands[g-1];

        // Subband minima
        const double Ei = isb.get_E();
        const double Ej = jsb.get_E();
        const double Ef = fsb.get_E();
        const double Eg = gsb.get_E();

        // Output form-factors if desired
        if(ff_flag)
            output_ff(W,subbands,i,j,f,g);

        // Calculate Delta k0^2 [QWWAD3, Eq. 10.228]
        //   twice the change in KE, see Smet (55)
        double Deltak0sqr = 0;
        if(i+j != f+g)
            Deltak0sqr=4*m*(Ei + Ej - Ef - Eg)/(hBar*hBar);	

        gsl_spline *FF = 0;
        double kimax = 0;
        double kjmax = 0;

        if(opt.vm.count("Ecutoff"))
        {
            const double Ecutoff = opt.get_numeric_option("Ecutoff")*e/1000;
            FF = FF_table(Deltak0sqr, epsilon, isb, jsb, fsb, gsb, T,nq,S_flag,Ecutoff); // Form factor table
            kimax = isb.k(Ecutoff);
            kjmax = jsb.k(Ecutoff);
        }
        else
        {
            FF = FF_table(Deltak0sqr, epsilon, isb, jsb, fsb, gsb, T,nq,S_flag); // Form factor table
            kimax=isb.get_k_max(T);
            kjmax=jsb.get_k_max(T);
        }

        gsl_interp_accel *acc = gsl_interp_accel_alloc(); // Creates accelerator for interpolation of FF

        /* calculate maximum value of ki & kj and hence kj step length	*/
        const double dki=kimax/((float)nki - 1); // step length for loop over ki
        const double dkj=kjmax/((float)nkj - 1); // step length for kj integration

        std::valarray<double> Wbar_integrand_ki(nki); // initialise integral for average scattering rate
        std::valarray<double> Wijfg(nki);             // Scattering rate for a given initial wave vector
        std::valarray<double> Ei_t(nki);              // Total energy of initial state (for output file) [meV]

        // calculate c-c rate for all ki
        for(unsigned int iki=0;iki<nki;iki++)
        {
            const double ki=dki*(float)iki; // carrier momentum

            // integrate over |kj|
            std::valarray<double> Wijfg_integrand_kj(nkj);

            for(unsigned int ikj=0;ikj<nkj;ikj++)
            {
                const double kj=dkj*(float)ikj; // carrier momentum

                // Find Fermi-Dirac occupation at kj
                const double P=jsb.f_FD_k(kj,T);

                // Integral over alpha
                std::valarray<double> Wijfg_integrand_alpha(nalpha);

                for(unsigned int ialpha=0;ialpha<nalpha;ialpha++)
                {
                    const double alpha=dalpha*(float)ialpha; // angle between ki and kj

                    // Compute (vector)kj-(vector)(ki) [QWWAD3, 10.221]
                    const double kij_sqr = ki*ki+kj*kj-2*ki*kj*cos(alpha);
                    const double kij = sqrt(kij_sqr);

                    // Can also pre-calculate a few of the terms needed inside the following loop
                    // to save time
                    const double kij_sqr_plus_Deltak0sqr = kij_sqr + Deltak0sqr;
                    const double two_kij_sqr_plus_Deltak0sqr = kij_sqr + kij_sqr_plus_Deltak0sqr;
                    const double two_kij_times_sqrt__kij_sqr_plus_Deltak0sqr__
                                     = 2 * kij * sqrt(kij_sqr_plus_Deltak0sqr);

                    // Now perform innermost integral (over theta)
                    std::valarray<double> Wijfg_integrand_theta(ntheta);

                    for(unsigned int itheta=0;itheta<ntheta;itheta++)
                    {
                        /* calculate argument of sqrt function=4*q_perp*q_perp,
                         * see [QWWAD3, 10.231],
                           to check for imaginary q_perp, if argument is positive, q_perp is
                           real and hence calculate scattering rate, otherwise ignore and move
                           onto next q_perp */
                        // Note that most of the terms here are computed before the loop, so we
                        // only need to look up the cos(theta)
                        const double q_perpsqr4 = two_kij_sqr_plus_Deltak0sqr - two_kij_times_sqrt__kij_sqr_plus_Deltak0sqr__ * cos_theta[itheta];

                        if(q_perpsqr4>=0) 
                        {
                            const double q_perp=sqrt(q_perpsqr4)/2; // in-plane momentum, |ki-kf|

                            // Find the form-factor at this wave-vector by looking it up in the
                            // spline we created earlier
                            Wijfg_integrand_theta[itheta] = gsl_spline_eval(FF, q_perp, acc);
                        }
                    } /* end theta */

                    Wijfg_integrand_alpha[ialpha] = integral(Wijfg_integrand_theta, dtheta);
                } /* end alpha */

                Wijfg_integrand_kj[ikj] = integral(Wijfg_integrand_alpha, dalpha) * P * kj;
            } /* end kj   */

            Wijfg[iki] = integral(Wijfg_integrand_kj,dkj);

            // Multiply by pre-factor [QWWAD3, 10.233]
            Wijfg[iki] *= m*e*e*e*e / (4*pi*hBar*hBar*hBar*(4*4*pi*pi*epsilon*epsilon));
            Ei_t[iki] = isb.E_total(ki) * 1000/e;

            /* calculate Fermi-Dirac weighted mean of scattering rates over the 
               initial carrier states, note that the integral step length 
               dE=2*sqr(hBar)*ki*dki/(2m)					*/
            Wbar_integrand_ki[iki] = Wijfg[iki]*ki*isb.f_FD_k(ki, T);
        } /* end ki	*/

        /* output scattering rate versus carrier energy=subband minima+in-plane
           kinetic energy						*/
        char	filename[9];	/* character string for output filename		*/
        sprintf(filename,"cc%i%i%i%i.r",i,j,f,g);
        write_table_xy(filename, Ei_t, Wijfg);

        const double Wbar = integral(Wbar_integrand_ki, dki)/(pi*isb.get_pop());

        fprintf(FccABCD,"%i %i %i %i %20.17le\n", i,j,f,g,Wbar);

        gsl_spline_free(FF);
        gsl_interp_accel_free(acc);
} /* end while over states */

fclose(FccABCD);	/* close weighted mean output file	*/

return EXIT_SUCCESS;
} /* end main */

/** Tabulate the matrix element defined as 
 *    C_if⁺(q,z') = ∫_{z'}^∞ dz ψ_i(z) ψ_f(z)/exp(qz)]
 *  for a given wavevector, with respect to position
 */
std::valarray<double> find_Cif_p(const std::valarray<double>& psi_if, 
                                 const std::valarray<double>& exp_qz,
                                 const std::valarray<double>& z)
{
    const size_t nz = z.size();
    std::valarray<double> Cif_p(nz);
    const double dz=z[1]-z[0];

    Cif_p[nz-1] = psi_if[nz-1] / exp_qz[nz-1] * dz;

    for(int iz = nz-2; iz >=0; iz--)
        Cif_p[iz] = Cif_p[iz+1] + psi_if[iz] / exp_qz[iz] * dz;

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
std::valarray<double> find_Cif_m(const std::valarray<double>& psi_if, 
                                 const std::valarray<double>& exp_qz,
                                 const std::valarray<double>& z)
{
    const size_t nz = z.size();
    std::valarray<double> Cif_m(nz);
    const double dz = z[1]-z[0];

    // Seed the first value as zero
    Cif_m[0] = 0;

    // Now, perform a block integration by summing on top of the previous
    // value in the array
    for(unsigned int iz = 1; iz < nz; iz++)
        Cif_m[iz] = Cif_m[iz-1] + psi_if[iz-1] * exp_qz[iz-1] * dz;

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
std::valarray<double> find_exp_qz(const double q, const std::valarray<double>& z)
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
double Iif(const unsigned int iz0,
           const std::valarray<double>& Cif_p,
           const std::valarray<double>& Cif_m, 
           const std::valarray<double>& exp_qz)
{
    return Cif_m[iz0]/exp_qz[iz0] + Cif_p[iz0]*exp_qz[iz0];
}

/* This function calculates the overlap integral over all four carrier
   states		*/
double A(const double   q_perp,
         const Subband &isb,
         const Subband &jsb,
         const Subband &fsb,
         const Subband &gsb)
{
 const std::valarray<double> z = isb.z_array();
 const size_t nz = z.size();
 const double dz = z[1] - z[0];

 // Convenience labels for wave-functions in each subband
 const std::valarray<double> psi_i = isb.psi_array();
 const std::valarray<double> psi_j = jsb.psi_array();
 const std::valarray<double> psi_f = fsb.psi_array();
 const std::valarray<double> psi_g = gsb.psi_array();

 // Products of wavefunctions can be computed in advance
 const std::valarray<double> psi_if = psi_i * psi_f;
 const std::valarray<double> psi_jg = psi_j * psi_g;

 const std::valarray<double> expTerm = find_exp_qz(q_perp, z);
 const std::valarray<double> Cjg_plus  = find_Cif_p(psi_jg, expTerm, z);
 const std::valarray<double> Cjg_minus = find_Cif_m(psi_jg, expTerm, z);

 std::valarray<double> Aijfg_integrand(nz);

 // Integral of i(=0) and f(=2) over z
 for(unsigned int iz=0;iz<nz;iz++)
 {
     const double Ijg = Iif(iz, Cjg_plus, Cjg_minus, expTerm);
     Aijfg_integrand[iz] = psi_if[iz] * Ijg;
 }

 const double Aijfg = integral(Aijfg_integrand, dz);

 return Aijfg;
}

/**
 * \brief returns the screening factor, referred to by Smet as e_sc
 */
double PI(const Subband &isb,
          const double   q_perp,
          const double   T)
{
    const double m    = isb.get_md_0();    // Effective mass at band-edge [kg]

    // Now perform the integration, equation 44 of Smet [QWWAD3, 10.238]
    double mu=isb.get_E();
    const double dmu=1e-3*e;
    double integral=0;

    double dI; // Intervals in I

    do
    {
        const double ki = isb.k(mu);

        // Find low-temperature polarizability *at this wave-vector*
        //
        // Equation 43 of Smet, QWWAD3, 10.236
        double P0 = m/(pi*hBar*hBar);

        if(q_perp>2*ki)
            P0 -= m/(pi*hBar*hBar)*sqrt(1-gsl_pow_2(2*ki/q_perp));

        dI = P0/(4*kB*T*gsl_pow_2(cosh((isb.get_Ef() - mu)/(2*kB*T))));
        integral += dI*dmu;
        mu+=dmu;
    }while(dI>integral/100); // continue until integral converges

    return integral;
}

/**
 *  \brief Compute the form factor [Aijfg/(esc q)]^2
 */
gsl_spline * FF_table(const double                 Deltak0sqr,
                      const double                 epsilon,
                      const Subband               &isb,
                      const Subband               &jsb,
                      const Subband               &fsb,
                      const Subband               &gsb,
                      const double                 T,
                      const size_t                 nq,
                      const bool                   S_flag,
                      const double                 E_cutoff)
{
    // Find maximum wave-vectors for calculation if not specified
    double kimax = 0.0; // Max value of ki [1/m]
    double kjmax = 0.0; // Max value of kj [1/m]

    if(E_cutoff > 0)
    {
        kimax = isb.k(E_cutoff*1.1);
        kjmax = jsb.k(E_cutoff*1.1);
    }
    else
    {
        kimax = isb.get_k_max(T*1.1);
        kjmax = jsb.get_k_max(T*1.1);
    }

    // maximum in-plane wave vector
    const double q_perp_max=sqrt(2*gsl_pow_2(kimax+kjmax)+Deltak0sqr+2*(kimax+kjmax)*
                 sqrt(gsl_pow_2(kimax+kjmax)+Deltak0sqr))/2;

    const double dq=q_perp_max/((float)(nq-1));	// interval in q_perp

    std::valarray<double> q_perp(nq);
    std::valarray<double> FF(nq);

    for(unsigned int iq=0;iq<nq;iq++)
    {
        q_perp[iq] = iq*dq;

        // Scattering matrix element (all 4 states)
        const double _Aijfg = A(q_perp[iq], isb, jsb, fsb, gsb);

        double _PI    = 0.0; // Polarizability
        double _Aiiii = 0.0; // Matrix element for lowest subband

        // Allow screening to be turned off
        if(S_flag)
        {
            _PI    = PI(isb, q_perp[iq], T);
            _Aiiii = A(q_perp[iq], isb, isb, isb, isb);
        }

        // Screening permittivity * wave vector
        // Note that the pole at q_perp=0 is avoided as long as screening is included
        const double esc_q = q_perp[iq] + 2*pi*e*e/(4*pi*epsilon) * _PI * _Aiiii;
        FF[iq] = _Aijfg*_Aijfg / (esc_q * esc_q);
    }

    // Fix singularity by "clipping" the top off it:
    if(!S_flag)
        FF[0] = FF[1];

    // Pack the table of FF vs q into a cubic spline
    gsl_spline *q_FF = gsl_spline_alloc(gsl_interp_cspline, nq);
    gsl_spline_init(q_FF, &(q_perp[0]), &(FF[0]), nq);

    return q_FF;
}

/* This function outputs the formfactors into files	*/
static void output_ff(const double        W, // Arbitrary well width to generate q
                      const std::vector<Subband> &subbands,
                      const unsigned int  i,
                      const unsigned int  j,
                      const unsigned int  f,
                      const unsigned int  g)
{
 char	filename[9];	/* output filename				*/
 FILE	*FA;		/* output file for form factors versus q_perp	*/

 /* First generate filename and then open file	*/

 sprintf(filename,"A%i%i%i%i.r", i, j, f, g);	
 if((FA=fopen(filename,"w"))==0)
 {
     std::cerr << "Error: Cannot open input file '" << filename << "'." << std::endl;
     exit(EXIT_FAILURE);
 }

 // Convenience labels for each subband (NB., these are indexed from 0)
 const Subband isb = subbands[i-1];
 const Subband jsb = subbands[j-1];
 const Subband fsb = subbands[f-1];
 const Subband gsb = subbands[g-1];

 for(unsigned int iq=0;iq<100;iq++)
 {
  const double q_perp=6*iq/(100*W); // In-plane scattering vector
  const double Aijfg=A(q_perp,isb,jsb,fsb,gsb);
  fprintf(FA,"%le %le\n",q_perp*W,gsl_pow_2(Aijfg));
 }

 fclose(FA);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
