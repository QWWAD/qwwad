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
#include <cmath>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qclsim-maths.h"
#include "qwwad-options.h"
#include "qwwad-schroedinger-donor-2D.h"

using namespace Leeds;
using namespace constants;

static bool repeat_lambda(const double                  lambda,
                          double                       &lambda_0,
                          const double                  E,
                          double                       &E0);

struct lambda_search_params
{
    const std::valarray<double>  &z;
    const double                  epsilon;
    const double                  mstar;
    const double                  r_d;
    const std::valarray<double>  &Vp;
    const size_t                  N_w;
    const double                  delta_E;
};

/**
 * \brief Find the energy of a carrier using a given Bohr radius
 */
static double find_E_at_lambda(double  lambda,
                               void   *params)
{
    const lambda_search_params *p = reinterpret_cast<lambda_search_params *>(params);
    SchroedingerSolverDonor2D se(p->mstar, p->Vp, p->z, p->epsilon, p->r_d, lambda, p->delta_E, 1);
    std::vector<State> solutions = se.get_solutions();

    return solutions[0].get_E();
}

/**
 * Find the minimum carrier energy, and corresponding Bohr radius using a fast search algorithm
 */
static void find_E_min_fast(double                      &E0,
                            double                      &lambda0,
                            const double                 lambda_start,
                            const double                 lambda_stop,
                            const std::valarray<double> &z,
                            const double                 epsilon,
                            const double                 mstar,
                            const double                 r_d,
                            const std::valarray<double> &V,
                            const size_t                 N_w,
                            const double                 delta_E)
{
    double _lambda_start = lambda_start;

    // Set up the numerical solver using GSL
    lambda_search_params params = {z, epsilon, mstar, r_d, V, N_w, delta_E};
    gsl_function f;
    f.function = &find_E_at_lambda;
    f.params   = &params;

    // First perform a very coarse search for a suitable estimate of a starting point
    const double Elo = GSL_FN_EVAL(&f, _lambda_start);
    const double Ehi = GSL_FN_EVAL(&f, lambda_stop);

    E0 = Elo + Ehi; // Set initial estimate as being higher than Elo and Ehi
    lambda0 = lambda_start;

    const double dlambda = (lambda_stop - lambda_start)/4; // Separation between endpoints [m]

    // Search for a suitable lambda value until we find which quadrant the mimimum lies in
    do
    {
        lambda0 += dlambda; // Increment the Bohr radius
        if(lambda0 >= lambda_stop)
        {
            std::cerr << "Can't find a minimum in this range" << std::endl;
            exit(EXIT_FAILURE);
        }
        E0 = GSL_FN_EVAL(&f, lambda0);
    }
    while((E0 > Elo) || (E0 > Ehi));
    _lambda_start = lambda0 - dlambda;
    
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(s, &f, lambda0, lambda_start, lambda_stop);

    size_t max_iter = 100; // Maximum number of iterations before giving up
    int status = 0;        // Error flag for GSL
    unsigned int iter=0;   // The number of iterations attempted so far

    // Variational calculation (search over lambda)
    do
    {
        ++iter;
        status  = gsl_min_fminimizer_iterate(s);
        lambda0 = gsl_min_fminimizer_x_minimum(s);
        E0      = gsl_min_fminimizer_f_minimum(s);
        const double lambda_lo = gsl_min_fminimizer_x_lower(s);
        const double lambda_hi = gsl_min_fminimizer_x_upper(s);
        status  = gsl_min_test_interval(lambda_lo, lambda_hi, 0.1e-10, 0.0);
        printf("r_d %le lambda %le energy %le meV\n",r_d,lambda0,E0/(1e-3*e));
    }while((status == GSL_CONTINUE) && (iter < max_iter));

    gsl_min_fminimizer_free(s);
}

/**
 * Find the minimum carrier energy using a linear search
 */
static void find_E_min_linear(double &E0,
                              double &lambda0,
                              const double lambda_start,
                              const double lambda_step,
                              const double lambda_stop,
                              const std::valarray<double> &z,
                              const double epsilon,
                              const double mstar,
                              const double r_d,
                              const std::valarray<double> &V,
                              const size_t N_w,
                              const double delta_E)
{
    double lambda=lambda_start; // Initial Bohr radius value [m]
    E0 = e;                     // Minimum energy of single donor 1eV

    lambda_search_params params = {z, epsilon, mstar, r_d, V, N_w, delta_E};
    bool solution_not_found = true; // True if we haven't found the solution yet

    // Variational calculation (search over lambda)
    do
    {
        const double E = find_E_at_lambda(lambda, &params);
        printf("r_d %le lambda %le energy %le meV\n",r_d,lambda,E/(1e-3*e));
        solution_not_found = repeat_lambda(lambda,lambda0,E,E0);
        lambda+=lambda_step; // increments Bohr radius
    }while((solution_not_found && lambda_stop < 0) || lambda<lambda_stop);
}

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
    opt.add_string_option ("lambdasearch",  "fast", "Method to use for locating Bohr radius (\"fast\" or \"linear\")");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    /* computational default values */
    const size_t N_w=99;             // Number of strips in w integration

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
        if(opt.get_string_option("lambdasearch") == "linear")
            find_E_min_linear(E0[i_d], lambda_0[i_d], lambda_start, lambda_step, lambda_stop, z, epsilon, mstar, r_d[i_d], V, N_w, delta_E);
        else if (opt.get_string_option("lambdasearch") == "fast")
        {
            if(lambda_stop < 0)
            {
                std::cerr << "Upper limit on Bohr radius must be set to a positive value using --lambdastop" << std::endl;
                exit(EXIT_FAILURE);
            }

            find_E_min_fast(E0[i_d], lambda_0[i_d], lambda_start, lambda_stop, z, epsilon, mstar, r_d[i_d], V, N_w, delta_E);
        }
        else
        {
            std::cerr << "Unrecognised search type: " << opt.get_string_option("lambdasearch") << std::endl;
            exit(EXIT_FAILURE);
        }

        SchroedingerSolverDonor2D se(mstar, V, z, epsilon, r_d[i_d], lambda_0[i_d], delta_E, 1);
        std::vector<State> solutions     = se.get_solutions();
        std::vector<State> solutions_chi = se.get_solutions_chi();

        std::valarray<double> psi(solutions[0].psi_array());
        std::valarray<double> chi(solutions_chi[0].psi_array());

        /* generate output filename (and open file for writing) using 
           the basis wf%i.r where the integer %i is the donor index i_d  */
        char   filename[9];     /* character string for wavefunction filename  */
        sprintf(filename,"wf%i.r",i_d);

        write_table_xyz(filename, z, psi, chi);
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

/**
 * \brief compares minimum value of energy for this lambda
 *
 * \param[in]     lambda   The new Bohr radius to be tested [m]
 * \param[in,out] lambda_0 The previous estimate of the Bohr radius [m]
 * \param[in]     E        The new binding energy [J]
 * \param[in,out] E0       The previous estimate of the binding energy [J]
 *
 * \returns True if the binding energy is lower than the previous estimate
 *
 * \details The new binding energy is compared with the previous value. If it
 *          is lower, then the Bohr radius and binding energy are overwritten
 *          and the function returns "true" to indicate that we haven't yet
 *          found the absolute minimum.
 */
static bool repeat_lambda(const double  lambda,
                          double       &lambda_0,
                          const double  E,
                          double       &E0)
{
    bool flag = false; // True if we find a new minimum 

    // Set new minimum if required
    if(E < E0)
    {
        E0       = E;
        lambda_0 = lambda;
        flag=true;
    }

    return flag;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
