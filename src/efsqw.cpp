/**
 * \file   efsqw.cpp
 * \brief  Calculate the 1-particle energies and wavefunctions in a single quantum well
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \details One particle energies and wavefunctions in a single
   quantum well.  Three basic theoretical approaches are
   contained within this sourcecode

   (i) Constant mass

   -hbar^2 d^2 Y + V Y = E Y,   -lw < z < +lw,  V=0
   ------  ----                  --        --
     2m*   dz^2                   2         2

   (ii) Different masses in well and barrier

   -hbar^2 d^2 Y + V Y = E Y,   -lw < z < +lw,  V=0
   ------  ----                  --        --
     2m*   dz^2                   2         2

   with the additional constraint of the boundary conditions

   Y and 1 dY , continuous
         - --
         m dz
 
   this represents the Hamiltonian P 1 P + V
                                     -
                                     m

   (iii) Different masses in well and barrier

   -hbar^2 d^2 Y + V Y = E Y,   -lw < z < +lw,  V=0
   ------  ----                  --        --
     2m*   dz^2                   2         2

   with the additional constraint of the boundary conditions

   Y and dY , continuous
         --
         dz

   this represents the Hamiltonian 1    P P    1     + V
                                   -           -
                                 sqrt(m)     sqrt(m)

   The code is based around point (ii). Point (i) is 
   implemented by the user selecting m_b=m_w.  Point
   (iii) is implemented by allowing different m_b and m_w
   in the evaluation of k and K, but m_b is artificially
   forced to equal m_w for the boundary conditions.

   The system is solved by expressing the
   standard condition as a function f(x)=0 and
   implementing a Newton-Raphson iterative method
   where the independent variable x is the energy.
   The first estimate of the energy is found by 
   following the function along the x axis until it
   changes sign then using a midpoint rule.
*/

#include <cstdio>
#include <cstdlib>
#include <strings.h>
#include <valarray>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

#define   N  1001      /* number of wavefunction sampling
                          points, must be odd to provide an
                          even number of strips for numerical 
                          integration.                            */

/**
 * Handler for command-line options
 */
class EFSQWOptions : public Options
{
    public:
        EFSQWOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("well-width,a", po::value<double>()->default_value(100),
                     "Width of quantum well [angstrom].")
                    
                    ("barrier-width,b", po::value<double>()->default_value(200),
                     "Width of barrier [angstrom]. Note that this is only used "
                     "for the purposes of outputting the data. The calculation here "
                     "assumes that the barriers are infinitely thick.  As such, the "
                     "wavefunctions do not decay to precisely zero at the boundaries.")

                    ("alt-KE,k", po::bool_switch()->default_value(false),
                     "Use alternative kinetic energy operator (1/m)PP")

                    ("well-mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass in well (relative to that of a free electron)")

                    ("barrier-mass,n", po::value<double>()->default_value(0.067),
                     "Effective mass in barrier (relative to that of a free electron)")

                    ("particle,p", po::value<char>()->default_value('e'),
                     "Particle to be used: 'e', 'h' or 'l'")

                    ("states,s", po::value<size_t>()->default_value(1),
                     "Number of states to find")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")
                    ;

                std::string doc("Find the eigenstates of a single finite quantum well. "
                                "The energies are written to the file \"E*.r\", and the "
                                "wavefunctions are written to \"wf_*i.r\" where the '*' "
                                "is replaced by the particle ID in each case and the "
                                "'i' is replaced by the number of the state");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the width of the quantum well [m]
        double get_well_width() const {return vm["well-width"].as<double>()*1e-10;}
        
        /// \returns the width of the barriers [m]
        double get_barrier_width() const {return vm["barrier-width"].as<double>()*1e-10;}

        /// \returns true if normal kinetic energy operator is to be used
        bool get_T_flag() const {return !vm["alt-KE"].as<bool>();}

        /// \returns the effective mass in the quantum well [kg]
        double get_well_mass() const {return vm["well-mass"].as<double>()*me;}

        /// \returns the effective mass in the quantum well [kg]
        double get_barrier_mass() const {return vm["barrier-mass"].as<double>()*me;}

        /// \returns the particle ID
        char get_particle() const {return vm["particle"].as<char>();}

        /// \returns the number of spatial points
        size_t get_n_states() const {return vm["states"].as<size_t>();}

        /// \returns the effective mass in the quantum well [J]
        double get_potential() const {return vm["potential"].as<double>()*1e-3*e;}
};

/// Parameters needed for computing square-well eigenstates
struct sqw_params
{
    double a;           ///< Well width [m]
    double m_B;         ///< Boundary mass [kg]
    double m_b;         ///< Barrier mass [kg]
    double m_w;         ///< Well mass [kg]
    double V;           ///< Barrier potential [J]
    bool   parity_flag; ///< True for odd states
};

double f(double  energy,
         void   *params);

void wavef(const double a,
           const double b,
           const double E,
           const double m_b,
           const double m_w,
           const double V,
           const int    i_state,
           const char   p,
           const bool   parity_flag);

int main(int argc,char *argv[])
{
    EFSQWOptions opt(argc, argv);
    const double a      = opt.get_well_width();
    const double b      = opt.get_barrier_width();
    const double m_w    = opt.get_well_mass();
    const double m_b    = opt.get_barrier_mass();
    const char   p      = opt.get_particle();
    const double V      = opt.get_potential();   
    const size_t state  = opt.get_n_states();
    const bool   T_flag = opt.get_T_flag();

    const double dx = 1e-4*e; // arbitrarily small energy increment [0.1meV]
    double x = dx; // first energy estimate

    /* Setting the barrier mass equal to the well mass within the boundary 
       conditions produces T=(1/m)PP */
    const double m_B = T_flag?m_b:m_w;

    std::vector<double> E; // Energies of states

    for(unsigned int i_state=1; i_state<=state;i_state++)
    {
        // deduce parity: false if even parity
        const bool parity_flag = (i_state%2 != 1);

        sqw_params params = {a, m_B, m_b, m_w, V, parity_flag};
        gsl_function F;
        F.function = &f;
        F.params   = &params;

        /// Look for a solution in the range (x,x+dx)
        gsl_root_fsolver *solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        int status;
        bool root_found = false;
        double x_end = x;
        bool V_limit_hit = false;

        do
        {
            x = x_end;
            x_end+=dx;

            // Set a coarse bracket for search
            double y1=f(x,     &params);
            double y2=f(x_end, &params);
            root_found = (y1*y2 < 0);
            V_limit_hit = std::isunordered(y1,y2);

            if(root_found)
            {
                E.push_back(0);
                double x_lo = x;
                double x_hi = x_end;
                gsl_root_fsolver_set(solver, &F, x_lo, x_hi);

                do
                {
                    status = gsl_root_fsolver_iterate(solver);
                    E[i_state-1] = gsl_root_fsolver_root(solver);
                    double x_lo = gsl_root_fsolver_x_lower(solver);
                    double x_hi = gsl_root_fsolver_x_upper(solver);
                    status = gsl_root_test_interval(x_lo, x_hi, 1e-15*e, 0);
                }while(status == GSL_CONTINUE);
                wavef(a,b,E[i_state-1],m_b,m_w,V,i_state,p,parity_flag);
            }
        }while(!(root_found or V_limit_hit));
        x+=dx;

        gsl_root_fsolver_free(solver);
    }
    
    char filename[9]; // output filename
    sprintf(filename,"E%c.r",p);
    std::valarray<double> _E(E.size());
    for(unsigned int i=0; i<E.size(); ++i)
        _E[i] = E[i] /(0.001*e);
    write_table_x(filename, _E, true);

    return EXIT_SUCCESS;
}

/**
 * \brief standard fw result =0
 *
 * \param[in] energy local energy
 * \param[in] params sqw params object
 */
double f(double  energy,
         void   *params)
{
    const sqw_params *p = reinterpret_cast<sqw_params *>(params);
    const double a   = p->a;
    const double m_B = p->m_B;
    const double m_b = p->m_b;
    const double m_w = p->m_w;
    const double V   = p->V;
    const bool   parity_flag = p->parity_flag;

    const double k=sqrt(2*m_w*energy)/hBar; // electron wave vector
    const double K=sqrt(2*m_b*(V-energy))/hBar; // wavefunction decay constant

    double result = 0.0;

    if(parity_flag)
        result = k*cot(k*a/2)/m_w+K/m_B;
    else
        result = k*tan(k*a/2)/m_w-K/m_B;

    return result;
}

/**
 * \brief calculates the uncorrelated one particle wavefunctions for the electron and hole and writes to an external file.
 *
 * \param[in] a           well width
 * \param[in] b           barrier width
 * \param[in] E           local energy
 * \param[in] m_b         barrier mass
 * \param[in] m_w         well mass
 * \param[in] V           barrier potential
 * \param[in] i_state     state index
 * \param[in] p           particle ID
 * \param[in] parity_flag true for odd states, false for even
 */
void wavef(const double a,
           const double b,
           const double E,
           const double m_b,
           const double m_w,
           const double V,
           const int    i_state,
           const char   p,
           const bool   parity_flag)
{
    double A;         /* In the well the wavefunction psi=Acoskz */
    const double B=1; /* and in the barrier  psi=Bexp(-Kz)       */
    double norm_int;  /* integral over all space of psi*psi      */

    // Generate z co-ordinates
    std::valarray<double> z(N);

    for (unsigned int i_z=0;i_z<N;i_z++)
        z[i_z]=i_z*(a+2*b)/(N-1)-(b+a/2);

    // Define k and K
    const double k=sqrt(2*m_w/hBar*E/hBar); // wave vector in the well
    const double K=sqrt(2*m_b/hBar*(V-E)/hBar); // decay constant in barrier
    
    std::valarray<double> psi(N); // wavefunction

    if(parity_flag) // odd parity wavefunction
    {
        A=B*exp(-K*a/2)/sin(k*a/2);

        for (unsigned int i_z=0;i_z<N;i_z++) // calculate wavefunctions
        {
            if (z[i_z]<(-a/2)) // Left barrier
            {
                psi[i_z]=-B*exp(-K*fabs(z[i_z]));
            }
            if ((z[i_z]>=(-a/2))&&(z[i_z]<(a/2))) // Well
            {
                psi[i_z]=A*sin(k*z[i_z]);
            }
            if (z[i_z]>=(a/2)) // Right barrier
            {
                psi[i_z]=B*exp(-K*z[i_z]);
            }
        }

        // normalisation integral for odd parity type I
        norm_int=gsl_pow_2(A)*(a/2-sin(k*a)/(2*k))-
            gsl_pow_2(B)*exp(-K*a)*(exp(-2*K*b)-1)/K;
    }
    else // even parity wavefunction
    {
        A=B*exp(-K*a/2)/cos(k*a/2);

        for (unsigned int i_z=0;i_z<N;i_z++) // calculate wavefunctions
        {
            if (z[i_z]<(-a/2)) // Left barrier
            {
                psi[i_z]=B*exp(-K*fabs(z[i_z]));
            }
            if ((z[i_z]>=(-a/2))&&(z[i_z]<(a/2))) // Well
            {
                psi[i_z]=A*cos(k*z[i_z]);
            }
            if (z[i_z]>=(a/2)) // Right barrier
            {
                psi[i_z]=B*exp(-K*z[i_z]);
            }
        }

        // normalisation integral for even parity type I
        norm_int=gsl_pow_2(A)*(a/2+sin(k*a)/(2*k))+
            gsl_pow_2(B)*exp(-K*a)*(1-exp(-2*K*b))/K;
    }

    // normalise wavefunction
    psi/=sqrt(norm_int);

    char filename[9]; // output filename
    sprintf(filename,"wf_%c%i.r",p,i_state);
    write_table_xy(filename, z, psi);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
