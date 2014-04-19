/**
 * \file   tdb.cpp
 * \brief  Calculate transmission coefficient for double barrier structure
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 *
 * \details calculates the transmission coefficient T(E)
 *          for a double barrier structure as a function of the 
 *          incident energy E of a particle.  The barriers can be of 
 *          unequal width and the particle can have a different mass in
 *          the barrier material to the `well', the barrier heights are
 *          however equal.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "struct.h"
#include "maths.h"
#include "qclsim-constants.h"
#include "qwwad-options.h"

using namespace Leeds;
using namespace constants;

/**
 * Handler for command-line options
 */
class TDBOptions : public Options
{
    public:
        TDBOptions(int argc, char* argv[])
        {
            try
            {
                program_specific_options->add_options()
                    ("left-barrier-width,a", po::value<double>()->default_value(100),
                     "Width of left barrier [angstrom].")

                    ("well-width,b", po::value<double>()->default_value(100),
                     "Width of well [angstrom].")

                    ("right-barrier-width,c", po::value<double>()->default_value(100),
                     "Width of right barrier [angstrom].")

                    ("well-mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass in well (relative to that of a free electron).")

                    ("barrier-mass,m", po::value<double>()->default_value(0.067),
                     "Effective mass in both barriers (relative to that of a free electron).")

                    ("potential", po::value<double>()->default_value(100),
                     "Barrier potential [meV]")

                    ("energy-step,d", po::value<double>()->default_value(0.01),
                     "Energy step [meV]")
                    ;

                std::string doc("Find the transmission coefficient through a double tunnelling barrier. "
                                "The values are written as a function of energy to the file T.r.");

                add_prog_specific_options_and_parse(argc, argv, doc);	
            }
            catch(std::exception &e)
            {
                std::cerr << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        /// \returns the step in energy required for the output file [J]
        double get_energy_step() const {return vm["energy-step"].as<double>()*1e-3*e;}

        /// \returns the width of the left barrier [m]
        double get_left_barrier_width() const {return vm["left-barrier-width"].as<double>()*1e-10;}

        /// \returns the width of the well [m]
        double get_well_width() const {return vm["well-width"].as<double>()*1e-10;}

        /// \returns the width of the right barrier [m]
        double get_right_barrier_width() const {return vm["right-barrier-width"].as<double>()*1e-10;}

        /// \returns the effective mass in the well [kg]
        double get_well_mass() const {return vm["well-mass"].as<double>()*me;}

        /// \returns the effective mass in the barriers [kg]
        double get_barrier_mass() const {return vm["barrier-mass"].as<double>()*me;}

        /// \returns the effective mass in the quantum well [J]
        double get_potential() const {return vm["potential"].as<double>()*1e-3*e;}
};

typedef std::complex<double> cdouble;

static cmat2x2 invcmat2x2 (const cmat2x2 M);
static cdouble detcmat2x2 (const cmat2x2 M);
static cmat2x2 cmat2x2mult(const cmat2x2 M1,
                           const cmat2x2 M2);

int main(int argc,char *argv[])
{
    const TDBOptions opt(argc, argv);

    const double dE  = opt.get_energy_step();         // [J]
    const double L1  = opt.get_left_barrier_width();  // [m]
    const double L2  = opt.get_well_width();          // [m]
    const double L3  = opt.get_right_barrier_width(); // [m]
    const double m_w = opt.get_well_mass();           // [kg]
    const double m_b = opt.get_barrier_mass();        // [kg]
    const double V   = opt.get_potential();           // [J]

    // Calculate interfaces
    const double I2=L1;
    const double I3=L1+L2;
    const double I4=L1+L2+L3;

    // Initialise energy

    double E=dE; 

    // Open output file
    FILE *FT=fopen("T.r","w");

    // Loop over energy
    do
    {
        const double k=sqrt(2*m_w*E)/hBar;
        const double K=sqrt(2*m_b*(V-E))/hBar;

        // Define transfer matrices
        cmat2x2 M1;
        M1.M[0][0] = 1;
        M1.M[0][1] = 1;
        M1.M[1][0] = cdouble(0.0,  k/m_w);
        M1.M[1][1] = cdouble(0.0, -k/m_w);

        cmat2x2 M2;
        M2.M[0][0] = 1;
        M2.M[0][1] = 1;
        M2.M[1][0] = +K/m_b;
        M2.M[1][1] = -K/m_b;

        cmat2x2 M3;
        M3.M[0][0] = exp(+K*I2);
        M3.M[0][1] = exp(-K*I2);
        M3.M[1][0] =  K*exp(+K*I2)/m_b;
        M3.M[1][1] = -K*exp(-K*I2)/m_b;

        cmat2x2 M4;
        M4.M[0][0] = cdouble(cos(k*I2),         +sin(k*I2));
        M4.M[0][1] = cdouble(cos(k*I2),         -sin(k*I2));
        M4.M[1][0] = cdouble(-k*sin(+k*I2)/m_w, +k*cos(k*I2)/m_w);
        M4.M[1][1] = cdouble(+k*sin(-k*I2)/m_w, -k*cos(k*I2)/m_w);

        cmat2x2 M5;
        M5.M[0][0] = cdouble(cos(k*I3),         +sin(k*I3));
        M5.M[0][1] = cdouble(cos(k*I3),         -sin(k*I3));
        M5.M[1][0] = cdouble(-k*sin(+k*I3)/m_w, +k*cos(k*I3)/m_w);
        M5.M[1][1] = cdouble(+k*sin(-k*I3)/m_w, -k*cos(k*I3)/m_w);

        cmat2x2 M6;
        M6.M[0][0] = exp(+K*I3);
        M6.M[0][1] = exp(-K*I3);
        M6.M[1][0] =  K*exp(+K*I3)/m_b;
        M6.M[1][1] = -K*exp(-K*I3)/m_b;

        cmat2x2 M7;
        M7.M[0][0] = exp(+K*I4);
        M7.M[0][1] = exp(-K*I4);
        M7.M[1][0] =  K*exp(+K*I4)/m_b;
        M7.M[1][1] = -K*exp(-K*I4)/m_b;

        cmat2x2 M8;
        M8.M[0][0] = cdouble(cos(k*I4),         +sin(k*I4));
        M8.M[0][1] = cdouble(cos(k*I4),         -sin(k*I4));
        M8.M[1][0] = cdouble(-k*sin(+k*I4)/m_w, +k*cos(k*I4)/m_w);
        M8.M[1][1] = cdouble(+k*sin(-k*I4)/m_w, -k*cos(k*I4)/m_w);

        cmat2x2 M=cmat2x2mult(invcmat2x2(M1),
                cmat2x2mult(M2,
                    cmat2x2mult(invcmat2x2(M3),
                        cmat2x2mult(M4,
                            cmat2x2mult(invcmat2x2(M5),
                                cmat2x2mult(M6,
                                    cmat2x2mult(invcmat2x2(M7),M8)
                                    )
                                )
                            )
                        )
                    )
                );

        const double T=1/(norm(M.M[0][0])); // Transission coeff

        fprintf(FT,"%20.17le %20.17le\n",E/(1e-3*e),T);
        E+=dE;
    }while(E<V);

    fclose(FT);

    return EXIT_SUCCESS;
}

/**
 * \brief Multiplies two complex 2x2 matrices together
 */
static cmat2x2 cmat2x2mult(const cmat2x2 M1,
                           const cmat2x2 M2)
{
    cmat2x2 M;

    M.M[0][0] = M1.M[0][0] * M2.M[0][0] + M1.M[0][1] * M2.M[1][0];
    M.M[0][1] = M1.M[0][0] * M2.M[0][1] + M1.M[0][1] * M2.M[1][1];
    M.M[1][0] = M1.M[1][0] * M2.M[0][0] + M1.M[1][1] * M2.M[1][0];
    M.M[1][1] = M1.M[1][0] * M2.M[0][1] + M1.M[1][1] * M2.M[1][1];

    return M;
}

/**
 * \brief Calculates the inverse of a complex 2x2
 */
static cmat2x2 invcmat2x2(const cmat2x2 M)
{
 cmat2x2 Minv;

 cdouble determinant = detcmat2x2(M);

 Minv.M[0][0] =  M.M[1][1] / determinant;
 Minv.M[0][1] = -M.M[0][1] / determinant;
 Minv.M[1][0] = -M.M[1][0] / determinant;
 Minv.M[1][1] =  M.M[0][0] / determinant;

 return Minv;
}

/**
 * \brief Calculates the determinant of a complex 2x2
 */
static cdouble detcmat2x2(const cmat2x2 M)
{
 return M.M[0][0] * M.M[1][1] - M.M[0][1] * M.M[1][0];
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
