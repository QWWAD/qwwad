/**
 * \file   qwwad_ef_donor_generic.cpp
 * \brief  Use variational technique to minimise Hamiltonian of donor state
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <utility>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include "qwwad/constants.h"
#include "qwwad/file-io.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

/**
 * \brief Identifiers for different states
 */
enum StateID
{
    STATE_1S,
    STATE_2S,
    STATE_2PX,
    STATE_2PZ
};

class Wavefunction3D
{
private:
    double    _lambda; ///< Bohr radius
    arma::vec _wf;     ///< Wavefunction
    arma::vec _V;      ///< Potential profile
    arma::vec _z;      ///< Spatial samples in z

public:
    double    epsilon; ///< Permittivity
    double    m;       ///< Effective mass
    double    r_i;     ///< z-position of dopant
    StateID   S;       ///< State ID

    Wavefunction3D(decltype(_wf) wf,
                   decltype(_V)  V,
                   decltype(_z)  z,
                   double        epsilon,
                   double        m,
                   double        r_i,
                   StateID       S)
        :
            _wf(std::move(wf)),
            _V(std::move(V)),
            _z(std::move(z)),
            epsilon(epsilon),
            m(m),
            r_i(r_i),
            S(S)
    {}

private:
    double get_Laplacian(double x,
                         double y,
                         unsigned int iz) const;
public:
    /**
     * \brief Set the Bohr radius
     */
    void set_lambda(const double lambda) {_lambda = lambda;}

    double get_psi(double x, double y, unsigned int iz) const;

    double get_energy_integrand(double       x,
                                double       y,
                                unsigned int iz) const;

    static double get_energy_integrand_y(double  y,
                                         void   *params);

    static double get_PD_integrand_y(double  y,
                                     void   *params);

    static double get_energy_integrand_x(double  x,
                                         void   *params);
    
    static double get_PD_integrand_x(double  x,
                                     void   *params);

    static double get_energy(double  lambda,
                             void   *params);

    double get_energy() const;
};

/**
 * \brief Configure command-line options
 */
Options configure_options(int argc, char** argv)
{
    Options opt;
    std::string doc("Find state of electron attached to a donor in a 2D system using a generic search");

    opt.add_option<double>     ("dcpermittivity,e",   13.18, "Bulk relative permittivity");
    opt.add_option<double>     ("mass,m",             0.067, "Bulk effective mass (relative to free electron)");
    opt.add_option<char>       ("particle,p",           'e', "ID of particle to be used: 'e', 'h' or 'l' for "
                                                             "electrons, heavy holes or light holes respectively.");
    opt.add_option<size_t>     ("subband",                1, "Principal quantum number of subband for which to find impurity state.");
    opt.add_option<std::string>("impuritystate",       "1s", "Symmetry of impurity state (1s, 2s, 2px, or 2pz)");
    opt.add_option<double>     ("donorposition,r",           "Location of donor ion [Angstrom]");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    return opt;
}

int main(int argc,char *argv[])
{
    const auto opt = configure_options(argc, argv);

    const auto epsilon    = opt.get_option<double>("dcpermittivity") * eps0; // Permittivity [F/m]
    const auto m          = opt.get_option<double>("mass") * me;             // Effective mass [kg]
    const auto p          = opt.get_option<char>  ("particle");              // Particle ID (e, h, or l)
    const auto subband    = opt.get_option<size_t>("subband");               // Principal quantum number of state to find

    StateID S = STATE_1S;

    const auto impuritystate_string = opt.get_option<std::string>("impuritystate");

    if (impuritystate_string == "1s")
        S = STATE_1S;
    else if (impuritystate_string == "2s")
        S = STATE_2S;
    else if (impuritystate_string == "2px")
        S = STATE_2PX;
    else if (impuritystate_string == "2pz")
        S = STATE_2PZ;
    else
    {
        std::cerr << "Unknown impurity state ID: " << impuritystate_string << std::endl;
        exit(EXIT_FAILURE);
    }

    arma::vec z; // Spatial location [m]
    arma::vec V; // Confining potential [J]
    read_table("v.r", z, V);

    std::ostringstream filename; // input filename
    filename << "wf_" << p << subband << ".r";

    arma::vec z_tmp; // Dummy file for unused spatial locations
    arma::vec wf;    // Wave function samples at each point [m^{-1/2}]
    read_table(filename.str(), z_tmp, wf);

    const auto lambda_0=4*pi*epsilon*(hBar/e)*(hBar/e)/m;/* Bohr	theory (1s)	*/

    // Get donor location [m].  If unspecified, assume it's in the middle
    const auto nz = z.size();
    auto r_d = (z[nz-1] + z[0])/2.0;

    if (opt.get_argument_known("donorposition") > 0)
        r_d = opt.get_option<double>("donorposition") * 1e-10;

    // Perform variational calculation for each donor/acceptor position
    double lambda=lambda_0;	// initial lambda guess

    // Double the estimate of Bohr radius if we're in a second orbital
    // This isn't correct for 2pz, but it's still better than the 1s
    // estimate!
    if((S==STATE_2S)||(S==STATE_2PX)||(S==STATE_2PZ))lambda*=2;

    // The 3D wavefunction calculator for the system
    Wavefunction3D wf3d(wf, V, z, epsilon, m, r_d, S);

    // Set up the numerical solver using GSL
    gsl_function f;
    f.function = &Wavefunction3D::get_energy;
    f.params   = &wf3d;

    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(s, &f, lambda*10, lambda/100, lambda*100);

    size_t max_iter = 100; // Maximum number of iterations before giving up
    int status = 0;        // Error flag for GSL
    unsigned int iter=0;   // The number of iterations attempted so far

    double E = 1000*e; // Minimum energy of carrier [J]

    // Vectors to store the history of the search
    std::vector<double> lambda_history;
    std::vector<double> E_history;

    // Variational calculation (search over lambda)
    do
    {
        ++iter;
        status  = gsl_min_fminimizer_iterate(s);

        if(status) {
            std::cerr << "GSL error in qwwad_ef_donor_generic: " << std::endl;
        }

        const auto lambda_lo = gsl_min_fminimizer_x_lower(s);
        const auto lambda_hi = gsl_min_fminimizer_x_upper(s);
        lambda = gsl_min_fminimizer_x_minimum(s);
        E      = gsl_min_fminimizer_f_minimum(s);

        // Update the search history
        lambda_history.push_back(lambda);
        E_history.push_back(E);

        status = gsl_min_test_interval(lambda_lo, lambda_hi, 0.1e-10, 0.0);
    }while((status == GSL_CONTINUE) && (iter < max_iter));

    gsl_min_fminimizer_free(s);

    // Save the search log
    write_table("searchlog.r",
                lambda_history,
                E_history);

    // Output neutral dopant binding energies (E) and 
    // Bohr radii (lambda) in meV and Angstrom respectively
    std::vector<unsigned int> indices(1, subband);
    std::vector<double> E_out(1, E*1000.0/e);
    std::vector<double> lambda_out(1,lambda*1.0e10);
    write_table("Ee.r", indices, E_out);
    write_table("l.r",  indices, lambda_out);


    // Get wavefunction at (x,y=0)
//    std::valarray<double> wf_out(nz);

//    for(unsigned int iz = 0; iz < nz; ++iz)
//        wf_out[iz] = Psi(wf[iz], lambda, 0, 0, z[iz], S);

//    std::ostringstream wf_file;
//    wf_file << "wf_" << p << subband << "_out.r";
//    write_table(wf_file.str(), z, wf_out);

    return EXIT_SUCCESS;
}

/// The Laplacian of Psi
double Wavefunction3D::get_Laplacian(double x,
                                     double y,
                                     unsigned int iz) const
{
    // For the in-plane derivative, use a very small step, for accuracy
    // For the growth direction, use the sample spacing
    const auto dxy = _lambda/100;
    const auto dz = _z[1] - _z[0];

    const auto psi_111 = get_psi(x, y, iz);
    const auto psi_110 = get_psi(x, y, iz-1);
    const auto psi_112 = get_psi(x, y, iz+1);
    const auto psi_101 = get_psi(x, y-dxy, iz);
    const auto psi_121 = get_psi(x, y+dxy, iz);
    const auto psi_011 = get_psi(x-dxy, y, iz);
    const auto psi_211 = get_psi(x+dxy, y, iz);

    const auto d2Pdx2 = (psi_011 - 2*psi_111 + psi_211)/(dxy*dxy);
    const auto d2Pdy2 = (psi_101 - 2*psi_111 + psi_121)/(dxy*dxy);
    const auto d2Pdz2 = (psi_110 - 2*psi_111 + psi_112)/(dz*dz);
    const auto laplace_Psi = d2Pdx2 + d2Pdy2 + d2Pdz2;

    return laplace_Psi;
}

/// Get the energy integrand
double Wavefunction3D::get_energy_integrand(double       x,
                                            double       y,
                                            unsigned int iz) const
{
    // Pre-calculate a couple of params
    const auto hBar_sq_by_2m  = hBar*hBar/(2.0*m);
    const auto e_sq_by_4pieps = e*e/(4.0*pi*epsilon);

    const auto z_dash = _z[iz] - r_i; // Separation from donor in z-direction [m]
    const auto r_xz   = hypot(x, z_dash);
    const auto r      = hypot(y, r_xz);

    const auto laplace_Psi = get_Laplacian(x, y, iz);

    // The integrand for the Hamiltonian expectation value
    // QWWAD 3, Eq. 5.142
    const auto Psixyz = get_psi(x, y, iz);

    return Psixyz*(-hBar_sq_by_2m*laplace_Psi
                   + (_V[iz]-e_sq_by_4pieps/r)*Psixyz);
}

struct Integrand_y_params
{
    const Wavefunction3D *wf3d;
    double                x;
    unsigned int          iz;
};

double Wavefunction3D::get_energy_integrand_y(double  y,
                                              void   *params)
{
    auto *p = reinterpret_cast<Integrand_y_params *>(params);

    return p->wf3d->get_energy_integrand(p->x, y, p->iz);
}

double Wavefunction3D::get_PD_integrand_y(double  y,
                                          void   *params)
{
    auto *p = reinterpret_cast<Integrand_y_params *>(params);
    auto psi = p->wf3d->get_psi(p->x, y, p->iz);

    return psi*psi;
}

struct Integrand_x_params
{
    const Wavefunction3D *wf3d;
    unsigned int          iz;
};

double Wavefunction3D::get_energy_integrand_x(double  x,
                                              void   *params)
{
    auto *p = reinterpret_cast<Integrand_x_params *>(params);
    auto wf3d = p->wf3d;

    // Perform integration over y, noting that a factor of 2 is included
    // to account for even symmetry
    gsl_function F;
    Integrand_y_params integrand_y_params = {wf3d, x, p->iz};
    auto w = gsl_integration_workspace_alloc(1000);
    F.function = &Wavefunction3D::get_energy_integrand_y;
    F.params   = &integrand_y_params;
    double result, error;
    gsl_integration_qags(&F, 0, 5*wf3d->_lambda, 0, 1e-3, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return 2.0*result;
}

double Wavefunction3D::get_PD_integrand_x(double  x,
                                          void   *params)
{
    auto *p = reinterpret_cast<Integrand_x_params *>(params);
    auto wf3d = p->wf3d;

    // Perform integration over y, noting that a factor of 2 is included
    // to account for even symmetry
    gsl_function F;
    Integrand_y_params integrand_y_params = {wf3d, x, p->iz};
    auto w = gsl_integration_workspace_alloc(1000);
    F.function = &Wavefunction3D::get_PD_integrand_y;
    F.params   = &integrand_y_params;
    double result, error;
    gsl_integration_qags(&F, 0, 5*wf3d->_lambda, 0, 1e-3, 1000, w, &result, &error);
    gsl_integration_workspace_free(w);

    return 2.0*result;
}

/**
 * \brief Calculates the expectation value (the energy) of the Hamiltonian operator
 */
double Wavefunction3D::get_energy(double  lambda,
                                  void   *params)
{
    auto *p = reinterpret_cast<Wavefunction3D *>(params);
    p->set_lambda(lambda);

    return p->get_energy();
}

double Wavefunction3D::get_energy() const
{
    const auto dz  = _z[1] - _z[0]; // z- (growth) direction step length [m]
    const auto nz  = _z.size();    // Number of spatial samples in z direction

    const size_t nslice = 1000; // Maximum number of slices for in-plane integration
    const double relerr = 1e-3; // Maximum relative error for in-plane integration

    // Integrands wrt z for calculating wavefunction overlap
    // and Hamiltonian
    arma::vec PD_integrand_z(nz);
    arma::vec H_integrand_z(nz);
    auto w = gsl_integration_workspace_alloc(nslice);

    // Compute integrand over the z-axis, skipping both end-points since we
    // need the 2nd derivatives
    for(unsigned int iz=1; iz < nz-1; ++iz)
    {
        Integrand_x_params integrand_x_params = {this, iz};
        double result, error;
        
        gsl_function F_energy;
        F_energy.function = &Wavefunction3D::get_energy_integrand_x;
        F_energy.params   = &integrand_x_params;
        gsl_integration_qags(&F_energy, 0, 5*_lambda, 0, relerr, nslice, w, &result, &error);
        H_integrand_z[iz]  = 2*result;

        gsl_function F_PD;
        F_PD.function = &Wavefunction3D::get_PD_integrand_x;
        F_PD.params   = &integrand_x_params;
        gsl_integration_qags(&F_PD, 0, 5*_lambda, 0, relerr, nslice, w, &result, &error);
        PD_integrand_z[iz] = 2*result; //simps(PD_integrand_xz, dxy);
    }

    gsl_integration_workspace_free(w);

    // Note that endpoints of the integral can keep their default value of zero, since
    // psi decays to zero at infinity

    // Compute the final value of the energy using Eq. 5.141, QWWAD3
    const auto H_exp = integral(H_integrand_z, dz);
    const auto norm  = integral(PD_integrand_z, dz);
    const auto E = H_exp/norm;

    return E;
}

/**
 * \brief The wave function psi(z)phi(r)
 */
double Wavefunction3D::get_psi(const double       x,
                               const double       y,
                               const unsigned int iz) const
{
    const double z_dash = std::abs(_z[iz] - r_i);
    const auto r_xy = hypot(x,y);
    const auto r    = hypot(r_xy,z_dash);

    double result = 0.0;

    switch(S)
    {
        case STATE_1S:
            result = _wf[iz]*exp(-r/_lambda);
            break;
        case STATE_2S:
            result = _wf[iz]*(1.0-r/_lambda)*exp(-r/_lambda);
            break;
        case STATE_2PX:
            result = _wf[iz]*fabs(x)*exp(-r/_lambda);
            break;
        case STATE_2PZ:
            result = _wf[iz]*fabs(z_dash)*exp(-r/_lambda);
            break;
        default:
            throw std::runtime_error("Unrecognised orbital");
    }

    return result;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
