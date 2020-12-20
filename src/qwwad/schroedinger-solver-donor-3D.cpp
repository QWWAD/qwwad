/**
 * \file   schroedinger-solver-donor-3D.cpp
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Schrodinger solver functions for donor with 3D symmetrical hydrogenic wavefunction
 */

#include "schroedinger-solver-donor-3D.h"

#include <cmath>
#include "constants.h"

namespace QWWAD
{
using namespace constants;

/**
 * \brief Initialise a Schroedinger solver for donors with spherical symmetry
 *
 * \param[in] m      Effective mass [kg]
 * \param[in] V      Potential profile [J]
 * \param[in] z      Spatial sampling points [m]
 * \param[in] eps    Permittivity of material [F/m]
 * \param[in] r_d    Donor location [m]
 * \param[in] lambda Bohr radius [m]
 * \param[in] dE     Energy separation to use in eigenvalue search [J]
 */
SchroedingerSolverDonor3D::SchroedingerSolverDonor3D(const double     m,
                                                     const arma::vec &V,
                                                     const arma::vec &z,
                                                     const double     eps,
                                                     const double     r_d,
                                                     const double     lambda,
                                                     const double     dE) :
    SchroedingerSolverDonor(m, V, z, eps, r_d, lambda, dE)
{}

/**
 * \brief Computes the binding energy integral \f$I_1\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [m^2]
 *
 * \details See Eq. 5.69, QWWAD3. The integral is solved analytically as
 *          \f[
 *             I_{1}=2\pi\left(
 *             \frac{\lambda\vert z^{\prime}\vert}{2}+\frac{\lambda^{2}}{4}\right)
 *             e^{-\frac{2\vert z^{\prime}\vert}{\lambda}}\label{I13D2}
 *          \f]
 */
auto SchroedingerSolverDonor3D::I_1(const double z_dash) const -> double
{
    return 2*pi*(fabs(z_dash)*_lambda/2 + _lambda*_lambda/4) * exp(-2*fabs(z_dash)/_lambda);
}

/**
 * \brief Computes the binding energy integral \f$I_2\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.77, QWWAD3. The integral evaluates to
 *          \f[
 *             I_2=2\pi\left(-\frac{z^\prime}{2}\mbox{e}^{-\frac{2\vert z^{\prime}\vert}
 *             {\lambda}}\right)
 *          \f]
 */
auto SchroedingerSolverDonor3D::I_2(const double z_dash) const -> double
{
    return 2*pi*(-z_dash/2)*exp(-2*fabs(z_dash)/_lambda);
}

/**
 * \brief Computes the binding energy integral \f$I_3\f$ for a 3D trial wavefunction
 *
 * \returns Binding energy integral [dimensionless]
 *
 * \details See Eq. 5.89, QWWAD3. The integral is solved analytically as
 *          \f[
 *            I_3=2\pi\left(\frac{\vert z^\prime\vert}{2\lambda}-\frac{3}{4}\right)
 *            \mbox{e}^{-\frac{2\vert z^\prime\vert}{\lambda}}\label{I33D3}
 *          \f]
 */
auto SchroedingerSolverDonor3D::I_3(const double z_dash) const -> double
{
    return 2*pi*(fabs(z_dash)/(2*_lambda)-0.75)*exp(-2*fabs(z_dash)/_lambda);
}

/**
 * \brief Computes the binding energy integral \f$I_4\f$ for a 3D trial wavefunction
 *
 * \param[in] z_dash displacement between electron and donor in z-direction [m]
 *
 * \returns Binding energy integral [m]
 *
 * \details See Eq. 5.92, QWWAD3. The integral is given by
 *          \f[
 *             I_4=2\pi\left(\frac{\lambda}{2}\mbox{e}^{-\frac{2\vert z^\prime\vert}{\lambda}}
 *          \f]
 */
auto SchroedingerSolverDonor3D::I_4(const double z_dash) const -> double
{
    return 2*pi*(_lambda/2)*exp(-2*fabs(z_dash)/_lambda);
}

auto
SchroedingerSolverDonor3D::calculate_psi_from_chi() -> std::vector<Eigenstate>
{
    std::vector<Eigenstate> solutions;

    auto z = get_z();

    for (auto ist : _solutions_chi)
    {
        const auto E   = ist.get_energy();
        const auto chi = ist.get_wavefunction_samples();
        const auto psi = exp(-abs(z - _r_d)/_lambda) * chi;

        const auto psi_state = Eigenstate(E,z,psi);

        solutions.push_back(psi_state);
    }

    return solutions;
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
