#include <gtest/gtest.h>
#include "qwwad-schroedinger-tridiagonal.h"
#include "qwwad/constants.h"

using namespace Leeds;
using namespace QWWAD;
using namespace constants;

TEST(SchroedingerSolverTridiag, parabolicInfTest)
{
    // Using slightly bonkers mass/length!
    const double L     = 1;
    const size_t nz    = 500;
    const double nst   = 50;

    std::valarray<double> m(nz);
    std::valarray<double> V(nz);
    std::valarray<double> z(nz);

    // Note that the well is actually bounded by the
    // box boundary conditions, so the walls lie at the points OUTSIDE the
    // list of z-values
    const double dz = L/(nz+1);

    // Create an infinite well (just rely on the box boundary)
    for(unsigned int iz = 0; iz < nz; ++iz)
    {
        z[iz] = dz*(iz+1);
        m[iz] = 1;
        V[iz] = 0.0;
    }
    EXPECT_DOUBLE_EQ(L-dz,z[z.size()-1]);

    SchroedingerSolverTridiag se(m, V, z, nst);

    // Check that only one state is found using default options
    const std::vector<State> solutions = se.get_solutions();
    EXPECT_EQ(solutions.size(), nst);

    // Check that all states are as expected
    for (unsigned int ist = 0; ist < nst; ++ist)
    {
        // Check energy of state (to within 1%)
        const double E = solutions.at(ist).get_E();
        const double E_expected = pow(hBar*pi*(ist+1),2.0)/(2.0*m[0]*L*L);
        EXPECT_NEAR(E_expected, E, E_expected/100);

        // Check normalisation of state
        const std::valarray<double> PD = solutions.at(ist).psi_squared();
        const double integral_norm = integral(PD,dz);
        EXPECT_NEAR(1.0, integral_norm, 1e-10);

        // Check expectation position (should be middle of well)
        const std::valarray<double> z_expt_dz = PD*z;
        const double z_expt = integral(z_expt_dz, dz);
        EXPECT_NEAR(L/2.0, z_expt, 0.001*L/2.0);

        // Check peak magnitude of probability density (allow 1% error)
        EXPECT_NEAR(2.0/L, PD.max(), 0.01*2.0/L);

        // Check that wavefunction is bounded
        // TODO: This is tricky because the samples at the edges
        //       aren't actually zero.  The box boundary conditions
        //       actually just force psi to zero at the points OUTSIDE
        //       the z-array.
//        EXPECT_NEAR(0.0, PD[0], PD.max()/100);
//        EXPECT_NEAR(0.0, PD[z.size()-1], PD.max()/100);
    }
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
