#include <gtest/gtest.h>
#include "../src/qwwad-schroedinger-infinite-well.h"
#include "../src/qclsim-constants.h"

using namespace Leeds;
using namespace constants;

TEST(SchroedingerSolverInfWell, defaultTest)
{
    const double m   = 1;
    const double L   = 1;
    const size_t nz  = 100;

    SchroedingerSolverInfWell se(m, L, nz);

    // Check that only one state is found using default options
    const std::vector<State> solutions = se.get_solutions();
    EXPECT_EQ(solutions.size(), 1);

    EXPECT_STREQ("infinite-square-well", se.get_name().c_str());

    // Check that its energy is as expected
    const double E = solutions.at(0).get_E();
    const double E_expected = hBar*hBar*pi*pi/2.0;
    EXPECT_DOUBLE_EQ(E_expected, E);
}

TEST(SchroedingerSolverInfWell, parabolicTest)
{
    // Using slightly bonkers mass/length!
    const double m     = 1;
    const double L     = 1;
    const size_t nz    = 500;
    const double alpha = 0;
    const double V     = 0;
    const double nst   = 100;

    SchroedingerSolverInfWell se(m, L, nz, alpha, V, nst);

    // Check that only one state is found using default options
    const std::vector<State> solutions = se.get_solutions();
    EXPECT_EQ(solutions.size(), nst);
    const std::valarray<double> z = se.get_z();
    const double dz = z[1] - z[0];

    // Check that all states are as expected
    for (unsigned int ist = 0; ist < nst; ++ist)
    {
        // Check energy of state
        const double E = solutions.at(ist).get_E();
        const double E_expected = pow(hBar*pi*(ist+1),2.0)/2.0;
        EXPECT_DOUBLE_EQ(E_expected, E);

        // Check normalisation of state
        const std::valarray<double> PD = solutions.at(ist).psi_squared();
        const double integral_norm = trapz(PD,dz);
        EXPECT_NEAR(1.0, integral_norm, 1e-10);

        // Check expectation position (should be middle of well)
        const double z_expt = trapz(PD*z, dz);
        EXPECT_NEAR(L/2.0, z_expt, 0.001*L/2.0);

        // Check peak magnitude of probability density (allow 1% error)
        EXPECT_NEAR(2.0/L, PD.max(), 0.01*2.0/L);

        // Check that wavefunction is bounded
        EXPECT_DOUBLE_EQ(0.0, PD[0]);
        EXPECT_NEAR(0.0, PD[z.size()-1], 1.0e-10);
    }
}

TEST(SchroedingerSolverInfWell, nonParabolicTest)
{
    // Using slightly bonkers mass/length!
    const double m     = 1;
    const double L     = 1;
    const size_t nz    = 500;
    const double alpha = 1;
    const double V     = 0;
    const double nst   = 100;

    SchroedingerSolverInfWell se(m, L, nz, alpha, V, nst);

    // Check that only one state is found using default options
    const std::vector<State> solutions = se.get_solutions();
    EXPECT_EQ(solutions.size(), nst);
    const std::valarray<double> z = se.get_z();
    const double dz = z[1] - z[0];

    // Check that all states are as expected
    for (unsigned int ist = 0; ist < nst; ++ist)
    {
        // Check energy of state
        const double E = solutions.at(ist).get_E();
        const double E_expected = 0.5 * (sqrt(pow(hBar*pi*(ist+1),2.0) + 1.0) - 1.0);
        EXPECT_DOUBLE_EQ(E_expected, E);

        // Check normalisation of state
        const std::valarray<double> PD = solutions.at(ist).psi_squared();
        const double integral_norm = trapz(PD,dz);
        EXPECT_NEAR(1.0, integral_norm, 1e-10);

        // Check expectation position (should be middle of well)
        const double z_expt = trapz(PD*z, dz);
        EXPECT_NEAR(L/2.0, z_expt, 0.001*L/2.0);

        // Check peak magnitude of probability density (allow 1% error)
        EXPECT_NEAR(2.0/L, PD.max(), 0.01*2.0/L);

        // Check that wavefunction is bounded
        EXPECT_DOUBLE_EQ(0.0, PD[0]);
        EXPECT_NEAR(0.0, PD[z.size()-1], 1.0e-10);
    }
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
