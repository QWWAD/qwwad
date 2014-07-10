#include <gtest/gtest.h>

#include "testable_abstract_schroedinger_solver.h"

class AbstractSchroedingerSolverTest : public ::testing::Test
{
protected:
    Leeds::TestableAbstractSchroedingerSolver * schroedinger_sovler;

    void SetUp()
    {
        std::valarray<double> temp = std::valarray<double>(); 
        schroedinger_sovler = new Leeds::TestableAbstractSchroedingerSolver( temp, temp );
    }

    void TearDown()
    {
        delete schroedinger_sovler;
    }

    std::vector<Leeds::State> make_some_solutions()
    {
        std::valarray<double> some_wavefunction = std::valarray<double>( 2 );
        some_wavefunction[0] = 1.0;
        some_wavefunction[1] = 2.0;

        std::vector<Leeds::State> solutions = std::vector<Leeds::State>();
        solutions.push_back( Leeds::State( 1.0, some_wavefunction ) );
        solutions.push_back( Leeds::State( 2.0, some_wavefunction ) );
        solutions.push_back( Leeds::State( 3.0, some_wavefunction ) );
        return solutions;
    }

    void expect_solutions_equal( std::vector<Leeds::State> const& expected_solutions,
                                 std::vector<Leeds::State> const& actual_solutions )
    {
        ASSERT_EQ( expected_solutions.size(), actual_solutions.size() );
        for( std::vector<Leeds::State>::const_iterator expected_state = expected_solutions.begin(),
                                                       end_exected_state = expected_solutions.end(),
                                                       actual_state = actual_solutions.begin(),
                                                       end_actual_state = actual_solutions.end();
             expected_state != end_exected_state && actual_state != end_actual_state;
             ++expected_state, ++actual_state )
        {
            EXPECT_DOUBLE_EQ( expected_state->get_E(), actual_state->get_E() );
            ASSERT_EQ( expected_state->size(), actual_state->size() );
            for( size_t iPsi=0; iPsi<expected_state->size(); ++iPsi )
            {
                EXPECT_DOUBLE_EQ( expected_state->psi( iPsi ), actual_state->psi( iPsi ) );
            }
        }
    }
};

TEST_F( AbstractSchroedingerSolverTest, get_solutions_when_solutions_found_and_not_converting_to_mev )
{
    std::vector<Leeds::State> expected_solutions = make_some_solutions();
    schroedinger_sovler->set_current_solutions( expected_solutions );
    std::vector<Leeds::State> returned_solutions = schroedinger_sovler->get_solutions();
    expect_solutions_equal( expected_solutions, returned_solutions );
}
