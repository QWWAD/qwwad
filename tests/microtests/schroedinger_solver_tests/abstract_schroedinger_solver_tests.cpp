#include <gtest/gtest.h>

#include "testable_abstract_schroedinger_solver.h"

using namespace QWWAD;

class AbstractSchroedingerSolverTest : public ::testing::Test
{
protected:
    TestableAbstractSchroedingerSolver *schroedinger_solver;

    void SetUp()
    {
        std::valarray<double> temp; 
        schroedinger_solver = new TestableAbstractSchroedingerSolver( temp, temp );
    }

    void TearDown()
    {
        delete schroedinger_solver;
    }

    std::vector<State> make_some_solutions()
    {
        std::valarray<double> some_wavefunction(2);
        some_wavefunction[0] = 1.0;
        some_wavefunction[1] = 2.0;

        auto solutions = std::vector<State>();
        solutions.push_back( State( 1.0, some_wavefunction ) );
        solutions.push_back( State( 2.0, some_wavefunction ) );
        solutions.push_back( State( 3.0, some_wavefunction ) );
        return solutions;
    }

    void expect_solutions_equal( std::vector<State> const& expected_solutions,
                                 std::vector<State> const& actual_solutions )
    {
        ASSERT_EQ( expected_solutions.size(), actual_solutions.size() );
        for( std::vector<State>::const_iterator expected_state    = expected_solutions.begin(),
                                                end_exected_state = expected_solutions.end(),
                                                actual_state      = actual_solutions.begin(),
                                                end_actual_state  = actual_solutions.end();
             expected_state != end_exected_state && actual_state != end_actual_state;
             ++expected_state,
             ++actual_state )
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
    auto expected_solutions = make_some_solutions();
    schroedinger_solver->set_current_solutions( expected_solutions );
    auto returned_solutions = schroedinger_solver->get_solutions();
    expect_solutions_equal( expected_solutions, returned_solutions );
}
