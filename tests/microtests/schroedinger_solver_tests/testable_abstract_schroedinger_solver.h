#ifndef LEEDS_TESTABLE_ABSTRACT_SCHROEDINGER_SOLVER_H
#define LEEDS_TESTABLE_ABSTRACT_SCHROEDINGER_SOLVER_H

#include "qwwad/schroedinger-solver.h"

#include <string>

namespace QWWAD
{
class TestableAbstractSchroedingerSolver : public SchroedingerSolver
{
public:
    TestableAbstractSchroedingerSolver( std::valarray<double> const &V,
                                        std::valarray<double> const &z,
                                        unsigned int const           nst_max=0 ) :
        SchroedingerSolver( V, z, nst_max )
    {}

    void set_current_solutions( std::vector<State> const solutions )
    {
        _solutions = solutions;
    }


    virtual std::string get_name()
    {
        return std::string();
    }

    virtual void calculate()
    {
    }
};
} // end namespace
#endif
