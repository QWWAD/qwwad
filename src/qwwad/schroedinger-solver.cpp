/**
 *  \file     schroedinger-solver.cpp
 *  \author   Jonathan Cooper 
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions for two-dimensional systems
 */

#include "schroedinger-solver.h"

#include <sstream>
#include <stdexcept>
#include <utility>

#include "constants.h"

namespace QWWAD
{
using namespace constants;

/**
 * Get the solutions to the Schroedinger equation.
 *
 * \details The solutions are computed on the first call to this function, but
 *          subsequent calls just recall the values and are hence much faster.
 */
auto SchroedingerSolver::get_solutions(const bool convert_to_meV) -> std::vector<Eigenstate>
{
    // Only calculate if we haven't done so yet
    if(_solutions.empty())
    {
        calculate();

        // Delete any states that are out of the desired energy range
        // Ideally, sub-classes should never compute anything outside this
        // range!
        for(auto it = _solutions.begin(); it != _solutions.end(); ++it)
        {
            if (_E_max_set && gsl_fcmp(it->get_energy(), _E_max, e*1e-12) == 1)
            {
                _solutions.erase(it);
            }
            else if (_E_min_set && gsl_fcmp(it->get_energy(), _E_min, e*1e-12) == -1)
            {
                _solutions.erase(it);
            }
        }
    }

    if(convert_to_meV)
    {
        std::vector<Eigenstate> sol_meV;

        for(auto sol_J : _solutions)
        {
            const auto E   = sol_J.get_energy();
            const auto z   = sol_J.get_position_samples();
            const auto psi = sol_J.get_wavefunction_samples();

            sol_meV.emplace_back(E*1000/e, z, psi);
        }

        return sol_meV;
    }
    else
        return _solutions;
}

/**
 * Assign arrays of potential and position
 *
 * \param[in] nst_max Maximum number of states to find
 *
 * \details If nst_max=0 (the default), all states will be found
 *          that lie within the range of the input potential profile
 */
SchroedingerSolver::SchroedingerSolver(decltype (_V)             V,
                                       decltype (_z)             z,
                                       const decltype(_nst_max)  nst_max) :
    _V(std::move(V)),
    _z(std::move(z)),
    _nst_max(nst_max),
    _E_min(0.0),
    _E_max(0.0),
    _E_min_set(false),
    _E_max_set(false),
    _solutions()
{}

/**
 * \brief Set the lower cut-off energy
 *
 * \param[in] E_min The new lower cut-off energy
 */
void SchroedingerSolver::set_E_min(const double E_min)
{
    // If the upper limit has been set, check that we're not
    // above it
    if(_E_max_set)
    {
        if(gsl_fcmp(_E_max, E_min, _E_max/1e6) != 1)
        {
            std::ostringstream oss;
            oss << "Desired lower cut-off energy: " << E_min * 1000/e << " meV is greater than upper cut-off energy: " << _E_max * 1000/e << " meV.";
            throw std::domain_error(oss.str());
        }
    }

    _E_min = E_min;
    _E_min_set = true;
}

/**
 * \brief Set the upper cut-off energy
 *
 * \param[in] E_max The new upper cut-off energy
 */
void SchroedingerSolver::set_E_max(const double E_max)
{
    // If the lower limit has been set, check that we're not
    // below that
    if(_E_min_set)
    {
        if(gsl_fcmp(_E_min, E_max, _E_min/1e6) != -1)
        {
            std::ostringstream oss;
            oss << "Desired upper cut-off energy: " << E_max * 1000/e << " meV is less than lower cut-off energy: " << _E_min * 1000/e << " meV.";
            throw std::domain_error(oss.str());
        }
    }

    _E_max = E_max;
    _E_max_set = true;
}
} // namespace QWWAD
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
