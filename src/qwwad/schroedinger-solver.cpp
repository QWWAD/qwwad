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
 * \brief Delete all existing solutions and recalculate
 */
void
SchroedingerSolver::refresh_solutions()
{
    _solutions.clear();
    _solutions = calculate();

    // Delete any states that are out of the desired energy range
    // Ideally, sub-classes should never compute anything outside this
    // range!
    for(auto it = _solutions.begin(); it != _solutions.end(); ++it)
    {
        auto E = it->get_energy();
        if (energy_above_range(E) || energy_below_range(E)) {
            _solutions.erase(it);
        }
    }
}

/**
 * \brief Get the solutions to the Schroedinger equation.
 *
 * \details The solutions are computed on the first call to this function, but
 *          subsequent calls just recall the values and are hence much faster.
 *
 * \param[in] convert_to_meV If true, answers are returned with energy in meV rather than J.
 *
 * \returns The set of eigenstates for this Hamiltonian
 */
auto SchroedingerSolver::get_solutions(const bool convert_to_meV) -> std::vector<Eigenstate>
{
    // Only calculate if we haven't done so yet
    if(_solutions.empty()) {
        refresh_solutions();
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
 * \brief Return the potential profile array
 *
 * \returns The potential at each point [J]
 */
auto
SchroedingerSolver::get_V() const -> decltype(V_)
{
    if(V_.empty()) {
        std::cerr << "Potential profile has not been set" << std::endl;
    }

    return V_;
}

/**
 * \brief Return the position at each point
 *
 * \returns The position at each point [m]
 */
auto
SchroedingerSolver::get_z() const -> decltype(_z)
{
    if(_z.empty()) {
        std::cerr << "Position profile has not been set" << std::endl;
    }

    return _z;
}

/**
 * \brief Set the lower cut-off energy
 *
 * \param[in] E_min The new lower cut-off energy
 */
void SchroedingerSolver::set_E_min(const double E_min)
{
    // If the upper limit has been set, check that we're not
    // above it
    if(E_max_set_)
    {
        if(gsl_fcmp(E_max_, E_min, E_max_/1e6) != 1)
        {
            std::ostringstream oss;
            oss << "Desired lower cut-off energy: " << E_min * 1000/e << " meV is greater than upper cut-off energy: " << E_max_ * 1000/e << " meV.";
            throw std::domain_error(oss.str());
        }
    }

    E_min_ = E_min;
    E_min_set_ = true;
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
    if(E_min_set_)
    {
        if(gsl_fcmp(E_min_, E_max, E_min_/1e6) != -1)
        {
            std::ostringstream oss;
            oss << "Desired upper cut-off energy: " << E_max * 1000/e << " meV is less than lower cut-off energy: " << E_min_ * 1000/e << " meV.";
            throw std::domain_error(oss.str());
        }
    }

    E_max_ = E_max;
    E_max_set_ = true;
}

/**
 * \brief Return true if the specified energy lies above the upper cut-off limit for the solver
 *
 * \param[in] E The energy to test [J]
 *
 * \return true if the upper cut-off is set, and the energy is above it.
 */
auto
SchroedingerSolver::energy_above_range(const double E) const -> bool
{
    bool result = false;

    if (E_max_set_ && gsl_fcmp(E, E_max_, e*1e-12) == 1) {
        result = true;
    }

    return result;
}

/**
 * \brief Return true if the specified energy lies below the lower cut-off limit for the solver
 *
 * \param[in] E The energy to test [J]
 *
 * \return true if the lower cut-off is set and the energy is below it.
 */
auto
SchroedingerSolver::energy_below_range(const double E) const -> bool
{
    bool result = false;

    if (E_min_set_ && gsl_fcmp(E, E_min_, e*1e-12) == -1) {
        result = true;
    }

    return result;
}

/**
 * \brief Find the minimum energy in which to search for solutions
 *
 * \details This either returns the value of E_min_ if it has been
 *          set, or the minimum value of the potential profile if not.
 *          Use set_E_min() if you want to force this to a particular
 *          value.
 *
 * \return The minimum energy for the search [J]
 */
auto
SchroedingerSolver::get_E_search_min() const -> double
{
    double E_search_min;

    if(E_min_set_) E_search_min = E_min_;
    else           E_search_min = get_V().min();

    return E_search_min;
}

/**
 * \brief Find the maximum energy in which to search for solutions
 *
 * \details This either returns the value of E_max_ if it has been
 *          set, or the maximum value of the potential profile if not.
 *          Use set_E_max() if you want to force this to a particular
 *          value.
 *
 * \return The maximum energy for the search [J]
 */
auto
SchroedingerSolver::get_E_search_max() const -> double
{
    double E_search_max;

    if(E_max_set_) E_search_max = E_max_;
    else           E_search_max = get_V().max();

    return E_search_max;
}

} // namespace QWWAD
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
