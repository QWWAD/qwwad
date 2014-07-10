/**
 *  \file     qwwad-schroedinger.cpp
 *  \author   Jonathan Cooper 
 *  \author   Alex Valavanis <a.valavanis@leeds.ac.uk>
 *  \brief    Implementatation of Schrodinger solver functions for two-dimensional systems
 */

#include "qwwad-schroedinger.h"

#include "qclsim-constants.h"

namespace Leeds
{

using namespace constants;

/**
 * Get the solutions to the Schroedinger equation.
 *
 * \details The solutions are computed on the first call to this function, but
 *          subsequent calls just recall the values and are hence much faster.
 */
std::vector<State> SchroedingerSolver::get_solutions(const bool convert_to_meV)
{
    // Only calculate if we haven't done so yet
    if(_solutions.empty())
    {
        calculate();

        // Delete any states that are out of the desired energy range
        // Ideally, sub-classes should never compute anything outside this
        // range!
        while(_E_cutoff_set && !_solutions.empty() && gsl_fcmp(_solutions.back().get_E(), _E_cutoff, e*1e-12) == 1)
            _solutions.pop_back();

        // Normalise wavefunctions
        for(std::vector<State>::iterator ist = _solutions.begin(); ist != _solutions.end(); ++ist)
            ist->normalise(_z);
    }

    if(convert_to_meV)
    {
        std::vector<State> sol_meV;

        for(std::vector<State>::iterator sol_J = _solutions.begin(); sol_J != _solutions.end(); ++sol_J)
            sol_meV.push_back(State(sol_J->get_E()*1000/e, sol_J->psi_array()));

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
SchroedingerSolver::SchroedingerSolver(const std::valarray<double> &V,
                                       const std::valarray<double> &z,
                                       const unsigned int           nst_max) :
    _V(V),
    _z(z),
    _nst_max(nst_max),
    _E_cutoff(0.0),
    _E_cutoff_set(false),
    _solutions()
{}

/**
 * \brief Set the cut-off energy
 *
 * \param[in] E The new cut-off energy
 */
void SchroedingerSolver::set_E_cutoff(const double E)
{
    if(gsl_fcmp(_V.min(), E, _V.min()/1e6) != -1)
    {
        std::ostringstream oss;
        oss << "Invalid cut-off energy: " << E * 1000/e << " meV is lower than band-edge potential: " << _V.min() * 1000/e << " meV.";
        throw std::domain_error(oss.str());
    }

    _E_cutoff = E;
    _E_cutoff_set = true;
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
