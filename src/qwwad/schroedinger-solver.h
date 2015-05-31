/**
 * \file   schroedinger-solver.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver
 */

#ifndef QWWAD_SCHROEDINGER_SOLVER_H
#define QWWAD_SCHROEDINGER_SOLVER_H

#include "eigenstate.h"

namespace QWWAD
{
/**
 * Abstract base class for any Schroedinger-equation solver
 *
 * \todo Cache the solutions and only calculate on first request
 */
class SchroedingerSolver
{
public:
    SchroedingerSolver(const std::valarray<double> &V,
                       const std::valarray<double> &z,
                       const unsigned int           nst_max=0);

    std::vector<Eigenstate> get_solutions(const bool convert_to_meV=false);

    /**
     * \returns the array of spatial positions [m]
     */
    std::valarray<double> get_z() const {return _z;}

    /**
     * \returns the potential profile [J]
     */
    std::valarray<double> get_V() const {return _V;}

    virtual std::string get_name() = 0;
    virtual ~SchroedingerSolver() {};

    void set_E_cutoff(const double E);

    /**
     * \brief Turn off filtering of solutions by energy
     */
    inline void unset_E_cutoff()
    {
        _E_cutoff_set = false;
    }

protected:
    virtual void calculate() = 0;

    std::valarray<double> _V; ///< Confining potential [J]
    std::valarray<double> _z; ///< Spatial points [m]
    unsigned int    _nst_max; ///< Maximum number of states to find

    // Options for specifying cut-off energy
    double _E_cutoff;     ///< Cut-off energy for solutions
    bool   _E_cutoff_set; ///< True if a cut-off energy has been set

    ///< Set of solutions to the Schroedinger equation
    std::vector<Eigenstate> _solutions;
};
} // namespace QWWAD

#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
