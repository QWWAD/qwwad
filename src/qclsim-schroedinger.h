/**
 * \file   qclsim-schroedinger.h
 * \author Jonathan Cooper
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Declarations for Schroedinger solver
 */

#ifndef QCLSIM_SCHROEDINGER_H
#define QCLSIM_SCHROEDINGER_H

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qclsim-linalg.h"

namespace Leeds
{

/** 
 * \brief The type of effective mass description to use in simulations.
 */
enum EffectiveMassType {
    CONSTANT_ME, /**< Constant effective mass across entire structure */
    VARIABLE_ME, /**< Spatially-varying effective mass through structure */
    NONPARABOLIC, /**< Energy-dependent effective mass (nonparabolic dispersion) */ 

    /**
     * \brief	Approximate method for solving nonparabolic Schroedinger equation
     *
     * \details	This approximate method uses a Taylor expansion of the nonparabolic effective
     * 		mass to simplify the eigenvalue problem obtained from the Schroedinger equation.
     * 		(Details of which can be found in Alharbi, Opt. Quant. Electron., 40, 551-559 (2008).
     *
     * 		This approximation is only valid for states which are energetically close the
     * 		conduction band edge. For states whose energy above the conduction band edge is
     * 		comparable with the band gap the effective mass is overestimated causing the energy of
     * 		the state to be lower than expected. Furthermore, the effective mass is increasingly
     * 		overestimated with increasing energy to a point where, at an energy approximately
     * 		one band gap above the conduction band edge, many states bunch together and become
     * 		degenerate.
     *
     * 		Users that are simulating deep quantum wells, where the conduction band offset is
     * 		comparable to that of the band gap, are advised to use the 'nonparabolic' option with
     * 		'fwf' in order to stop this situation occurring.
     */
    TAYLOR
};

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

    std::vector<State> get_solutions(const bool convert_to_meV=false);

    /**
     * \returns the array of spatial positions [m]
     */
    std::valarray<double> get_z() const {return _z;}

    virtual std::string get_name() = 0;
    virtual ~SchroedingerSolver() {};

protected:
    virtual void calculate() = 0;

    std::valarray<double> _V; ///< Confining potential [J]
    std::valarray<double> _z; ///< Spatial points [m]
    unsigned int    _nst_max; ///< Maximum number of states to find

    ///< Set of solutions to the Schroedinger equation
    std::vector<State> _solutions;
};

/**
 * Schroedinger solver that uses a full generalised matrix
 */
class SchroedingerSolverFull : public SchroedingerSolver
{
public:
    SchroedingerSolverFull(const std::valarray<double>& me,
                           const std::valarray<double>& alpha,
                           const std::valarray<double>& V,
                           const std::valarray<double>& z,
                           const unsigned int           nst_max=0);

    std::string get_name() {return "full";}

private:
    std::valarray<double> A;
    void calculate();
};

/**
 * Schroedinger equation solver (using Taylor approximation)
 */
class SchroedingerSolverTaylor : public SchroedingerSolver
{
public:
    SchroedingerSolverTaylor(const std::valarray<double> &me,
                             const std::valarray<double> &alpha,
                             const std::valarray<double> &V,
                             const std::valarray<double> &z,
                             const unsigned int           nst_max=0);
    
    std::string get_name() {return "Taylor";}
private:
    void calculate();
    std::valarray<double> AB; ///< Upper triangle of Hamiltonian matrix
    std::valarray<double> BB; ///< Lower triangle of Hamiltonian matrix
};

/**
 * Solver for Schroedinger's equation using a tridiagonal Hamiltonian matrix
 */
class SchroedingerSolverTridiag : public SchroedingerSolver
{
public:
    SchroedingerSolverTridiag(const std::valarray<double>& me,
                              const std::valarray<double>& V,
                              const std::valarray<double>& z,
                              const unsigned int           nst_max=0);

    std::string get_name() {return "tridiagonal";}
private:
    void calculate();
    std::valarray<double> diag; ///< Diagonal elements of matrix
    std::valarray<double> sub; ///< Sub-diagonal elements of matrix
};
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
