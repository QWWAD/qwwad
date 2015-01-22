/**
 * \file qwwad-fileio.cpp
 *
 * \brief Convenience functions for reading common data from files
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#include "qclsim-constants.h"
#include "qclsim-fileio.h"
#include "qwwad-fileio.h"

/**
 * \brief Reads subband minima from file
 *
 * \param[in] p A character to represent the particle being considered
 *
 * \returns A set of subband minima for each state in the form
 *
 * \details Data is read from a file called 'E*.r', where * is the
 *          particle ID.  The data is expected in the form:
 *          COLUMN 1: State index
 *          COLUMN 2: Energy [meV]
 */
std::valarray<double> read_E(char p)
{
    char filename[9];	// filename string

    std::valarray<double> indices;
    std::valarray<double> E;

    sprintf(filename,"E%c.r",p);
    Leeds::read_table(filename, indices, E);

    E *= 1e-3*Leeds::constants::e; // convert meV->J

    return E;
}

/**
 * \brief reads subband populations from N.r
 *
 * \details The input file is expected to take the form:
 *          COLUMN 1: Subband index
 *          COLUMN 2: Population [10^10 cm^{-2}]
 */
std::valarray<double> read_populations(int n)
{
    std::valarray<unsigned int> indices(n);
    std::valarray<double>       N(n);
    Leeds::read_table("N.r", indices, N, n);

    N *= 1e+10*1e+4; // convert from units of 10^10cm^{-2}->m^{-2}

    return N;
}

/**
 * \brief Scans the file v.r and returns the maximum value of the potential
 *
 * \returns Maximum potential [J]
 */
double Vmax()
{
    std::valarray<double> z;
    std::valarray<double> V;
    Leeds::read_table("v.r", z, V);
    return V.max();
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
