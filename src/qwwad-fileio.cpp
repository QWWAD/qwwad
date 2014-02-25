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
    Leeds::read_table_xy(filename, indices, E);

    E *= 1e-3*Leeds::constants::e; // convert meV->J

    return E;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
