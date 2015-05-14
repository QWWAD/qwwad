/**
 * \file file-io-deprecated.h
 *
 * \brief Convenience functions for reading common data from files
 *
 * \author Paul Harrison  <p.harrison@shu.ac.uk>
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_FILE_IO_DEPRECATED_H
#define QWWAD_FILE_IO_DEPRECATED_H
#include <valarray>

namespace QWWAD
{
std::valarray<double> read_E(char p);
std::valarray<double> read_populations(int n);
double Vmax();
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
