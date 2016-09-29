/**
 * \file   file-io.cpp
 * \brief  Functions for reading and writing data from standard input
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 */

#include "file-io.h"

namespace QWWAD
{
FileLinesNotAsExpected::FileLinesNotAsExpected(const std::string &fname,
                                               const size_t       nexpected,
                                               const size_t       nread) :
    filename(fname),
    nlines_expected(nexpected),
    nlines_read(nread)
    {}

FileLinesNotAsExpected::FileLinesNotAsExpected(const FileLinesNotAsExpected &other) :
    filename(other.filename),
    nlines_expected(other.nlines_expected),
    nlines_read(other.nlines_read)
{}

#if HAVE_NOEXCEPT
const char * FileLinesNotAsExpected::what() const noexcept
#else
const char * FileLinesNotAsExpected::what() const throw()
#endif
{
    std::ostringstream oss;
    oss << filename << " contains " << nlines_read
        << " lines of data. Expected " << nlines_expected;
    const std::string  str  = oss.str();
    const char        *cstr = str.c_str();
    return cstr;


}

void parse_items(std::istream &stream)
{
    stream.clear();
}
} // namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
