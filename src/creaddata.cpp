/**
 * \file   creaddata.cpp
 * \brief  Functions for reading data from standard input
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 */

#include "creaddata.h"

namespace Leeds
{
FileLinesExceedBufferSize::FileLinesExceedBufferSize(const char   *fname,
                                                     const size_t  buffer_size) :
    filename(fname),
    _buffer_size(buffer_size)
{}

FileLinesExceedBufferSize::FileLinesExceedBufferSize(const FileLinesExceedBufferSize &other) :
    filename(other.filename),
    _buffer_size(other._buffer_size)
{}

/** Message to describe the error */
#if HAVE_NOEXCEPT
const char* FileLinesExceedBufferSize::what() const noexcept
#else
const char* FileLinesExceedBufferSize::what() const throw()
#endif
{
    std::ostringstream oss;
    oss << "Number of lines in " << filename
        << " exceeds the maximum permitted value (" << _buffer_size << ")";
    const std::string  str  = oss.str();
    const char        *cstr = str.c_str();
    return cstr;
}

FileLinesNotAsExpected::FileLinesNotAsExpected(const char   *fname,
                                               const size_t  nexpected,
                                               const size_t  nread) :
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

/** Checks that a property lies in the closed interval [0,1] */
void check_c_interval_0_1(double* px)
{
    if(*px < 0.0 or *px > 1.0)
    {
        std::ostringstream oss;
        oss << "Value (" << *px << ") lies outside the closed interval [0,1].";
        throw std::domain_error(oss.str());
    }
}

/** Checks that a property is positive and nonzero */
void check_positive(double* pW)
{
    if(*pW <= 0)
    {
        std::ostringstream oss;
        oss << "Nonpositive value (" << *pW << ") detected.";
        throw std::domain_error(oss.str());
    }
}

/** Checks that a property is not negative */
void check_not_negative(double* pW)
{
    if(*pW < 0)
    {
        std::ostringstream oss;
        oss << "Negative value (" << *pW << ") detected.";
        throw std::domain_error(oss.str());
    }
}
} // namespace Leeds
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
