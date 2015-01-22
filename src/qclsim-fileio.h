/**
 * \file   qclsim-fileio.h
 * \brief  Functions for reading data from standard input 
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \todo   Use C++11 variadic templates to make this much cleaner and more generic!
 * \todo   Replace string tokenisation with C++ istringstream handling
 */

#ifndef CREADDATA_H
#define CREADDATA_H

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstring>
#include <cstdlib>
#include <valarray>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace Leeds
{

/**
 * Exception that occurs when a file contains too many lines to fit in a buffer
 */
class FileLinesExceedBufferSize : public std::exception {
public:
FileLinesExceedBufferSize(const char   *fname,
                          const size_t  buffer_size);

FileLinesExceedBufferSize(const FileLinesExceedBufferSize &other);

#if HAVE_NOEXCEPT
virtual const char* what() const noexcept;
~FileLinesExceedBufferSize () noexcept {}
#else
virtual const char* what() const throw();
~FileLinesExceedBufferSize () throw() {}
#endif

private:
std::string filename;     ///< The name of the file with too many lines
size_t      _buffer_size; ///< The buffer size that has been exceeded
};

/** Exception that occurs when a file contains wrong number of lines */
class FileLinesNotAsExpected : public std::exception
{
public:
FileLinesNotAsExpected(const char   *fname,
                       const size_t  nexpected,
                       const size_t  nread);

FileLinesNotAsExpected(const FileLinesNotAsExpected &other);

#if HAVE_NOEXCEPT
virtual const char * what() const noexcept;
~FileLinesNotAsExpected () noexcept {}
#else
virtual const char * what() const throw();
~FileLinesNotAsExpected () throw() {}
#endif

private:
std::string filename;        ///< The name of the file with wrong number of lines
size_t      nlines_expected; ///< Number of lines that were expected
size_t      nlines_read;     ///< Number of lines that were read
};

// TODO: Make this configurable
/// Maximum number of lines allowed in input files
const size_t nlines_max = 10000;

/** Functions to check data */
void check_c_interval_0_1(double*);
void check_positive(double*);
void check_not_negative(double*);

/**
 * Read an array of size n from a single line
 *
 * \param[out] dest   Destination array for input data
 * \param[in]  n      Number of data items to read from line
 * \param[in]  stream Stream from which to read data
 *
 * \return 0 if scan was successful, 1 if there was an error
 */
    template <class T>
int read_line_array(std::valarray<T> &dest, const size_t n, std::istream& stream)
{
    std::streamsize nbytes = 100; // Initial size of buffer
    int scan_result=1; // Flag to show whether scan was successful [1=error]

    if(!stream.good())
        throw std::runtime_error("Could not read stream");

    // Buffer for line data
    char* linebuffer = new char[nbytes+1];

    if(stream.getline(linebuffer, nbytes) && linebuffer[0] != '\0')
    {
        unsigned short i=0; // Loop counter for data items

        // Pointer to a token on the line
        char* pch=strtok(linebuffer, "\t ");

        /* Loop over all expected items on the line and read them to
         * array one by one */
        for(i=0; i<n; i++){
            if(pch == NULL)
                throw std::runtime_error("Data missing on at least one line");

            dest[i]=atof(pch); // Copy data to array
            pch=strtok(NULL, "\t "); // Read next data item
        }

        scan_result=0;
    }

    delete[] linebuffer;
    return scan_result;
}

/**
 * Read an array of unknown size from a single line
 *
 * \param[in]  stream - Stream from which to read data
 */
    template <class T>
void read_line_array_u(std::valarray<T>& dest, std::istream& stream)
{
    std::streamsize nbytes=10000; // Initial size of buffer

    std::vector<T> dest_tmp; // Temp storage for output data

    if(!stream)
        throw std::runtime_error("Could not read stream");
    
    // Buffer for line data
    char* linebuffer = new char[nbytes+1];

    // Read line from stream into input buffer
    if(!stream.getline(linebuffer, nbytes) or linebuffer[0] == '\0')
    {
        delete[] linebuffer;
        throw std::runtime_error("Blank input line detected");
    }

    // Get a pointer to first token on line
    char* pch=strtok(linebuffer, "\t ");

    while(pch != NULL)
    {
        if(dest_tmp.size() > nlines_max)
            throw std::length_error("Buffer overflow. Too many lines in input file");

        dest_tmp.push_back(atof(pch)); // Copy data to array
        pch = strtok(NULL, "\t "); // Try to read next data from line
    } 

    delete[] linebuffer;

    dest.resize(dest_tmp.size());
    std::copy(dest_tmp.begin(), dest_tmp.end(), &dest[0]);
}


/** 
 * Read a single numerical value from a line of input
 *
 * \param[out] dest    The variable to which the data will be written
 * \param[in]  stream  The stream from which to read the data
 *
 * \return "0" if read was successful, "1" if not
 */
template <class T>
int read_line(T& dest, std::istream& stream)
{
    int scan_result = 1; // Flag showing whether scan successful

    if(!stream)
        throw std::runtime_error("Could not read stream");

    std::string linebuffer; // Buffer for line data

    if(getline(stream, linebuffer) and !linebuffer.empty())
    {
        std::istringstream oss(linebuffer);

        if(oss >> dest)
            scan_result = 0; // Mark scan as successful
        else
            throw std::runtime_error("Some data missing on line");
    }

    return scan_result;
}

/** 
 * Read two data values from a line of input
 *
 * \param[out] destx   The variable to which the 1st data value will be written
 * \param[out] desty   The variable to which the 2nd data value will be written
 * \param[in]  stream  The stream from which to read the data
 *
 * \return "0" if read was successful, "1" if not
 */
template <class Tx, class Ty>
int read_line(Tx& destx, Ty& desty, std::istream& stream)
{
    int scan_result = 1; // Flag showing whether scan successful

    if(!stream)
        throw std::runtime_error("Could not read stream");

    std::string linebuffer; // Buffer for line data

    if(getline(stream, linebuffer) and !linebuffer.empty())
    {
        std::istringstream oss(linebuffer);

        if(oss >> destx >> desty)
            scan_result = 0; // Mark scan as successful
        else
            throw std::runtime_error("Some data missing on line");
    }

    return scan_result;
}

/** 
 * \brief Read 3 data values from a line of input
 *
 * \param[out] destx   The destination for the first data item
 * \param[out] desty   The destination for the second data item
 * \param[out] destz   The destination for the third data item
 * \param[in]  stream  The input stream from which to read data.
 *
 * \details The first three whitespace-delimited values on a line are stored in
 *          the output variables \c dest1, \c dest2 and \c dest3.
 *
 * \returns 0 if successful, 1 if not.
 */
template <class Tx, class Ty, class Tz>
int read_line(Tx &destx, Ty &desty, Tz &destz, std::istream& stream)
{
    int scan_result = 1; // Flag showing whether scan successful

    if(!stream)
        throw std::runtime_error("Could not read stream");

    std::string linebuffer; // Buffer for line data

    if(getline(stream, linebuffer) and !linebuffer.empty())
    {
        std::istringstream oss(linebuffer);

        if(oss >> destx >> desty >> destz)
            scan_result = 0; // Mark scan as successful
        else
            throw std::runtime_error("Some data missing on line");
    }

    return scan_result;
}

/** 
 * \brief Read 4 data values from a line of input
 *
 * \param[out] destx   The destination for the 1st data item
 * \param[out] desty   The destination for the 2nd data item
 * \param[out] destz   The destination for the 3rd data item
 * \param[out] destu   The destination for the 4th data item
 * \param[in]  stream  The input stream from which to read data.
 *
 * \details The first 4 whitespace-delimited values on a line are stored in
 *          the output variables \c destx, \c desty, \c destz and \c destu.
 *
 * \returns 0 if successful, 1 if not.
 */
template <class Tx, class Ty, class Tz, class Tu>
int read_line(Tx &destx, Ty &desty, Tz &destz, Tu &destu, std::ifstream& stream)
{
    int scan_result = 1; // Flag showing whether scan successful

    if(!stream)
        throw std::runtime_error("Could not read stream");

    std::string linebuffer; // Buffer for line data

    if(getline(stream, linebuffer) and !linebuffer.empty())
    {
        std::istringstream oss(linebuffer);

        if(oss >> destx >> desty >> destz >> destu)
            scan_result = 0; // Mark scan as successful
        else
            throw std::runtime_error("Some data missing on line");
    }

    return scan_result;
}

/**
 * Read numerical data from a file containing data in a single column
 *
 * \param[in]  fname Filename from which to read data
 * \param[out] x     Value array into which data will be written
 */
    template <class T>
void read_table(const char* fname, std::valarray<T>& x)
{
    std::ifstream stream(fname);

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    std::vector<T> x_temp;
    unsigned int nlines=0;

    while(!stream.eof())
    {
        if(nlines >= nlines_max)
            throw FileLinesExceedBufferSize(fname, nlines_max);

        T buffer = 0; // Buffer for input data

        // If data is valid, stick it into temp vector
        if(!read_line(buffer, stream))
            x_temp.push_back(buffer);
    }
    
    stream.close();	
    
    // Copy data into output array
    x.resize(x_temp.size());
    std::copy(x_temp.begin(), x_temp.end(), &x[0]);
}


/**
 * Write a single array of numerical data to a file
 *
 * \param[in] fname     Filename to which to read data
 * \param[in] x         Value array containing data
 * \param[in] with_num  Add an initial column containing the line number
 * \param[in] precision The number of decimal places to use in output
 */
    template <class T>
void write_table_x(const char* fname,
        const std::valarray<T>& x,
        const bool with_num = false,
        const int precision = 12)
{
    std::ofstream stream(fname);
    const size_t nx = x.size();

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    stream << std::setprecision(precision) << std::scientific;
    for(unsigned int i=0; i<nx; i++)
    {
        if(with_num)
            stream << i+1 << std::setprecision(precision) << std::scientific << "\t" << x[i] << std::endl;
        else
            stream << std::setprecision(precision) << std::scientific << x[i] << std::endl;
    }

    stream.close();	
}


/**
 * Read numerical data from a file containing data in two columns
 *
 * \param[in]  fname      Filename from which to read data
 * \param[out] x          Value array into which data from 1st column will be
 *                        written
 * \param[out] y          Value array into which data from 2nd column will be
 *                        written
 * \param[in]  n_expected The number of lines of data expected in the file. If
 *                        you don't know the number, just omit this parameter
 *                        or set it to zero.
 *
 * \todo At the moment, the n_expected value is just used for checking the
 *       size of the data.  It might be sensible to allow it to be used for
 *       sizing the output arrays.  Probably a bit more efficient.
 */
template <class Tx, class Ty>
void read_table(const char* fname,
                std::valarray<Tx>& x,
                std::valarray<Ty>& y,
                size_t n_expected = 0)
{
    std::ifstream stream(fname);

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    std::vector<Tx> x_temp;
    std::vector<Ty> y_temp;
    unsigned int nlines=0;

    while(!stream.eof()){
        if(nlines >= nlines_max)
            throw FileLinesExceedBufferSize(fname, nlines_max);

        Tx buffer_x = 0; // Buffer for x input data
        Ty buffer_y = 0; // Buffer for y input data

        // If data is valid, stick it into temp vector
        if(!read_line(buffer_x, buffer_y, stream))
        {
            x_temp.push_back(buffer_x);
            y_temp.push_back(buffer_y);
        }
    }

    const size_t nx = x_temp.size();
    const size_t ny = y_temp.size();
    x.resize(nx);
    y.resize(ny);

    if(nx != ny)
    {
        std::ostringstream oss;
        oss << "Columns in " << fname << " are different lengths.";
        throw std::runtime_error(oss.str());
    }

    if(n_expected != 0 and nx != n_expected)
        throw FileLinesNotAsExpected(fname, nx, n_expected);

    // Copy data into output array
    std::copy(x_temp.begin(), x_temp.end(), &x[0]);
    std::copy(y_temp.begin(), y_temp.end(), &y[0]);

    stream.close();	
}


/**
 * Write two arrays of numerical data to columns in a file
 *
 * \param[in] fname     Filename to which to read data
 * \param[in] x         Value array containing x data
 * \param[in] y         Value array containing y data
 * \param[in] with_num  Add an initial column containing the line number
 * \param[in] precision Precision with which to output floating point numbers
 */
    template <class Tx, class Ty>
void write_table_xy(const char              *fname,
                    const std::valarray<Tx> &x,
                    const std::valarray<Ty> &y,
                    const bool               with_num = false,
                    const size_t             precision = 12)
{
    std::ofstream stream(fname);
    const size_t nx = x.size();
    const size_t ny = y.size();

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    if(nx != ny)
    {
        std::ostringstream oss;
        oss << "x and y data have different sizes: nx = " << nx << ", ny = " << ny << ".";
        throw std::runtime_error(oss.str());
    }

    for(unsigned int i=0; i<nx; i++)
    {
        if(with_num)
            stream << i+1 << "\t" << std::setprecision(precision)
                            << std::scientific << x[i] << "\t" << y[i] << std::endl;
        else
            stream << std::setprecision(precision) << std::scientific << x[i] << "\t" << y[i] << std::endl;
    }

    stream.close();	
}


/**
 * Read numerical data from a file containing data in three columns
 *
 * \param[in]  fname Filename from which to read data
 * \param[out] x     Value array into which data from 1st column will be written
 * \param[out] y     Value array into which data from 2nd column will be written
 * \param[out] z     Value array into which data from 3rd column will be written
 */
    template <class Tx, class Ty, class Tz>
void read_table(const char* fname,
                std::valarray<Tx>& x,
                std::valarray<Ty>& y,
                std::valarray<Tz>& z)
{
    std::ifstream stream(fname);

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    std::vector<Tx> x_temp;
    std::vector<Ty> y_temp;
    std::vector<Tz> z_temp;
    unsigned int nlines=0;

    while(!stream.eof()){
        if(nlines >= nlines_max)
            throw FileLinesExceedBufferSize(fname, nlines_max);

        Tx buffer_x = 0; // Buffer for x input data
        Ty buffer_y = 0; // Buffer for y input data
        Tz buffer_z = 0; // Buffer for z input data

        // If data is valid, stick it into temp vector
        if(!read_line(buffer_x, buffer_y, buffer_z, stream))
        {
            x_temp.push_back(buffer_x);
            y_temp.push_back(buffer_y);
            z_temp.push_back(buffer_z);
        }
    }

    const size_t nx = x_temp.size();
    const size_t ny = y_temp.size();
    const size_t nz = z_temp.size();
    x.resize(nx);
    y.resize(ny);
    z.resize(nz);

    if(nx != ny or nx != nz or ny != nz)
    {
        std::ostringstream oss;
        oss << "Columns in " << fname << " have different lengths: nx = " << nx
            << ", ny = " << ny << ", nz = " << nz;
        throw std::runtime_error(oss.str());
    }

    // Copy data into output array
    std::copy(x_temp.begin(), x_temp.end(), &x[0]);
    std::copy(y_temp.begin(), y_temp.end(), &y[0]);
    std::copy(z_temp.begin(), z_temp.end(), &z[0]);

    stream.close();	
}


/**
 * Read numerical data from a file containing data in four columns
 *
 * \param[in]  fname Filename from which to read data
 * \param[out] x     Value array into which data from 1st column will be written
 * \param[out] y     Value array into which data from 2nd column will be written
 * \param[out] z     Value array into which data from 3rd column will be written
 * \param[out] u     Value array into which data from 4th column will be written
 */
    template <class Tx, class Ty, class Tz, class Tu>
void read_table(const char* fname,
                std::valarray<Tx>& x,
                std::valarray<Ty>& y,
                std::valarray<Tz>& z,
                std::valarray<Tu>& u)
{
    std::ifstream stream(fname);

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    std::vector<Tx> x_temp;
    std::vector<Ty> y_temp;
    std::vector<Tz> z_temp;
    std::vector<Tz> u_temp;
    unsigned int nlines=0;

    while(!stream.eof()){
        if(nlines >= nlines_max)
            throw FileLinesExceedBufferSize(fname, nlines_max);

        Tx buffer_x = 0; // Buffer for x input data
        Ty buffer_y = 0; // Buffer for y input data
        Tz buffer_z = 0; // Buffer for z input data
        Tz buffer_u = 0; // Buffer for u input data

        // If data is valid, stick it into temp vector
        if(!read_line(buffer_x, buffer_y, buffer_z, buffer_u, stream))
        {
            x_temp.push_back(buffer_x);
            y_temp.push_back(buffer_y);
            z_temp.push_back(buffer_z);
            u_temp.push_back(buffer_u);
        }
    }

    const size_t nx = x_temp.size();
    const size_t ny = y_temp.size();
    const size_t nz = z_temp.size();
    const size_t nu = u_temp.size();
    x.resize(nx);
    y.resize(ny);
    z.resize(nz);
    u.resize(nu);

    if(nx != ny || nx != nz || ny != nz || nu != nz)
    {
        std::ostringstream oss;
        oss << "Columns in " << fname << " have different lengths: nx = " << nx
            << ", ny = " << ny << ", nz = " << nz << ", nu = " << nu;
        throw std::runtime_error(oss.str());
    }

    // Copy data into output array
    std::copy(x_temp.begin(), x_temp.end(), &x[0]);
    std::copy(y_temp.begin(), y_temp.end(), &y[0]);
    std::copy(z_temp.begin(), z_temp.end(), &z[0]);
    std::copy(u_temp.begin(), u_temp.end(), &u[0]);

    stream.close();
}

/**
 * Write three arrays of numerical data to columns in a file
 *
 * \param[in] fname     Filename to which to read data
 * \param[in] x         Value array containing x data
 * \param[in] y         Value array containing y data
 * \param[in] z         Value array containing z data
 * \param[in] with_num  Add an initial column containing the line number
 */
    template <class Tx, class Ty, class Tz>
void write_table_xyz(const char* fname,
        const std::valarray<Tx>& x,
        const std::valarray<Ty>& y,
        const std::valarray<Tz>& z,
        const bool with_num = false)
{
    std::ofstream stream(fname);
    const size_t nx = x.size();
    const size_t ny = y.size();
    const size_t nz = z.size();

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    if(nx != ny or nx != nz or ny != nz)
    {
        std::ostringstream oss;
        oss << "x, y and z data have different sizes: nx = " << nx << ", ny = " << ny << ", nz = " << nz << ".";
        throw std::runtime_error(oss.str());
    }

    for(unsigned int i=0; i<nx; i++)
    {
        if(with_num)
            stream << i+1 << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << std::endl;
        else
            stream << x[i] << "\t" << y[i] << "\t" << z[i] << std::endl;
    }

    stream.close();	
}

/**
 * Write four arrays of numerical data to columns in a file
 *
 * \param[in] fname     Filename to which to read data
 * \param[in] x         Value array containing x data
 * \param[in] y         Value array containing y data
 * \param[in] z         Value array containing z data
 * \param[in] u         Value array containing u data
 * \param[in] with_num  Add an initial column containing the line number
 */
    template <class Tx, class Ty, class Tz, class Tu>
void write_table_xyzu(const char* fname,
        const std::valarray<Tx>& x,
        const std::valarray<Ty>& y,
        const std::valarray<Tz>& z,
        const std::valarray<Tu>& u,
        const bool with_num = false)
{
    std::ofstream stream(fname);
    const size_t nx = x.size();
    const size_t ny = y.size();
    const size_t nz = z.size();
    const size_t nu = u.size();

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    if(nx != ny || nx != nz || ny != nz || nu != nz)
    {
        std::ostringstream oss;
        oss << "x, y, z and u data have different sizes: nx = " << nx << ", ny = " << ny << ", nz = " << nz << ", nu = " << nu << ".";
        throw std::runtime_error(oss.str());
    }

    for(unsigned int i=0; i<nx; i++)
    {
        if(with_num)
            stream << i+1 << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << u[i] << std::endl;
        else
            stream << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << u[i] << std::endl;
    }

    stream.close();	
}

/** 
 * Read 3 numerical data items and a character string from a line of input
 */
    template <class T1, class T2, class T3>
int read_line_xyz_char(T1& dest1, T2& dest2, T3& dest3, char*& dest4, std::ifstream& stream)
{
    std::streamsize nbytes=100; // Initial size of buffer
    int scan_result = 1; // Flag showing whether scan successful

    if(!stream)
        throw std::runtime_error("Could not read stream");

    char* linebuffer = new char[nbytes+1];

    if(stream.getline(linebuffer, nbytes) and linebuffer[0] != '\0')
    {
        // Pointer to a token on the line
        char* pch=strtok(linebuffer, "\t ");

        if(pch == NULL)
        {
            delete[] linebuffer;
            throw std::runtime_error("Some data missing on at least one line");
        }

        dest1=static_cast<T1>(atof(pch));
        pch=strtok(NULL, "\t ");

        if(pch == NULL)
        {
            delete[] linebuffer;
            throw std::runtime_error("Some data missing on at least one line");
        }

        dest2=static_cast<T2>(atof(pch));
        pch=strtok(NULL, "\t ");

        if(pch == NULL)
        {
            delete[] linebuffer;
            throw std::runtime_error("Some data missing on at least one line");
        }

        dest3=static_cast<T3>(atof(pch));
        pch=strtok(NULL, "\t ");

        if(pch == NULL)
        {
            delete[] linebuffer;
            throw std::runtime_error("Some data missing on at least one line");
        }

        dest4=strdup(pch);
        scan_result = 0; // Mark scan as successful
    }

    delete[] linebuffer;
    return scan_result;
}
} // namespace Leeds
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
