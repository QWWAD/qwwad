/**
 * \file   file-io.h
 * \brief  Functions for reading and writing data from standard input 
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \todo   Use C++11 variadic templates to make this much cleaner and more generic!
 * \todo   Replace string tokenisation with C++ istringstream handling
 */

#ifndef QWWAD_FILE_IO_H
#define QWWAD_FILE_IO_H

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstring>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <iostream>

namespace QWWAD
{
/** Exception that occurs when a file contains wrong number of lines */
class FileLinesNotAsExpected : public std::exception
{
public:
FileLinesNotAsExpected(const std::string &fname,
                       const size_t       nexpected,
                       const size_t       nread);

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
template <template<typename, typename...> class Tcontainer,
          class T>
int read_line_array(Tcontainer<T> &dest, const size_t n, std::istream& stream)
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
template <template<typename, typename...> class Tcontainer,
          class T>
void read_line_array_u(Tcontainer<T>& dest, std::istream& stream)
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
        dest_tmp.push_back(atof(pch)); // Copy data to array
        pch = strtok(NULL, "\t "); // Try to read next data from line
    } 

    delete[] linebuffer;

    dest.resize(dest_tmp.size());
    std::copy(dest_tmp.begin(), dest_tmp.end(), &dest[0]);
}

void parse_items(std::istream &stream);

/**
 * \brief Recursively read all the data items in a stream into destination variables
 *
 * \param[in]  stream    The stream from which to read
 * \param[out] destx     The destination for the next data item
 * \param[out] remainder A set of destinations for all remaining data items
 *
 * \details You probably don't want to call this directly.  Use read_line
 *          where possible to get data from a single line
 */
template <class Tnext, class... Tremainder>
void parse_items(std::istream  &stream,
                 Tnext         &destx,
                 Tremainder    &...remainder)
{
    if(!stream)
    {
        throw std::runtime_error("Could not read stream");
    }

    // Try to read a single item
    if(!(stream >> destx))
    {
        std::cout << destx << std::endl;
        throw std::runtime_error("Could not read item");
    }

    // Recursively read the remaining items
    parse_items(stream, remainder...);
}

/**
 * \brief Read data items from a line in a stream
 *
 * \param[in]  stream       The stream from which to read
 * \param[out] destinations A set of locations to store the data items
 *
 * \return 0 if successful, 1 if not
 */
template <class... Targs>
int read_line(std::istream &stream,
              Targs        &...destinations)
{
    int scan_result = 1; // Flag showing whether scan successful

    if(!stream)
    {
        throw std::runtime_error("Could not read stream");
    }

    std::string linebuffer; // Buffer for line data

    if(getline(stream, linebuffer) and !linebuffer.empty())
    {
        std::istringstream oss(linebuffer);

        try
        {
            parse_items(oss, destinations...);
            scan_result = 0;
        }
        catch(std::runtime_error &e)
        {
            std::ostringstream err_ss;
            err_ss << "Data missing on line: '" << linebuffer << "'";
            throw std::runtime_error(err_ss.str());
        }
    }

    return scan_result;
}

/**
 * Read numerical data from a file containing data in a single column
 *
 * \param[in]  fname Filename from which to read data
 * \param[out] x     Value array into which data will be written
 */
template <class Tstring,
          template<typename, typename...> class Tcontainer,
          class T>
void read_table(const Tstring fname, Tcontainer<T>& x)
{
    std::ifstream stream(fname);

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << fname;
        throw std::runtime_error(oss.str());
    }

    std::vector<T> x_temp;

    while(!stream.eof())
    {
        T buffer = 0; // Buffer for input data

        // If data is valid, stick it into temp vector
        if(!read_line(stream, buffer))
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
template <class Tstring,
          template<typename, typename...> class Tcontainer,
          class T>
void write_table(const Tstring        fname,
                 const Tcontainer<T> &x,
                 const bool           with_num = false,
                 const int            precision = 12)
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
template <class Tstring,
          template<typename, typename...> class Tcontainerx,
          template<typename, typename...> class Tcontainery,
          class Tx,
          class Ty>
void read_table(const Tstring    fname,
                Tcontainerx<Tx> &x,
                Tcontainery<Ty> &y,
                const size_t     n_expected = 0)
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

    while(!stream.eof()){
        Tx buffer_x = 0; // Buffer for x input data
        Ty buffer_y = 0; // Buffer for y input data

        // If data is valid, stick it into temp vector
        try
        {
            if(!read_line(stream, buffer_x, buffer_y))
            {
                x_temp.push_back(buffer_x);
                y_temp.push_back(buffer_y);
            }
        }
        catch (std::runtime_error &e)
        {
            std::ostringstream err_st;
            err_st << "Error reading " << fname << std::endl
                   << e.what();
            throw std::runtime_error(err_st.str());
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
template<class Tstring,
         template<typename, typename...> class Tcontainerx,
         template<typename, typename...> class Tcontainery,
         class Tx,
         class Ty
        >
void write_table(const Tstring          fname,
                 const Tcontainerx<Tx> &x,
                 const Tcontainery<Ty> &y,
                 const bool             with_num = false,
                 const size_t           precision = 12)
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
        {
            stream << i+1 << "\t";
        }

        stream << std::setprecision(precision)
               << std::scientific
               << x[i] << "\t" << y[i] << std::endl;
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
template<template<typename, typename...> class Tcontainerx,
         template<typename, typename...> class Tcontainery,
         template<typename, typename...> class Tcontainerz,
         class Tx,
         class Ty,
         class Tz
        >
void read_table(const char      *fname,
                Tcontainerx<Tx> &x,
                Tcontainery<Ty> &y,
                Tcontainerz<Tz> &z)
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

    while(!stream.eof()){
        Tx buffer_x = 0; // Buffer for x input data
        Ty buffer_y = 0; // Buffer for y input data
        Tz buffer_z = 0; // Buffer for z input data

        // If data is valid, stick it into temp vector
        if(!read_line(stream, buffer_x, buffer_y, buffer_z))
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

    // Copy data into output container
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
template<class Tstring,
         template<typename, typename...> class Tcontainerx,
         template<typename, typename...> class Tcontainery,
         template<typename, typename...> class Tcontainerz,
         template<typename, typename...> class Tcontaineru,
         class Tx,
         class Ty,
         class Tz,
         class Tu>
void read_table(const Tstring    fname,
                Tcontainerx<Tx>& x,
                Tcontainery<Ty>& y,
                Tcontainerz<Tz>& z,
                Tcontaineru<Tu>& u)
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
    std::vector<Tu> u_temp;

    while(!stream.eof()){
        Tx buffer_x; // Buffer for x input data
        Ty buffer_y; // Buffer for y input data
        Tz buffer_z; // Buffer for z input data
        Tu buffer_u; // Buffer for u input data

        // If data is valid, stick it into temp vector
        if(!read_line(stream, buffer_x, buffer_y, buffer_z, buffer_u))
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
template <class Tstring,
          template<typename, typename...> class Tcontainerx,
          template<typename, typename...> class Tcontainery,
          template<typename, typename...> class Tcontainerz,
          class Tx,
          class Ty,
          class Tz>
void write_table(const Tstring          fname,
                 const Tcontainerx<Tx> &x,
                 const Tcontainery<Ty> &y,
                 const Tcontainerz<Tz> &z,
                 const bool             with_num = false)
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
template<class Tstring,
         template<typename, typename...> class Tcontainerx,
         template<typename, typename...> class Tcontainery,
         template<typename, typename...> class Tcontainerz,
         template<typename, typename...> class Tcontaineru,
         class Tx,
         class Ty,
         class Tz,
         class Tu>
void write_table(const Tstring          fname,
                 const Tcontainerx<Tx> &x,
                 const Tcontainery<Ty> &y,
                 const Tcontainerz<Tz> &z,
                 const Tcontaineru<Tu> &u,
                 const bool             with_num = false)
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
} // namespace
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
