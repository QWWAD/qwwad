/**
 * \file   qclsim-linalg.h
 * \brief  Linear algebra utilities
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \date   2013-01-10
 */

#ifndef QCLSIM_LINALG_H
#define QCLSIM_LINALG_H

#if HAVE_CONFIG_H
# include "config.h"
#endif //HAVE_CONFIG_H

#include <valarray>
#include <vector>
#include <sstream>
#include "qclsim-fileio.h"
#include "qclsim-maths.h"

// If we don't have the official LAPACK C bindings, we need to declare
// the necessary external functions
//
// TODO: Kill all this once all supported platforms have LAPACK >= 3.4
#if !HAVE_LAPACKE

// Declaration of external LAPACK functions
extern "C" {
    /** Determine double-precision machine parameters */
    double dlamch_(const char* returnValue);

    /** Solve general eigenproblem */
    void dgeev_(const char *JOBVL,
                const char *JOBVR,
                const int  *N,
                double      A[],
                const int  *LDA,
                double      WR[],
                double      WI[],
                double      VL[],
                const int  *LDVL,
                double      VR[],
                const int  *LDVR,
                double      WORK[],
                int        *LWORK,
                int        *INFO);

    /** Solve symmetric-definite banded eigenproblem*/
    void dsbgvx_(const char* JOBZ,
                 const char* RANGE,
                 const char* UPLO,
                 const int* N, const int* KA, const int* KB, double AB[], const int* LDAB,
                 double BB[], const int* LDBB, double* Q, const int* LDQ, const double* VL,
                 const double* VU, const int* IL, const int* IU, const double* ABSTOL, int* M,
                 double W[], double Z[], const int* LDZ, double WORK[], int IWORK[],
                 int IFAIL[], int* INFO);

    /** Solve real symmetric tridiagonal eigenproblem */
    void dstevx_(const char* JOBZ, const char* RANGE, const int* N, double D[],
                 double E[], const double* VL, const double* VU,
                 const int* IL, const int* IU, const double* ABSTOL, int* M, 
                 double W[], double Z[], const int* LDZ, double WORK[], int IWORK[], 
                 int IFAIL[], int* INFO);

    /**
     * Solve a tridiagonal inversion problem: \f$ x := A^{-1} x\f$
     */
    void dgtsv_(const int* N, const int* NRHS, double DL[], double D[], 
                double DU[], double B[], int* LDB, int* INFO);

    /**
     * Factorise real symmetric positive definate tridagonal matrix
     */
    void dpttrf_(const int* N, double D[], double E[], int* INFO);

    /**
     * Solve real symmetric positive definite tridagonal matrix
     */
    void dpttrs_(const int* N, const int* NRHS, double D[], double E[],
                 double B[], const int* LDB, int* INFO);

    /**
     * Factorise general matrix
     */
    void dgetrf_(const int* M, const int* N, double A[], const int* LDA,
                 int IPIV[], int* INFO);

    /**
     * Solve general matrix
     */
    void dgetrs_(const char* TRANS, const int* N, const int* NRHS, double A[],
                 const int* LDA, int IPIV[], double B[], const int* LDB, int* INFO);
} // extern "C"
#endif // !HAVE_LAPACKE

extern "C" {
    /**
     * Tridiagonal matrix multiplication: \f$B := \alpha A X + \beta B\f$
     */
    void dlagtm_(const char* TRANS, int* N, const int* NRHS, double* ALPHA,
                 double DL[], double D[], double DU[], double X[], int* LDX,
                 double* BETA, double B[], int* LDB);
}

namespace Leeds
{

/**
 * \brief A 3D vector
 */
template <class T>
class vector3D {
public:
    T x; ///< The 'x' component
    T y; ///< The 'y' component
    T z; ///< The 'z' component

    vector3D() :
        x(0), y(0), z(0)
    {}

    vector3D(T x, T y, T z) :
        x(x),
        y(y),
        z(z)
    {}

    vector3D& operator+=(const vector3D &rhs)
    {
        x+=rhs.x;
        y+=rhs.y;
        z+=rhs.z;

        return *this;
    }

    vector3D& operator*=(const T &rhs)
    {
        x*=rhs;
        y*=rhs;
        z*=rhs;

        return *this;
    }
};

template <class T>
inline vector3D<T> operator+(vector3D<T> lhs, const vector3D<T> &rhs)
{
    lhs += rhs;
    return lhs;
}

template <class T>
inline vector3D<T> operator*(vector3D<T> lhs, const T &rhs)
{
    lhs *= rhs;
    return lhs;
}

/**
 * \brief A solution to an eigenvalue problem
 */
template <class T>
class EVP_solution {
    protected:
        T _E; ///< Eigenvalue
        std::valarray<T> _psi; ///< Eigenvector

    public:
        /**
         * Initialise an EVP solution by copying known eigenvalue/
         * eigenvector pair
         */
        EVP_solution(T E, std::valarray<T> psi)
            : _E(E), _psi(psi)
        {}

        EVP_solution(const size_t n) : 
            _E(0),
            _psi(0.0,n)
        {}

        /**
         * Return the value of element i of an eigenvector
         */
        T psi(const unsigned int i) const
        {
            if(i >= _psi.size())
            {
                std::ostringstream oss;
                oss << "Requested element " << i << " in eigenvector of size " << _psi.size() << ".";
                throw std::domain_error(oss.str());
            }

            return _psi[i];
        }

        /**
         * Return the entire eigenvector as an array
         */
        std::valarray<T> psi_array() const
        {
            return _psi;
        }

        size_t size() const
        {
            return _psi.size();
        }

        /** Return the eigenvalue */
        T get_E() const
        {
            return _E;
        }

        /** 
         * Scale the eigenvector by a constant factor
         *
         * \param[in] a A constant scaling factor
         */
        void scale(const T a)
        {
            _psi*=a;
        }

        /**
         * Normalise the square-integral of the eigenvector
         * 
         * \param[in] x An array of samples over which the eigenvector
         *              is mapped
         *
         * This is useful for things such as normalising wavefunctions
         * in quantum mechanics.  In this case, 'x' would be the
         * spatial grid over which the wavefunction is defined.
         */
        void normalise(const std::valarray<T>& x)
        {
            if(x.size() != _psi.size())
            {
                std::ostringstream oss;
                oss << "Cannot normalise eigenvector.  Grid contains " << x.size()
                    << " samples but eigenvector contains " << _psi.size()
                    << ".";
                throw std::runtime_error(oss.str());
            }

            const T dx = x[1] - x[0]; // Grid step size

            // Normalisation factor
            const T A = sqrt(trapz(psi_squared(), dx));

            scale(1.0/A);
        }

        /**
         * Get the squared value of the eigenvector at a given point
         *
         * \param[in] i The index of the point for which to find eigenvector squared
         */
        T psi_squared(const unsigned int i) const
        {
            if(i >= _psi.size())
            {
                std::ostringstream oss;
                oss << "Requested element " << i << " in an eigenvector of size " << _psi.size();
                throw std::domain_error(oss.str());
            }

            return psi(i)*psi(i);
        }

        /**
         * Get an array of squared eigenvector values
         */
        std::valarray<T> psi_squared() const
        {
            return _psi*_psi;
        }

        /**
         * Find the largest square-magnitude of any element in a set of eigenvectors
         */
        static T psi_squared_max(const std::vector< EVP_solution<T> > &EVP)
        {
            double PDmax = 0.0;

            // Loop through all EVP solutions
            for(typename std::vector< EVP_solution<T> >::const_iterator st = EVP.begin(); st != EVP.end(); ++st)
            {
                // If this EVP solution has highest square-magnitude so far
                // store its value
                std::valarray<T> PD = st->psi_squared();
                PDmax = GSL_MAX_DBL(PDmax, PD.max());
            }

            return PDmax;
        }

        /** 
         * \brief Read a set of eigenvalue problem solutions from file. (Eigenvalues and eigenvectors.)
         *
         * \param[in]   Eigenval_name       The name of the file which holds the eigenvalues in a single column
         * \param[in]   Eigenvect_prefix    Prefix of the files holding the eigenvectors
         * \param[in]   Eigenvect_ext       Extension of the files holding the eigenvectors
         * \param[in]   eigenvalue_scale    Value by which all eigenvalues will be divided upon
         *                                  read
         * \param[in]   ignore_first_column True if first column of eigenvalue file should be
         *                                  ignored
         * 
         * \returns  A vector containing the eigen-solutions
         *
         * \details Reads in eigenvectors and eigenvealues from files into a vector with elements of type
         *          <EVP_solution>.
         */
        static std::vector< EVP_solution<T> >
            read_from_file(const std::string &Eigenval_name,
                           const std::string &Eigenvect_prefix,
                           const std::string &Eigenvect_ext,
                           const T            eigenvalue_scale    = 1.0,
                           const bool         ignore_first_column = false)
            {
                std::vector < EVP_solution<T> > solutions;

                // Read eigenvalues into tempory memory
                std::valarray<T> E_temp;

                if(ignore_first_column)
                {
                    std::valarray<unsigned int> indices;
                    read_table_xy(Eigenval_name.c_str(), indices, E_temp);
                }
                else
                    read_table_x(Eigenval_name.c_str(), E_temp);

                E_temp /= eigenvalue_scale;

                // Set number of states
                const size_t nst = E_temp.size();
                if(nst==0)
                {
                    std::ostringstream oss;
                    oss << Eigenval_name << " appears to be empty. Is this the correct eigenvalue input file?.";
                    throw std::runtime_error(oss.str());
                }

                // Read first eigenvector into tempory memory to get size of vectors
                std::valarray<T> z_temp;
                std::valarray<T> psi_temp;
                std::string Eigenvect_name = Eigenvect_prefix + "1" + Eigenvect_ext;
                read_table_xy(Eigenvect_name.c_str(), z_temp, psi_temp);

                if(z_temp.size() == 0)
                {
                    std::ostringstream oss;
                    oss << "No data found in " << Eigenvect_name << ". Is this the correct eigenvector input file?";
                    throw std::runtime_error(oss.str());
                }

                // Resize permanent store of eigen-solutions to correct size and copy in first eigen-solution
                const size_t psi_size = z_temp.size();
                solutions.resize(nst, EVP_solution<T>(psi_size));
                solutions[0] = EVP_solution<T>(E_temp[0], psi_temp);

                // Read in remaining eigenvectors and copy into permanent store
                for(unsigned int ist=1; ist<nst; ist++){
                    std::stringstream Eigenvect_name_sstream;
                    Eigenvect_name_sstream << Eigenvect_prefix << ist+1 << Eigenvect_ext;
                    Eigenvect_name = Eigenvect_name_sstream.str();
                    read_table_xy(Eigenvect_name.c_str(), z_temp, psi_temp, psi_size);
                    solutions[ist] = EVP_solution<T>(E_temp[ist], psi_temp);
                }
                
                return solutions;
            }
        
        /** 
         * \brief Write a set of eigenvalue problem solutions to file. (Eigenvalues and eigenvectors.)
         *
         * \param[in]  Eigenval_name       The name of the file to which the eigenvalues will be written
         * \param[in]  Eigenvect_prefix    Prefix of the files holding the eigenvectors
         * \param[in]  Eigenvect_ext       Extension of the files holding the eigenvectors
         * \param[in]  solutions           Vector from which to read the eigen-solutions
         * \param[in]  z                   Array of coordinates to which eigenvectors reference
         *
         */
        template <class Tz>
            static void write_to_file(const std::string& Eigenval_name,
                    const std::string& Eigenvect_prefix,
                    const std::string& Eigenvect_ext,
                    const std::vector< EVP_solution<T> >& solutions,
                    const std::valarray<Tz>& z,
                    const bool               with_num=false)
            {
                // TODO Should the eigen-states be output with a higher precision?
                //      (Like fermi energy is...)
                // Output eigenvalues
                std::valarray<T> E_temp(solutions.size());
                for(unsigned int ist=0; ist < solutions.size(); ist++)
                    E_temp[ist] = solutions[ist].get_E();
                write_table_x(Eigenval_name.c_str(), E_temp, with_num, 17);

                // Output eigenvectors
                for(unsigned int ist=0; ist < solutions.size(); ist++){
                    std::stringstream Eigenvect_name_sstream;
                    Eigenvect_name_sstream << Eigenvect_prefix << ist+1 << Eigenvect_ext;
                    std::string Eigenvect_name = Eigenvect_name_sstream.str();
                    write_table_xy(Eigenvect_name.c_str(), z, solutions[ist].psi_array(), false, 17);
                }
            }
};

std::vector< EVP_solution<double> >
eigen_general(double       A[],
              const double VL,
              const double VU,
              int          N,
              unsigned int n_max=0);

std::vector< EVP_solution<double> >
eigen_banded(double       AB[],
             double       BB[],
             const double VL,
             const double VU,
             int          n,
             unsigned int n_max = 0);

std::vector< EVP_solution<double> >
eigen_tridiag(double       D[],
              double       E[],
              const double VL,
              const double VU,
              int          n,
              unsigned int n_max = 0);

std::valarray<double>
solve_cyclic_matrix(std::valarray<double> A_sub,
                    std::valarray<double> A_diag,
                    double                cyclic,
                    std::valarray<double> b);

void matrixProduct(double*      pB,
                   double*      pA,
                   const size_t N);

typedef EVP_solution<double> State;

double z_av(const State                 &i,
            const std::valarray<double> &z);

double mij(const State                 &i,
           const State                 &j,
           const std::valarray<double> &z);

} // namespace Leeds
#endif //QCLSIM_LINALG_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
