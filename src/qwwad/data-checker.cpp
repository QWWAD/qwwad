/**
* \file  data-checker.h
* \brief Class for testing data integrity
*
* \author Alex Valavanis <a.valavanis@leeds.ac.uk>
*/

#include "data-checker.h"

namespace QWWAD {
/**
 * \brief Check that all data items are positive and nonzero
 */
void DataChecker::check_positive() const
{
    for (unsigned int i = 0; i < _data.size(); ++i)
    {
        if(_data(i) <= 0)
        {
            std::ostringstream oss;
            oss << "Nonpositive value (" << _data(i) << ") detected.";
            throw std::domain_error(oss.str());
        }
    }
}

/**
 * \brief Check that all data items are not negative
 */
void DataChecker::check_not_negative() const
{
    for (unsigned int i = 0; i < _data.size(); ++i)
    {
        if(_data(i) < 0)
        {
            std::ostringstream oss;
            oss << "Negative value (" << _data(i) << ") detected.";
            throw std::domain_error(oss.str());
        }
    }
}

/**
 * \brief Convenience wrapper for testing if a data set is positive
 *
 * \param[in] data The data set to check
 */
void DataChecker::check_positive(const arma::vec &data)
{
    DataChecker checker(data);
    checker.check_positive();
}

/**
 * \brief Convenience wrapper for testing if a data set is non-negative
 *
 * \param[in] data The data set to check
 */
void DataChecker::check_not_negative(const arma::vec &data)
{
    DataChecker checker(data);
    checker.check_not_negative();
}

/**
 * Checks that a property lies in the closed interval [0,1]
 *
 * \todo Move this into the class
 */
void check_c_interval_0_1(const double *px)
{
    if(*px < 0.0 or *px > 1.0)
    {
        std::ostringstream oss;
        oss << "Value (" << *px << ") lies outside the closed interval [0,1].";
        throw std::domain_error(oss.str());
    }
}

} // namespace QWWAD 
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
