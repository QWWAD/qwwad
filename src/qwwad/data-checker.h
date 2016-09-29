/**
 * \file  data-checker.h
 * \brief Class for testing data integrity
 *
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */

#ifndef QWWAD_DATA_CHECK_H
#define QWWAD_DATA_CHECK_H

#include <armadillo>

namespace QWWAD {

class DataChecker
{
private:
    arma::vec _data; // The data to test

public:
    DataChecker(const arma::vec &data)
        : _data(data)
    {}

    void check_positive() const;
    void check_not_negative() const;

    static void check_positive(const arma::vec &data);
    static void check_not_negative(const arma::vec &data);
};

void check_c_interval_0_1(double*);
} // namespace QWWAD 
#endif
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
