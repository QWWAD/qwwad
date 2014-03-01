/**
 * \file   run-all-unit-tests.cpp
 * \brief  Invokes all unit tests for the QWWAD package
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 */
#include <gtest/gtest.h>

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
