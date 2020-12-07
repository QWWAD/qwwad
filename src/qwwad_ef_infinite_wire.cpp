/**
 * \file  qwwad_ef_infinite_wire.cpp
 * \brief Envelope Function Infinite Wire
 *
 * \details This program calculates the eigenfunctions and eigenenergies of
 *          an infinitely deep rectangular cross-section quantum wire. The 
 *          relevant parameters are passed via command line arguments.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <gsl/gsl_math.h>

#include "qwwad/file-io.h"
#include "qwwad/constants.h"
#include "qwwad/options.h"

using namespace QWWAD;
using namespace constants;

auto main(int argc, char **argv) -> int
{
    Options opt;

    std::string doc("Find eigenstates of an infinite rectangular quantum wire.");

    opt.add_option<double>("ywidth,y",   100,   "Width of quantum wire in y-direction [angstrom].");
    opt.add_option<double>("zwidth,z",   100,   "Width of quantum wire in z-direction [angstrom].");
    opt.add_option<double>("mass,m",     0.067, "Effective mass (relative to free electron).");
    opt.add_option<size_t>("nz,N",       100,   "Number of spatial points for output file.");
    opt.add_option<size_t>("nst,s",      1,     "Number of states to find.");
    opt.add_option<char>  ("particle,p", 'e',   "ID of particle to be used: 'e', 'h' or 'l', for "
            "electrons, heavy holes or light holes respectively.");

    opt.add_prog_specific_options_and_parse(argc, argv, doc);

    const auto Ly = opt.get_option<double>("ywidth") * 1e-10; // wire width in y-direction [m]
    const auto Lz = opt.get_option<double>("zwidth") * 1e-10; // wire width in z-direction [m]
    const auto p  = opt.get_option<char>("particle");         // particle ID (e, h or l)
    const auto m  = opt.get_option<double>("mass") * me;      // effective mass [kg]
    const auto N  = opt.get_option<size_t>("nz");             // number of spatial steps
    const auto s  = opt.get_option<size_t>("nst");            // number of states

    std::vector<unsigned int> y_index(s*s); // y-index of each state
    std::vector<unsigned int> z_index(s*s); // z-index of each state
    std::vector<double>       E(s*s);       // Energy of each state [meV]

    // Loop over all y and z state indices
    unsigned int ist = 0;

    for(unsigned int in_y=1;in_y<=s;in_y++)
    {
        for(unsigned int in_z=1;in_z<=s;in_z++)
        {
            y_index[ist] = in_y;
            z_index[ist] = in_z;
            E[ist]       = gsl_pow_2(pi*hBar)/(2*m)*(gsl_pow_2(in_y/Ly)+gsl_pow_2(in_z/Lz)) /(1e-3*e);

            std::vector<double> y(N*N);
            std::vector<double> z(N*N);
            std::vector<double> psi(N*N);

            const double dy = Ly/(N-1);
            const double dz = Ly/(N-1);

            unsigned int izy = 0;

            // Loop over 2D space
            for(unsigned int iy=0;iy<N;iy++)
            {
                y[izy] = (float)iy*dy;
                const double psi_y=sqrt(2/Ly)*sin(in_y*pi*y[izy]/Ly);

                for(unsigned int iz=0;iz<N;iz++)
                {
                    z[izy] = (float)iz*dz;
                    const double psi_z=sqrt(2/Lz)*sin(in_z*pi*z[izy]/Lz);
                    psi[izy] = psi_y * psi_z;

                    ++izy;
                }
            }

            std::ostringstream cd_filename;
            cd_filename << "cd" << in_y << in_z << ".r";
            write_table(cd_filename.str(), y, z, psi);

            ++ist; // Increment the overall state index
        } // end in_z
    } // end in_y

    std::ostringstream E_filename;
    E_filename << "E" << p << ".r";
    write_table(E_filename.str(), y_index, z_index, E);

    return EXIT_SUCCESS;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
