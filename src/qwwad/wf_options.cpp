/**
 * \file   wf_options.cpp
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 * \brief  Implementations of common methods for programmes that access
 * 	   wavefunction and energy files.
 */

#include "wf_options.h"

namespace QWWAD {
/**
 * \brief Default constructor
 */
WfOptions::WfOptions()
{
    add_option<std::string>("wffileprefix","wf_e", "Prefix of wavefunction filenames.");
    add_option<std::string>("wffileext",   ".r",   "File extension for wavefunction files.");
    add_option<std::string>("energyfile",  "Ee.r", "Filename of energy file.");
}

/*
 * \brief Constructs the filename for a wavefunction with a given number
 *
 * \param[in] i The index of the state
 *
 * \return The filename
 *
 * \details The filename is in the form: <wf_out_prefix>i<wf_out_ext>,
 *          where i is the index of the state
 */
auto WfOptions::get_wf_filename(const int ist) const -> std::string
{
    if(ist<1) {
        throw std::runtime_error("Trying to get a wf filename with an index that is less than 1!");
    }

    const auto prefix = get_wf_prefix();
    const auto ext    = get_wf_ext();

    std::ostringstream oss;
    oss << prefix << ist << ext;

    return oss.str();
}

auto WfOptions::get_wf_prefix() const -> std::string
{
    const auto prefix = get_option<std::string>("wffileprefix");
    return prefix;
}
        
auto WfOptions::get_wf_ext() const -> std::string
{
    const auto ext = get_option<std::string>("wffileext");
    return ext;
}

auto WfOptions::get_energy_filename() const -> std::string{
    const auto filename = get_option<std::string>("energyfile");
    return filename;
}
} // end namespace
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
