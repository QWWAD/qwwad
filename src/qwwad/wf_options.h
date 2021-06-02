/**
 * \file   wf_options.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 *
 * \brief  Declaration of WfOptions class. Used to extend the Options class
 * 	   to include options for defing wanefunction and energy i/o files.
 */

#ifndef QWWAD_WF_OPTIONS_H
#define QWWAD_WF_OPTIONS_H

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include "options.h"

namespace QWWAD {

class WfOptions : public Options
{
public:
    WfOptions();

    [[nodiscard]] auto get_wf_filename(int ist) const -> std::string;
    [[nodiscard]] auto get_wf_prefix() const -> std::string;
    [[nodiscard]] auto get_wf_ext() const -> std::string;
    [[nodiscard]] auto get_energy_filename() const -> std::string;
};
} // end namespace
#endif // QWWAD_WF_OPTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
