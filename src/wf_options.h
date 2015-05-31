/**
 * \file   wf_options.h
 * \author Jonathan Cooper <jdc.tas@gmail.com>
 * \author Alex Valavanis  <a.valavanis@leeds.ac.uk>
 *
 * \brief  Declaration of WfOptions class. Used to extend the Options class
 * 	   to include options for defing wanefunction and energy i/o files.
 */

#ifndef WF_OPTIONS_H
#define WF_OPTIONS_H

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "qwwad/options.h"

namespace QWWAD {
/**
 * \brief Control flags for input, output or in/out options
 */
enum WfOptionMode
{
    WF_OPTION_MODE_IN,  ///< Only include options for reading wf data
    WF_OPTION_MODE_OUT, ///< Only include options for writing wf data
    WF_OPTION_MODE_IO   ///< Include options both for reading and writing
};

class WfOptions : public Options
{
    public:
        WfOptions(const WfOptionMode mode = WF_OPTION_MODE_IO);

        /* Return paths for input files */
        std::string get_wf_input_path(const int ist) const;
        std::string get_wf_input_prefix() const;
        std::string get_wf_input_ext() const;
        std::string get_energy_input_path() const;
        std::string get_potential_input_path() const;

        /* Return paths for output files */
        std::string get_wf_output_path(const int ist) const;
        std::string get_wf_output_prefix() const;
        std::string get_wf_output_ext() const;
        std::string get_energy_output_path() const;
        std::string get_potential_output_path() const;

private:
        WfOptionMode _mode; ///< Whether we're using input/output mode
};
} // end namespace
#endif // WFOPTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
