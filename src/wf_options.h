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

#include "qwwad-options.h"

class WfOptions : public Options
{
    public:
        /* Constructor prototype */
        WfOptions();

        /* Return paths for input files */
        std::string get_wf_input_path(const int ist) const;
        std::string get_wf_input_prefix() const;
        std::string get_wf_input_ext() const{ return vm["wf-input-ext"].as<std::string>(); }
        std::string get_energy_input_path() const;
        std::string get_potential_input_path() const;

        /* Return paths for output files */
        std::string get_wf_output_path(const int ist) const;
        std::string get_wf_output_prefix() const;
        std::string get_wf_output_ext() const;
        std::string get_energy_output_path() const;
        std::string get_potential_output_path() const;
};

#endif // WFOPTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
