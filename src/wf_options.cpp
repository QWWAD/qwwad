/**
 * \file   wf_options.cpp
 * \author Jonathan Cooper
 * \date   2012-01-31
 * \brief  Implementations of common methods for programms that access
 * 	   wavefunction and energy files.
 */

#include "wf_options.h"
#include <error.h>

WfOptions::WfOptions()
{
    program_specific_options->add_options()
        ("wf-input-prefix",
         po::value<std::string>()->default_value("wf_e"),
         "Set prefix of wavefunction input files.")

        ("wf-input-ext",
         po::value<std::string>()->default_value(".dat"),
         "Set file extension of wavefunction input files.")

        ("potential-input",
         po::value<std::string>()->default_value("Vcb.dat"),
         "Set filename of potential input file.")

        ("energy-input",
         po::value<std::string>()->default_value("Ee.dat"),
         "Set filename of energy input file.")

        ("input-dir",
         po::value<std::string>(),
         "Select directory containing input files.")

        ("wf-output-prefix",
         po::value<std::string>(),
         "Set prefix of wavefunction output files.")

        ("wf-output-ext",
         po::value<std::string>(),
         "Set file extension of wavefunction output files.")

        ("potential-output",
         po::value<std::string>(),
         "Set filename of potential output file.")

        ("energy-output",
         po::value<std::string>(),
         "Set filename of energy output file.")

        ("output-dir",
         po::value<std::string>(),
         "Select directory containing output files.")

        ("output-ext",
         po::value<std::string>()->default_value(".out"),
         "Set default output extension.")
        ;
}


/*
 * \brief Constructs the output filenames with paths according to the details below
 *
 * \details
 * 	  Ouputs new the selected energys and wfs into a new directory 
 *	  <ouput_dir> with the filenames <energy_output_file> and
 *	  <wf_out_prefix>i<wf_out_ext>; for i= 1,2,3...
 *	  If no output directory is specified then the files are
 *	  written to the current directory with the new output filenames.
 *	  If no output filenames are specified then the files are
 *	  written to the output directory with the same input filenames.
 *	  If niether are specified but the input file directory is specified
 *	  (and is not '.') then the files are written to the cwd.
 *	  If none of the above are specified then the file are written to the
 *	  cwd with the same filenames only with '.out' appened to the end.
 *
 * \TODO Add check to make sure output directory exist and if not create it!
 */
std::string WfOptions::get_wf_input_path(const int ist) const{
    if(ist<1)
        error(EXIT_FAILURE, 0, "Trying to get a wf filename with an index that is less than 1!");
    
    if(vm.count("input-dir"))
        return vm["input-dir"].as<std::string>() + vm["wf-input-prefix"].as<std::string>() +
            boost::lexical_cast<std::string>(ist) + vm["wf-input-ext"].as<std::string>();
    else
        return vm["wf-input-prefix"].as<std::string>() +
            boost::lexical_cast<std::string>(ist) + vm["wf-input-ext"].as<std::string>();
}

std::string WfOptions::get_wf_input_prefix() const{
    if(vm.count("input-dir"))
        return vm["input-dir"].as<std::string>() + vm["wf-input-prefix"].as<std::string>();
    else
        return vm["wf-input-prefix"].as<std::string>();
}

std::string WfOptions::get_energy_input_path() const{
    if(vm.count("input-dir"))
        return vm["input-dir"].as<std::string>() + vm["energy-input"].as<std::string>();
    else
        return vm["energy-input"].as<std::string>();
}

std::string WfOptions::get_potential_input_path() const{
    if(vm.count("input-dir"))
        return vm["input-dir"].as<std::string>() + vm["potential-input"].as<std::string>();
    else
        return vm["potential-input"].as<std::string>();
}

std::string WfOptions::get_wf_output_path(const int ist) const{

    std::string wf_output_path;
    std::string default_wf_output_ext = vm["output-ext"].as<std::string>();

    /* If output dir specified add to start of output path*/
    if(vm.count("output-dir")){
        wf_output_path = vm["output-dir"].as<std::string>();
        default_wf_output_ext = "";
    }
    else if(vm.count("input-dir") &&
            (vm["input-dir"].as<std::string>().compare(".")!=0 ||
             vm["input-dir"].as<std::string>().compare("./")!=0)){
        wf_output_path = "";
        default_wf_output_ext = "";
    }

    // Set wf output prefix
    if(vm.count("wf-output-prefix")) wf_output_path += vm["wf-output-prefix"].as<std::string>();
    else wf_output_path += vm["wf-input-prefix"].as<std::string>();


    // Add number into wavefunction path
    wf_output_path += boost::lexical_cast<std::string>(ist);

    // Add output extention
    if(vm.count("wf-output-ext")) wf_output_path += vm["wf-output-ext"].as<std::string>();
    else wf_output_path += vm["wf-input-ext"].as<std::string>();

    // If new paths set remove default extension
    if(vm.count("energy-output") && 
            vm.count("potential-ouput") &&
            (vm.count("wf-output-prefix") ||
             vm.count("wf-output-ext")))
        default_wf_output_ext = "";
    // Add the default extention to the paths. This will be blank unless nothing specified!
    wf_output_path += default_wf_output_ext;

    return wf_output_path;
}

std::string WfOptions::get_wf_output_prefix() const{

    std::string wf_output_prefix;

    /* If output dir specified add to start of output path*/
    if(vm.count("output-dir"))
        wf_output_prefix = vm["output-dir"].as<std::string>();
    else
        wf_output_prefix = "";

    // Set wf output prefix
    if(vm.count("wf-output-prefix")) wf_output_prefix += vm["wf-output-prefix"].as<std::string>();
    else wf_output_prefix += vm["wf-input-prefix"].as<std::string>();

    return wf_output_prefix;
}

std::string WfOptions::get_wf_output_ext() const{

    std::string wf_output_ext;
    std::string default_wf_output_ext = vm["output-ext"].as<std::string>();

    /* If output dir specified add to start of output path*/
    if(vm.count("output-dir")){
        default_wf_output_ext = "";
    }
    else if(vm.count("input-dir") &&
            (vm["input-dir"].as<std::string>().compare(".")!=0 ||
             vm["input-dir"].as<std::string>().compare("./")!=0)){
        default_wf_output_ext = "";
    }

    // Add output extention
    if(vm.count("wf-output-ext")) wf_output_ext += vm["wf-output-ext"].as<std::string>();
    else wf_output_ext += vm["wf-input-ext"].as<std::string>();

    // If new paths set remove default extension
    if(vm.count("energy-output") && 
            vm.count("potential-ouput") &&
            (vm.count("wf-output-prefix") ||
             vm.count("wf-output-ext")))
        default_wf_output_ext = "";
    // Add the default extention to the paths. This will be blank unless nothing specified!
    wf_output_ext += default_wf_output_ext;

    return wf_output_ext;
}

std::string WfOptions::get_energy_output_path() const{

    std::string energy_output_path;
    std::string default_energy_output_ext = vm["output-ext"].as<std::string>();

    // If output dir specified add to start of output path
    if(vm.count("output-dir")){
        energy_output_path = vm["output-dir"].as<std::string>();
        default_energy_output_ext = "";
    }
    else if(vm.count("input-dir") &&
            (vm["input-dir"].as<std::string>().compare(".")!=0 ||
             vm["input-dir"].as<std::string>().compare("./")!=0)){
        energy_output_path = "";
        default_energy_output_ext = "";
    }

    // Set energy output file if specified, else use input filename
    if(vm.count("energy-output")) energy_output_path += vm["energy-output"].as<std::string>();
    else energy_output_path += vm["energy-input"].as<std::string>();

    // If new paths set remove default extension
    if(vm.count("energy-output") && 
            vm.count("potential-ouput") &&
            (vm.count("wf-output-prefix") ||
             vm.count("wf-output-ext")))
        default_energy_output_ext = "";
    // Add the default extention to the paths. This will be blank unless nothing specified!
    energy_output_path += default_energy_output_ext;

    return energy_output_path;
}

std::string WfOptions::get_potential_output_path() const{

    std::string potential_output_path;
    std::string default_potential_output_ext = vm["output-ext"].as<std::string>();

    // If output dir specified add to start of output path
    if(vm.count("output-dir")){
        potential_output_path = vm["output-dir"].as<std::string>();
        default_potential_output_ext = "";
    }
    else if(vm.count("input-dir") &&
            (vm["input-dir"].as<std::string>().compare(".")!=0 ||
             vm["input-dir"].as<std::string>().compare("./")!=0)){
        default_potential_output_ext = "";
    }

    // Set potential output file if specified. else use input filename
    if(vm.count("potential-output")) potential_output_path += vm["potential-output"].as<std::string>();
    else potential_output_path += vm["potential-input"].as<std::string>();

    // If new paths set remove default extension
    if(vm.count("energy-output") && 
            vm.count("potential-ouput") &&
            (vm.count("wf-output-prefix") ||
             vm.count("wf-output-ext")))
        default_potential_output_ext = "";
    // Add the default extention to the paths. This will be blank unless nothing specified!
    potential_output_path += default_potential_output_ext;

    return potential_output_path;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
