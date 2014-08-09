/**
 * \file  qwwad-options.cpp
 * \brief Implementaions of common methods for program options
 */

#if HAVE_CONFIG_H
# include "config.h"
#endif

#include "qwwad-options.h"
#include <iostream>
#include <fstream>

Options::Options() :
    generic_options_commandline(new po::options_description("Generic options")),
    generic_options_any(new po::options_description("Configuration options")),
    config_filename("qwwad.cfg"),
    vm(),
    program_specific_options(new po::options_description("Program-specific options"))
{
    generic_options_commandline->add_options()
        ("help,h",    "display a help message")
        ("version,v", "display the version of the program")

        ("config,c", 
         po::value(&config_filename)->default_value("qwwad.cfg"),
         "name of the configuration file to be used")
        ;

    generic_options_any->add_options()
        ("verbose,V", po::bool_switch(),
         "display lots of information about calculation")
        ;
}

/** Copy constructor */
Options::Options(const Options &opt) :
    generic_options_commandline(new po::options_description()),
    generic_options_any(new po::options_description()),
    config_filename(opt.config_filename),
    vm(opt.vm),
    program_specific_options(new po::options_description())
{
    generic_options_commandline->add(*opt.generic_options_commandline);
    generic_options_any->add(*opt.generic_options_any);
    program_specific_options->add(*opt.program_specific_options);
}

Options & Options::operator=(const Options &opt)
{
    // Create a temporary copy of other object
    Options tmp(opt);

    // Swap data with the temporary object
    std::swap(generic_options_commandline, tmp.generic_options_commandline);
    std::swap(generic_options_any,         tmp.generic_options_any);
    std::swap(config_filename,             tmp.config_filename);
    std::swap(vm,                          tmp.vm);
    std::swap(program_specific_options,    tmp.program_specific_options);

    return *this;
}

Options::~Options()
{
    delete generic_options_any;
    delete generic_options_commandline;
    delete program_specific_options;
}

/**
 * \brief Print the version of the program to screen and then quit
 *
 * \param prog_name The name of the program
 */
void Options::print_version_then_exit(char* prog_name) const
{
    std::cout << prog_name << " " << PACKAGE_VERSION << std::endl
              << "Copyright (c) 2014 Paul Harrison and Alex Valavanis." << std::endl
              << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>." << std::endl
              << "This is free software: you are free to change and redistribute it." << std::endl
              << "There is NO WARRANTY, to the extent permitted by law." << std::endl
              << std::endl
              << "Any use of this software in published work must be accompanied by a " << std::endl
              << "citation of the textbook: \"Quantum Wells, Wires and Dots\" (4th edition), " << std::endl
              << "Paul Harrison and Alexander Valavanis, Wiley, Chichester (2015). " << std::endl
              << "in addition to any works cited in the source code." << std::endl;

    exit (EXIT_SUCCESS);
}

/**
 * \brief Add program-specific options to the list of options and then parse
 * \param[in] argc    The number of command-line arguments
 * \param[in] argv    The list of command-line arguments
 * \param[in] summary A short (1 line) documentation string describing the purpose of the program.
 *                    This is displayed at the top when the user gives the "--help" option.
 * \param[in] details More detailed notes to display below the program options in "--help" output
 *
 * \details Configuration is read in the following (decreasing) order of
 *          preference:
 *          - Command-line arguments
 *          - Configuration file items
 *          - Environment variables
 *
 *          Command-line arguments are given after the name of the program,
 *          and may use either the short form e.g.
 *          \code
 *            efiw -L5
 *          \endcode
 *          or the long form e.g.
 *          \code
 *            efiw --length=5
 *          \endcode
 *
 *          The equivalent configuration file items may be specified using the 
 *          format
 *          <tt>key = \<value\></tt>, for example
 *          \code
 *            length = 5
 *          \endcode
 *          The default configuration file is called \c qwwad.cfg but the name
 *          can be overridden using the \c --config command-line option.
 *
 *          The environment variables must be in block-capitals, and must begin 
 *          with the prefix \c QWWAD_ for example
 *          \code
 *            export QWWAD_LENGTH=5
 *          \endcode
 *
 * \todo The environment-parser will currently accept \a all known options, 
 *       some of which don't make sense.  For example,
 *       \code
 *         export QWWAD_HELP=1
 *       \endcode
 *       will force all programs to just display the help message and exit.
 *
 * \todo It would be better to use a smarter parser that ignores undesired options
 *       and unknown options.
 */
void Options::add_prog_specific_options_and_parse(int          argc,
                                                  char        *argv[],
                                                  std::string  summary,
                                                  std::string  details)
{
    try {
        // Allow all options to be given on the command-line
        po::options_description command_line_options;
        command_line_options.add(*generic_options_commandline);
        command_line_options.add(*program_specific_options);
        command_line_options.add(*generic_options_any);

        // First read everything specified on the command-line
        po::store(po::parse_command_line(argc, argv, command_line_options), vm);
        po::notify(vm);

        // Now, read from config file (if available)
        std::ifstream config_filestream(config_filename.c_str(), std::ifstream::in);

        po::options_description config_options;
        config_options.add(*generic_options_any);
        config_options.add(*program_specific_options);

        if (config_filestream)
        {
            if (get_verbose())
            {
                std::cout << "Reading configuration from " 
                          << config_filename << std::endl;
            }

            po::store(po::parse_config_file(config_filestream, config_options, true), vm);
            po::notify(vm);
        }
        else
        {
            if (get_verbose())
                std::cout << "No configuration file found" << std::endl;
        }

        // Finally, look for any suitable-looking environment variables
        po::store(po::parse_environment(config_options, "QWWAD_"), vm);
        po::notify(vm);

        // Post-processing for default options...
        if(vm.count("help"))
        {
            std::cout << summary << std::endl
                      << command_line_options << std::endl
                      << details << std::endl
                      << std::endl
                      << "Report bugs to " << PACKAGE_BUGREPORT << std::endl;

            exit(EXIT_SUCCESS);
        }
        // Display the version number and copyright notice
        if (vm.count ("version")) 
            print_version_then_exit(argv[0]);
    }
    catch(std::exception& e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * \brief Get the value of a floating-point numerical option
 *
 * \param[in] name The long name of the option
 *
 * \return The value of the option as a floating-point number
 */
double Options::get_numeric_option(const std::string &name) const
{
    const std::string val_str = vm[name].as<std::string>();
    std::istringstream iss(val_str);
    double val;

    if (!(iss >> val))
    {
        std::ostringstream oss("Can't read ");
        oss << name;
        throw oss.str();
    }
    else
        return val;
}

/**
 * \brief Adds a floating-point option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] description   A short description of what the option does
 */
void Options::add_numeric_option(const std::string &name,
                                 const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::value<std::string>(),
         description.c_str());
}

/**
 * \brief Adds a floating-point option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] default_value The default value of the option
 * \param[in] description   A short description of what the option does
 */
void Options::add_numeric_option(const std::string &name,
                                 const double       default_value,
                                 const std::string &description)
{
    std::ostringstream oss;
    oss << default_value;

    program_specific_options->add_options()
        (name.c_str(),
         po::value<std::string>()->default_value(oss.str()),
         description.c_str());
}

/**
 * \brief Get the value of a size option
 *
 * \param[in] name The long name of the option
 *
 * \return The value of the option as a positive integer
 */
size_t Options::get_size_option(const std::string &name) const
{
    const size_t val = vm[name].as<size_t>();
    return val;
}

/**
 * \brief Adds a size option to the program (i.e., a +ve integer)
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] default_value The default value of the option
 * \param[in] description   A short description of what the option does
 */
void Options::add_size_option(const std::string &name,
                              const size_t       default_value,
                              const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::value<size_t>()->default_value(default_value),
         description.c_str());
}

/**
 * \brief Get the value of a character option
 *
 * \param[in] name The long name of the option
 *
 * \return The value of the option as a character
 */
char Options::get_char_option(const std::string &name) const
{
    const char val = vm[name].as<char>();
    return val;
}

/**
 * \brief Adds a character option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] default_value The default value of the option
 * \param[in] description   A short description of what the option does
 */
void Options::add_char_option(const std::string &name,
                              const char         default_value,
                              const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::value<char>()->default_value(default_value),
         description.c_str());
}

/**
 * \brief Get the value of a string option
 *
 * \param[in] name The long name of the option
 *
 * \return The value of the option as a string
 */
std::string Options::get_string_option(const std::string &name) const
{
    const std::string val = vm[name].as<std::string>();
    return val;
}

/**
 * \brief Adds a string option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] default_value The default value of the option
 * \param[in] description   A short description of what the option does
 */
void Options::add_string_option(const std::string &name,
                                const std::string &default_value,
                                const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::value<std::string>()->default_value(default_value),
         description.c_str());
}

/**
 * \brief Adds a string option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] description   A short description of what the option does
 */
void Options::add_string_option(const std::string &name,
                                const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::value<std::string>(),
         description.c_str());
}

/**
 * \brief Get the value of a switch
 *
 * \param[in] name The long name of the option
 *
 * \return True if the option was specified
 */
bool Options::get_switch(const std::string &name) const
{
    const bool val = vm[name].as<bool>();
    return val;

}

/**
 * \brief Adds a switch option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] description   A short description of what the option does
 *
 * \details The switch is FALSE unless the option is specified
 */
void Options::add_switch(const std::string &name,
                         const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::bool_switch()->default_value(false),
         description.c_str());
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
