/**
 * \file   qwwad-options.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Basic options for all programs in QWWAD simulation suite
 */
#ifndef QWWAD_OPTIONS_H
#define QWWAD_OPTIONS_H

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace QWWAD {
/**
 * \brief Common options for all QCLsim programs
 *
 * \details This is an abstract base class that cannot be instantiated 
 *          directly.  The intention is that the a subclass should be
 *          created, which adds any extra options needed for a program.
 *          The virtual \c print() method also needs to be defined for all
 *          subclasses.
 */
class Options
{
    private:
        /**
         * Generic options such as "--help" that are common to all programs
         * but don't have any effect on the configuration of programs.
         * These are only suitable for specification on the command-line.
         */
        po::options_description* generic_options_commandline;

        /**
         * Generic options that are common to all (or most) programs.
         * These normally relate to some sorts of global setting, like
         * precision parameters.  They can be given on the 
         * command-line, in config files or as environment variables.
         */
        po::options_description *generic_options_any;

        po::options_description config_options;

        std::string              config_filename; ///< Configuration filename
        
        void print_version_then_exit(char* prog_name) const;

        std::string name_mapper(std::string in) const;

    protected:
        /**
         * Storage for the (raw) values entered on the command-line
         */
        po::variables_map vm;

    public:
        bool get_argument_known(const std::string &name) const;

        /**
         * \brief Adds an option to the program, with a default argument specified
         *
         * \param[in] name          The name of the option ("<long form>,<short form>")
         * \param[in] default_value The default value of the option
         * \param[in] description   A short description of what the option does
         */
        template <typename T>
        void add_option(const std::string &name,
                        const T            default_value,
                        const std::string &description)
        {
            program_specific_options->add_options()
                (name.c_str(),
                 po::value<T>()->default_value(default_value),
                 description.c_str());
        }

        /**
         * \brief Adds an option to the program
         *
         * \param[in] name          The name of the option ("<long form>,<short form>")
         * \param[in] description   A short description of what the option does
         */
        template <typename T>
        void add_option(const std::string &name,
                        const std::string &description)
        {
            program_specific_options->add_options()
                (name.c_str(),
                 po::value<T>(),
                 description.c_str());
        }

        /**
         * \brief Get the value of an option
         *
         * \param[in] name The long name of the option
         *
         * \returns The value of the option
         */
        template <typename T>
        T get_option(const std::string &name) const
        {
            return vm[name].as<T>();
        }

        void add_prog_specific_options_and_parse(const int     argc,
                                                 char ** const argv,
                                                 std::string   summary);
    protected:
        /**
         * \brief The additional options for a specific program
         *
         * \details This should be defined for each subclass of 
         *          \c Options.  It should then be appended to the
         *          generic options and run through the parser using
         *          \c add_prog_specific_options_and_parse(...)
         */
        po::options_description* program_specific_options;

    public:
        Options();
        Options(const Options &options);
        Options & operator=(const Options &options);
        virtual ~Options();

        /**
         * \brief Return whether or not verbose output is desired
         *
         * \returns \c true if verbose output is wanted
         */
        bool get_verbose() const {return vm["verbose"].as<bool>();}
};

/**
 * \brief Adds a floating-point option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] description   A short description of what the option does
 *
 * \details For floating-point numbers, the value is actually stored internally as a string.
 *          This is so that the --help message doesn't display the value to arbitrary
 *          precision.
 */
template <>
inline void Options::add_option<double>(const std::string &name,
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
 *
 * \details For floating-point numbers, the value is actually stored internally as a string.
 *          This is so that the --help message doesn't display the value to arbitrary
 *          precision.
 */
template <>
inline void Options::add_option<double>(const std::string &name,
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
 * \brief Adds a switch option to the program
 *
 * \param[in] name          The name of the option ("<long form>,<short form>")
 * \param[in] description   A short description of what the option does
 */
template <>
inline void Options::add_option<bool>(const std::string &name,
                                      const std::string &description)
{
    program_specific_options->add_options()
        (name.c_str(),
         po::bool_switch()->default_value(false),
         description.c_str());
}

/**
 * \brief Get the value of a floating-point numerical option
 *
 * \param[in] name The long name of the option
 *
 * \returns The value of the option as a floating-point number.
 *
 * \details In the case of floating-point arguments, the values are actually stored
 *          as strings in the variable map so that they don't appear to arbitrary
 *          precision in the --help text.  This specialisation of the get_options()
 *          function converts the variable to a floating-point number on output.
 */
template <>
inline double Options::get_option<double>(const std::string &name) const
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
} // end namespace
#endif // QWWAD_OPTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
