/**
 * \file   qwwad-options.h
 * \author Alex Valavanis <a.valavanis@leeds.ac.uk>
 * \brief  Basic options for all programs in QWWAD simulation suite
 */
#ifndef QWWAD_OPTIONS_H
#define QWWAD_OPTIONS_H

#include <boost/program_options.hpp>

namespace po = boost::program_options;

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

        std::string              config_filename; ///< Configuration filename
        
        void print_version_then_exit(char* prog_name) const;
    public:
        /**
         * Storage for the (raw) values entered on the command-line
         */
        po::variables_map vm;

        double      get_numeric_option(const std::string &name) const;
        size_t      get_size_option(const std::string &name) const;
        char        get_char_option(const std::string &name) const;
        std::string get_string_option(const std::string &name) const;
        bool        get_switch(const std::string &name) const;

        void add_numeric_option(const std::string &name,
                                const std::string &description);

        void add_numeric_option(const std::string &name,
                                const double       default_value,
                                const std::string &description);

        void add_size_option(const std::string &name,
                             const size_t       default_value,
                             const std::string &description);

        void add_char_option(const std::string &name,
                             const char         default_value,
                             const std::string &description);

        void add_string_option(const std::string &name,
                               const std::string &default_value,
                               const std::string &description);

        void add_switch(const std::string &name,
                        const std::string &description);

        // Common options for all programs
        void add_prog_specific_options_and_parse(int          argc,
                                                 char        *argv[],
                                                 std::string  summary,
                                                 std::string  details="");
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
#endif // QWWAD_OPTIONS_H
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
