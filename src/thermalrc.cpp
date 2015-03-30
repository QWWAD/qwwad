#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <valarray>
#include "qwwad-options.h"
#include "qclsim-fileio.h"

#if HAVE_LAPACKE
#include <lapacke.h>
#endif

using namespace Leeds;

class ThermalRCOptions: public Options
{
    double dc; // Duty cycle
    double f; // Pulse repetition rate [Hz]

    public:
    ThermalRCOptions(int argc, char* argv[]);

    /** 
     * Electrical power dissipation while QCL is switched on [W]
     * It is assumed that all power is dissipated in the active
     * region heterostructure (i.e. contacts have zero resistance!)
     */
    double get_power() const {return vm["power"].as<double>();}

    /// Return fractional duty cycle (i.e. 0 to 1)
    double get_duty_cycle() const {return dc;}

    /// Return pulse repetition rate [Hz]
    double get_f_rep() const {return f;}

    /// Return heatsink temperature [K]
    double get_heatsink_temperature() const {return vm["Tsink"].as<double>();}

    /// Return number of pulses to simulate
    size_t get_n_rep() const {return vm["nrep"].as<size_t>();}

    /// Return ridge area [m^2]
    double get_area() const {return vm["area"].as<double>()*1e-6;}
    double get_R() const {return vm["resistance"].as<double>();}
    double get_C() const {return vm["capacitance"].as<double>();}

    std::string get_infile() const {return vm["infile"].as<std::string>();}

    void print() const {}
};


// Define and parse all user options and return them in a structure
ThermalRCOptions::ThermalRCOptions(int argc, char *argv[]) :
    dc(0.02),
    f(10e3)
{
    program_specific_options->add_options()
        ("area", po::value<double>()->default_value(0.85*0.14),
         "QCL ridge area [mm^2]")

        ("infile", po::value<std::string>()->default_value("thermal_layers.dat"),
         "Waveguide layers data file")

        ("Tsink,T", po::value<double>()->default_value(80.0), 
         "set heatsink temperature [K]")

        ("dc,d", po::value<double>()->default_value(2), 
         "set duty cycle for pulse train [%]")

        ("frequency,f", po::value<double>()->default_value(10), 
         "set pulse repetition rate [kHz]")

        ("power,P", po::value<double>()->default_value(17.65),
         "set pulse power [W]")

        ("nrep", po::value<size_t>()->default_value(1), 
         "set number of pulse periods to simulate")

	("resistance,R", po::value<double>()->default_value(18),
	 "set thermal resistance [K/W]")

	("capacitance,C", po::value<double>()->default_value(1e-6),
	 "set thermal capacitance [J/K]")
        ;

    std::string doc = "Calculate temperature variation in active region over time";
    add_prog_specific_options_and_parse(argc,argv,doc);

    // Check that heatsink temperature is positive
    if (vm["Tsink"].as<double>() <= 0.0)
    {
        std::ostringstream oss;
        oss << "Heatsink temperature, "
            << vm["Tsink"].as<double>()
            << " is not positive.";
        throw std::domain_error(oss.str());
    }

    // Check that duty cycle is positive and
    // rescale to a decimal value
    if(vm.count("dc"))
    {
        dc = vm["dc"].as<double>() * 0.01;

        if(dc <= 0.0 or dc >= 1.0)
        {
            std::ostringstream oss;
            oss << "Specified duty cycle, " << dc << ", is invalid.";
            throw std::domain_error(oss.str());
        }
    }

    // Check that frequency is positive and
    // rescale to Hz
    if(vm.count("frequency"))
    {
        f = vm["frequency"].as<double>() * 1.0e3;

        if(f <= 0)
            throw std::domain_error ("Pulse repetition rate must "
                    "be positive.");
    }

    // Check that power is positive
    if(vm["power"].as<double>() <= 0.0)
    {
        std::ostringstream oss;
        oss << "Electrical power dissipation, "
            << vm["power"].as<double>()
            << "is not positive.";
        throw std::domain_error(oss.str());
    }

    if(vm["area"].as<double>() <= 0.0)
    {
        std::ostringstream oss;
        oss << "Ridge area, "
            << vm["area"].as<double>()
            << "is not positive.";
        throw std::domain_error(oss.str());
    }
    
    if(vm["resistance"].as<double>() <= 0.0)
    {
        std::ostringstream oss;
        oss << "Thermal resistance, "
            << vm["resistance"].as<double>()
            << "is not positive.";
        throw std::domain_error(oss.str());
    }
    
    if(vm["capacitance"].as<double>() <= 0.0)
    {
        std::ostringstream oss;
        oss << "Thermal capacitance, "
            << vm["capacitance"].as<double>()
            << "is not positive.";
        throw std::domain_error(oss.str());
    }
}

int main(int argc, char *argv[])
{
    // Grab user preferences
    ThermalRCOptions opt(argc, argv);

    double dt_max=1.0/(1000*opt.get_f_rep());
    double pw = opt.get_duty_cycle()/opt.get_f_rep();
    
    double _f_rep = opt.get_f_rep(); // Pulse repetition rate [Hz]
    double t_period = 1.0/_f_rep; // Length of a period [s]
    size_t nt_per = ceil(t_period/dt_max); // Divide period into steps
    double dt = t_period/float(nt_per); // Time-increment to use
    
    if(opt.get_verbose())
        printf("dt=%.4f ns.\n",dt*1e9);

    size_t _n_rep = opt.get_n_rep(); // Number of pulse repetitions
    unsigned int i=0;

    double R = opt.get_R(); // K/W
    double C = opt.get_C();

    std::vector<double> _q;
    std::vector<double> _t;

    // Indices of the time samples in the middle of each pulse
    std::valarray<unsigned int> _it_mid(_n_rep);
    std::valarray<double> t_mid(_n_rep); // Time at middle of each pulse
    std::valarray<double> T_mid(_n_rep); // Temperature at middle of each pulse
   
    // Generate power pulse train
    for(unsigned int iper=0; iper<_n_rep; iper++)
    {
        const double t_start = t_period*(float)iper; // Time at start of period [s]
	
	// Step through time...
        for(unsigned int it=0; it<nt_per; it++)
        {
            _t.push_back(t_start+dt*it);
            i++;

	    // If this time-step is within the pulse, then "switch on" the
	    // electrical power
	    if (dt*it <= pw)
		    _q.push_back(opt.get_power());
	    else _q.push_back(0);

	    // Find the middle of the pulse
	    if(dt*it <= pw/2.0)
		    _it_mid[iper] = i;

	}
    }

    size_t nt=_t.size();

    std::valarray<double> DL(nt-1); // Lower subdiagonal
    std::valarray<double> D(nt);    // Diagonal
    std::valarray<double> DU(nt-1); // Upper subdiagonal

    std::valarray<double> B(nt); // RHS vector

    for(unsigned int it = 0; it < nt; ++it)
    {
	    if(it != nt-1)
	    {
		    DL[it] = -C/dt;
		    DU[it] = 0;
	    }

	    D[it] = 1/R + C/dt;
	    B[it] = _q[it];
    }

#if HAVE_LAPACKE
    LAPACKE_dgtsv(LAPACK_COL_MAJOR,
                  nt,
		  1,
		  &DL[0], &D[0], &DU[0],
		  &B[0], nt);
#else
    int INFO = 1;
    int N    = nt;
    int NRHS = 1;
    dgtsv_(&N, &NRHS, &DL[0], &D[0], &DU[0], &B[0], &N, &INFO);
#endif

    std::valarray<double> t(&_t[0], nt);
    std::valarray<double> q(&_q[0], nt);
    std::valarray<double> T = B + opt.get_heatsink_temperature();

    write_table("T-t.dat", std::valarray<double>(1e6*t), T);

    for(unsigned int irep = 0; irep < _n_rep; irep++)
    {
	    t_mid[irep] = t[_it_mid[irep]];
	    T_mid[irep] = T[_it_mid[irep]];
    }

    write_table("T-mid_t.dat", t_mid, T_mid);
    write_table("q-t.dat",     t, q);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
