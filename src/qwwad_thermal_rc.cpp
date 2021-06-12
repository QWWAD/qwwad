#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <valarray>

#include "qwwad/options.h"
#include "qwwad/file-io.h"
#include "qwwad/linear-algebra.h"

#include <lapack.h>

using namespace QWWAD;

class ThermalRCOptions: public Options
{
    double dc = 0.02; // Duty cycle
    double f  = 10e3; // Pulse repetition rate [Hz]

    public:
    ThermalRCOptions(int argc, char** argv);

    /// Return fractional duty cycle (i.e. 0 to 1)
    [[nodiscard]] auto get_duty_cycle() const -> double {return dc;}

    /// Return pulse repetition rate [Hz]
    [[nodiscard]] auto get_f_rep() const -> double {return f;}

    void print() const {}
};


// Define and parse all user options and return them in a structure
ThermalRCOptions::ThermalRCOptions(int argc, char **argv)
{
    add_option<double>("Tsink,T",       80.0,   "set heatsink temperature [K]");
    add_option<double>("dc,d",           2.9,   "set duty cycle for pulse train [%]");
    add_option<double>("frequency,f",   10.0,   "set pulse repetition rate [kHz]");
    add_option<double>("power,P",       17.65,  "set pulse power [W]");
    add_option<size_t>("nrep",           1,     "set number of pulse periods to simulate");
    add_option<double>("resistance,R",  18.0,   "set thermal resistance [K/W]");
    add_option<double>("capacitance,C", 1.0e-6, "set thermal capacitance [J/K]");

    std::string doc = "Calculate temperature variation in active region over time";
    add_prog_specific_options_and_parse(argc,argv,doc);

    // Check that heatsink temperature is positive
    auto Tsink = get_option<double>("Tsink");
    if (Tsink <= 0.0)
    {
        std::ostringstream oss;
        oss << "Heatsink temperature, "
            << Tsink
            << " is not positive.";
        throw std::domain_error(oss.str());
    }

    // Check that duty cycle is positive and
    // rescale to a decimal value
    if(get_argument_known("dc"))
    {
        dc = get_option<double>("dc") * 0.01;

        if(dc <= 0.0 or dc >= 1.0)
        {
            std::ostringstream oss;
            oss << "Specified duty cycle, " << dc << ", is invalid.";
            throw std::domain_error(oss.str());
        }
    }

    // Check that frequency is positive and
    // rescale to Hz
    if(get_argument_known("frequency"))
    {
        f = get_option<double>("frequency") * 1.0e3;

        if(f <= 0) {
            throw std::domain_error ("Pulse repetition rate must "
                    "be positive.");
        }
    }

    auto power = get_option<double>("power");

    // Check that power is positive
    if(power <= 0.0)
    {
        std::ostringstream oss;
        oss << "Electrical power dissipation, "
            << power
            << "is not positive.";
        throw std::domain_error(oss.str());
    }

    auto resistance = get_option<double>("resistance");

    if(resistance <= 0.0)
    {
        std::ostringstream oss;
        oss << "Thermal resistance, "
            << resistance
            << "is not positive.";
        throw std::domain_error(oss.str());
    }

    auto capacitance = get_option<double>("capacitance");

    if(capacitance <= 0.0)
    {
        std::ostringstream oss;
        oss << "Thermal capacitance, "
            << capacitance
            << "is not positive.";
        throw std::domain_error(oss.str());
    }
}

auto main(int argc, char *argv[]) -> int
{
    // Grab user preferences
    ThermalRCOptions opt(argc, argv);

    double dt_max=1.0/(1000*opt.get_f_rep());
    double pw = opt.get_duty_cycle()/opt.get_f_rep();
    
    double _f_rep = opt.get_f_rep(); // Pulse repetition rate [Hz]
    double t_period = 1.0/_f_rep; // Length of a period [s]
    size_t nt_per = ceil(t_period/dt_max); // Divide period into steps
    double dt = t_period/float(nt_per); // Time-increment to use
    
    if(opt.get_verbose()) {
        printf("dt=%.4f ns.\n",dt*1e9);
    }

    auto _n_rep = opt.get_option<size_t>("nrep"); // Number of pulse repetitions
    unsigned int i=0;

    auto R = opt.get_option<double>("resistance"); // K/W
    auto C = opt.get_option<double>("capacitance");

    std::vector<double> _q;
    std::vector<double> _t;

    // Indices of the time samples in the middle of each pulse
    std::valarray<unsigned int> _it_mid(_n_rep);
    std::valarray<double> t_mid(_n_rep); // Time at middle of each pulse
    std::valarray<double> T_mid(_n_rep); // Temperature at middle of each pulse
   
    /**
     * Electrical power dissipation while QCL is switched on [W]
     * It is assumed that all power is dissipated in the active
     * region heterostructure (i.e. contacts have zero resistance!)
     */
    auto power = opt.get_option<double>("power");

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
	    if (dt*it <= pw) {
		    _q.push_back(power);
            } else {
                _q.push_back(0);
            }

	    // Find the middle of the pulse
	    if(dt*it <= pw/2.0) {
                _it_mid[iper] = i;
            }

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

    int INFO = 1;
    int N    = nt;
    int NRHS = 1;
    dgtsv_(&N, &NRHS, &DL[0], &D[0], &DU[0], &B[0], &N, &INFO);

    std::valarray<double> t(&_t[0], nt);
    std::valarray<double> q(&_q[0], nt);
    auto T_H = opt.get_option<double>("Tsink"); // Heat-sink temperature [K]
    std::valarray<double> T = B + T_H;

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
