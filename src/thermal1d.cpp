#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_math.h>
#include "qwwad-options.h"
#include "qclsim-material.h"
#include "material_library.h"
#include "qclsim-material-specification.h"
#include "qclsim-material-property.h"
#include "qclsim-maths.h"
#include "qclsim-fileio.h"
#include "qclsim-linalg.h"
#include <glibmm/ustring.h>

#if HAVE_LAPACKE
#include <lapacke.h>
#endif

using namespace Leeds;

/**
 *  Find the thermal conductivity of a material [W/m/K]
 *  \param[in] T   temperature [K]
 *
 *  \todo Figure out where all these values come from!
 *  \todo These values only work for a limited range of
 *        temperatures. Restrict the domain accordingly?
 */
static double thermal_cond(MaterialSpecification &mat, const double T)
{
    double k = 0.0;

    try
    {
        k = mat.get_prop_val_x("thermal-conductivity-vs-alloy");
    }
    catch(std::exception &e)
    {
        try
        {
            double k0_1  = mat.get_prop_val_0("thermal-conductivity-0K-1");
            double k0_2  = mat.get_prop_val_0("thermal-conductivity-0K-2");
            double tau_1 = mat.get_prop_val_0("thermal-conductivity-decay-index-1");
            double tau_2 = mat.get_prop_val_0("thermal-conductivity-decay-index-2");
            double k_1 = k0_1*pow(T,tau_1);
            double k_2 = k0_2*pow(T,tau_2);
            k = lin_interp(k_1,k_2,mat.alloy);
        }
        catch(std::exception &e)
        {
            try
            {
                const double k0  = mat.get_prop_val_0("thermal-conductivity-0K");
                const double tau = mat.get_prop_val_0("thermal-conductivity-decay-index");
                k = k0 * pow(T,tau);
            }
            catch(std::exception &e)
            {
                try
                {
                    const double k_hi    = mat.get_prop_val_0("thermal-conductivity-high-T");
                    const double k_decay = mat.get_prop_val_0("thermal-conductivity-inverse-T");
                    k = k_hi + k_decay/T;
                }
                catch(std::exception &e)
                {
                    k = mat.xml->get_property("thermal-conductivity-T")->get_val(T);
                }
            }
        }
    }

    return k;
}

class Thermal1DOptions: public Options
{
    double dc; // Duty cycle
    double f; // Pulse repetition rate [Hz]

public:
    Thermal1DOptions(int argc, char* argv[]);

    /// Return fractional duty cycle (i.e. 0 to 1)
    double get_duty_cycle() const {return dc;}

    /// Return pulse repetition rate [Hz]
    double get_f_rep() const {return f;}

    void print() const;
};


// Define and parse all user options and return them in a structure
Thermal1DOptions::Thermal1DOptions(int argc, char *argv[]) :
    dc(0.02),
    f(10e3)
{
    std::string doc = "Calculate temperature variation in active region over time";

    add_size_option   ("active,a",                   2, "Index of active-region layer");
    add_numeric_option("area",                   0.119, "QCL ridge area [mm^2]");
    add_string_option ("infile",  "thermal_layers.dat", "Waveguide layers data file");
    add_numeric_option("Tsink,T",                 80.0, "Heatsink temperature [K]");
    add_numeric_option("dy,y",                 1.00e-7, "Spatial resolution [m]");
    add_numeric_option("dc,d",                       2, "Duty cycle for pulse train [%]");
    add_numeric_option("frequency,f",               10, "Pulse repetition rate [kHz]");
    add_numeric_option("power,P",                17.65, "Pulse power [W]");
    add_size_option   ("nrep",                       1, "Number of pulse periods to simulate");

    add_prog_specific_options_and_parse(argc,argv,doc);

    // Check that heatsink temperature is positive
    const double Tsink = get_numeric_option("Tsink");
    if (Tsink <= 0.0)
    {
        std::ostringstream oss;
        oss << "Heatsink temperature, " << Tsink << " is not positive.";
        throw std::domain_error(oss.str());
    }

    // Check that spatial step is positive
    const double dy = get_numeric_option("dy");
    if(dy <= 0.0)
        throw std::domain_error ("Spatial resolution must be positive");

    // Check that duty cycle is positive and
    // rescale to a decimal value
    dc = get_numeric_option("dc") * 0.01;

    if(dc <= 0.0 or dc >= 1.0)
    {
        std::ostringstream oss;
        oss << "Specified duty cycle, " << dc << " is invalid.";
        throw std::domain_error(oss.str());
    }

    // Check that frequency is positive and
    // rescale to Hz
    f = get_numeric_option("frequency") * 1.0e3;

    if(f <= 0)
        throw std::domain_error ("Pulse repetition rate must "
                "be positive.");

    // Check that power is positive
    const double power = get_numeric_option("power");
    if(power <= 0.0)
    {
        std::ostringstream oss;
        oss << "Electrical power dissipation, " << power << " is not positive.";
        throw std::domain_error(oss.str());
    }

    // Check that area is positive
    const double area = get_numeric_option("area");
    if(area <= 0.0)
    {
        std::ostringstream oss;
        oss << "Ridge area, " << area << " is not positive.";
        throw std::domain_error(oss.str());
    }

    if(get_verbose()) print();
}


void Thermal1DOptions::print() const
{
    std::cout << "Heat sink temperature = " << get_numeric_option("Tsink")  << " K"           << std::endl
              << "Power                 = " << get_numeric_option("power")  << " W"           << std::endl
              << "Frequency             = " << f/1.0e3                      << " kHz"         << std::endl
              << "Duty cycle            = " << get_duty_cycle()*100         << "%"            << std::endl
              << "Period length         = " << 1e6/f                        << " microsecond" << std::endl
              << "Number of periods     = " << get_size_option("nrep")                        << std::endl
              << "Spatial resolution    = " << get_numeric_option("dy")*1e6 << " micron"      << std::endl
              << "Ridge area            = " << get_numeric_option("area")   << " mm^2"        << std::endl;
}


class Thermal1DData {
public:
    Thermal1DData(const Thermal1DOptions& opt);
    std::vector<MaterialSpecification> mat_layer; ///< Material in each layer
    std::valarray<double> d;                      ///< Layer thickness [m]
};

Thermal1DData::Thermal1DData(const Thermal1DOptions& opt) :
    mat_layer(0),
    d(0)
{
    std::vector<double> d_tmp; // Temp storage for layer thickness
    const size_t nbuff = 10000;

    std::string infile(opt.get_string_option("infile"));
    std::ifstream stream(infile.c_str());

    if(!stream.is_open())
    {
        std::ostringstream oss;
        oss << "Could not open " << infile << ".";
        throw std::runtime_error(oss.str());
    }
    
    while(!stream.eof())
    {
        // Check for buffer overflow
        if(d_tmp.size() >= nbuff)
            throw std::length_error("Buffer overflow.  Increase size of QCLSIM_BUFFER_SIZE variable");

        double thickness_buffer = 0.0;
        double alloy_buffer     = 0.0;
        double doping_buffer    = 0.0;
        char* mat_string;

        // Read four columns from input file:
        // - thickness
        // - alloy fraction
        // - doping
        // - material name
        if(!read_line_xyz_char(thickness_buffer,
                               alloy_buffer,
                               doping_buffer,
                               mat_string,
                               stream))
        {
            check_positive(&thickness_buffer);
            check_c_interval_0_1(&alloy_buffer);
            check_not_negative(&doping_buffer);
            
            const double n3d = doping_buffer * gsl_pow_3(100); // Rescale doping to 1/m^3
            mat_layer.push_back(MaterialSpecification(mat_string,
                                                      alloy_buffer,
                                                      n3d));

            free(mat_string);

            // Copy thickness to array and scale to metres
            d_tmp.push_back(thickness_buffer*1e-6);
        }
    }

    if(opt.get_verbose())
        std::cout << "Read " << d_tmp.size() << " layers from " << infile << std::endl;

    if(d_tmp.empty())
    {
        std::ostringstream oss;
        oss << "Could not read any layers from " << infile;
        throw std::runtime_error(oss.str());
    }

    d.resize(d_tmp.size());

    for(unsigned int iL = 0; iL < d_tmp.size(); ++iL)
        d[iL] = d_tmp[iL];
}

static double calctave(const std::valarray<double> &g,
                       const std::valarray<double> &T);

static void calctemp(double dt,
                     double *Told,
                     double q_old[],
                     double q_new[], 
                     std::vector<MaterialSpecification> &mat,
                     std::valarray<double>& T,
                     Thermal1DOptions& opt);

/**
 * \brief Find the specific heat capacity of a material [J/kg/K]
 *
 * \param[in] mat The material to use
 * \param[in] T   Temperature [K]
 *
 * \returns The specific heat capacity in J/kg/K.
 *
 * \details If parameters are available, then the temperature-dependent
 *          specific heat capacity is computed.  If not, then a constant
 *          value is found
 */
static double cp(MaterialSpecification &mat, const double T)
{
    double cp = 0.0;
    try
    {
       cp = mat.xml->get_property("specific-heat-capacity-T")->get_val(T);
    }
    catch (std::exception &e)
    {
       cp = mat.get_prop_val_0("specific-heat-capacity");
    }

    if (cp <= 0)
    {
        std::ostringstream oss;
        oss << "Negative specific heat capacity " << cp <<
               " was calculated for " << mat.xml->get_description() << " at T= " << T << " K";
        std::cerr << oss.str() << std::endl;
    }

    return cp;
}

int main(int argc, char *argv[])
{
    // Grab user preferences
    Thermal1DOptions opt(argc, argv);
    Thermal1DData data(opt);

    const double dy = opt.get_numeric_option("dy");
    const double L  = data.d.sum(); // Length of structure [m]
    const size_t ny = ceil(L/dy);   // Find number of points in structure

    if(opt.get_verbose())
    {
        std::cout << "ny = " << ny << std::endl;
    }

    // Thickness of active region [m]
    const unsigned int iAR    = opt.get_size_option("active");
    const double       L_AR   = data.d[iAR];
    const double       area   = opt.get_numeric_option("area")*1e-6; // [m^2]
    const double       volume = L_AR * area;  // Active region volume (L x h x w) [m^3]

    // Power density in active region [W/m^3]
    const double power_density = opt.get_numeric_option("power")/volume;
    const double pw = opt.get_duty_cycle()/opt.get_f_rep();

    if(opt.get_verbose())
    {
        printf("Power density = %5.2e W/m3.\n", power_density);
        printf("Pulse width = %5.1f ns.\n",pw*1e9);
    }

    std::vector<MaterialSpecification> mat(ny); // Material at each point
    std::valarray<double> y(ny);                // Spatial coordinates [m]
    std::valarray<double> g(ny);                // Power density profile [W/m^3]
    mat[0] = data.mat_layer[0];

    double bottom_of_layer=0;
    unsigned int iy=1;

    for(unsigned int iL=0; iL < data.d.size(); iL++){
        bottom_of_layer += data.d[iL];

        // Check that we haven't finished filling the array and that
        // the next point is still in this layer
        while(iy<ny and y[iy-1]+dy < bottom_of_layer)
        {
            y[iy]   = dy*iy;
            mat[iy] = data.mat_layer[iL];

            // Assume that only the active-region is heated
            if (iL == iAR)
                g[iy] = power_density;

            iy++;
        }
    }

    double _Tsink = opt.get_numeric_option("Tsink");
    
    // Spatial temperature profile through structure [K].
    // Assume that initially all points are in thermal equilibrium with heat sink.
    std::valarray<double> T(_Tsink, ny); 
    
    double* Told = new double[ny];

    for(unsigned int iy=0; iy<ny; iy++) // sink
        Told[iy]=_Tsink;

    FILE* FT=fopen("struct.dat","w");

    for(unsigned int iy=0; iy<ny; iy++)
        fprintf(FT,"%20.8e  %s\n", iy*dy*1e6, mat[iy].xml->get_description().c_str());

    fclose(FT);

    // A null vector, just for convenience in linear algebra
    double* null_vec = new double[ny]();

    // TODO: The Crank-Nicolson method allows quite large timesteps to be 
    // used... however, it's only going to give a sane result if the 
    // timestep is much shorter than the smallest "feature" in the time
    // period
    //double dt_max=timestep(data, T);
    double dt_max=1.0/(1000*opt.get_f_rep());

    double _f_rep = opt.get_f_rep(); // Pulse repetition rate [Hz]
    double time_period = 1.0/_f_rep; // Length of a period [s]
    size_t nt_per = ceil(time_period/dt_max); // Divide period into steps
    double dt = time_period/float(nt_per); // Time-increment to use

    if(opt.get_verbose())
        printf("dt=%.4f ns.\n",dt*1e9);

    const size_t _n_rep = opt.get_size_option("nrep"); // Number of pulse repetitions
   
    // Samples of average T_AR at each time-step 
    std::valarray<double> t(nt_per * _n_rep);
    std::valarray<double> T_avg(nt_per * _n_rep);
    
    // Samples of average T_AR in the middle of each pulse
    std::valarray<double> t_mid(_n_rep); // Time at middle of each pulse
    std::valarray<double> T_mid(_n_rep); // Average AR temperature at middle of each pulse

    // Samples of maximum T_AR in each pulse
    std::valarray<double> t_max(_n_rep); // Time at which peak temp occurs in each pulse
    std::valarray<double> T_max(_n_rep); // Peak AR temperature of each pulse

    // Samples of minimum T_AR in each pulse
    std::valarray<double> t_min(_n_rep); // Time at which minimum temp occurs in each pulse
    std::valarray<double> T_min(1e9, _n_rep); // Minimum AR temperature of each pulse

    // Samples of average T_AR across final pulse
    std::valarray<double> t_period(nt_per);
    std::valarray<double> T_period(nt_per);

    // Rising and falling edge of final pulse in temperature profile
    std::vector<double> _t_rise;
    std::vector<double> _T_rise;
    std::vector<double> _t_fall;
    std::vector<double> _T_fall;

    for(unsigned int iper=0; iper<_n_rep; iper++)
    {
        double t_start = time_period*iper; // Time at start of period [s]
        t_min[iper] = t_start;
        t_max[iper] = t_start;

        // Step through time...
        for(unsigned int it=0; it < nt_per; it++)
        {
            // Index of time_sample relative to start of pulse train
            const unsigned int it_total = it + nt_per*iper;
            t[it_total] = t_start+dt*it;

            // Heating term at this time-step and at the last
            // timestep
            double* q_old = null_vec;
            double* q_now = null_vec;

            // If this time-step is within the pulse, then
            // "switch on" the electrical power
            if (dt*it <= pw)
                q_now = &(g[0]);

            // Likewise for the previous time-step
            if (it > 0 and dt*(it-1) <= pw)
                q_old = &(g[0]);

            // Calculate the spatial temperature profile at this 
            // timestep
            calctemp(dt, Told, q_old, q_now, mat, T, opt);

            // Find spatial average of T_AR
            T_avg[it_total] = calctave(g, T);
            
            // Find T_AR at middle of the pulse
	    if(dt*it <= pw/2.0)
            {
		    t_mid[iper] = t[it_total];
                    T_mid[iper] = T_avg[it_total];
            }

            // Copy temperatures to old time step 
            for(unsigned int iy=0; iy<ny; iy++)
                Told[iy] = T[iy];

            // Find maximum AR temperature
            if(T_avg[it_total] > T_max[iper])
            {
                T_max[iper] = T_avg[it_total];
                t_max[iper] = t[it_total];
            }

            // Find minimum AR temperature
            if(T_avg[it_total] < T_min[iper])
            {
                T_min[iper] = T_avg[it_total];
                t_min[iper] = t[it_total];
            }

            // Record rising and falling edges if this is the last pulse
            if(iper == _n_rep-1)
            {
                if(fmod(t[it_total], time_period) <= pw) // if during pulse
                {
                    _t_rise.push_back(t[it_total]*1e6);
                    _T_rise.push_back(T_avg[it_total]);
                }
                else
                {
                    _t_fall.push_back(t[it_total]*1e6);
                    _T_fall.push_back(T_avg[it_total]);
                }

                t_period[it] = t[it_total]*1e6;
                T_period[it] = T_avg[it_total];
            }
        }// end time loop

        // Print progress to screen
        if(opt.get_verbose())
        {
            printf("Period=%u Tmax= %.4f K at t=%.2f microseconds\n",
                   iper+1,
                   T_max[iper],
                   t_max[iper]*1e6);
        }
    }// end period loop

    std::valarray<double> t_rise(&_t_rise[0], _t_rise.size());
    std::valarray<double> T_rise(&_T_rise[0], _T_rise.size());
    std::valarray<double> t_fall(&_t_fall[0], _t_fall.size());
    std::valarray<double> T_fall(&_T_fall[0], _T_fall.size());

    write_table("T_t.dat", std::valarray<double>(1e6*t), T_avg);
    write_table("T-mid_t.dat",std::valarray<double>(1e6*t_mid), T_mid);
    write_table("Tmax_t.dat", std::valarray<double>(1e6*t_max), T_max);
    write_table("Tmin_t.dat", std::valarray<double>(1e6*t_min), T_min);
    write_table("Trise_t.dat", t_rise, T_rise);
    write_table("Tfall_t.dat", t_fall, T_fall);
    write_table("T_y.dat", y, T);
    write_table("T-period_t.dat", t_period, T_period);

    delete[] Told;
    delete[] null_vec;

    return EXIT_SUCCESS;
}


// Calculate the spatial temperature profile across the device at the
// next time step in the sequence.
static void calctemp(double dt,
                     double *Told,
                     double q_old[],
                     double q_new[], 
                     std::vector<MaterialSpecification> &mat,
                     std::valarray<double>& T,
                     Thermal1DOptions& opt)
{
    const size_t ny = mat.size();
    const double dy = opt.get_numeric_option("dy");
    double dy_sq = dy*dy;
    double* RHS_subdiag   = new double[ny-1];
    double* RHS_diag      = new double[ny];
    double* RHS_superdiag = new double[ny-1];
    double* LHS_subdiag   = new double[ny-1];
    double* LHS_diag      = new double[ny];
    double* LHS_superdiag = new double[ny-1];

    // Note that the bottom of the device is not calculated.  We leave it
    // set to the heatsink temperature (Dirichlet boundary)
    LHS_diag[0] = 1.0;
    LHS_superdiag[0] = 0.0;

    RHS_diag[0] = 1.0;
    RHS_superdiag[0] = 0.0;

    T[0] = 0.0;
    double k_last = thermal_cond(mat[0],Told[0]);
    double k      = thermal_cond(mat[1],Told[1]);
    double k_next = thermal_cond(mat[2],Told[2]);

    double rho_cp = 0;

    for(unsigned int iy=1; iy<ny-1; iy++)
    {
        const double rho = mat[iy].get_prop_val_x("density");

        // Product of density and spec. heat cap [J/(m^3.K)]
        const double _cp = cp(mat[iy], Told[iy]);
        rho_cp = rho * _cp;

        // Find interface values of the thermal conductivity using
        // Eq. 3.25 in Craig's thesis.
        double kn=(2*k*k_next)/(k+k_next);
        double ks=(2*k*k_last)/(k+k_last);
        double r = dt/(2.0*rho_cp);
        double alpha = r*ks/dy_sq;
        double gamma = r*kn/dy_sq;

        RHS_subdiag[iy-1] = alpha;
        RHS_superdiag[iy] = gamma;
        RHS_diag[iy] = 1.0 - (alpha+gamma);

        LHS_subdiag[iy-1] = -alpha;
        LHS_superdiag[iy] = -gamma;
        LHS_diag[iy] = 1.0 + (alpha+gamma);

        T[iy] = r*(q_old[iy]+q_new[iy]);

        k_last = k;
        k = k_next;
        k_next = thermal_cond(mat[iy+1],Told[iy+1]);
    }

    // At last point, use Neumann boundary, i.e. dT/dy=0, which gives
    // T[n] = T[n-2] in the finite-difference approximation
    double kns=(2*k_next*k)/(k_next+k);
    const double rho = mat[ny-1].get_prop_val_x("density");
    rho_cp = rho * cp(mat[ny-1], Told[ny-1]);
    double r = dt/(2.0*rho_cp);
    double alpha_gamma = r*kns/dy_sq;
    RHS_subdiag[ny-2] = 2.0*alpha_gamma;
    RHS_diag[ny-1] = 1.0 - 2.0*alpha_gamma;
    LHS_subdiag[ny-2] = -2.0*alpha_gamma;
    LHS_diag[ny-1] = 1.0 + 2.0*alpha_gamma;
    T[ny-1] = r*(q_old[ny-1]+q_new[ny-1]);

    int N = ny;
    double scale = 1;
    char TRANS='N';
    int NRHS=1;

    // Perform matrix multiplication using LAPACK
    // result:-> result + RHS*Told
    dlagtm_(&TRANS, &N, &NRHS, &scale, RHS_subdiag, RHS_diag, RHS_superdiag,
            Told, &N, &scale, &T[0], &N);

#if HAVE_LAPACKE
    int INFO=LAPACKE_dgtsv(LAPACK_COL_MAJOR, // Use column-major ordering
            ny, // Order of matrix
            1,  // Number of right-hand-side
            LHS_subdiag,
            LHS_diag,
            LHS_superdiag,
            &T[0],
            ny);
#else
    int INFO=1;
    dgtsv_(&N, &NRHS, LHS_subdiag, LHS_diag, LHS_superdiag, &T[0], &N, &INFO);
#endif
    
    if(INFO!=0)
        throw std::runtime_error("Could not solve matrix inversion problem. Check all input parameters!");

    delete[] RHS_diag;
    delete[] RHS_subdiag;
    delete[] RHS_superdiag;
    delete[] LHS_diag;
    delete[] LHS_subdiag;
    delete[] LHS_superdiag;
}

/**
 * Find average temperature inside active region
 *
 * \param[in] g  The power density spatial profile [W/m^3]
 * \param[in] T  Temperature profile [K]
 *
 * \returns The spatial average of temperature inside active region [K]
 *
 * \details The active region is identified as being the section in which
 *          heat is injected.  There may be a more general way to do this,
 *          such that heating in the contacts etc can be modelled without
 *          affecting this function.
 */
static double calctave(const std::valarray<double> &g,
                       const std::valarray<double> &T)
{
    double T_AR_cumulative=0;
    unsigned int n_AR=0;
    const double epsilon = g.max()/1e20;

    for(unsigned int iy=0; iy < T.size(); iy++)
    {
        // Consider this point if the power density is nonzero
        // Note that this assumes that only the active-region is heated
        if(gsl_fcmp(g[iy],0.0,epsilon) == 1)
        {
            n_AR++;
            T_AR_cumulative += T[iy];
        }
    }

    return T_AR_cumulative/float(n_AR);
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
