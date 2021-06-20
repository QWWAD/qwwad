#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_math.h>
#include "qwwad/options.h"
#include "qwwad/material.h"
#include "qwwad/material-library.h"
#include "qwwad/debye.h"
#include "qwwad/material-property-numeric.h"
#include "qwwad/maths-helpers.h"
#include "qwwad/file-io.h"
#include "qwwad/linear-algebra.h"
#include <glibmm/ustring.h>

using namespace QWWAD;

/**
 * Find the thermal conductivity of a material [W/m/K]
 *
 * \param[in] mat The material system
 * \param[in] x   Alloy fraction (if applicable)
 * \param[in] T   Temperature [K]
 *
 * \todo Figure out where all these values come from!
 * \todo These values only work for a limited range of
 *       temperatures. Restrict the domain accordingly?
 */
static auto thermal_cond(const Material &mat,
                           const double    x,
                           const double    T) -> double
{
    double k = 0.0;

    try
    {
        k = mat.get_property_value("thermal-conductivity-vs-alloy", x);
    }
    catch(std::exception &e)
    {
        try
        {
            double k0_1  = mat.get_property_value("thermal-conductivity-0K-1");
            double k0_2  = mat.get_property_value("thermal-conductivity-0K-2");
            double tau_1 = mat.get_property_value("thermal-conductivity-decay-index-1");
            double tau_2 = mat.get_property_value("thermal-conductivity-decay-index-2");
            double k_1 = k0_1*pow(T,tau_1);
            double k_2 = k0_2*pow(T,tau_2);
            k = lin_interp(k_1,k_2, x);
        }
        catch(std::exception &e)
        {
            try
            {
                const double k0  = mat.get_property_value("thermal-conductivity-0K");
                const double tau = mat.get_property_value("thermal-conductivity-decay-index");
                k = k0 * pow(T,tau);
            }
            catch(std::exception &e)
            {
                try
                {
                    const double k_hi    = mat.get_property_value("thermal-conductivity-high-T");
                    const double k_decay = mat.get_property_value("thermal-conductivity-inverse-T");
                    k = k_hi + k_decay/T;
                }
                catch(std::exception &e)
                {
                    k = mat.get_property_value("thermal-conductivity-T", T);
                }
            }
        }
    }

    return k;
}

class Thermal1DOptions: public Options
{
    double dc = 0.02; // Duty cycle
    double f  = 10e3; // Pulse repetition rate [Hz]

public:
    Thermal1DOptions(int argc, char** argv);

    /// Return fractional duty cycle (i.e. 0 to 1)
    [[nodiscard]] auto get_duty_cycle() const -> double {return dc;}

    /// Return pulse repetition rate [Hz]
    [[nodiscard]] auto get_f_rep() const -> double {return f;}

    void print() const;
};


// Define and parse all user options and return them in a structure
Thermal1DOptions::Thermal1DOptions(int argc, char ** argv)
{
    std::string doc = "Calculate temperature variation in active region over time";

    add_option<size_t>     ("active,a",                   2, "Index of active-region layer");
    add_option<double>     ("area",                   0.119, "QCL ridge area [mm^2]");
    add_option<std::string>("infile",  "thermal_layers.dat", "Waveguide layers data file");
    add_option<double>     ("Tsink,T",                 80.0, "Heatsink temperature [K]");
    add_option<double>     ("dy,y",                 1.00e-7, "Spatial resolution [m]");
    add_option<double>     ("dc,d",                       2, "Duty cycle for pulse train [%]");
    add_option<double>     ("frequency,f",               10, "Pulse repetition rate [kHz]");
    add_option<double>     ("power,P",                17.65, "Pulse power [W]");
    add_option<size_t>     ("nrep",                       1, "Number of pulse periods to simulate");

    add_prog_specific_options_and_parse(argc,argv,doc);

    // Check that heatsink temperature is positive
    const auto Tsink = get_option<double>("Tsink");
    if (Tsink <= 0.0)
    {
        std::ostringstream oss;
        oss << "Heatsink temperature, " << Tsink << " is not positive.";
        throw std::domain_error(oss.str());
    }

    // Check that spatial step is positive
    const auto dy = get_option<double>("dy");
    if(dy <= 0.0) {
        throw std::domain_error ("Spatial resolution must be positive");
    }

    // Check that duty cycle is positive and
    // rescale to a decimal value
    dc = get_option<double>("dc") * 0.01;

    if(dc <= 0.0 or dc >= 1.0)
    {
        std::ostringstream oss;
        oss << "Specified duty cycle, " << dc << " is invalid.";
        throw std::domain_error(oss.str());
    }

    // Check that frequency is positive and
    // rescale to Hz
    f = get_option<double>("frequency") * 1.0e3;

    if(f <= 0) {
        throw std::domain_error ("Pulse repetition rate must "
                                 "be positive.");
    }

    // Check that power is positive
    const auto power = get_option<double>("power");
    if(power <= 0.0)
    {
        std::ostringstream oss;
        oss << "Electrical power dissipation, " << power << " is not positive.";
        throw std::domain_error(oss.str());
    }

    // Check that area is positive
    const auto area = get_option<double>("area");
    if(area <= 0.0)
    {
        std::ostringstream oss;
        oss << "Ridge area, " << area << " is not positive.";
        throw std::domain_error(oss.str());
    }

    if(get_verbose()) {
        print();
    }
}


void Thermal1DOptions::print() const
{
    std::cout << "Heat sink temperature = " << get_option<double>("Tsink")  << " K"           << std::endl
              << "Power                 = " << get_option<double>("power")  << " W"           << std::endl
              << "Frequency             = " << f/1.0e3                      << " kHz"         << std::endl
              << "Duty cycle            = " << get_duty_cycle()*100         << "%"            << std::endl
              << "Period length         = " << 1e6/f                        << " microsecond" << std::endl
              << "Number of periods     = " << get_option<size_t>("nrep")                     << std::endl
              << "Spatial resolution    = " << get_option<double>("dy")*1e6 << " micron"      << std::endl
              << "Ridge area            = " << get_option<double>("area")   << " mm^2"        << std::endl;
}

class Thermal1DData {
public:
    Thermal1DData(const Thermal1DOptions &opt,
                  const MaterialLibrary  &material_library);
    std::vector<Material> mat_layer; ///< Material in each layer
    arma::vec x;         ///< Alloy composition in each layer
    arma::vec d;         ///< Layer thickness [m]
};

Thermal1DData::Thermal1DData(const Thermal1DOptions &opt,
                             const MaterialLibrary  &material_library)
{
    arma::vec doping; // Unused doping data
    std::vector<std::string> mat_name; // Materialname
    auto infile = opt.get_option<std::string>("infile");
    read_table(infile, d, x, doping, mat_name);

    d *= 1e-6; // Rescale thickness to metres

    for(auto const &name : mat_name) {
        const auto *material = material_library.get_material(name);
        mat_layer.push_back(*material);
    }

    if(opt.get_verbose()) {
        std::cout << "Read " << d.size() << " layers from " << infile << std::endl;
    }

    if(d.empty()) {
        std::ostringstream oss;
        oss << "Could not read any layers from " << infile;
        throw std::runtime_error(oss.str());
    }
}

static auto calctave(const arma::vec &g,
                       const arma::vec &T) -> double;

static auto calctemp(double dt,
                          arma::vec        &Told,
                          arma::vec  const &q_old,
                          arma::vec  const &q_new,
                          arma::uvec const &iLayer,
                          const std::vector<Material> &mat,
                          const arma::vec   &x,
                          const std::vector<DebyeModel> &dm_layer,
                          const arma::vec   &rho_layer,
                          Thermal1DOptions& opt) -> arma::vec;

auto main(int argc, char *argv[]) -> int
{
    constexpr float CM2_TO_M2 = 1e-6;
    constexpr float S_TO_NS   = 1e9;

    // Grab user preferences
    MaterialLibrary material_library("");
    Thermal1DOptions opt(argc, argv);
    const Thermal1DData data(opt, material_library);

    const auto dy = opt.get_option<double>("dy");
    const auto L  = sum(data.d); // Length of structure [m]
    const size_t ny = ceil(L/dy);   // Find number of points in structure

    if(opt.get_verbose()) {
        std::cout << "ny = " << ny << std::endl;
    }

    // Thickness of active region [m]
    const auto iAR    = opt.get_option<size_t>("active");
    const auto L_AR   = data.d[iAR];
    const auto area   = opt.get_option<double>("area")*CM2_TO_M2;
    const auto volume = L_AR * area;  // Active region volume (L x h x w) [m^3]

    // Power density in active region [W/m^3]
    const auto power_density = opt.get_option<double>("power")/volume;
    const auto pw = opt.get_duty_cycle()/opt.get_f_rep();

    if(opt.get_verbose())
    {
        printf("Power density = %5.2e W/m3.\n", power_density);
        printf("Pulse width = %5.1f ns.\n",pw*S_TO_NS);
    }

    auto const y = arma::linspace(0, L, ny);                // Spatial coordinates [m]

    arma::vec  g = arma::zeros(ny); // Power density profile [W/m^3]
    arma::uvec iLayer = arma::zeros<arma::uvec>(ny);     // Index of layer containing each point

    double bottom_of_layer=0;
    unsigned int iy=1;

    std::vector<DebyeModel> dm_layer;
    const size_t nL = data.d.size();
    arma::vec rho_layer = arma::zeros(nL);

    // Loop through each layer and figure out which points it contains
    for(unsigned int iL=0; iL < nL; iL++){
        bottom_of_layer += data.d[iL];

        // Check that we haven't finished filling the array and that
        // the next point is still in this layer
        while(iy<ny and y(iy-1)+dy < bottom_of_layer)
        {
            iLayer(iy) = iL;    // Note the layer containing this point

            // Assume that only the active-region is heated
            if (iL == iAR) {
                g(iy) = power_density;
            }

            iy++;
        }

        double T_D = 0.0;
        double M   = 0.0;
        unsigned int natoms = 0;

        try {
            T_D = data.mat_layer[iL].get_property_value("debye-temperature", data.x[iL]);
            M   = data.mat_layer[iL].get_property_value("molar-mass", data.x[iL]);
            natoms = data.mat_layer[iL].get_property_value("natoms");
            rho_layer[iL] = data.mat_layer[iL].get_property_value("density", data.x[iL]);
        } catch (std::exception &e) {
            std::cerr << "Could not find material parameters for "
                      << data.mat_layer[iL].get_description() << std::endl;
            exit(EXIT_FAILURE);
        }

        dm_layer.emplace_back(T_D, M, natoms);
    }

    const auto _Tsink = opt.get_option<double>("Tsink");
    
    // Spatial temperature profile through structure [K].
    // Assume that initially all points are in thermal equilibrium with heat sink.
    arma::vec T = arma::zeros(ny);
    T.fill(_Tsink);
    
    arma::vec Told = arma::zeros(ny);
    Told.fill(_Tsink);

    std::ofstream FT("struct.dat");

    for(unsigned int iy=0; iy<ny; iy++)
    {
        auto const iL = iLayer(iy); // Look up layer containing this point
        auto const mat = data.mat_layer[iL]; // Get the material in the layer

        // Now save the material to file
        FT << iy*dy*1e6 << "\t" << mat.get_description() << std::endl;
    }

    FT.close();

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

    if(opt.get_verbose()) {
        printf("dt=%.4f ns.\n",dt*1e9);
    }

    const auto _n_rep = opt.get_option<size_t>("nrep"); // Number of pulse repetitions
   
    // Samples of average T_AR at each time-step 
    arma::vec t = arma::zeros(nt_per * _n_rep);
    arma::vec T_avg = arma::zeros(nt_per * _n_rep);
    
    // Samples of average T_AR in the middle of each pulse
    arma::vec t_mid = arma::zeros(_n_rep); // Time at middle of each pulse
    arma::vec T_mid = arma::zeros(_n_rep); // Average AR temperature at middle of each pulse

    // Samples of maximum T_AR in each pulse
    arma::vec t_max = arma::zeros(_n_rep); // Time at which peak temp occurs in each pulse
    arma::vec T_max = arma::zeros(_n_rep); // Peak AR temperature of each pulse

    // Samples of minimum T_AR in each pulse
    arma::vec t_min = arma::zeros(_n_rep); // Time at which minimum temp occurs in each pulse
    arma::vec T_min = arma::zeros(_n_rep); // Minimum AR temperature of each pulse
    T_min.fill(1e9);

    // Samples of average T_AR across final pulse
    arma::vec t_period = arma::zeros(nt_per);
    arma::vec T_period = arma::zeros(nt_per);

    // Temperature profile at end of final pulse
    arma::vec T_y_max = arma::zeros(ny);

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
            arma::vec q_old = arma::zeros(ny);
            arma::vec q_now = arma::zeros(ny);

            // If this time-step is within the pulse, then
            // "switch on" the electrical power
            if (dt*it <= pw) {
                q_now = g;
            }

            // Likewise for the previous time-step
            if (it > 0 and dt*(it-1) <= pw) {
                q_old = g;
            }

            // Calculate the spatial temperature profile at this 
            // timestep
            T = calctemp(dt, Told, q_old, q_now, iLayer, data.mat_layer, data.x, dm_layer, rho_layer, opt);

            // Find spatial average of T_AR
            T_avg(it_total) = calctave(g, T);
            
            // Find T_AR at middle of the pulse
	    if(dt*it <= pw/2.0) {
                t_mid(iper) = t(it_total);
                T_mid(iper) = T_avg(it_total);
            }

            // Copy temperatures to old time step 
            for(unsigned int iy=0; iy<ny; iy++) {
                Told(iy) = T(iy);
            }

            // Find maximum AR temperature
            if(T_avg(it_total) > T_max(iper)) {
                T_max(iper) = T_avg(it_total);
                t_max(iper) = t(it_total);

                for(unsigned int iy = 0; iy < ny; ++iy) {
                    T_y_max(iy) = T(iy);
                }
            }

            // Find minimum AR temperature
            if(T_avg(it_total) < T_min(iper))
            {
                T_min(iper) = T_avg(it_total);
                t_min(iper) = t(it_total);
            }

            // Record rising and falling edges if this is the last pulse
            if(iper == _n_rep-1)
            {
                if(fmod(t(it_total), time_period) <= pw) // if during pulse
                {
                    _t_rise.push_back(t(it_total)*1e6);
                    _T_rise.push_back(T_avg(it_total));
                }
                else
                {
                    _t_fall.push_back(t(it_total)*1e6);
                    _T_fall.push_back(T_avg(it_total));
                }

                t_period(it) = t(it_total)*1e6;
                T_period(it) = T_avg(it_total);
            }
        }// end time loop

        // Print progress to screen
        if(opt.get_verbose())
        {
            printf("Period=%u Tmax= %.4f K at t=%.2f microseconds\n",
                   iper+1,
                   T_max(iper),
                   t_max(iper)*1e6);
        }
    }// end period loop

    write_table("T_t.dat",    arma::vec(1e6*t), T_avg);
    write_table("T-mid_t.dat",arma::vec(1e6*t_mid), T_mid);
    write_table("Tmax_t.dat", arma::vec(1e6*t_max), T_max);
    write_table("Tmin_t.dat", arma::vec(1e6*t_min), T_min);
    write_table("Trise_t.dat", _t_rise, _T_rise);
    write_table("Tfall_t.dat", _t_fall, _T_fall);
    write_table("T_y.dat", y, T);
    write_table("T_y_max.dat", y, T_y_max);
    write_table("T-period_t.dat", t_period, T_period);

    return EXIT_SUCCESS;
}

// Calculate the spatial temperature profile across the device at the
// next time step in the sequence.
static auto calctemp(double dt,
                          arma::vec  &Told,
                          arma::vec  const &q_old,
                          arma::vec  const &q_new,
                          arma::uvec const &iLayer,
                          const std::vector<Material> &mat_layer,
                          const arma::vec &x,
                          const std::vector<DebyeModel> &dm_layer,
                          const arma::vec &rho_layer,
                          Thermal1DOptions& opt) -> arma::vec
{
    const auto ny = iLayer.size();
    const auto dy = opt.get_option<double>("dy");
    const auto dy_sq = dy*dy;

    // Note that the bottom of the device is not calculated.  We leave it
    // set to the heatsink temperature (Dirichlet boundary)
    arma::vec LHS_diag      = arma::ones(ny);
    arma::vec LHS_subdiag   = arma::zeros(ny-1);
    arma::vec LHS_superdiag = arma::zeros(ny-1);

    // Material parameter matrix for RHS of Crank-Nicolson solver
    arma::vec B_diag      = arma::ones(ny);
    arma::vec B_superdiag = arma::zeros(ny-1);
    arma::vec B_subdiag   = arma::zeros(ny-1);

    // Indices of layers containing the current, previous and next points
    auto iL_prev = iLayer(0);
    auto iL_this = iLayer(1);
    auto iL_next = iLayer(2);

    double k_prev = thermal_cond(mat_layer[iL_prev], x(iL_prev), Told(0));
    double k_this = thermal_cond(mat_layer[iL_this], x(iL_this), Told(1));
    double k_next = thermal_cond(mat_layer[iL_next], x(iL_next), Told(2));

    double rho_cp = 0;

    arma::vec q = arma::zeros(ny); // heating vector for RHS of Crank-Nicolson solver

    for(unsigned int iy=1; iy<ny-1; iy++)
    {
        iL_this = iLayer(iy); // Update the current layer index
        iL_next = iLayer(iy+1);

        // Product of density and spec. heat cap [J/(m^3.K)]
        const auto _cp = dm_layer[iL_this].get_cp(Told(iy));
        rho_cp = rho_layer(iL_this) * _cp;

        // Find interface values of the thermal conductivity using
        // Eq. 3.25 in Craig's thesis.
        const auto kn=(2*k_this*k_next)/(k_this+k_next);
        const auto ks=(2*k_this*k_prev)/(k_this+k_prev);
        const auto r = dt/(2.0*rho_cp);
        const auto alpha = r*ks/dy_sq;
        const auto gamma = r*kn/dy_sq;

        B_subdiag(iy-1) = alpha;
        B_superdiag(iy) = gamma;
        B_diag(iy)      = 1.0 - (alpha+gamma);

        LHS_subdiag(iy-1) = -alpha;
        LHS_superdiag(iy) = -gamma;
        LHS_diag(iy)      = 1.0 + (alpha+gamma);

        q(iy) = r*(q_old(iy)+q_new(iy));

        k_prev = k_this;
        k_this = k_next;
        k_next = thermal_cond(mat_layer[iL_next], x(iL_next), Told(iy+1));
    }

    // At last point, use Neumann boundary, i.e. dT/dy=0, which gives
    // T[n] = T[n-2] in the finite-difference approximation
    double kns=(2*k_next*k_this)/(k_next+k_this);
    const double rho = rho_layer(iLayer(ny-1));
    rho_cp = rho * dm_layer[iLayer(ny-1)].get_cp(Told(ny-1));
    double r = dt/(2.0*rho_cp);
    double alpha_gamma = r*kns/dy_sq;
    B_subdiag(ny-2) = 2.0*alpha_gamma;
    B_diag(ny-1)    = 1.0 - 2.0*alpha_gamma;

    LHS_subdiag(ny-2) = -2.0*alpha_gamma;
    LHS_diag(ny-1) = 1.0 + 2.0*alpha_gamma;

    q(ny-1) = r*(q_old(ny-1)+q_new(ny-1));

    // Perform matrix multiplication to get the RHS vector of the
    // Crank-Nicolson solver
    auto RHS = multiply_vec_tridiag(B_subdiag,
                                    B_diag,
                                    B_superdiag,
                                    Told,
                                    q);

    auto T = solve_tridiag(LHS_subdiag,
                           LHS_diag,
                           LHS_superdiag,
                           RHS);

    return T;
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
static auto calctave(const arma::vec &g,
                       const arma::vec &T) -> double
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
