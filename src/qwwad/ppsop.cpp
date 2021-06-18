#include "ppsop.h"

#include <cstring>
#include <sstream>
#include <stdexcept>

#include "constants.h"

using namespace QWWAD;
using namespace constants;
/**
 * \brief Returns the spin-orbit parameter lambda
 *
 * \param[in] type Atomic species
 *
 * \return The lambda parameter
 */
auto
lambda(const char *type) -> double
{
    double val = 0.0;

    // TODO: Make this a switch-case over an enum?
    if(strcmp(type,"SI") == 0) {
        val = 0.000106*h*c*Rinf;
    } else if(strcmp(type,"GE") == 0) {
        val = 0.00058*h*c*Rinf;
    } else if(strcmp(type,"GAASmz") == 0) {
        val = 0.000402*h*c*Rinf;
    } else if(strcmp(type,"ASGAmz") == 0) {
        val = 0.000402*1.38*h*c*Rinf;
    } else if(strcmp(type,"CDTE") == 0) {
        val = 0.055*(0.0250-0.0090)*h*c*Rinf;
    } else if(strcmp(type,"TECD") == 0) {
        val = 0.055*(0.0250+0.0090)*h*c*Rinf;
    } else if(strcmp(type,"CDTEcb") == 0) {
        val = 0.343*(0.002-0.0002)*h*c*Rinf;
    } else if(strcmp(type,"TECDcb") == 0) {
        val = 0.343*(0.002+0.0002)*h*c*Rinf;
    } else {
        std::ostringstream oss;
        oss << "Error atom type '" << type << "' undefined in spin-orbit parameter set!";
        throw std::runtime_error(oss.str());
    }

    return val;
}
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99 :
