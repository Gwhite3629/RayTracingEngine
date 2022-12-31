#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "constants.h"
#include "photonmap.h"
#include "kdtree.h"

#include <vector>
#include <complex>
#include <functional>
#include <algorithm>

struct properties {
    double c; // Speed of light
    double permittivity;
    double permeability;
    double n; // Refractive index
    double A;
    double B;
    double absorption;
    double scattering;
    std::complex<double> Cn;
};

struct point {
    std::vector<double> normal;
};

class Surface {
public:
    void surface(properties material, properties atmosphere)
    {
        atm = atmosphere;
        mat = material;
    }

    void populate(void)
    {

    }

    std::complex<double> refractive_index(photon light)
    {
        std::complex<double> RI;
        RI.real(mat.A + (mat.B / (light.p_data.wavelength*light.p_data.wavelength)));
        RI.imag((mat.absorption * light.p_data.wavelength)/(4*pi()));
        return RI;
    }

    // Perpendicular component
    double reflected_S(photon light, point p)
    {
        std::complex<double> RI = refractive_index(light);
        set_atm_CRI(light);
        double theta_r = p.normal[0]*light.direction[0] + p.normal[1]*light.direction[1] + p.normal[1]*light.direction[1];
        double theta_i = theta_r;
        double theta_t = std::__complex_asin((RI * sin(theta_i))/atm.Cn).real();

        std::complex<double> C1 = RI * cos(theta_i);
        std::complex<double> C2 = atm.Cn * cos(theta_t);
        return pow(std::__complex_abs((C1-C2)/(C1+C2)), 2);
    }

    // Parallel component
    double reflected_P(photon light, point p)
    {
        std::complex<double> RI = refractive_index(light);
        set_atm_CRI(light);
        double theta_r = p.normal[0]*light.direction[0] + p.normal[1]*light.direction[1] + p.normal[1]*light.direction[1];
        double theta_i = theta_r;
        double theta_t = std::__complex_asin((RI * sin(theta_i))/atm.Cn).real();

        std::complex<double> C1 = RI * cos(theta_t);
        std::complex<double> C2 = atm.Cn * cos(theta_i);
        return pow(std::__complex_abs((C1-C2)/(C1+C2)), 2);
    }

    double transmitted_S(photon light, point p)
    {
        return 1 - reflected_S(light, p);
    }

    double transmitted_P(photon light, point p)
    {
        return 1 - reflected_P(light, p);
    }

    double reflected(photon light, point p)
    {
        double S_C = reflected_S(light, p);
        double P_C = reflected_P(light, p);
        return 0.5*(S_C + P_C);
    }

    double transmitted(photon light, point p)
    {
        return 1 - reflected(light, p);
    }

private:
    std::complex<double> Z;
    properties atm;
    properties mat;

    void set_atm_CRI(photon light)
    {
        atm.Cn.real(atm.A + (atm.B / (light.p_data.wavelength * light.p_data.wavelength)));
        atm.Cn.imag((atm.absorption * light.p_data.wavelength)/(4*pi()));
    }
};

#endif // _SURFACE_H_