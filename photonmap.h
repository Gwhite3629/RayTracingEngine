#ifndef _PHOTONMAP_H_
#define _PHOTONMAP_H_

#include "constants.h"

#include <vector>
#include <complex>
#include <functional>
#include <algorithm>

struct photon_data {
    double speed;
    double frequency;
};

struct derived_data {
    double wavelength;
    double wavenumber;
    std::vector<double> wavevector;
    std::vector<double> momentum;
    double energy;
};

class photon {
public:

    photon_data data;
    derived_data p_data;
    std::vector<double> position;
    std::vector<double> direction;

    void derive (void) {
        p_data = {
            wavelength(),
            wavenumber(),
            wavevector(),
            momentum(),
            energy()
        };
    }

private:

    inline double wavelength(void) {
        return data.speed/data.frequency;
    }

    inline double wavenumber(void) {
        return 2.0*pi()/wavelength();
    }

    inline std::vector<double> wavevector(void) {
        std::vector<double> out = direction;
        std::transform(out.begin(), out.end(), out.begin(), [w = wavenumber()](double e){return e *= w;});
        return out;
    }

    inline std::vector<double> momentum(void) {
        std::vector<double> out = wavevector();
        std::transform(out.begin(), out.end(), out.begin(), [p = reduced_plank](double e){return e *= p;});
        return out;
    }

    inline double energy(void) {
        return plank * data.frequency;
    }
};

#endif // _PHOTONMAP_H_