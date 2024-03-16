#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "GalSim.h"
#include <omp.h>
#include <vector>
#include <complex>

namespace py = pybind11;

py::array_t<double> getFluxVec(
    const py::object& profileObj,
    const py::array_t<double>& xyCoords
){
    if (xyCoords.ndim() != 2 || xyCoords.shape(0) != 2) {
        throw std::runtime_error("xyCoords must be a 2D array with shape (2, n)");
    }

    auto xy = xyCoords.unchecked<2>(); // Access without bounds checking for efficiency
    size_t n_points = xyCoords.shape(1);
    std::vector<double> fluxes(n_points);

    galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(profileObj);
    #pragma omp parallel for
    for(size_t i = 0; i < n_points; ++i) {
        fluxes[i] = profile.xValue(galsim::Position<double>(xy(0, i), xy(1, i)));
    }

    size_t dim = std::sqrt(n_points);
    return py::array_t<double>({dim, dim}, fluxes.data());
}

// Utility function to generate rfftfreq
std::vector<double> rfftfreq(int n, double scale) {
    std::vector<double> result(n/2 + 1);
    for(int i = 0; i <= n / 2; ++i) {
        result[i] = i / (scale * n);
    }
    return result;
}

// Utility function to generate fftfreq
std::vector<double> fftfreq(int n, double scale) {
    std::vector<double> result(n);
    double val = 1.0 / (n * scale);
    for(int i = 0; i < n; ++i) {
        if (i < (n + 1) / 2) {
            result[i] = i * val;
        } else {
            result[i] = (i - n) * val;
        }
    }
    return result;
}

py::array_t<std::complex<double>> mulFourier(
    double scale,
    const py::object& profileObj,
    const py::array_t<std::complex<double>>& gal_kprof
){

    // Prepare two kx, ky grids
    size_t dim_y = gal_kprof.shape(0);
    size_t dim_x = dim_y / 2 + 1;
    size_t n_points = dim_x * dim_y;
    auto x = rfftfreq(dim_y, scale / M_PI / 2.0);
    auto y = fftfreq(dim_y, scale / M_PI / 2.0);

    // Prepare the output array, galaxy array and psf gsobj
    auto gal = gal_kprof.unchecked<2>();
    std::vector<std::complex<double>> fluxes(n_points);
    galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(profileObj);

    #pragma omp parallel for
    for(size_t i = 0; i < n_points; ++i) {
        fluxes[i] = gal(i / dim_x, i % dim_x) * profile.kValue(
            galsim::Position<double>(
                x[i % dim_x],
                y[i / dim_x]
            )
        );
    }

    return py::array_t<std::complex<double>>({dim_y, dim_x}, fluxes.data());
}

PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Pybind11 interface for GalSim flux and Fourier computations";
    m.def(
        "getFluxVec", &getFluxVec,
            "Get flux values at multiple x,y coordinates",
        py::arg("profileObj"),
        py::arg("xyCoords")
    );
    m.def(
        "mulFourier", &mulFourier,
            "Get Fourier transform values at multiple kx, ky coordinates",
        py::arg("scale"),
        py::arg("profileObj"),
        py::arg("gal_kprof")
    );
}
