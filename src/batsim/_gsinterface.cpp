#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "GalSim.h"
#include <omp.h>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <cmath>


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



py::array_t<double> convolvePsf(
    double scale,
    const py::object& profileObj,
    const py::array_t<double>& gal_prof
){
    auto input_buf = gal_prof.request();
    size_t dim = input_buf.shape[0];

    // Allocate FFTW arrays
    double* in = fftw_alloc_real(dim * dim);
    fftw_complex* out = fftw_alloc_complex(dim * (dim / 2 + 1));

    // Copy the input data to FFTW input
    std::memcpy(in, input_buf.ptr, sizeof(double) * dim * dim);

    // Plan and execute forward FFT
    fftw_plan p_forward = fftw_plan_dft_r2c_2d(dim, dim, in, out, FFTW_ESTIMATE);
    fftw_execute(p_forward);

    // Frequency grids
    auto x_freqs = rfftfreq(dim, scale / M_PI / 2.0);
    auto y_freqs = fftfreq(dim, scale / M_PI / 2.0);
    // Galsim object
    galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(profileObj);

    // Process FFT result using profileObj
    #pragma omp parallel for
    for (size_t y = 0; y < dim; ++y) {
        for (size_t x = 0; x < (dim / 2 + 1); ++x) {
            size_t index = y * (dim / 2 + 1) + x;
            std::complex<double> fft_val(out[index][0], out[index][1]);
            std::complex<double> result = fft_val * profile.kValue(
                galsim::Position<double>(
                    x_freqs[x],
                    y_freqs[y]
                )
            );
            out[index][0] = result.real();
            out[index][1] = result.imag();
        }
    }

    // Allocate array for inverse FFT result
    double* ifft_out = fftw_alloc_real(dim * dim);

    // Plan and execute inverse FFT
    fftw_plan p_backward = fftw_plan_dft_c2r_2d(dim, dim, out, ifft_out, FFTW_ESTIMATE);
    fftw_execute(p_backward);

    // Normalize the inverse FFT result
    for (size_t i = 0; i < dim * dim; ++i) {
        ifft_out[i] /= (dim * dim);
    }

    // Wrap the result in a numpy array
    auto result = py::array_t<double>(input_buf.size, ifft_out);
    result.resize({dim, dim});

    // Cleanup
    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(in);
    fftw_free(out);
    fftw_free(ifft_out);

    return result;
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
        "convolvePsf", &convolvePsf,
            "Convolve galaxy profile with PSF",
        py::arg("scale"),
        py::arg("profileObj"),
        py::arg("gal_prof")
    );
}
