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
    const double scale,
    const py::object& gsobj,
    const py::array_t<double>& xy_coords
){
    if (xy_coords.ndim() != 2 || xy_coords.shape(0) != 2) {
        throw std::runtime_error("xy_coords must be a 2D array with shape (2, n)");
    }

    auto xy = xy_coords.unchecked<2>();
    const int n_points = xy_coords.shape(1);
    std::vector<double> fluxes(n_points);

    double area = scale * scale;
    galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(gsobj);
    #pragma omp parallel for
    for(int i = 0; i < n_points; ++i) {
        fluxes[i] = profile.xValue(
            galsim::Position<double>(xy(0, i), xy(1, i))
        ) * area;
    }

    int dim = std::sqrt(n_points);
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

bool is_c_contiguous(const py::array& arr) {
    return arr.flags() & py::array::c_style;
}

// Convolve gal_prof (defined in configuration space) with Galsim PSF object
// Then perform a down sampling
py::array_t<double> convolvePsf(
    const double scale,
    const py::object& gsobj,
    const py::array_t<double>& gal_prof,
    const int downsample_ratio,
    const int ngrid
){

    bool test = is_c_contiguous(gal_prof);
    if (! test) {
        throw std::runtime_error(
            "Input galaxy array is not continuous in memory"
        );
    }
    auto info = gal_prof.request();
    int dim = info.shape[0];

    // down sampled scale and dimension
    double scale2 = scale * downsample_ratio;
    int dim2 = dim / downsample_ratio;

    // Frequency grids for the down sampled signal
    const auto x_freqs2 = rfftfreq(dim2, scale2 / M_PI / 2.0);
    const auto y_freqs2 = fftfreq(dim2, scale2 / M_PI / 2.0);

    // Allocate FFTW arrays with pointers
    double* in = static_cast<double*>(info.ptr);
    fftw_complex* out = fftw_alloc_complex(dim * (dim / 2 + 1));
    fftw_complex* out2 = fftw_alloc_complex(dim2 * (dim2 / 2 + 1));

    // Plan and execute forward FFT
    fftw_plan p_forward = fftw_plan_dft_r2c_2d(dim, dim, in, out, FFTW_ESTIMATE);
    fftw_execute(p_forward);

    // Galsim object
    const galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(gsobj);

    // Process FFT result using gsobj
    #pragma omp parallel for
    for (int y2 = 0; y2 < dim2; ++y2) {
        int y = (y2 >= dim2 / 2) ? (dim - dim2 + y2) : y2;
        for (int x2 = 0; x2 < (dim2 / 2 + 1); ++x2) {
            int x = x2;
            int index = y * (dim / 2 + 1) + x;
            int index2 = y2 * (dim2 / 2 + 1) + x2;
            std::complex<double> fft_val(out[index][0], out[index][1]);
            std::complex<double> result = fft_val * profile.kValue(
                galsim::Position<double>(
                    x_freqs2[x2],
                    y_freqs2[y2]
                )
            );
            out2[index2][0] = result.real();
            out2[index2][1] = result.imag();
        }
    }
    // Cleanup fftw
    fftw_destroy_plan(p_forward);
    fftw_free(out);

    // Allocate array for inverse FFT result
    double* ifft_out = fftw_alloc_real(dim2 * dim2);
    // Plan and execute inverse FFT
    fftw_plan p_backward = fftw_plan_dft_c2r_2d(dim2, dim2, out2, ifft_out, FFTW_ESTIMATE);
    fftw_execute(p_backward);
    // Cleanup fftw
    fftw_destroy_plan(p_backward);
    fftw_free(out2);

    // Wrap the result in a numpy array
    // prevent memory leakage
    auto result = py::array_t<double>({ngrid, ngrid});
    // Use unchecked for faster access
    auto r = result.mutable_unchecked<2>();
    // Fill with 0
    std::fill(result.mutable_data(), result.mutable_data() + ngrid * ngrid, 0.0);

    int dim2_center = dim2 / 2;
    int res_center = ngrid / 2;

    // Normalize the inverse FFT result
    const int norm_factor = dim2 * dim2;
    if (dim2_center >= res_center) {
        // shrinking
        int start = 0;
        int end = ngrid;
        for (int y = start; y < end; ++y) {
            int yf = y - res_center + dim2_center;
            for (int x = start; x < end; ++x) {
                // Calculate the corresponding source index in ifft_out
                int xf = x - res_center + dim2_center;
                // Copy and normalize
                r(y, x) = ifft_out[yf * dim2 + xf] / norm_factor;
            }
        }
    }
    else {
        // padding zeros
        int start = res_center - dim2_center;
        int end = res_center + dim2_center;
        for (int y = start; y < end; ++y) {
            int yf = y - res_center + dim2_center;
            for (int x = start; x < end; ++x) {
                // Calculate the corresponding source index in ifft_out
                int xf = x - res_center + dim2_center;
                // Copy and normalize
                r(y, x) = ifft_out[yf * dim2 + xf] / norm_factor;
            }
        }

    }

    // Cleanup fftw
    fftw_free(ifft_out);
    return result;
}

PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Pybind11 interface for GalSim flux and Fourier computations";
    m.def(
        "getFluxVec", &getFluxVec,
            "Get flux values at multiple x,y coordinates",
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("xy_coords")
    );
    m.def(
        "convolvePsf", &convolvePsf,
            "Convolve galaxy profile with PSF",
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("gal_prof"),
        py::arg("downsample_ratio"),
        py::arg("ngrid")
    );
}
