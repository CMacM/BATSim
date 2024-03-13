#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstring>
#include "GalSim.h"
#include "fftw3.h"
#include <omp.h>
#include <Eigen/Dense>

namespace py = pybind11;

template <typename T>
std::vector<std::vector<T>> reshapeTo2D(const std::vector<T>& arr) {
    int length = arr.size();
    int dim = std::sqrt(length);
    if (dim * dim != length) {
        throw std::runtime_error("The length of the array is not a perfect square.");
    }
    std::vector<std::vector<T>> result(dim, std::vector<T>(dim));

    for (int i = 0; i < length; ++i) {
        int row = i / dim;
        int col = i % dim;
        result[row][col] = arr[i];
    }

    return result;
}

// Adjust the function to accept a py::object instead of a direct reference
// to galsim::SBProfile. This allows passing Python objects.
std::vector<double> getFluxVec(py::object profileObj, 
                               const std::vector<double>& xCoords, 
                               const std::vector<double>& yCoords) {
    galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(profileObj);

    

    // Pre-construct a vector of Position objects
    std::vector<galsim::Position<double>> positions;
    positions.reserve(xCoords.size());
    for(size_t i = 0; i < xCoords.size(); ++i) {
        positions.emplace_back(xCoords[i], yCoords[i]);
    }

    // Get flux values at multiple x,y coordinates using the pre-constructed positions
    std::vector<double> fluxes(xCoords.size());
    #pragma omp parallel for
    for(size_t i = 0; i < positions.size(); ++i) {
        double flux = profile.xValue(positions[i]);
        fluxes[i] = flux;
    }
    return fluxes;
}

std::vector<double> convolveImage(py::object gal_obj,
                                  py::object psf_obj,
                                  const std::vector<double>& galx,
                                  const std::vector<double>& galy,
                                  const std::vector<double>& psfx,
                                  const std::vector<double>& psfy,
                                  const double& rescale) {

    // Get flux values at multiple x,y coordinates using the pre-constructed positions
    std::vector<double> gal_fluxes = getFluxVec(gal_obj, galx, galy);  
    std::vector<double> psf_fluxes = getFluxVec(psf_obj, psfx, psfy);

    // Define sizes and required padding
    const int gal_nn = sqrt(galx.size());
    const int psf_nn = sqrt(psfx.size());

    int pad_width = 0;
    if (gal_nn > psf_nn) {
        pad_width = (gal_nn - psf_nn) / 2;
    }
    else if (psf_nn > gal_nn) {
        pad_width = (psf_nn - gal_nn) / 2;
    }

    // Compute sizes for cropping in Fourier space
    const int nn_cut = int(gal_nn / rescale);
    const int maxN = gal_nn / 2 + nn_cut / 2;
    const int minN = gal_nn / 2 - nn_cut / 2;

    // Allocate FFTW complex arrays for forward and backward transforms
    fftw_complex *gal_data, *psf_data, *conv_data, *cropped_data;
    fftw_plan p_fwd_gal, p_fwd_psf, p_bwd_conv;

    gal_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gal_nn * gal_nn);
    psf_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gal_nn * gal_nn); // Padding to galaxy size
    conv_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * gal_nn * gal_nn);
    cropped_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nn_cut * nn_cut);

    // Plans for forward and inverse FFTs
    p_fwd_gal = fftw_plan_dft_2d(gal_nn, gal_nn, gal_data, gal_data, FFTW_FORWARD, FFTW_ESTIMATE);
    p_fwd_psf = fftw_plan_dft_2d(gal_nn, gal_nn, psf_data, psf_data, FFTW_FORWARD, FFTW_ESTIMATE);
    p_bwd_conv = fftw_plan_dft_2d(nn_cut, nn_cut, cropped_data, cropped_data, FFTW_BACKWARD, FFTW_ESTIMATE);

    // Copy galaxy data into the FFTW input
    #pragma omp parallel for
    for (int i = 0; i < gal_nn; ++i) {
        for (int j = 0; j < gal_nn; ++j) {
            int index = i * gal_nn + j;
            gal_data[index][0] = gal_fluxes[index]; // Real part
            gal_data[index][1] = 0; // Imaginary part is zero
        }
    }

    // Zero-pad and copy PSF data into the FFTW input
    memset(psf_data, 0, sizeof(fftw_complex) * gal_nn * gal_nn); // Initialize with zeros
    #pragma omp parallel for
    for (int i = 0; i < psf_nn; ++i) {
        for (int j = 0; j < psf_nn; ++j) {
            int index = (i + pad_width) * gal_nn + (j + pad_width);
            psf_data[index][0] = psf_fluxes[i * psf_nn + j]; // Real part
            psf_data[index][1] = 0; // Imaginary part is zero
        }
    }

    // Execute forward FFTs
    fftw_execute(p_fwd_gal);
    fftw_execute(p_fwd_psf);

    // Multiply in frequency domain (convolution)
    #pragma omp parallel for
    for (int i = 0; i < gal_nn * gal_nn; ++i) {
        conv_data[i][0] = (gal_data[i][0] * psf_data[i][0] - gal_data[i][1] * psf_data[i][1]) / (gal_nn * gal_nn);
        conv_data[i][1] = (gal_data[i][0] * psf_data[i][1] + gal_data[i][1] * psf_data[i][0]) / (gal_nn * gal_nn);
    }

    // Crop in Fourier space (select a smaller region)
    #pragma omp parallel for
    for (int i = minN; i < maxN; ++i) {
        for (int j = minN; j < maxN; ++j) {
            int index_src = i * gal_nn + j;
            int index_dst = (i - minN) * nn_cut + (j - minN);
            cropped_data[index_dst][0] = conv_data[index_src][0];
            cropped_data[index_dst][1] = conv_data[index_src][1];
        }
    }

    // Execute inverse FFT on the cropped data
    fftw_execute(p_bwd_conv);

    // Output the real part of the result, normalized
    std::vector<double> real_prof(nn_cut * nn_cut);
    #pragma omp parallel for
    for (int i = 0; i < nn_cut; ++i) {
        for (int j = 0; j < nn_cut; ++j) {
            int index = i * nn_cut + j;
            // Normalize by the size of the transform to account for FFTW's unnormalized IFFT
            real_prof[index] = cropped_data[index][0] / (nn_cut * nn_cut);
        }
    }

    // Cleanup
    fftw_destroy_plan(p_fwd_gal);
    fftw_destroy_plan(p_fwd_psf);
    fftw_destroy_plan(p_bwd_conv);
    fftw_free(gal_data);
    fftw_free(psf_data);
    fftw_free(conv_data);
    fftw_free(cropped_data);

    return real_prof;
}

std::string version() {
    return galsim::version();
}

PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Get flux values at multiple x,y coordinates";
    m.def("getFluxVec", &getFluxVec, "Get flux values at multiple x,y coordinates",
          py::arg("profileObj"), py::arg("xCoords"), py::arg("yCoords"));
    m.def("version", &version, "Get version of the module");
    m.def("convolveImage", &convolveImage, "Convolve two images",
          py::arg("gal_obj"), py::arg("psf_obj"), py::arg("galx"), py::arg("galy"), 
          py::arg("psfx"), py::arg("psfy"), py::arg("rescale"));
}