#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "GalSim.h"
#include <omp.h>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace py = pybind11;

// ---------- Utility helpers ----------

std::vector<double> rfftfreq(int n, double scale) {
    std::vector<double> result(n / 2 + 1);
    for (int i = 0; i <= n / 2; ++i) {
        result[i] = static_cast<double>(i) / (scale * n);
    }
    return result;
}

std::vector<double> fftfreq(int n, double scale) {
    std::vector<double> result(n);
    const double val = 1.0 / (n * scale);
    for (int i = 0; i < n; ++i) {
        result[i] = (i < (n + 1) / 2) ? i * val : (i - n) * val;
    }
    return result;
}

bool is_c_contiguous(const py::array& arr) {
    return (arr.flags() & py::array::c_style) != 0;
}

int round_up_multiple(int n, int m = 16) {
    if (m <= 0) throw std::runtime_error("round_up_multiple requires m > 0");
    return ((n + m - 1) / m) * m;
}

void copy_centered_into(
    const double* src,
    int src_dim,
    double* dst,
    int dst_dim
) {
    std::fill(dst, dst + dst_dim * dst_dim, 0.0);

    const int src_c = src_dim / 2;
    const int dst_c = dst_dim / 2;

    for (int y = 0; y < src_dim; ++y) {
        const int yy = y - src_c + dst_c;
        if (yy < 0 || yy >= dst_dim) continue;

        for (int x = 0; x < src_dim; ++x) {
            const int xx = x - src_c + dst_c;
            if (xx < 0 || xx >= dst_dim) continue;

            dst[yy * dst_dim + xx] = src[y * src_dim + x];
        }
    }
}

void ifftshift2d_inplace(double* arr, int dim) {
    std::vector<double> tmp(dim * dim);
    const int shift = dim / 2;  // matches your Stamp convention

    for (int y = 0; y < dim; ++y) {
        const int yy = (y + shift) % dim;
        for (int x = 0; x < dim; ++x) {
            const int xx = (x + shift) % dim;
            tmp[yy * dim + xx] = arr[y * dim + x];
        }
    }
    std::copy(tmp.begin(), tmp.end(), arr);
}

void fftshift2d_inplace(double* arr, int dim) {
    // For even dim this is identical to ifftshift.
    // For odd dim your Stamp convention still wants shift = dim/2.
    std::vector<double> tmp(dim * dim);
    const int shift = dim / 2;

    for (int y = 0; y < dim; ++y) {
        const int yy = (y + dim - shift) % dim;
        for (int x = 0; x < dim; ++x) {
            const int xx = (x + dim - shift) % dim;
            tmp[yy * dim + xx] = arr[y * dim + x];
        }
    }
    std::copy(tmp.begin(), tmp.end(), arr);
}

py::array_t<double> crop_centered_raw(
    const double* src,
    int src_dim,
    int out_dim,
    double norm = 1.0
) {
    auto result = py::array_t<double>({out_dim, out_dim});
    auto out = result.mutable_unchecked<2>();

    const int src_c = src_dim / 2;
    const int out_c = out_dim / 2;

    for (int y = 0; y < out_dim; ++y) {
        const int sy = y - out_c + src_c;
        for (int x = 0; x < out_dim; ++x) {
            const int sx = x - out_c + src_c;
            if (sx >= 0 && sx < src_dim && sy >= 0 && sy < src_dim) {
                out(y, x) = src[sy * src_dim + sx] * norm;
            } else {
                out(y, x) = 0.0;
            }
        }
    }

    return result;
}

inline double cubic_weight(double t) {
    t = std::abs(t);
    if (t < 1.0) {
        return 1.5 * t * t * t - 2.5 * t * t + 1.0;
    } else if (t < 2.0) {
        return -0.5 * t * t * t + 2.5 * t * t - 4.0 * t + 2.0;
    } else {
        return 0.0;
    }
}

inline double bicubic_sample_zero_padded(
    const double* src,
    int dim,
    double y,
    double x
) {
    const int x0 = static_cast<int>(std::floor(x));
    const int y0 = static_cast<int>(std::floor(y));

    auto at = [&](int yy, int xx) -> double {
        if (xx < 0 || xx >= dim || yy < 0 || yy >= dim) return 0.0;
        return src[yy * dim + xx];
    };

    double sum = 0.0;
    for (int j = -1; j <= 2; ++j) {
        const double wy = cubic_weight(y - (y0 + j));
        for (int i = -1; i <= 2; ++i) {
            const double wx = cubic_weight(x - (x0 + i));
            sum += at(y0 + j, x0 + i) * wy * wx;
        }
    }
    return sum;
}

// ---------- Flux sampling ----------

py::array_t<double> getFluxVec(
    const double scale,
    const galsim::SBProfile& gsobj,
    const py::array_t<double>& xy_coords
) {
    auto info = xy_coords.request();
    if (info.ndim != 2 || info.shape[0] != 2) {
        throw std::runtime_error("xy_coords must have shape (2, N)");
    }

    const int n_points = static_cast<int>(info.shape[1]);
    const int dim = static_cast<int>(std::llround(std::sqrt(static_cast<double>(n_points))));
    if (dim * dim != n_points) {
        throw std::runtime_error("xy_coords second dimension must be a perfect square");
    }

    auto xy = xy_coords.unchecked<2>();
    auto result = py::array_t<double>({dim, dim});
    auto out = result.mutable_data();

    const double pixel_area = scale * scale;

    gsobj.xValue(galsim::Position<double>(xy(0, 0), xy(1, 0)));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_points; ++i) {
        out[i] = gsobj.xValue(
            galsim::Position<double>(xy(0, i), xy(1, i))
        ) * pixel_area;
    }

    return result;
}

// ---------- Fine-grid convolution with internal additive padding ----------

py::array_t<double> convolvePsfFine(
    const double scale,
    const galsim::SBProfile& psf,
    const py::array_t<double>& gal_prof,
    const int pad_pixels = 16
) {
    if (!is_c_contiguous(gal_prof)) {
        throw std::runtime_error("gal_prof must be C-contiguous");
    }

    auto info = gal_prof.request();
    if (info.ndim != 2 || info.shape[0] != info.shape[1]) {
        throw std::runtime_error("gal_prof must be a square 2D array");
    }

    const int dim = static_cast<int>(info.shape[0]);
    const double* gal_in = static_cast<double*>(info.ptr);

    if (pad_pixels < 0) {
        throw std::runtime_error("pad_pixels must be >= 0");
    }

    const int padded_dim = round_up_multiple(dim + 2 * pad_pixels, 16);

    double* padded_real = fftw_alloc_real(padded_dim * padded_dim);
    fftw_complex* padded_k = fftw_alloc_complex(padded_dim * (padded_dim / 2 + 1));

    if (!padded_real || !padded_k) {
        if (padded_real) fftw_free(padded_real);
        if (padded_k) fftw_free(padded_k);
        throw std::runtime_error("FFTW allocation failed");
    }

    copy_centered_into(gal_in, dim, padded_real, padded_dim);

    ifftshift2d_inplace(padded_real, padded_dim);

    fftw_plan p_forward = fftw_plan_dft_r2c_2d(
        padded_dim, padded_dim, padded_real, padded_k, FFTW_ESTIMATE
    );
    fftw_plan p_backward = fftw_plan_dft_c2r_2d(
        padded_dim, padded_dim, padded_k, padded_real, FFTW_ESTIMATE
    );

    if (!p_forward || !p_backward) {
        if (p_forward) fftw_destroy_plan(p_forward);
        if (p_backward) fftw_destroy_plan(p_backward);
        fftw_free(padded_k);
        fftw_free(padded_real);
        throw std::runtime_error("FFTW plan creation failed");
    }

    fftw_execute(p_forward);

    const auto kx = rfftfreq(padded_dim, scale);
    const auto ky = fftfreq(padded_dim, scale);

    psf.kValue(galsim::Position<double>(2.0 * M_PI * kx[0], 2.0 * M_PI * ky[0]));

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < padded_dim; ++y) {
        for (int x = 0; x < (padded_dim / 2 + 1); ++x) {
            const int idx = y * (padded_dim / 2 + 1) + x;

            const std::complex<double> img_k(
                padded_k[idx][0], padded_k[idx][1]
            );

            const std::complex<double> psf_k = psf.kValue(
                galsim::Position<double>(
                    2.0 * M_PI * kx[x],
                    2.0 * M_PI * ky[y]
                )
            );

            const std::complex<double> prod = img_k * psf_k;
            padded_k[idx][0] = prod.real();
            padded_k[idx][1] = prod.imag();
        }
    }

    fftw_execute(p_backward);

    fftshift2d_inplace(padded_real, padded_dim);

    const double inv_norm =
        1.0 / static_cast<double>(padded_dim) / static_cast<double>(padded_dim);

    auto result = crop_centered_raw(padded_real, padded_dim, dim, inv_norm);

    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(padded_k);
    fftw_free(padded_real);

    return result;
}

// ---------- Output-grid resampling ----------

py::array_t<double> resampleToGrid(
    const py::array_t<double>& image,
    const double in_scale,
    const double out_scale,
    const int out_dim
) {
    if (!is_c_contiguous(image)) {
        throw std::runtime_error("image must be C-contiguous");
    }

    auto info = image.request();
    if (info.ndim != 2 || info.shape[0] != info.shape[1]) {
        throw std::runtime_error("image must be a square 2D array");
    }
    if (in_scale <= 0.0 || out_scale <= 0.0) {
        throw std::runtime_error("in_scale and out_scale must be > 0");
    }
    if (out_dim <= 0) {
        throw std::runtime_error("out_dim must be > 0");
    }

    const int in_dim = static_cast<int>(info.shape[0]);
    const double* src = static_cast<double*>(info.ptr);

    auto result = py::array_t<double>({out_dim, out_dim});
    auto out = result.mutable_unchecked<2>();

    const double ratio = out_scale / in_scale;
    const double src_c = 0.5 * static_cast<double>(in_dim - 1);
    const double dst_c = 0.5 * static_cast<double>(out_dim - 1);

    const double area_ratio = (out_scale * out_scale) / (in_scale * in_scale);

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < out_dim; ++y) {
        const double sy = src_c + (static_cast<double>(y) - dst_c) * ratio;
        for (int x = 0; x < out_dim; ++x) {
            const double sx = src_c + (static_cast<double>(x) - dst_c) * ratio;
            out(y, x) = bicubic_sample_zero_padded(src, in_dim, sy, sx) * area_ratio;
        }
    }

    return result;
}

// ---------- Coarse-grid integration to detector pixels ----------
// This is equivalent to convolving with a pixel response

py::array_t<double> integrateToDetector(
    const py::array_t<double>& fine_image,
    const int downsample_ratio
) {
    if (!is_c_contiguous(fine_image)) {
        throw std::runtime_error("fine_image must be C-contiguous");
    }

    auto info = fine_image.request();
    if (info.ndim != 2 || info.shape[0] != info.shape[1]) {
        throw std::runtime_error("fine_image must be a square 2D array");
    }
    if (downsample_ratio <= 0) {
        throw std::runtime_error("downsample_ratio must be > 0");
    }

    const int fine_dim = static_cast<int>(info.shape[0]);
    if (fine_dim % downsample_ratio != 0) {
        throw std::runtime_error(
            "fine_image dimension must be divisible by downsample_ratio"
        );
    }

    const int coarse_dim = fine_dim / downsample_ratio;
    const double* src = static_cast<double*>(info.ptr);

    auto result = py::array_t<double>({coarse_dim, coarse_dim});
    auto out = result.mutable_unchecked<2>();

    #pragma omp parallel for schedule(static)
    for (int cy = 0; cy < coarse_dim; ++cy) {
        for (int cx = 0; cx < coarse_dim; ++cx) {
            double sum = 0.0;
            const int y0 = cy * downsample_ratio;
            const int x0 = cx * downsample_ratio;

            for (int fy = 0; fy < downsample_ratio; ++fy) {
                const int y = y0 + fy;
                const int row = y * fine_dim;
                for (int fx = 0; fx < downsample_ratio; ++fx) {
                    const int x = x0 + fx;
                    sum += src[row + x];
                }
            }
            out(cy, cx) = sum;
        }
    }

    return result;
}

// ---------- Final crop/pad ----------

py::array_t<double> centerCropOrPad(
    const py::array_t<double>& image,
    const int out_dim
) {
    if (!is_c_contiguous(image)) {
        throw std::runtime_error("image must be C-contiguous");
    }

    auto info = image.request();
    if (info.ndim != 2 || info.shape[0] != info.shape[1]) {
        throw std::runtime_error("image must be a square 2D array");
    }
    if (out_dim <= 0) {
        throw std::runtime_error("out_dim must be > 0");
    }

    const int in_dim = static_cast<int>(info.shape[0]);
    const double* src = static_cast<double*>(info.ptr);

    return crop_centered_raw(src, in_dim, out_dim, 1.0);
}

PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Pybind11 interface for GalSim fine-grid sampling, convolution, and resampling";

    m.def(
        "getFluxVec",
        &getFluxVec,
        "Sample flux values on a provided coordinate grid",
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("xy_coords")
    );

    m.def(
        "convolvePsfFine",
        &convolvePsfFine,
        "Convolve a fine-grid image with a GalSim profile on the same fine grid",
        py::arg("scale"),
        py::arg("psf"),
        py::arg("gal_prof"),
        py::arg("pad_pixels") = 16
    );

    m.def(
        "resampleToGrid",
        &resampleToGrid,
        "Resample an image from in_scale to out_scale on a centered output grid",
        py::arg("image"),
        py::arg("in_scale"),
        py::arg("out_scale"),
        py::arg("out_dim")
    );

    m.def(
        "integrateToDetector",
        &integrateToDetector,
        "Integrate a fine-grid image onto a detector grid by block summation",
        py::arg("fine_image"),
        py::arg("downsample_ratio")
    );

    m.def(
        "centerCropOrPad",
        &centerCropOrPad,
        "Center-crop or zero-pad an image to out_dim x out_dim",
        py::arg("image"),
        py::arg("out_dim")
    );
}