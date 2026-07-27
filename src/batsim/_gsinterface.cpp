#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/conduit/pybind11_conduit_v1.h>
#include "GalSim.h"
#include <omp.h>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace py = pybind11;

template <typename CoordT, typename FluxT>
py::array_t<FluxT> getFluxVecTyped(
    const double scale,
    const galsim::SBProfile& gsobj,
    const py::array_t<CoordT, py::array::c_style | py::array::forcecast>& xy_coords
) {
    auto xy = xy_coords.template unchecked<2>();

    if (xy_coords.ndim() != 2 || xy_coords.shape(0) != 2) {
        throw std::runtime_error("xy_coords must have shape (2, n_points)");
    }

    const py::ssize_t n_points = xy_coords.shape(1);
    const int dim = static_cast<int>(std::sqrt(static_cast<double>(n_points)));
    const int n_used = dim * dim;

    if (n_used != n_points) {
        throw std::runtime_error("xy_coords.shape[1] must be a perfect square");
    }

    auto result = py::array_t<FluxT, py::array::c_style>({dim, dim});
    auto out = result.mutable_data();

    const double area = scale * scale;

    // Pre-warm GalSim's internal cache with a single serial call.
    gsobj.xValue(
        galsim::Position<double>(
            static_cast<double>(xy(0, 0)),
            static_cast<double>(xy(1, 0))
        )
    );

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n_used; ++i) {
        const double x = static_cast<double>(xy(0, i));
        const double y = static_cast<double>(xy(1, i));

        const double flux = gsobj.xValue(
            galsim::Position<double>(x, y)
        ) * area;
        out[i] = static_cast<FluxT>(flux);
    }

    return result;
}

template <typename ComplexT>
py::array_t<ComplexT> getPsfKValueTyped(
    const double scale,
    const galsim::SBProfile& gsobj,
    const int n
) {
    if (n <= 0) {
        throw std::runtime_error("n must be positive");
    }

    const int nkx = n / 2 + 1;
    auto result = py::array_t<ComplexT>({n, nkx});
    auto out = result.mutable_data();

    const double dk = 2.0 * M_PI / (static_cast<double>(n) * scale);

    // Pre-warm GalSim's internal cache with a single serial call.
    gsobj.kValue(galsim::Position<double>(0.0, 0.0));

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < n; ++y) {
        const int ky_index = (y < (n + 1) / 2) ? y : y - n;
        const double ky = dk * static_cast<double>(ky_index);

        for (int x = 0; x < nkx; ++x) {
            const double kx = dk * static_cast<double>(x);
            const std::complex<double> value = gsobj.kValue(
                galsim::Position<double>(kx, ky)
            );

            out[y * nkx + x] = ComplexT(
                static_cast<typename ComplexT::value_type>(value.real()),
                static_cast<typename ComplexT::value_type>(value.imag())
            );
        }
    }

    return result;
}

PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Pybind11 interface for GalSim flux sampling";
    m.def(
        "getFluxVec",
        &getFluxVecTyped<double, double>,
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("xy_coords"),
        "Sample galaxy flux using float64 coordinates and output."
    );
    m.def(
        "getFluxVec64",
        &getFluxVecTyped<double, double>,
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("xy_coords"),
        "Sample galaxy flux using float64 coordinates and output."
    );
    m.def(
        "getFluxVec32",
        &getFluxVecTyped<float, float>,
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("xy_coords"),
        "Sample galaxy flux using float32 coordinates and output."
    );
    m.def(
        "getPsfKValue64",
        &getPsfKValueTyped<std::complex<double>>,
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("n"),
        "Sample PSF kValue on a float64/complex128 rFFT frequency grid."
    );
    m.def(
        "getPsfKValue32",
        &getPsfKValueTyped<std::complex<float>>,
        py::arg("scale"),
        py::arg("gsobj"),
        py::arg("n"),
        "Sample PSF kValue on a float32/complex64 rFFT frequency grid."
    );
}
