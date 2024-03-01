#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <thread>
#include <vector>
#include "GalSim.h"

namespace py = pybind11;

// Adjust the function to accept a py::object instead of a direct reference
// to galsim::SBProfile. This allows passing Python objects.
std::vector<double> getFluxVec(py::object profileObj, 
                               const std::vector<double>& xCoords, 
                               const std::vector<double>& yCoords) {
    // Use py::cast to attempt to convert the Python object to a C++ galsim::SBProfile reference
    // This assumes that the Python object is a wrapper around a galsim::SBProfile or its subclass instance
    galsim::SBProfile& profile = py::cast<galsim::SBProfile&>(profileObj);

    // Get flux values at multiple x,y coordinates
    std::vector<double> fluxes;
    fluxes.reserve(xCoords.size());
    for(size_t i = 0; i < xCoords.size(); ++i) {
        double flux = profile.xValue(galsim::Position<double>(xCoords[i], yCoords[i]));
        fluxes.push_back(flux);
    }
    return fluxes;
}



std::string version() {
    return galsim::version();
}

// Binding code
PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Get flux values at multiple x,y coordinates";
    m.def("getFluxVec", &getFluxVec, "Get flux values at multiple x,y coordinates",
          py::arg("profileObj"), py::arg("xCoords"), py::arg("yCoords"));
    m.def("version", &version, "Get version of the module");
}