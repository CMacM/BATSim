// Test file
// Aim is to interface directly with GalSim C++ layer to get flux values
// at multiple x,y coordinates and pass them back to python as an array
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "GalSim.h"

namespace py = pybind11;

// Function to get flux values at multiple x,y coordinates
std::vector<double> getFluxVec(const galsim::SBProfile& prof, const std::vector<double>& x, 
                                const std::vector<double>& y)
{
    std::vector<double> fluxVec;
    for (int i = 0; i < x.size(); i++)
    {   
        fluxVec.push_back(prof.xValue(galsim::Position<double>(x[i], y[i])));
    }
    return fluxVec;
}

// Binding code
PYBIND11_MODULE(_gsinterface, m) {
    m.doc() = "Get flux values at multiple x,y coordinates";
    m.def("getFluxVec", &getFluxVec, "Get flux values at multiple x,y coordinates");
}