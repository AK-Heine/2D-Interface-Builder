#include "Python.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "backend.h"

namespace py = pybind11;
using float_2d_vec = std::vector<std::vector<float>>;

PYBIND11_MODULE(backend, m)
{
    m.doc() = "backend c++ implementation"; // optional module docstring
    py::bind_vector<std::vector<float>>(m, "FloatVector");
    py::bind_vector<float_2d_vec>(m, "FloatVector2D");
    m.def("find_coincidence", &find_coincidence, "A function that finds coincidence points.");
    m.def("backend_routine", &backend_routine, "Full C++ routine to generate coincidence matrices.");
}
