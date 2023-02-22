#include <pybind11/pybind11.h>
#include <pybind11/stl.h>   
#include <vector>
#include <tuple>
#include "graph.h"
#include "fasthare.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace pybind11::literals;


FastHare::fasthare_output fasthare_reduction(
    FastHare::ski sk_ising, 
    // If not empty, load the SK Hamiltonian graph from file
    string file = "", 
    double alpha = 1.0)
{
    if (!file.empty()) {
        FastHare fh(file, alpha);
        fh.run();
        return fh.get_output();
    } else  {
        FastHare fh(sk_ising, alpha);
        fh.run();
        return fh.get_output();
    }
}

namespace py = pybind11;

PYBIND11_MODULE(fasthare, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: fasthare

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("fasthare_reduction", &fasthare_reduction,
        py::arg("sk_ising") = FastHare::ski(),
        py::arg("file") = string(""), 
        py::arg("alpha") = 1.0, 
         R"pbdoc(
        Reduce/compress an SK Hamiltonian to obtain an equilvalent
        albeit smaller SK Hamiltonian
    )pbdoc");

    // m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
    //     Subtract two numbers

    //     Some other explanation about the subtract function.
    // )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
