#include <pybind11/pybind11.h>

#include "evd.hpp"

#define EVD_READWRITE(pyname, cppname) \
    .def_readwrite(#pyname, &evdOptions::cppname)

#define EVD_TRIVIAL_READWRITE(name) EVD_READWRITE(name, name)

int evd_process(evdOptions *opts);

namespace py = pybind11;

PYBIND11_MODULE(evdlib, m) {

    py::class_<evdOptions>(m, "Evd")
        .def(py::init<>())
        EVD_TRIVIAL_READWRITE(inputDS)
        EVD_TRIVIAL_READWRITE(outputFolder)
        EVD_TRIVIAL_READWRITE(outputCompressedSlcFolder)
        EVD_TRIVIAL_READWRITE(compSlc)

        EVD_READWRITE(weightsDS, wtsDS)
        EVD_READWRITE(minimumNeighbors, minNeighbors)

        EVD_TRIVIAL_READWRITE(miniStackCount)
        EVD_TRIVIAL_READWRITE(blocksize)
        EVD_TRIVIAL_READWRITE(memsize)

        EVD_READWRITE(halfWindowX, Nx)
        EVD_READWRITE(halfWindowY, Ny)

        EVD_TRIVIAL_READWRITE(method)
        EVD_TRIVIAL_READWRITE(bandWidth)

        .def("print", &evdOptions::print)
        .def("run", [](evdOptions& self) {
            evd_process(&self);
        })
        ;
}
