#include <pybind11/pybind11.h>

#include "vect2rast/vect2rast.hpp"

PYBIND11_MODULE(vect2rast, m) {
  m.doc() = "Python bindings to the vect2rast library";
  m.attr("__author__") = pybind11::cast(vect2rast::author());
  m.attr("__version__") = pybind11::cast(vect2rast::version());
  m.def("return_one", &vect2rast::return_one, "Return one in all circumstances.");
}
