#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include <numpy/arrayobject.h>
#include <string>
#include "PROconfig.h"

// expected by PROfit for printing stuff
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

template<class Obj>
class PyWrapper {
  public:
    Obj *o;

    PyWrapper():
      o(NULL) {}

    PyWrapper(Obj *o):
      o(o) {}

    ~PyWrapper() {
      if (o) {
        delete o;
      }
      o = NULL;
    };
};

PyWrapper<PROfit::PROconfig> *_init_PROconfig(const std::string &xml) {
  return new PyWrapper(new PROfit::PROconfig(xml));
}

PyWrapper<PROfit::PROpeller> *_init_PROpeller() {
  return new PyWrapper(new PROfit::PROpeller());
}

int _add(int a, int b) {
  return a + b;
}

namespace py = pybind11;

PYBIND11_MODULE(profit, m) {
    m.doc() =  "Python interface for core functionality of PROfit fitting library.";
    m.def("add", &_add, "Add two numbers");
    py::class_<PyWrapper<PROfit::PROconfig>>(m, "PyPROconfig")
        .def(py::init<>());

    m.def("init_PROconfig", &_init_PROconfig, "Initialize PROconfig object");
    m.def("init_PROpeller", &_init_PROpeller, "Initialize PROpeller object");
}
