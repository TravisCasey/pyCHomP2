/// Homology.h
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#pragma once

#include <memory>

#include "Chain.h"
#include "Complex.h"
#include "Integer.h"
#include "MorseComplex.h"
#include "MorseMatching.h"

/// Homology
inline std::shared_ptr<Complex> Homology(std::shared_ptr<Complex> base,
                                         Integer match_dim = -1,
                                         bool verbose = false) {
  std::shared_ptr<Complex> next = base;
  do {
    base = next;
    next.reset(new MorseComplex(base, match_dim));
  } while (next->size() != base->size());
  return base;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void HomologyBinding(py::module &m) {
  m.def("Homology", &Homology, py::arg("base"), py::arg("match_dim") = -1,
        py::arg("verbose") = false);
}
