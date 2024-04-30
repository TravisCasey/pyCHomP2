/// ConnectionMatrix.h
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#pragma once

#include <memory>
#include <vector>

#include "GradedComplex.h"
#include "Integer.h"
#include "MorseGradedComplex.h"

/// ConnectionMatrix
///   Computes the minimal morse complex obtained through a series of matchings.
///   Setting `truncate` flag disables matching of cells with grade exceeding
///   `max_grade`.
///   Setting match_dim stops matching at specified dimension; homology up to
///   dimension match_dim - 1 will be correct. Default value -1 uses the full
///   complex.
inline std::shared_ptr<GradedComplex> ConnectionMatrix(
    std::shared_ptr<GradedComplex> base, Integer match_dim = -1,
    bool truncate = false, Integer max_grade = 0, bool verbose = false) {
  std::shared_ptr<GradedComplex> next = base;
  do {
    base = next;
    next = MorseGradedComplex(base, match_dim, truncate, max_grade, verbose);
  } while (next->complex()->size() != base->complex()->size());
  return base;
}

/// ConnectionMatrixTower
///   Computes the minimal morse complex obtained through a series of matchings.
///   Returns the entire tower of intermediate complexes as a vector.
///   Setting `truncate` flag disables matching of cells with grade exceeding
///   `max_grade`.
///   Setting match_dim stops matching at specified dimension; homology up to
///   dimension match_dim - 1 will be correct. Default value -1 uses the full
///   complex.
inline std::vector<std::shared_ptr<GradedComplex>> ConnectionMatrixTower(
    std::shared_ptr<GradedComplex> base, Integer match_dim = -1,
    bool truncate = false, Integer max_grade = 0, bool verbose = false) {
  std::vector<std::shared_ptr<GradedComplex>> tower;
  std::shared_ptr<GradedComplex> next = base;
  do {
    tower.push_back(next);
    base = next;
    next = MorseGradedComplex(base, match_dim, truncate, max_grade, verbose);
  } while (next->complex()->size() != base->complex()->size());
  return tower;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void ConnectionMatrixBinding(py::module &m) {
  m.def("ConnectionMatrix", &ConnectionMatrix, py::arg("base"),
        py::arg("match_dim") = -1, py::arg("truncate") = false,
        py::arg("max_grade") = 0, py::arg("verbose") = false);
  m.def("ConnectionMatrixTower", &ConnectionMatrixTower, py::arg("base"),
        py::arg("match_dim") = -1, py::arg("truncate") = false,
        py::arg("max_grade") = 0, py::arg("verbose") = false);
}
