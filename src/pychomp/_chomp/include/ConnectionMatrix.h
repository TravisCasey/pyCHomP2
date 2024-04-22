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
inline std::shared_ptr<GradedComplex> ConnectionMatrix(
    std::shared_ptr<GradedComplex> base, bool truncate = false,
    Integer max_grade = 0) {
  std::shared_ptr<GradedComplex> next = base;
  do {
    base = next;
    next = MorseGradedComplex(base, truncate, max_grade);
  } while (next->complex()->size() != base->complex()->size());
  return base;
}

/// ConnectionMatrixTower
///   Computes the minimal morse complex obtained through a series of matchings.
///   Returns the entire tower of intermediate complexes as a vector.
///   Setting `truncate` flag disables matching of cells with grade exceeding
///   `max_grade`.
inline std::vector<std::shared_ptr<GradedComplex>> ConnectionMatrixTower(
    std::shared_ptr<GradedComplex> base, bool truncate = false,
    Integer max_grade = 0) {
  std::vector<std::shared_ptr<GradedComplex>> tower;
  std::shared_ptr<GradedComplex> next = base;
  do {
    tower.push_back(next);
    base = next;
    next = MorseGradedComplex(base, truncate, max_grade);
  } while (next->complex()->size() != base->complex()->size());
  return tower;
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void ConnectionMatrixBinding(py::module &m) {
  m.def("ConnectionMatrix", &ConnectionMatrix, py::arg("base"),
        py::arg("truncate") = false, py::arg("max_grade") = 0);
  m.def("ConnectionMatrixTower", &ConnectionMatrixTower, py::arg("base"),
        py::arg("truncate") = false, py::arg("max_grade") = 0);
}
