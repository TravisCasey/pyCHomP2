/// Grading.h
/// Shaun Harker
/// 2018-03-09
/// MIT LICENSE

#pragma once

#include <unordered_set>

#include "Complex.h"
#include "CubicalComplex.h"
#include "Integer.h"
#include "common.h"

/// construct_grading
///   Define a grading on the complex `c` subject to values on the
///   top-dimensional cells.
std::function<Integer(Integer)> construct_grading(
    std::shared_ptr<Complex> c,
    std::function<Integer(Integer)> top_cell_grading) {
  // Copy top_cell_grading (with offset)
  std::vector<Integer> top_cell_grading_;
  top_cell_grading_.resize(c->size(c->dimension()));
  Integer num_nontop_cells_ = c->size() - c->size(c->dimension());
  for (auto v : (*c)(c->dimension())) {
    top_cell_grading_[v - num_nontop_cells_] = top_cell_grading(v);
  }

  // Returned grading function calculates minimum value in the top-dimensional
  // star of `x`.
  return [=](Integer x) {
    Integer min_value = -1;

    for (auto v : c->topstar(x)) {
      auto new_val = top_cell_grading_[v - num_nontop_cells_];

      if (min_value == -1) {
        min_value = new_val;
      } else {
        min_value = std::min(min_value, new_val);
      }
    }

    return min_value;
  };
}

/// inclusion_grading
///   Define a grading based on inclusion in `included`. All cells of `c` in
///   the closure of `included` are graded 0; others are graded 1.
std::function<Integer(Integer)> inclusion_grading(
    std::shared_ptr<Complex> c, std::unordered_set<Integer> const& included) {
  auto included_closure = c->closure(included);
  return [=](Integer x) {
    if (included_closure.count(x)) return 0;
    return 1;
  };
}

/// cubical_nerve
///   Define a grading on a cubical complex `c` selecting those cells which have
///   all their vertices' positions in `positions`, up to dimension `max_dim`.
std::function<Integer(Integer)> cubical_nerve(
    std::shared_ptr<CubicalComplex> c,
    std::unordered_set<Integer> const& positions, Integer max_dim = -1) {
  CubicalComplex complex = *c;
  if (max_dim == -1) max_dim = complex.dimension();
  Integer vertex_count = complex.size(0);

  return [=](Integer x) {
    if (complex.cell_dim(x) > max_dim) return 1;
    for (Integer y : complex.closure({x})) {
      // y < vertex_count iff y is a 0-cell
      if (y < vertex_count && positions.count(complex.cell_pos(y)) == 0) {
        return 1;
      }
    }
    return 0;
  };
}

/// Python Bindings

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

inline void GradingBinding(py::module& m) {
  m.def("construct_grading", &construct_grading);
  m.def("inclusion_grading", &inclusion_grading);
  m.def("cubical_nerve", &cubical_nerve, py::arg("complex"),
        py::arg("positions"), py::arg("max_dim") = -1);
}
