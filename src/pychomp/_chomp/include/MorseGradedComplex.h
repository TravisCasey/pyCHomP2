/// MorseGradedComplex.h
/// Shaun Harker
/// 2017-07-20
/// MIT LICENSE

#pragma once

#include <memory>
#include <vector>

#include "Chain.h"
#include "GradedComplex.h"
#include "Integer.h"
#include "MorseComplex.h"
#include "MorseMatching.h"

/// MorseGradedComplex
inline std::shared_ptr<GradedComplex> MorseGradedComplex(
    std::shared_ptr<GradedComplex> base_graded_complex,
    std::shared_ptr<MorseMatching> matching) {
  std::shared_ptr<MorseComplex> complex(
      new MorseComplex(base_graded_complex->complex(), matching));

  // Convert indices of cells to compute new graded_complex mapping
  std::vector<Integer> graded_complex_mapping(complex->size());
  for (auto x : *complex) {
    Chain included = complex->include({x});
    graded_complex_mapping[x] = base_graded_complex->value(*included.begin());
  }

  return std::shared_ptr<GradedComplex>(new GradedComplex(
      complex, [=](Integer x) { return graded_complex_mapping[x]; }));
}

/// MorseGradedComplex
inline std::shared_ptr<GradedComplex> MorseGradedComplex(
    std::shared_ptr<GradedComplex> base_graded_complex, Integer match_dim = -1,
    bool truncate = false, Integer max_grade = 0, bool verbose = false) {
  // Compute matching first, then pass to above
  std::shared_ptr<MorseMatching> matching(MorseMatching::compute_matching(
      base_graded_complex, match_dim, truncate, max_grade, verbose));

  return MorseGradedComplex(base_graded_complex, matching);
}

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void MorseGradedComplexBinding(py::module &m) {
  m.def("MorseGradedComplex",
        (std::shared_ptr<GradedComplex>(*)(std::shared_ptr<GradedComplex>,
                                           std::shared_ptr<MorseMatching>)) &
            MorseGradedComplex);
  m.def("MorseGradedComplex",
        (std::shared_ptr<GradedComplex>(*)(std::shared_ptr<GradedComplex>,
                                           Integer, bool, Integer, bool)) &
            MorseGradedComplex,
        py::arg("base"), py::arg("match_dim") = -1, py::arg("truncate") = false,
        py::arg("max_grade") = 0, py::arg("verbose") = false);
}
