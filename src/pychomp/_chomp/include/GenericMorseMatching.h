/// GenericMorseMatching.h
/// Shaun Harker
/// 2018-02-16
/// MIT LICENSE

#pragma once

#include <memory>
#include <ostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "Chain.h"
#include "Complex.h"
#include "GradedComplex.h"
#include "Integer.h"
#include "MorseMatching.h"

class GenericMorseMatching : public MorseMatching {
 public:
  /// GenericMorseMatching
  ///   Delegating constructor for non-graded complex
  ///   Assigns trivial grading (all 0)
  GenericMorseMatching(std::shared_ptr<Complex> complex_ptr,
                       bool verbose = false)
      : GenericMorseMatching(std::make_shared<GradedComplex>(
                                 complex_ptr, [](Integer i) { return 0; }),
                             verbose) {}

  /// GenericMorseMatching
  GenericMorseMatching(std::shared_ptr<GradedComplex> graded_complex_ptr,
                       bool truncate = false, Integer max_grade = 0,
                       bool verbose = false) {
    GradedComplex const& graded_complex = *graded_complex_ptr;
    Complex const& complex = *graded_complex.complex();
    Integer N = complex.size();

    if (verbose) {
      std::cout << "Generic Morse Matching on " << N << " cells.";
      std::cout << std::endl;
    }

    mate_.resize(N, -1);
    priority_.resize(N);
    Integer num_processed = 0;
    std::vector<Integer> boundary_count(N);
    std::unordered_set<Integer> coreducible;
    std::unordered_set<Integer> ace_candidates;

    auto bd = [&](Integer x) {
      Chain result;
      auto x_val = graded_complex.value(x);
      for (auto y : complex.boundary({x})) {
        auto y_val = graded_complex.value(y);
        if (y_val > x_val) {
          throw std::logic_error("graded_complex closure property failed.");
        }
        if (x_val == y_val) result += y;
      }
      return result;
    };

    auto cbd = [&](Integer x) {
      Chain result;
      auto x_val = graded_complex.value(x);
      for (auto y : complex.coboundary({x})) {
        if (x_val == graded_complex.value(y)) result += y;
      }
      return result;
    };

    if (verbose) {
      std::cout << "Initializing..." << std::endl;
      operations_ = N;
      bar_prev_ = -1;
      num_processed = 0;
    }

    Integer M = 0;
    for (auto x : complex) {
      if (not truncate || graded_complex.value(x) <= max_grade) {
        ++M;
        boundary_count[x] = bd(x).size();
        switch (boundary_count[x]) {
          case 0:
            ace_candidates.insert(x);
            break;
          case 1:
            coreducible.insert(x);
            break;
        }
      }
      if (verbose) progress_(++num_processed);
    }

    auto process = [&](Integer y) {
      priority_[y] = graded_complex.value(y) * complex.size() + num_processed++;
      coreducible.erase(y);
      ace_candidates.erase(y);
      for (auto x : cbd(y)) {
        boundary_count[x] -= 1;
        switch (boundary_count[x]) {
          case 0:
            coreducible.erase(x);
            ace_candidates.insert(x);
            break;
          case 1:
            coreducible.insert(x);
            break;
        }
      }
    };

    if (verbose) {
      std::cout << "Matching cells..." << std::endl;
      operations_ = M;
      bar_prev_ = -1;
      num_processed = 0;
    }

    while (num_processed < M) {
      if (not coreducible.empty()) {
        Integer K, Q;

        // Extract K
        auto it = coreducible.begin();
        K = *it;
        coreducible.erase(it);

        // Find mate Q
        for (auto x : bd(K))
          if (mate_[x] == -1) {
            Q = x;
            break;
          }

        mate_[K] = Q;
        mate_[Q] = K;
        process(Q);
        process(K);
      } else {
        Integer A;
        auto it = ace_candidates.begin();
        A = *it;
        ace_candidates.erase(it);
        mate_[A] = A;
        process(A);
      }
      if (verbose) progress_(num_processed);
    }

    if (verbose) {
      std::cout << "Computing critical cells..." << std::endl;
      operations_ = N;
      bar_prev_ = -1;
      num_processed = 0;
    }

    // Compute critical cells
    Integer D = complex.dimension();
    begin_.resize(D + 2);
    Integer idx = 0;
    for (Integer d = 0; d <= D; ++d) {
      begin_[d] = idx;
      for (auto v : complex(d)) {
        if (not truncate || graded_complex.value(v) <= max_grade) {
          if (mate(v) == v) {
            reindex_.push_back({v, idx});
            ++idx;
          }
        }
        if (verbose) progress_(++num_processed);
      }
    }
    begin_[D + 1] = idx;

    if (verbose && N != 0) {
      std::cout << "Reduced to " << idx << " critical cells, a reduction of ";
      std::cout << (100 - (100 * idx) / N) << "%." << std::endl << std::endl;
    }
  }

  /// critical_cells
  std::pair<BeginType const&, ReindexType const&> critical_cells(void) const {
    return {begin_, reindex_};
  }

  /// mate
  Integer mate(Integer x) const { return mate_[x]; }

  /// priority
  Integer priority(Integer x) const { return priority_[x]; }

 private:
  std::vector<Integer> mate_;
  std::vector<Integer> priority_;
  BeginType begin_;
  ReindexType reindex_;
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void GenericMorseMatchingBinding(py::module& m) {
  py::class_<GenericMorseMatching, std::shared_ptr<GenericMorseMatching>>(
      m, "GenericMorseMatching")
      .def(py::init<std::shared_ptr<Complex>, bool>(), py::arg("base"),
           py::arg("verbose") = false)
      .def(py::init<std::shared_ptr<GradedComplex>, bool, Integer, bool>(),
           py::arg("base"), py::arg("truncate") = false,
           py::arg("max_grade") = 0, py::arg("verbose") = false)
      .def("mate", &GenericMorseMatching::mate)
      .def("priority", &GenericMorseMatching::priority);
}
