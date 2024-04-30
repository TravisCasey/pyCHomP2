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
                       Integer match_dim = -1, bool verbose = false)
      : GenericMorseMatching(std::make_shared<GradedComplex>(
                                 complex_ptr, [](Integer i) { return 0; }),
                             match_dim, false, 0, verbose) {}

  /// GenericMorseMatching
  ///   Computes Morse matching based on the mate algorithm.
  ///   Setting `truncate` flag stops matching on cells with grade exceeding
  ///   `max_grade`.
  ///   Setting match_dim stops matching at specified dimension; homology up to
  ///   dimension match_dim - 1 will be correct. Default value -1 uses the full
  ///   complex.

  // TODO change complex ptr behavior to match that of CubicalMorseMatching
  // TODO current implementation matches down K -> Q, which prevents allowing
  //   cells of match_dim match with cells of match_dim + 1 without having to
  //   process those cells. Switching this should allow for that improvement;
  //   However, this does not matter if it is run through CubicalMorseMatching
  //   first as the higher cells are removed anyway.
  GenericMorseMatching(std::shared_ptr<GradedComplex> graded_complex_ptr,
                       Integer match_dim = -1, bool truncate = false,
                       Integer max_grade = 0, bool verbose = false) {
    GradedComplex const& graded_complex = *graded_complex_ptr;
    Complex const& complex = *graded_complex.complex();

    // If match_dim is invalid (also default value of -1), set D to
    // dimension of the complex instead.
    Integer D = (match_dim > 0 && match_dim <= complex.dimension())
                    ? match_dim
                    : complex.dimension();

    // N is the number of cells of dimension D and below
    // Cells >= top_begin are of top dimenstion w.r.t. D,
    // and should not be matched up.
    Range top_range = complex(D);
    Integer N = *(top_range.end());
    Integer top_begin = *(top_range.begin());

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
      Integer x_val = graded_complex.value(x);
      for (Integer y : complex.boundary({x})) {
        Integer y_val = graded_complex.value(y);
        if (y_val > x_val) {
          throw std::logic_error("graded_complex closure property failed.");
        }
        if (x_val == y_val) result += y;
      }
      return result;
    };

    auto cbd = [&](Integer x) {
      Chain result;
      // Compute coboundary only up to dimension D
      if (x >= top_begin) return result;

      Integer x_val = graded_complex.value(x);
      for (Integer y : complex.coboundary({x})) {
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
    for (Integer x = 0; x < N; ++x) {
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
      priority_[y] = graded_complex.value(y) * M + num_processed++;
      coreducible.erase(y);
      ace_candidates.erase(y);
      for (Integer x : cbd(y)) {
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

    Integer A, K, Q;
    while (num_processed < M) {
      if (not coreducible.empty()) {
        // Extract K
        auto it = coreducible.begin();
        K = *it;
        coreducible.erase(it);

        // Find mate Q
        for (Integer x : bd(K))
          if (mate_[x] == -1) {
            Q = x;
            break;
          }

        mate_[K] = Q;
        mate_[Q] = K;
        process(Q);
        process(K);
      } else {
        // Extract A
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
    begin_.resize(D + 2);
    Integer idx = 0;
    for (Integer d = 0; d <= D; ++d) {
      begin_[d] = idx;
      for (Integer v : complex(d)) {
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
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void GenericMorseMatchingBinding(py::module& m) {
  py::class_<GenericMorseMatching, std::shared_ptr<GenericMorseMatching>>(
      m, "GenericMorseMatching")
      .def(py::init<std::shared_ptr<Complex>, Integer, bool>(), py::arg("base"),
           py::arg("match_dim") = -1, py::arg("verbose") = false)
      .def(py::init<std::shared_ptr<GradedComplex>, Integer, bool, Integer,
                    bool>(),
           py::arg("base"), py::arg("match_dim") = -1,
           py::arg("truncate") = false, py::arg("max_grade") = 0,
           py::arg("verbose") = false)
      .def("mate", &GenericMorseMatching::mate)
      .def("priority", &GenericMorseMatching::priority);
}
