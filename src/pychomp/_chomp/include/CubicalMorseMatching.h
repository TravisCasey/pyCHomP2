/// CubicalMorseMatching.h
/// Shaun Harker
/// 2018-02-16
/// MIT LICENSE

#pragma once

#include <memory>
#include <ostream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

#include "GradedComplex.h"
#include "Integer.h"
#include "MorseMatching.h"

class CubicalMorseMatching : public MorseMatching {
 public:
  /// CubicalMorseMatching
  ///   Delegating constructor for non-graded complex
  ///   Assigns trivial grading (all 0)
  CubicalMorseMatching(std::shared_ptr<CubicalComplex> complex_ptr,
                       bool verbose = false)
      : CubicalMorseMatching(std::make_shared<GradedComplex>(
                                 complex_ptr, [](Integer i) { return 0; }),
                             verbose) {}

  /// CubicalMorseMatching
  ///   Computes Morse matching based on hypercube templates.
  ///   Setting `truncate` flag stops matching on cells with grade exceeding
  ///   `max_grade`.
  CubicalMorseMatching(std::shared_ptr<GradedComplex> graded_complex_ptr,
                       bool truncate = false, Integer max_grade = 0,
                       bool verbose = false)
      : graded_complex_(graded_complex_ptr) {
    complex_ =
        std::dynamic_pointer_cast<CubicalComplex>(graded_complex_->complex());
    if (not complex_)
      throw std::invalid_argument(
          "CubicalMorseMatching must be constructed with a CubicalComplex");

    type_size_ = complex_->type_size();
    Integer N = complex_->size();
    Integer D = complex_->dimension();
    Integer num_processed;

    if (verbose) {
      std::cout << "Cubical Morse Matching on " << N << " cells.";
      std::cout << std::endl;
      operations_ = N;
      bar_prev_ = -1;
      num_processed = 0;
    }

    // In matching construction, mates are only found from queens to kings.
    // Found kings are cached in the prev_kings and next_kings sets.
    Integer cell_mate;
    std::unique_ptr<std::unordered_set<Integer>> prev_kings;
    std::unique_ptr<std::unordered_set<Integer>> next_kings =
        std::make_unique<std::unordered_set<Integer>>();

    Integer idx = 0;
    begin_.resize(D + 2);
    for (Integer d = 0; d <= D; ++d) {
      begin_[d] = idx;

      // Transfer ownership of next_kings and initialize empty set
      prev_kings = std::move(next_kings);
      next_kings = std::make_unique<std::unordered_set<Integer>>();

      for (Integer v : (*complex_)(d)) {
        // Discard:
        //   1. Fringe cells
        //   2. Cells that are truncated by grade
        //   3. Kings from previous iteration
        if (!complex_->rightfringe(v) &&
            (!truncate || graded_complex_->value(v) <= max_grade) &&
            !prev_kings->count(v)) {
          // Find mate
          cell_mate = mate_(v, D, true);
          if (cell_mate == v) {
            // Ace
            reindex_.push_back({v, idx});
            ++idx;
          } else {
            // Queen - cahce result for next iteration
            next_kings->insert(cell_mate);
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
  std::pair<BeginType const &, ReindexType const &> critical_cells(void) const {
    return {begin_, reindex_};
  }

  /// mate
  Integer mate(Integer x) const { return mate_(x, complex_->dimension()); }

  /// priority
  Integer priority(Integer x) const { return type_size_ - x % type_size_; }

 private:
  Integer type_size_;
  std::shared_ptr<GradedComplex> graded_complex_;
  std::shared_ptr<CubicalComplex> complex_;

  Integer mate_(Integer cell, Integer D, bool initial = false) const {
    // The initial option is set to true when called within constructor.
    // This only proposes mates of higher dimensions (i.e. kings).

    if (complex_->rightfringe(cell)) return cell;
    Integer shape = complex_->cell_shape(cell);
    Integer type_offset;
    Integer proposed_mate;

    for (Integer d = 0, bit = 1; d < D; ++d, bit <<= 1L) {
      // If initial only search for kings.
      if (initial && shape & bit) continue;

      Integer type_offset = type_size_ * complex_->TS()[shape ^ bit];
      Integer proposed_mate = complex_->cell_pos(cell) + type_offset;

      // A proposed mate must:
      //   1. Be in the same grade
      //   2. Not be a fringe cell
      //   3. Not be matched with another cell
      if (!complex_->rightfringe(proposed_mate) &&
          graded_complex_->value(proposed_mate) ==
              graded_complex_->value(cell) &&
          proposed_mate == mate_(proposed_mate, d)) {
        return proposed_mate;  // found mate
      }
    }
    return cell;  // Did not mate => critical cell
  }
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void CubicalMorseMatchingBinding(py::module &m) {
  py::class_<CubicalMorseMatching, std::shared_ptr<CubicalMorseMatching>>(
      m, "CubicalMorseMatching")
      .def(py::init<std::shared_ptr<CubicalComplex>, bool>(), py::arg("base"),
           py::arg("verbose") = false)
      .def(py::init<std::shared_ptr<GradedComplex>, bool, Integer, bool>(),
           py::arg("base"), py::arg("truncate") = false,
           py::arg("max_grade") = 0, py::arg("verbose") = false)
      .def("mate", &CubicalMorseMatching::mate)
      .def("priority", &CubicalMorseMatching::priority);
}
