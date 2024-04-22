/// CubicalMorseMatching.h
/// Shaun Harker
/// 2018-02-16
/// MIT LICENSE

#pragma once

#include <memory>
#include <stdexcept>
#include <utility>

#include "GradedComplex.h"
#include "Integer.h"
#include "MorseMatching.h"

class CubicalMorseMatching : public MorseMatching {
 public:
  /// CubicalMorseMatching
  ///   Delegating constructor for non-graded complex
  ///   Assigns trivial grading (all 0)
  CubicalMorseMatching(std::shared_ptr<CubicalComplex> complex_ptr)
      : CubicalMorseMatching(std::make_shared<GradedComplex>(
            complex_ptr, [](Integer i) { return 0; })) {}

  /// CubicalMorseMatching
  ///   Computes Morse matching based on hypercube templates.
  ///   Setting `truncate` flag stops matching on cells with grade exceeding
  ///   `max_grade`.
  CubicalMorseMatching(std::shared_ptr<GradedComplex> graded_complex_ptr,
                       bool truncate = false, Integer max_grade = 0)
      : graded_complex_(graded_complex_ptr) {
    complex_ =
        std::dynamic_pointer_cast<CubicalComplex>(graded_complex_->complex());
    if (not complex_)
      throw std::invalid_argument(
          "CubicalMorseMatching must be constructed with a CubicalComplex");

    type_size_ = complex_->type_size();
    Integer D = complex_->dimension();

    Integer idx = 0;
    begin_.resize(D + 2);
    for (Integer d = 0; d <= D; ++d) {
      begin_[d] = idx;
      for (auto v : (*complex_)(d)) {
        if (truncate && graded_complex_->value(v) > max_grade) continue;
        if (!complex_->rightfringe(v) && mate(v) == v) {
          reindex_.push_back({v, idx});
          ++idx;
        }
      }
    }
    begin_[D + 1] = idx;
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
  BeginType begin_;
  ReindexType reindex_;

  Integer mate_(Integer cell, Integer D) const {
    if (complex_->rightfringe(cell)) return cell;  // Not strictly necessary
    Integer shape = complex_->cell_shape(cell);
    Integer position = complex_->cell_pos(cell);

    for (Integer d = 0, bit = 1; d < D; ++d, bit <<= 1L) {
      Integer type_offset = type_size_ * complex_->TS()[shape ^ bit];
      Integer proposed_mate = position + type_offset;

      // A proposed mate must:
      //   1. Be in the same grade
      //   2. Not be a fringe cell
      //   3. Not already be matched with another cell
      if (graded_complex_->value(proposed_mate) ==
              graded_complex_->value(cell) &&
          !complex_->rightfringe(proposed_mate) &&
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
      .def(py::init<std::shared_ptr<CubicalComplex>>(), py::arg("base"))
      .def(py::init<std::shared_ptr<GradedComplex>, bool, Integer>(),
           py::arg("base"), py::arg("truncate") = false,
           py::arg("max_grade") = 0)
      .def("mate", &CubicalMorseMatching::mate)
      .def("priority", &CubicalMorseMatching::priority);
}
