/// MorseMatching.h
/// Shaun Harker
/// 2017-07-19
/// MIT LICENSE

#pragma once

#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "Complex.h"
#include "GradedComplex.h"
#include "Integer.h"

class MorseMatching {
 public:
  typedef std::vector<Integer>
      BeginType;  // To store location of first cell of each dim
  typedef std::vector<std::pair<Integer, Integer>>
      ReindexType;  // To convert indexing

  /// mate
  virtual Integer mate(Integer x) const = 0;

  /// priority
  virtual Integer priority(Integer x) const = 0;

  /// critical_cells
  virtual std::pair<BeginType const &, ReindexType const &> critical_cells(
      void) const = 0;

  /// compute_matching
  static std::shared_ptr<MorseMatching> compute_matching(
      std::shared_ptr<Complex> complex, bool verbose);

  /// compute_matching
  static std::shared_ptr<MorseMatching> compute_matching(
      std::shared_ptr<GradedComplex> graded_complex, bool truncate,
      Integer max_grade, bool verbose);

 protected:
  BeginType begin_;
  ReindexType reindex_;
  Integer operations_;
  const Integer bar_total_ = 50;
  Integer bar_prev_;

  /// progress_
  ///   Prints progress bar to the command line.
  ///   Matching class must set operations_ as the max number for processed.
  ///   Before each new progress bar, set bar_prev_ to -1.
  void progress_(Integer processed) {
    // bars to display
    Integer bars;
    if (operations_ != 0) {
      bars = processed * bar_total_ / operations_;
    } else {
      bars = bar_total_;
    }

    if (bars != bar_prev_) {
      bar_prev_ = bars;
      std::cout << "[";
      for (int i = 0; i < bar_total_; ++i) {
        if (i < bars) {
          std::cout << "=";
        } else if (i == bars) {
          std::cout << ">";
        } else {
          std::cout << " ";
        }
      }

      // progress percent
      std::cout << "] " << ((100 * bars) / bar_total_) << "%\r";
      std::cout.flush();
    }

    // completed progress bar
    if (bars == bar_total_) std::cout << std::endl;
    return;
  }
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void MorseMatchingBinding(py::module &m) {
  py::class_<MorseMatching, std::shared_ptr<MorseMatching>>(m, "MorseMatching")
      .def("mate", &MorseMatching::mate)
      .def("priority", &MorseMatching::priority);
}
