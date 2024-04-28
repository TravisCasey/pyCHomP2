/// MorseComplex.h
/// Shaun Harker
/// 2017-07-19
/// MIT LICENSE

#pragma once

#include <memory>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "Chain.h"
#include "Complex.h"
#include "Integer.h"
#include "Iterator.h"
#include "MorseMatching.h"

class MorseComplex : public Complex {
 public:
  /// MorseComplex constructor
  MorseComplex(std::shared_ptr<Complex> arg_base,
               std::shared_ptr<MorseMatching> arg_matching)
      : base_(arg_base), matching_(arg_matching) {
    auto begin_reindex = matching_->critical_cells();

    // Cell iterators
    begin_.clear();
    for (auto i : begin_reindex.first) {
      begin_.push_back(Iterator(i));
    }
    dim_ = begin_.size() - 2;

    // Include and project
    auto const& reindex = begin_reindex.second;
    for (auto pair : reindex) {
      include_.push_back(pair.first);
    }
    project_ =
        std::unordered_map<Integer, Integer>(reindex.begin(), reindex.end());

    // boundary
    bd_.resize(size());
    for (auto ace : *this) {
      bd_[ace] = lower(base()->boundary(include({ace})));
    }

    // coboundary
    cbd_.resize(size());
    for (auto ace : *this) {
      for (auto bd_cell : bd_[ace]) cbd_[bd_cell] += ace;
    }
  }

  /// delegating constructor for unmatched complexes
  MorseComplex(std::shared_ptr<Complex> arg_base, bool verbose = false)
      : MorseComplex(arg_base,
                     MorseMatching::compute_matching(arg_base, verbose)) {}

  /// boundary
  virtual Chain boundary(Chain const& c) const final {
    Chain result;
    for (auto x : c) result += bd_[x];
    return result;
  }

  /// coboundary
  virtual Chain coboundary(Chain const& c) const final {
    Chain result;
    for (auto x : c) result += cbd_[x];
    return result;
  }

  /// column
  ///   Apply "callback" method to every element in ith column of
  ///   boundary matrix
  virtual void column(
      Integer i, std::function<void(Integer)> const& callback) const final {
    for (auto x : bd_[i]) callback(x);
  };

  /// row
  ///   Apply "callback" method to every element in ith row of
  ///   boundary matrix
  virtual void row(Integer i,
                   std::function<void(Integer)> const& callback) const final {
    for (auto x : cbd_[i]) callback(x);
  };

  // Feature

  /// base
  std::shared_ptr<Complex> base(void) const { return base_; }

  /// matching
  std::shared_ptr<MorseMatching> matching(void) const { return matching_; }

  /// include
  Chain include(Chain const& c) {
    Chain result;
    for (auto x : c) result += include_[x];
    return result;
  }

  /// project
  Chain project(Chain const& c) {
    Chain result;
    for (auto x : c) {
      if (project_.count(x) > 0) result += project_[x];
    }
    return result;
  }

  /// lift
  Chain lift(Chain const& c) {
    Chain included = include(c);
    Chain canonical;
    Chain gamma;
    std::tie(canonical, gamma) = flow(base()->boundary(included));
    return included + gamma;
  }

  /// lower
  Chain lower(Chain const& c) {
    Chain canonical;
    Chain gamma;
    std::tie(canonical, gamma) = flow(c);
    return project(canonical);
  }

  /// flow
  std::pair<Chain, Chain> flow(Chain const& input) const {
    Chain canonical, gamma;
    auto isQueen = [&](Integer x) { return x < matching_->mate(x); };

    // Partial ordering on queens handled by priority queue
    auto compare = [&](Integer x, Integer y) {
      return matching_->priority(x) < matching_->priority(y);
    };
    std::priority_queue<Integer, std::vector<Integer>, decltype(compare)>
        priority(compare);

    auto process = [&](Integer x) {
      if (isQueen(x)) priority.push(x);
      canonical += x;
    };

    for (auto x : input) process(x);

    while (not priority.empty()) {
      auto queen = priority.top();
      priority.pop();
      if (canonical.count(queen) == 0) continue;
      auto king = matching_->mate(queen);
      gamma += king;
      base()->column(king, process);
    }

    return {canonical, gamma};
  }

  /// colift
  Chain colift(Chain const& c) {
    Chain included = include(c);
    Chain cocanonical;
    Chain cogamma;
    std::tie(cocanonical, cogamma) = coflow(base()->coboundary(included));
    return included + cogamma;
  }

  /// colower
  Chain colower(Chain const& c) {
    Chain cocanonical;
    Chain cogamma;
    std::tie(cocanonical, cogamma) = coflow(c);
    return project(cocanonical);
  }

  /// coflow
  ///   Dualization of flow exchanges role of kings and queens for their
  ///   cocell counterparts.
  std::pair<Chain, Chain> coflow(Chain const& input) const {
    Chain cocanonical, cogamma;
    auto isKing = [&](Integer x) { return x > matching_->mate(x); };

    // Partial ordering on kings handled by priority queue
    auto compare = [&](Integer x, Integer y) {
      return matching_->priority(x) > matching_->priority(y);
    };
    std::priority_queue<Integer, std::vector<Integer>, decltype(compare)>
        priority(compare);

    auto process = [&](Integer x) {
      if (isKing(x)) priority.push(x);
      cocanonical += x;
    };

    for (auto x : input) process(x);

    while (not priority.empty()) {
      auto king = priority.top();
      priority.pop();
      if (cocanonical.count(king) == 0) continue;
      auto queen = matching_->mate(king);
      cogamma += queen;
      base()->row(queen, process);
    }

    return {cocanonical, cogamma};
  }

 private:
  std::shared_ptr<Complex> base_;
  std::shared_ptr<MorseMatching> matching_;
  std::vector<Integer> include_;
  std::unordered_map<Integer, Integer> project_;
  std::vector<Chain> bd_;
  std::vector<Chain> cbd_;
};

/// Python Bindings

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

inline void MorseComplexBinding(py::module& m) {
  py::class_<MorseComplex, std::shared_ptr<MorseComplex>, Complex>(
      m, "MorseComplex")
      .def(py::init<std::shared_ptr<Complex>, std::shared_ptr<MorseMatching>>(),
           py::arg("base"), py::arg("matching"))
      .def(py::init<std::shared_ptr<Complex>, bool>(), py::arg("base"),
           py::arg("verbose") = false)
      .def("include", &MorseComplex::include)
      .def("project", &MorseComplex::project)
      .def("lift", &MorseComplex::lift)
      .def("lower", &MorseComplex::lower)
      .def("flow", &MorseComplex::flow)
      .def("colift", &MorseComplex::colift)
      .def("colower", &MorseComplex::colower)
      .def("coflow", &MorseComplex::coflow)
      .def("base", &MorseComplex::base)
      .def("matching", &MorseComplex::matching);
}
