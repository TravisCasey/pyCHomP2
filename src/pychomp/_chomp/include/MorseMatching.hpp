/// MorseMatching.hpp
/// Shaun Harker
/// 2018-02-23
/// MIT LICENSE

#include <memory>

#include "Complex.h"
#include "CubicalComplex.h"
#include "CubicalMorseMatching.h"
#include "GenericMorseMatching.h"
#include "GradedComplex.h"
#include "MorseMatching.h"

inline std::shared_ptr<MorseMatching> MorseMatching::compute_matching(
    std::shared_ptr<Complex> complex) {
  if (std::dynamic_pointer_cast<CubicalComplex>(complex)) {
    return std::make_shared<CubicalMorseMatching>(
        std::dynamic_pointer_cast<CubicalComplex>(complex));
  } else {
    return std::make_shared<GenericMorseMatching>(complex);
  }
}

inline std::shared_ptr<MorseMatching> MorseMatching::compute_matching(
    std::shared_ptr<GradedComplex> graded_complex, bool truncate = false,
    Integer max_grade = 0) {
  if (std::dynamic_pointer_cast<CubicalComplex>(graded_complex->complex())) {
    return std::make_shared<CubicalMorseMatching>(graded_complex, truncate,
                                                  max_grade);
  } else {
    return std::make_shared<GenericMorseMatching>(graded_complex, truncate,
                                                  max_grade);
  }
}
