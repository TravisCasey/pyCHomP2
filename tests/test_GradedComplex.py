"""Testing script for GradedComplex.h"""

from pychomp import *

def test_cubical_graded_complex():
    X = CubicalComplex([3, 4])

    def top_grading(cell):
        if cell in {36, 37, 38}:
            return 0
        if cell in {39, 40, 41}:
            return 1
        if cell in {42, 43, 44}:
            return 2
        if cell in {45, 46, 47}:
            return 3
        return 0

    grading = construct_grading(X, top_grading)
    gradX = GradedComplex(X, grading)

    assert gradX.complex().size() == X.size()
    assert gradX.complex().dimension() == X.dimension()

    assert gradX.count()[0] == [7, 10, 3]
    assert gradX.count()[1] == [3, 6, 3]
    assert gradX.count()[2] == [2, 6, 3]
    assert gradX.count()[3] == [0, 2, 3]

    for cell in gradX.complex():
        assert gradX.value(cell) == grading(cell)
