"""Testing script for CubicalMorseMatching.h"""

from pychomp import *

def test_matching():
    # Initialized with default (trivial) grading
    X = CubicalComplex([5, 6, 2, 5])
    X_match = CubicalMorseMatching(X)

    # Check trichotomy
    for cell in range(X.size()):
        mate = X_match.mate(cell)
        assert mate < X.size()
        assert mate == cell or X_match.mate(mate) == cell

    def grading(cell):
        if cell in {0, 1, 3, 4, 12, 15, 24, 25}:
            return 0
        return 1

    Y = CubicalComplex([3, 4])
    Y_grad = GradedComplex(Y, grading)
    Y_match = CubicalMorseMatching(Y_grad)

    # Check trichotomy
    for cell in range(Y.size()):
        mate = Y_match.mate(cell)
        assert mate < Y.size()
        assert mate == cell or Y_match.mate(mate) == cell

    Z = CubicalComplex([3, 4])
    Z_grad = GradedComplex(Z, grading)
    Z_match = CubicalMorseMatching(Z_grad)

    # Check trichotomy
    for cell in range(Z.size()):
        mate = Z_match.mate(cell)
        assert mate < Z.size()
        assert mate == cell or Z_match.mate(mate) == cell

def test_truncated_matching():
    X = CubicalComplex([3, 4])
    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0
    grad_X = GradedComplex(X, grading)
    X_match = CubicalMorseMatching(grad_X, truncate=True, max_grade=0)

    # Check trichotomy
    for cell in range(X.size()):
        if grading(cell) == 0:
            mate = X_match.mate(cell)
            assert mate < X.size() and grading(mate) == grading(cell)
            assert mate == cell or X_match.mate(mate) == cell
