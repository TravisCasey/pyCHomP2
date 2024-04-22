"""Testing script for Grading.h"""

from pychomp import *

def test_top_grading():
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

    X = CubicalComplex([3, 4])
    grading = construct_grading(X, top_grading)

    for cell in range(36, 48):
        assert grading(cell) == top_grading(cell)

    assert grading(0) == 0
    assert grading(25) == 0
    assert grading(15) == 0
    assert grading(27) == 0
    assert grading(28) == 1
    assert grading(7) == 1
    assert grading(21) == 2
    assert grading(20) == 1
    assert grading(35) == 3

def test_inclusion_grading():
    include = {1, 6, 12, 18, 24, 25, 39}
    full = {0, 1, 3, 4, 6, 7, 12, 15, 18, 24, 25, 27, 28, 39}

    X = CubicalComplex([3, 4])
    grading = inclusion_grading(X, include)

    for cell in X:
        if cell in full:
            assert grading(cell) == 0
        else:
            assert grading(cell) == 1

def test_cubical_nerve():
    positions = {0, 1, 2, 3, 4, 9, 10, 11, 12, 13, 14, 18}
    edges = {27, 28, 29, 30, 36, 37, 38, 39, 40, 54, 55, 63, 64, 65,
             81, 82, 83, 84, 85, 90, 99}
    square_count = 9
    X = CubicalComplex([3, 3, 3])
    grading = cubical_nerve(X, positions, max_dim=2)

    squares = 0
    for cell in X:
        if cell in positions or cell in edges:
            assert grading(cell) == 0
        elif X.cell_dim(cell) == 2 and grading(cell) == 0:
            squares += 1
        else:
            assert grading(cell) == 1
    assert squares == square_count
