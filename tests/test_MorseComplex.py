"""Testing script for MorseComplex.h"""

from pychomp import *

"""
CubicalComplex([3, 4])
`grading` defined to keep 0-cells and 1-cells only.


        ^           ^           ^
        |           |           |
       33          34          35        >
        |           |           |       /
        |           |           |      /
      > 9----21----10----22----11----23
     /  |           |           |
    /   |           |           |
       30          31          32        >
        |           |           |       /
        |           |           |      /
      > 6----18---- 7----19---- 8----20
     /  |           |           |
    /   |           |           |
       27          28          29        >
        |           |           |       /
        |           |           |      /
      > 3----15---- 4----16---- 5----17
     /  |           |           |
    /   |           |           |
       24          25          26        >
        |           |           |       /
        |           |           |      /
      > 0----12---- 1----13---- 2----14
     /  ^           ^           ^
    /   |           |           |


After cubical morse matching, critical cells are:

0D: 11
1D: 24, 25, 27, 28, 30, 31

(note fringe cells are discarded)

"""

def test_instantiation():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    assert morse_X.dimension() == 2
    assert len(morse_X) == 7
    assert morse_X.size() == 7
    assert morse_X.size(0) == 1
    assert morse_X.size(1) == 6
    assert morse_X.size(2) == 0
    assert morse_X.count() == [1, 6, 0]

    Y = CubicalComplex([3, 4])
    morse_Y = MorseComplex(Y)

    assert morse_Y.dimension() == 2
    assert len(morse_Y) == 1
    assert morse_Y.size() == 1
    assert morse_Y.size(0) == 1
    assert morse_Y.size(1) == 0
    assert morse_Y.size(2) == 0
    assert morse_Y.count() == [1, 0, 0]

def test_include_project():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    included = []
    for cell in morse_X:
        included.append(morse_X.include({cell}))
    assert included == [{11}, {24}, {25}, {27}, {28}, {30}, {31}]

    projected = []
    for chain in included:
        projected.append(morse_X.project(chain))
    assert projected == [{0}, {1}, {2}, {3}, {4}, {5}, {6}]
    assert morse_X.project({13}) == set()

def test_boundary():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    # One iteration reduces down to minimal morse complex
    for cell in morse_X:
        assert morse_X.boundary({cell}) == set()

def test_coboundary():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    # One iteration reduces down to minimal morse complex
    for cell in morse_X:
        assert morse_X.coboundary({cell}) == set()

def test_flow():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    # Queens
    assert morse_X.flow({0}) == (
        {11}, {12, 13, 26, 29, 32}
    )
    assert morse_X.flow({4}) == (
        {11}, {16, 29, 32}
    )
    assert morse_X.flow({10}) == (
        {11}, {22}
    )

    # Ace
    assert morse_X.flow({11}) == (
        {11}, set()
    )

    # Top-dimensional cells are either kings or aces -> canonical
    for cell in X(1):
        assert morse_X.flow({cell}) == (
            {cell}, set()
        )

def test_lift_lower():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    assert morse_X.lift(morse_X.project({11})) == {
        11
    }
    assert morse_X.lift(morse_X.project({25})) == {
        13, 16, 25, 26
    }
    assert morse_X.lift(morse_X.project({30})) == {
        18, 19, 21, 22, 30, 32
    }

    # As this is the minimal morse complex, each cell in the morse
    # complex should lift to a homology generator, which is a cycle
    for cell in morse_X:
        assert X.boundary(morse_X.lift({cell})) == set()

    # The canonical "part" of all 0-cells is just {11},
    # which projects to morse_X.project({11}) == {0}
    for cell in X(0):
        assert morse_X.lower({cell}) == morse_X.project({11})

    # Kings are canonical but project to the empty set.
    # Aces are canonical and project to themselves.
    for cell in X(1):
        assert morse_X.lower({cell}) == morse_X.project({cell})

def test_coflow():
    # The coflow is the dual of the flow.
    # The kings of the coflow are the cocells corresponding to
    # queens of the flow, and vice versa.
    # Critical cells correspond to critical cocells
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    # 0-cocells are either queens or aces -> cocanonical
    for cocell in X(0):
        assert morse_X.coflow({cocell}) == (
            {cocell}, set()
        )

    # Aces
    for cocell in morse_X:
        assert morse_X.coflow(morse_X.include({cocell})) == (
            morse_X.include({cocell}), set()
        )

    # Kings
    assert morse_X.coflow({12}) == (
        {24, 33, 23}, {0}
    )
    assert morse_X.coflow({29}) == (
        {17, 23, 27, 28, 33, 34, 35}, {0, 1, 2, 3, 4, 5}
    )

def test_colift_colower():
    X = CubicalComplex([3, 4])

    def grading(cell):
        if X.cell_dim(cell) == 2:
            return 1
        return 0

    grad_X = GradedComplex(X, grading)
    morse_X = MorseGradedComplex(grad_X, truncate=True, max_grade=0).complex()

    assert morse_X.colift(morse_X.project({11})) == {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
    }
    for cocell in morse_X(1):
        assert morse_X.colift({cocell}) == morse_X.include({cocell})

    # As this is the minimal morse complex, each cocell in the morse
    # complex should colift to a cohomology generator, which is a cocycle
    # Hence has no coboundary in the 0-grade original complex.
    for cocell in morse_X:
        lift_cbd = X.coboundary(morse_X.colift({cocell}))
        for cell in lift_cbd:
            assert grading(cell) != 0

    # Queens are cocanonical but project to the empty set.
    # Aces are cocanonical and project to themselves.
    for cocell in X(0):
        assert morse_X.colower({cocell}) == morse_X.project({cocell})

    assert morse_X.colower({12}) == {1}
    assert morse_X.colower({29}) == {3, 4}
