"""Testing script for CubicalComplex.h"""

from pychomp import *

"""
CubicalComplex([3, 4])


        ^           ^           ^
        |           |           |
       33    45    34    46    35    47  >
        |           |           |       /
        |           |           |      /
      > 9----21----10----22----11----23
     /  |           |           |
    /   |           |           |
       30    42    31    43    32    44  >
        |           |           |       /
        |           |           |      /
      > 6----18---- 7----19---- 8----20
     /  |           |           |
    /   |           |           |
       27    39    28    40    29    41  >
        |           |           |       /
        |           |           |      /
      > 3----15---- 4----16---- 5----17
     /  |           |           |
    /   |           |           |
       24    36    25    37    26    38  >
        |           |           |       /
        |           |           |      /
      > 0----12---- 1----13---- 2----14
     /  ^           ^           ^
    /   |           |           |

"""

"""
CubicalComplex([3, 3, 3])
Numbered 0-cells and 1-cells


    Third floor

        ^           ^           ^
        |   ^       |   ^       |   ^
       78  /       79  /       80  /     >
        | 105       | 106       | 107   /
        |/          |/          |/     /
      >24----51----25----52----26----53
     /  |           |           |
    /   |   ^       |   ^       |   ^
       75  /       76  /       77  /     >
        | 102       | 103       | 104   /
        |/          |/          |/     /
      >21----48----22----49----23----50
     /  |           |           |
    /   |   ^       |   ^       |   ^
       72  /       73  /       74  /     >
        | 99        | 100       | 101   /
        |/          |/          |/     /
      >18----45----19----46----20----47
     /  |           |           |
    /   ^           ^           ^

---------------------------------------------------

    Second floor

        ^           ^           ^
        |   ^       |   ^       |   ^
       69  /       70  /       71  /     >
        | 96        | 97        | 98    /
        |/          |/          |/     /
      >15----42----16----43----17----44
     /  |           |           |
    /   |   ^       |   ^       |   ^
       66  /       67  /       68  /     >
        | 93        | 94        | 95    /
        |/          |/          |/     /
      >12----39----13----40----14----41
     /  |           |           |
    /   |   ^       |   ^       |   ^
       63  /       64  /       65  /     >
        | 90        | 91        | 92    /
        |/          |/          |/     /
      > 9----36----10----37----11----38
     /  |           |           |
    /   ^           ^           ^

---------------------------------------------------

    First floor

        ^           ^           ^
        |   ^       |   ^       |   ^
       60  /       61  /       62  /     >
        | 87        | 88        | 89    /
        |/          |/          |/     /
      > 6----33---- 7----34---- 8----35
     /  |           |           |
    /   |   ^       |   ^       |   ^
       57  /       58  /       59  /     >
        | 84        | 85        | 86    /
        |/          |/          |/     /
      > 3----30---- 4----31---- 5----32
     /  |           |           |
    /   |   ^       |   ^       |   ^
       54  /       55  /       56  /     >
        | 81        | 82        | 83    /
        |/          |/          |/     /
      > 0----27---- 1----28---- 2----29
     /
    /

"""


def test_cells():
    X = CubicalComplex([3, 4])

    assert X.dimension() == 2
    assert len(X) == 48
    assert X.size() == 48
    assert X.size(0) == 12
    assert X.size(1) == 24
    assert X.size(2) == 12
    assert X.count() == [12, 24, 12]
    assert X.boxes() == [3, 4]

    assert X.cell_type(1) == 0
    assert X.cell_shape(1) == 0b00
    assert X.cell_pos(1) == 1
    assert X.cell_dim(1) == 0
    assert X.cell_index([1, 0], 0b0) == 1
    assert X.coordinates(1) == [1, 0]
    assert X.barycenter(1) == [2, 0]

    assert X.cell_type(29) == 2
    assert X.cell_shape(29) == 0b10
    assert X.cell_pos(29) == 5
    assert X.cell_dim(29) == 1
    assert X.cell_index([2, 1], 0b10) == 29
    assert X.coordinates(29) == [2, 1]
    assert X.barycenter(29) == [4, 3]

    assert X.cell_type(42) == 3
    assert X.cell_shape(42) == 0b11
    assert X.cell_pos(42) == 6
    assert X.cell_dim(42) == 2
    assert X.cell_index([0, 2], 0b11) == 42
    assert X.coordinates(42) == [0, 2]
    assert X.barycenter(42) == [1, 5]

    Y = CubicalComplex([3, 3, 3])
    assert Y.dimension() == 3
    assert len(Y) == 216
    assert Y.size() == 216
    assert Y.size(-1) == 0
    assert Y.size(0) == 27
    assert Y.size(1) == 81
    assert Y.size(2) == 81
    assert Y.size(3) == 27
    assert Y.size(4) == 0
    assert Y.count() == [27, 81, 81, 27]
    assert Y.boxes() == [3, 3, 3]

    # Test iterator
    cell_list = []
    for cell in Y:
      cell_list.append(cell)
    assert cell_list == list(range(216))
    for cell in Y(0):
      assert cell == 0
      break
    for cell in Y(1):
      assert cell == 27
      break
    for cell in Y(2):
      assert cell == 108
      break

    assert Y.cell_type(17) == 0
    assert Y.cell_shape(17) == 0b000
    assert Y.cell_pos(17) == 17
    assert Y.cell_dim(17) == 0
    assert Y.cell_index([2, 2, 1], 0b000) == 17
    assert Y.coordinates(17) == [2, 2, 1]
    assert Y.barycenter(17) == [4, 4, 2]

    assert Y.cell_type(102) == 3
    assert Y.cell_shape(102) == 0b100
    assert Y.cell_pos(102) == 21
    assert Y.cell_dim(102) == 1
    assert Y.cell_index([0, 1, 2], 0b100) == 102
    assert Y.coordinates(102) == [0, 1, 2]
    assert Y.barycenter(102) == [0, 2, 5]


def test_boundary():
    X = CubicalComplex([3, 4])

    assert X.boundary({0}) == set()
    assert X.boundary({12}) == {0, 1}
    assert X.boundary({17}) == {5, 6}
    assert X.boundary({23}) == {0, 11}
    assert X.boundary({34}) == {1, 10}
    assert X.boundary({40}) == {16, 19, 28, 29}
    assert X.boundary({41}) == {17, 20, 29, 30}
    assert X.boundary({46}) == {13, 22, 34, 35}
    assert X.boundary({47}) == {14, 23, 24, 35}

    assert X.boundary({1, 4, 7}) == set()
    assert X.boundary({25, 28, 31}) == {1, 10}
    assert X.boundary({25, 28, 31, 34}) == set()
    assert X.boundary({16, 19}) == {4, 5, 7, 8}
    assert X.boundary({36, 37, 39, 40}) == {12, 13, 18, 19, 24, 26, 27, 29}

    Y = CubicalComplex([3, 3, 3])

    assert Y.boundary({9}) == set()
    assert Y.boundary({31}) == {4, 5}
    assert Y.boundary({35}) == {8, 9}
    assert Y.boundary({62}) == {8, 11}
    assert Y.boundary({101}) == {2, 20}
    for cell in Y(0): assert len(Y.boundary({cell})) == 0
    for cell in Y(1): assert len(Y.boundary({cell})) == 2
    for cell in Y(2): assert len(Y.boundary({cell})) == 4
    for cell in Y(3): assert len(Y.boundary({cell})) == 6

    assert Y.boundary({36, 39, 64}) == {9, 12}
    assert Y.boundary({36, 39, 63, 64}) == set()


def test_coboundary():
    X = CubicalComplex([3, 4])

    assert X.coboundary({0}) == {12, 23, 24, 33}
    assert X.coboundary({10}) == {21, 22, 31, 34}
    assert X.coboundary({23}) == {44, 47}
    assert X.coboundary({13}) == {37, 46}
    assert X.coboundary({27}) == {38, 39}
    assert X.coboundary({36}) == set()

    assert X.coboundary({4, 5, 7, 8}) == {15, 17, 18, 20, 25, 26, 31, 32}
    assert X.coboundary({26, 29, 32, 35}) == {37, 38, 40, 41, 43, 44, 46, 47}
    assert X.coboundary({12, 15, 18, 21}) == set()
    assert X.coboundary({36, 37, 38}) == set()

    Y = CubicalComplex([3, 3, 3])
    assert Y.coboundary({0}) == {27, 53, 54, 78, 81, 99}
    assert Y.coboundary({13}) == {39, 40, 64, 67, 85, 94}
    for cell in Y(0): assert len(Y.coboundary({cell})) == 6
    for cell in Y(1): assert len(Y.coboundary({cell})) == 4
    for cell in Y(2): assert len(Y.coboundary({cell})) == 2
    for cell in Y(3): assert len(Y.coboundary({cell})) == 0


def test_closure():
    X = CubicalComplex([3, 4])

    assert X.closure({0}) == {0}
    assert X.closure({16}) == {4, 5, 16}
    assert X.closure({38}) == {2, 3, 5, 6, 14, 17, 26, 27, 38}

    assert X.closure({0, 4, 5, 16, 40}) == {0, 4, 5, 7, 8, 16, 19, 28, 29, 40}

    Y = CubicalComplex([3, 3, 3])

    assert Y.closure({25}) == {25}
    assert Y.closure({51}) == {24, 25, 51}
    assert Y.closure({43, 97, 98}) == {16, 17, 25, 26, 43, 97, 98}


def test_star():
    X = CubicalComplex([3, 4])

    assert X.star({0}) == {0, 12, 23, 24, 33, 36, 44, 45, 47}
    assert X.star({19}) == {19, 40, 43}
    assert X.star({40}) == {40}

    assert X.star({0, 12, 13}) == X.star({0}) | X.star({12}) | X.star({13})

    # `topstar` takes in a single cell (integer) and outputs a list
    # very different behavior than `star`
    assert set(X.topstar(0)) == {36, 44, 45, 47}
    assert set(X.topstar(19)) == {40, 43}
    assert set(X.topstar(40)) == {40}

    Y = CubicalComplex([3, 3, 3])

    assert {0, 27, 53, 54, 78, 81, 99}.issubset(Y.star({0}))
    assert Y.star({215}) == {215}

    assert len(Y.topstar(0)) == 8
    assert set(Y.topstar(215)) == {215}


def test_leftright():
    X = CubicalComplex([3, 4])

    assert X.left(7, 0) == 18
    assert X.left(7, 1) == 28
    assert X.left(9, 0) == 20
    assert X.left(9, 1) == 30
    assert X.left(12, 0) == 0
    assert X.left(12, 1) == 45
    assert X.left(29, 0) == 40
    assert X.left(29, 1) == 5
    assert X.left(33, 0) == 44
    assert X.left(33, 1) == 9
    assert X.left(36, 0) == 24
    assert X.left(36, 1) == 12
    assert X.left(45, 0) == 33
    assert X.left(45, 1) == 21

    assert X.right(7, 0) == 19
    assert X.right(7, 1) == 31
    assert X.right(9, 0) == 21
    assert X.right(9, 1) == 33
    assert X.right(12, 0) == 1
    assert X.right(12, 1) == 36
    assert X.right(29, 0) == 41
    assert X.right(29, 1) == 8
    assert X.right(33, 0) == 45
    assert X.right(33, 1) == 0
    assert X.right(36, 0) == 25
    assert X.right(36, 1) == 15
    assert X.right(45, 0) == 34
    assert X.right(45, 1) == 12

    assert not X.leftfringe(7)
    assert X.leftfringe(9)
    assert X.leftfringe(12)
    assert not X.leftfringe(29)
    assert X.leftfringe(33)
    assert not X.leftfringe(36)
    assert not X.leftfringe(45)

    assert not X.rightfringe(7)
    assert not X.rightfringe(9)
    assert not X.rightfringe(12)
    assert not X.rightfringe(29)
    assert X.rightfringe(33)
    assert not X.rightfringe(36)
    assert X.rightfringe(45)

    Y = CubicalComplex([3, 3, 3])

    assert not Y.leftfringe(13)
    assert Y.leftfringe(6)
    assert Y.leftfringe(20)
    assert Y.leftfringe(4)
    assert Y.leftfringe(28)
    assert not Y.leftfringe(41)
    assert Y.leftfringe(60)
    assert Y.leftfringe(66)

    assert not Y.rightfringe(13)
    assert not Y.rightfringe(6)
    assert Y.rightfringe(105)
    assert not Y.rightfringe(77)
    assert Y.rightfringe(50)


def test_coords():
    X = CubicalComplex([3, 4])

    assert X.mincoords(0) == 0b11
    assert X.maxcoords(0) == 0b00

    assert X.mincoords(17) == 0b00
    assert X.maxcoords(17) == 0b01

    assert X.mincoords(33) == 0b01
    assert X.maxcoords(33) == 0b10

    assert X.mincoords(43) == 0b00
    assert X.maxcoords(43) == 0b00

    Y = CubicalComplex([3, 3, 3])

    assert Y.mincoords(21) == 0b001
    assert Y.maxcoords(21) == 0b100

    assert Y.mincoords(14) == 0b000
    assert Y.maxcoords(14) == 0b001

    assert Y.mincoords(91) == 0b010
    assert Y.maxcoords(91) == 0b000

    assert Y.mincoords(35) == 0b100
    assert Y.maxcoords(35) == 0b011


def test_parallelneighbors():
    X = CubicalComplex([3, 4])

    assert set(X.parallelneighbors(0)) == {0}
    assert set(X.parallelneighbors(12)) == {12, 13, 23}
    assert set(X.parallelneighbors(31)) == {28, 31, 34}
    assert set(X.parallelneighbors(40)) == {36, 37, 38, 39, 40, 41, 42, 43, 44}

    Y = CubicalComplex([3, 3, 3])

    assert set(Y.parallelneighbors(7)) == {7}
    assert set(Y.parallelneighbors(44)) == {43, 44, 45}
    assert set(Y.parallelneighbors(84)) == {84, 93, 102}
    assert set(Y.parallelneighbors(65)) == {62, 65, 68}
