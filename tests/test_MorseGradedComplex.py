"""Testing script for MorseGradedComplex.h"""

from pychomp import *

def test_cubical_morse():
    def S1_grading(cell):
        if cell in {0, 1, 3, 4, 12, 15, 24, 25}:
            return 0
        return 1

    X = CubicalComplex([3, 4])
    X_grad = GradedComplex(X, S1_grading)
    X_morse = MorseGradedComplex(X_grad)
    X_morse_trunc = MorseGradedComplex(X_grad, truncate = True, max_grade = 0)

    assert X_morse.count()[0] == X_morse_trunc.count()[0]
    assert X_morse.count()[0] == [1, 1, 0]
    assert X_morse.count()[1] == [1, 1, 1]
    assert list(X_morse_trunc.count().keys()) == [0]

    # 0 grade complex is fully reduced; boundary is the zero operator
    for cell in X_morse.complex():
        assert X_morse.value(cell) == 1 or X_morse.complex().boundary({cell}) == set()

    Y_morse = MorseGradedComplex(X_morse)
    Y_morse_trunc = MorseGradedComplex(Y_morse, truncate = True)
    assert Y_morse.count()[0] == X_morse.count()[0]
    assert Y_morse.count()[1] == [0, 0, 1]
    assert Y_morse_trunc.count() == X_morse_trunc.count()
