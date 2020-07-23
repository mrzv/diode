import numpy as np
import diode

def test_cube():
    points = np.array([[0.,0.,0.],
                       [0.,0.,1.],
                       [0.,1.,0.],
                       [1.,0.,0.]])
    f = diode.fill_alpha_shapes(points)

    expected_f = [([1], 0.0),
                  ([2], 0.0),
                  ([0], 0.0),
                  ([3], 0.0),
                  ([1, 0], 0.25),
                  ([0, 3], 0.25),
                  ([2, 0], 0.25),
                  ([1, 3], 0.5),
                  ([2, 1], 0.5),
                  ([2, 3], 0.5),
                  ([1, 0, 3], 0.5),
                  ([2, 0, 3], 0.5),
                  ([2, 1, 0], 0.5),
                  ([2, 1, 3], 0.75),
                  ([2, 1, 0, 3], 0.75)]

    for s1,s2 in zip(sorted(f), sorted(expected_f)):
        assert(s1 == s2)

def test_square():
    points = np.array([[0.,0.],[1.,0.],[0.,1.]])
    f = diode.fill_alpha_shapes(points)

    expected_f = [([0], 0.0),
                  ([1], 0.0),
                  ([2], 0.0),
                  ([0, 1], 0.25),
                  ([0, 2], 0.25),
                  ([1, 2], 0.5),
                  ([0, 1, 2], 0.5)]

    for s1,s2 in zip(sorted(f), sorted(expected_f)):
        assert(s1 == s2)

def is_sorted(lst, key = lambda x: x):
    return all(key(lst[i]) <= key(lst[i+1]) for i in range(len(lst) - 1))

def test_sorted():
    np.random.seed(42)
    for dim in [2,3]:
        points = np.random.random((1000,dim))
        f = diode.fill_alpha_shapes(points)

        assert(is_sorted(f, key = lambda x: (x[1], len(x[0]))))


def test_periodic():
    np.random.seed(42)
    for dim in [2,3]:
        points = np.random.random((1000,dim))
        f = diode.fill_periodic_alpha_shapes(points)

        assert(is_sorted(f, key = lambda x: (x[1], len(x[0]))))
    
