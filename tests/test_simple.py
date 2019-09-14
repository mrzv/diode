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
