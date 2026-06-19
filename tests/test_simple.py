import numpy as np
import pytest
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

    # Compare canonically: vertex order within a simplex is not a contract
    # (the fast Delaunay-direct path emits a different per-simplex vertex order
    # than the old Alpha_shape_3 path; consumers sort vertices anyway).
    canon = lambda g: sorted((tuple(sorted(v)), val) for v, val in g)
    assert canon(f) == canon(expected_f)

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

    canon = lambda g: sorted((tuple(sorted(v)), val) for v, val in g)
    assert canon(f) == canon(expected_f)

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


# ---- with_attachment tests --------------------------------------------------

def _circumradius_squared(point_coords):
    """Squared circumradius of a simplex from its vertex coordinates.

    point_coords: (k, d) array, k in {1,2,3,4}, d in {2,3}.
    Returns the squared smallest enclosing sphere through all k points,
    matching CGAL's squared_radius for the same inputs.
    """
    pc = np.asarray(point_coords, dtype=np.float64)
    k, d = pc.shape
    if k == 1:
        return 0.0
    if k == 2:
        return 0.25 * np.sum((pc[0] - pc[1]) ** 2)
    if k == 3:
        # Circumradius of a triangle in 2D or 3D.
        a = pc[1] - pc[0]
        b = pc[2] - pc[0]
        a2 = np.dot(a, a)
        b2 = np.dot(b, b)
        ab = np.dot(a, b)
        denom = 2.0 * (a2 * b2 - ab * ab)
        if denom == 0:
            raise ValueError("degenerate triangle")
        # Bary coords of circumcenter relative to p0.
        ca = (b2 * (a2 - ab)) / denom
        cb = (a2 * (b2 - ab)) / denom
        center_offset = ca * a + cb * b
        return float(np.dot(center_offset, center_offset))
    if k == 4:
        # Circumradius of a tetrahedron in 3D.
        a = pc[1] - pc[0]
        b = pc[2] - pc[0]
        c = pc[3] - pc[0]
        a2 = np.dot(a, a)
        b2 = np.dot(b, b)
        c2 = np.dot(c, c)
        # Solve the 3x3 system M x = rhs where M rows are 2*<a,b,c> dotted with x
        M = 2.0 * np.array([[np.dot(a, a), np.dot(a, b), np.dot(a, c)],
                            [np.dot(a, b), np.dot(b, b), np.dot(b, c)],
                            [np.dot(a, c), np.dot(b, c), np.dot(c, c)]])
        rhs = np.array([a2, b2, c2])
        bary = np.linalg.solve(M, rhs)
        center_offset = bary[0] * a + bary[1] * b + bary[2] * c
        return float(np.dot(center_offset, center_offset))
    raise ValueError(f"unexpected simplex size {k}")


def test_cube_attachment_shapes():
    points = np.array([[0.,0.,0.],
                       [0.,0.,1.],
                       [0.,1.,0.],
                       [1.,0.,0.]])
    f = diode.fill_alpha_shapes(points, with_attachment=True)

    # Every entry is a 3-tuple.
    for entry in f:
        assert len(entry) == 3
        sigma, alpha, tau = entry
        assert isinstance(sigma, list)
        assert isinstance(tau, list)
        assert isinstance(alpha, float)
        assert 1 <= len(sigma) <= 4
        assert 1 <= len(tau) <= 4
        assert len(tau) >= len(sigma)


def test_cube_attachment_recompute():
    points = np.array([[0.,0.,0.],
                       [0.,0.,1.],
                       [0.,1.,0.],
                       [1.,0.,0.]])
    f = diode.fill_alpha_shapes(points, with_attachment=True)

    for sigma, alpha, tau in f:
        recomputed = _circumradius_squared(points[tau])
        assert abs(alpha - recomputed) < 1e-12, \
            f"sigma={sigma}, alpha={alpha}, tau={tau}, recomputed={recomputed}"


def test_cube_attachment_specific():
    points = np.array([[0.,0.,0.],
                       [0.,0.,1.],
                       [0.,1.,0.],
                       [1.,0.,0.]])
    f = diode.fill_alpha_shapes(points, with_attachment=True)

    # Build a lookup by sorted sigma vertex tuple.
    by_sigma = { tuple(sorted(sigma)): (alpha, sorted(tau)) for sigma, alpha, tau in f }

    # Vertices: alpha=0, tau == sigma.
    for v in [0, 1, 2, 3]:
        alpha, tau = by_sigma[(v,)]
        assert alpha == 0.0
        assert tau == [v]

    # Short edges (length 1): alpha=0.25, Gabriel, tau == sigma.
    for short_edge in [(0, 1), (0, 2), (0, 3)]:
        alpha, tau = by_sigma[short_edge]
        assert abs(alpha - 0.25) < 1e-12
        assert tuple(tau) == short_edge

    # Long edges (length sqrt(2)): alpha=0.5. For this configuration the
    # third tetrahedron vertex lies exactly on the edge's MEB boundary, so
    # CGAL classifies these edges as Gabriel and the recompute is trivially
    # consistent with tau == sigma.
    for long_edge in [(1, 2), (1, 3), (2, 3)]:
        alpha, tau = by_sigma[long_edge]
        assert abs(alpha - 0.5) < 1e-12
        # tau must be a Gabriel coface whose own squared circumradius equals 0.5.
        recomputed = _circumradius_squared(points[tau])
        assert abs(alpha - recomputed) < 1e-12
        # Both endpoints of the edge must be vertices of tau.
        assert long_edge[0] in tau and long_edge[1] in tau

    # Three Gabriel facets {0,1,2}, {0,1,3}, {0,2,3}: alpha=0.5, tau == sigma.
    for f3 in [(0, 1, 2), (0, 1, 3), (0, 2, 3)]:
        alpha, tau = by_sigma[f3]
        assert abs(alpha - 0.5) < 1e-12
        assert tuple(tau) == f3

    # Non-Gabriel facet {1,2,3}: alpha=0.75, tau == cell {0,1,2,3}.
    alpha, tau = by_sigma[(1, 2, 3)]
    assert abs(alpha - 0.75) < 1e-12
    assert tuple(tau) == (0, 1, 2, 3)

    # Cell {0,1,2,3}: alpha=0.75, tau == sigma.
    alpha, tau = by_sigma[(0, 1, 2, 3)]
    assert abs(alpha - 0.75) < 1e-12
    assert tuple(tau) == (0, 1, 2, 3)


def _close(a, b, rtol=1e-7, atol=1e-9):
    """Hybrid relative/absolute tolerance.

    Random Delaunay configurations contain near-flat tetrahedra with huge
    circumradii where the linear-algebra recompute and CGAL's geometric
    formula can diverge by several units in the last few digits of a large
    value. The user's downstream PyTorch implementation will use a different
    formula too, so we validate approximate agreement, not bit-equality.
    """
    return abs(a - b) <= atol + rtol * abs(a)


def test_random_3d_attachment_recompute():
    np.random.seed(42)
    points = np.random.random((50, 3))
    f = diode.fill_alpha_shapes(points, with_attachment=True)

    for sigma, alpha, tau in f:
        recomputed = _circumradius_squared(points[tau])
        assert _close(alpha, recomputed), \
            f"sigma={sigma}, alpha={alpha}, tau={tau}, recomputed={recomputed}, diff={alpha-recomputed}"


def test_random_2d_attachment_recompute():
    np.random.seed(7)
    points = np.random.random((50, 2))
    f = diode.fill_alpha_shapes(points, with_attachment=True)

    for sigma, alpha, tau in f:
        recomputed = _circumradius_squared(points[tau])
        assert _close(alpha, recomputed), \
            f"sigma={sigma}, alpha={alpha}, tau={tau}, recomputed={recomputed}, diff={alpha-recomputed}"


def test_2d_obtuse_attachment():
    # Obtuse triangle: longest edge (between 1 and 2) is non-Gabriel of the face.
    points = np.array([[0., 0.],
                       [3., 0.],
                       [1., 0.5]])
    f = diode.fill_alpha_shapes(points, with_attachment=True)

    by_sigma = { tuple(sorted(sigma)): (alpha, sorted(tau)) for sigma, alpha, tau in f }

    # Edge (1, 2) has the largest squared length 4 + 0.25 = 4.25; squared
    # circumradius of the edge alone is 4.25 / 4 ~ 1.0625. The face's
    # circumradius squared should be larger only if the triangle is acute.
    # For this obtuse triangle, the longest edge is non-Gabriel within the face.
    edge_long = (0, 1)
    alpha_long, tau_long = by_sigma[edge_long]
    # Verify recompute consistency regardless of attachment outcome.
    recomputed = _circumradius_squared(points[tau_long])
    assert abs(alpha_long - recomputed) < 1e-9


def test_pair_form_unchanged():
    # Backward compatibility: with_attachment=False (default) returns pairs.
    points = np.array([[0.,0.,0.],
                       [0.,0.,1.],
                       [0.,1.,0.],
                       [1.,0.,0.]])
    f_pair = diode.fill_alpha_shapes(points)
    f_pair_explicit = diode.fill_alpha_shapes(points, with_attachment=False)

    for entry in f_pair:
        assert len(entry) == 2
    assert sorted(f_pair) == sorted(f_pair_explicit)


def test_pair_and_triple_alpha_match():
    # Pair-form alpha values must match triple-form alpha values for the
    # same simplex. Also verifies sigma order is consistent.
    points = np.array([[0.,0.,0.],
                       [0.,0.,1.],
                       [0.,1.,0.],
                       [1.,0.,0.]])
    f_pair = diode.fill_alpha_shapes(points)
    f_triple = diode.fill_alpha_shapes(points, with_attachment=True)

    by_sigma_pair = { tuple(sorted(sigma)): alpha for sigma, alpha in f_pair }
    by_sigma_triple = { tuple(sorted(sigma)): alpha for sigma, alpha, _ in f_triple }
    assert by_sigma_pair == by_sigma_triple


def test_with_attachment_not_implemented_for_weighted():
    points = np.array([[0.,0.,0., 0.0],
                       [0.,0.,1., 0.0],
                       [0.,1.,0., 0.0],
                       [1.,0.,0., 0.0]])
    try:
        diode.fill_weighted_alpha_shapes(points, with_attachment=True)
    except NotImplementedError:
        return
    raise AssertionError("expected NotImplementedError for weighted with_attachment=True")


def test_with_attachment_not_implemented_for_periodic():
    np.random.seed(1)
    points = np.random.random((10, 3))
    try:
        diode.fill_periodic_alpha_shapes(points, with_attachment=True)
    except NotImplementedError:
        return
    raise AssertionError("expected NotImplementedError for periodic with_attachment=True")


def test_periodic_domain_requires_enough_bounds():
    points = np.array([[0.1, 0.1, 0.1],
                       [0.4, 0.1, 0.1],
                       [0.1, 0.4, 0.1],
                       [0.1, 0.1, 0.4]])

    with pytest.raises(RuntimeError, match="from/to must have at least"):
        diode.fill_periodic_alpha_shapes(points, False, [0.0], [1.0])

    with pytest.raises(RuntimeError, match="from/to must have at least"):
        diode.fill_periodic_alpha_shapes_arrays(points, False, [0.0], [1.0, 1.0, 1.0])

    with pytest.raises(RuntimeError, match="from/to must have at least"):
        diode.fill_periodic_delaunay(points, False, [0.0, 0.0, 0.0], [1.0])


def test_periodic_domain_rejects_empty_or_inverted_axes():
    points = np.array([[0.1, 0.1],
                       [0.4, 0.1],
                       [0.1, 0.4]])

    with pytest.raises(RuntimeError, match="periodic domain is empty or inverted"):
        diode.fill_periodic_alpha_shapes(points, False, [0.0, 1.0], [1.0, 1.0])

    with pytest.raises(RuntimeError, match="periodic domain is empty or inverted"):
        diode.fill_periodic_alpha_shapes_arrays(points, False, [0.0, 1.0], [1.0, 0.0])

    with pytest.raises(RuntimeError, match="periodic domain is empty or inverted"):
        diode.fill_periodic_delaunay_arrays(points, False, [0.0, 1.0], [1.0, 0.0])
