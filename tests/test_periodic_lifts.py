import itertools

import diode
import numpy as np
import pytest


def periodic_cloud():
    # Seeded 3D cloud with a pair whose shortest edge crosses the boundary.
    points = np.random.default_rng(703).random((250, 3))
    points[0] = [0.01, 0.5, 0.5]
    points[1] = [0.99, 0.5, 0.5]
    return points


def simplex_set(arrays):
    # Compare Delaunay outputs independently of their per-dimension ordering.
    return {tuple(int(v) for v in row) for array in arrays for row in np.sort(array, axis=1)}


@pytest.mark.parametrize("exact", [False, True])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_periodic_delaunay_lifts_contract(exact, dtype):
    dim = 3
    points = periodic_cloud().astype(dtype)
    bbox_min = np.zeros(dim)
    bbox_max = np.ones(dim)
    vertices, offsets = diode.fill_periodic_delaunay_lifts_arrays(
        points, exact=exact, bbox_min=bbox_min, bbox_max=bbox_max
    )

    assert len(vertices) == len(offsets) == dim + 1
    for simplex_dim, (vertex_rows, offset_rows) in enumerate(zip(vertices, offsets)):
        assert vertex_rows.dtype == np.int64
        assert offset_rows.dtype == np.int64
        assert vertex_rows.flags.c_contiguous
        assert offset_rows.flags.c_contiguous
        assert vertex_rows.shape[1:] == (simplex_dim + 1,)
        assert offset_rows.shape == (*vertex_rows.shape, dim)
        assert np.all(vertex_rows[:, 1:] > vertex_rows[:, :-1])
        assert np.all(offset_rows[:, 0] == 0)
        assert len({tuple(row) for row in vertex_rows}) == len(vertex_rows)

    old_vertices = diode.fill_periodic_delaunay_arrays(
        points, exact, bbox_min.tolist(), bbox_max.tolist()
    )
    assert simplex_set(vertices) == simplex_set(old_vertices)
    counts = [len(rows) for rows in vertices]
    assert sum((-1) ** d * count for d, count in enumerate(counts)) == 0

    width = bbox_max - bbox_min
    edge_displacements = {}
    for vertex_row, offset_row in zip(vertices[1], offsets[1]):
        lifted = points[vertex_row] + offset_row * width
        edge_displacements[tuple(vertex_row)] = lifted[1] - lifted[0]

    for vertex_rows, offset_rows in zip(vertices[2:], offsets[2:]):
        for vertex_row, offset_row in zip(vertex_rows, offset_rows):
            lifted = points[vertex_row] + offset_row * width
            for i, j in itertools.combinations(range(len(vertex_row)), 2):
                key = (int(vertex_row[i]), int(vertex_row[j]))
                np.testing.assert_allclose(
                    lifted[j] - lifted[i], edge_displacements[key], rtol=0, atol=1e-6
                )


@pytest.mark.parametrize("exact", [False, True])
def test_periodic_boundary_edge_uses_short_lift(exact):
    points = periodic_cloud()
    vertices, offsets = diode.fill_periodic_delaunay_lifts_arrays(
        points, exact=exact, bbox_min=[0, 0, 0], bbox_max=[1, 1, 1]
    )
    row_index = np.flatnonzero(np.all(vertices[1] == [0, 1], axis=1))
    assert row_index.shape == (1,)
    lifted = points[vertices[1][row_index[0]]] + offsets[1][row_index[0]]
    np.testing.assert_allclose(lifted[1] - lifted[0], [-0.02, 0.0, 0.0], atol=1e-12)

@pytest.mark.parametrize("exact", [False, True])
def test_periodic_boundary_edge_has_short_lift_alpha(exact):
    dim = 3
    points = periodic_cloud()
    values = {
        tuple(sorted(int(vertex) for vertex in simplex)): alpha
        for simplex, alpha in diode.fill_periodic_alpha_shapes(
            points, exact, [0] * dim, [1] * dim
        )
    }
    assert values[(0, 1)] == pytest.approx(0.0001, rel=1e-10, abs=1e-15)


@pytest.mark.parametrize("exact", [False, True])
def test_rectangular_translated_lifts_match_empty_spheres_and_alpha(exact):
    width = np.array([1.0, 1.25, 1.5])
    origin = np.array([-3.0, 5.0, 2.0])
    points = origin + periodic_cloud() * width
    vertices, offsets = diode.fill_periodic_delaunay_lifts_arrays(
        points, exact=exact, bbox_min=origin, bbox_max=origin + width
    )
    alpha = {
        tuple(sorted(int(v) for v in simplex)): value
        for simplex, value in diode.fill_periodic_alpha_shapes(
            points, exact, origin.tolist(), (origin + width).tolist()
        )
    }
    assert simplex_set(vertices) == alpha.keys()
    assert set(vertices[0][:, 0]) == set(range(len(points)))

    # Recover each sphere from the exported physical geometry, independently
    # of the alpha computation. Tetrahedron alpha is its squared circumradius.
    lifted = points[vertices[3]] + offsets[3] * width
    edges = lifted[:, 1:] - lifted[:, :1]
    squared_lengths = np.sum(edges * edges, axis=2)
    relative_centers = np.linalg.solve(2 * edges, squared_lengths[..., None])[..., 0]
    radii_squared = np.sum(relative_centers * relative_centers, axis=1)
    np.testing.assert_allclose(
        [alpha[tuple(row)] for row in vertices[3]], radii_squared,
        rtol=1e-7, atol=1e-10,
    )

    # Wrong lattice signs, axis extents, or origin handling can give internally
    # consistent alpha/lifts but fail the physical periodic empty-ball property.
    centers = lifted[:, 0] + relative_centers
    displacement = (centers[:, None] - points + width / 2) % width - width / 2
    distances_squared = np.sum(displacement * displacement, axis=2)
    assert np.all(distances_squared >= radii_squared[:, None] - 1e-9)


@pytest.mark.parametrize(
    "points,bbox_min,bbox_max",
    [
        (np.array([[0.1, 0.2, 0.3], [np.nan, 0.3, 0.4]]), [0, 0, 0], [1, 1, 1]),
        (np.array([[0.1, 0.2, 0.3], [1.0, 0.3, 0.4]]), [0, 0, 0], [1, 1, 1]),
        (np.array([[0.1, 0.2, 0.3], [0.3, 0.4, 0.5]]), [1, 0, 0], [0, 1, 1]),
        (np.array([[0.1, 0.2, 0.3], [0.3, 0.4, 0.5]]), [-np.inf, 0, 0], [1, 1, 1]),
        (np.array([[0.1, 0.2, 0.3], [0.3, 0.4, 0.5]]), [-1e308, 0, 0], [1e308, 1, 1]),
    ],
)
def test_periodic_delaunay_lifts_validate_input(points, bbox_min, bbox_max):
    with pytest.raises(RuntimeError):
        diode.fill_periodic_delaunay_lifts_arrays(
            points, bbox_min=bbox_min, bbox_max=bbox_max
        )


@pytest.mark.parametrize(
    "fill",
    [
        diode.fill_periodic_alpha_shapes,
        diode.fill_periodic_alpha_shapes_slow,
        diode.fill_periodic_alpha_shapes_arrays,
        diode.fill_periodic_delaunay,
        diode.fill_periodic_delaunay_arrays,
        diode.fill_periodic_delaunay_lifts_arrays,
    ],
    ids=lambda fill: fill.__name__,
)
@pytest.mark.parametrize("exact", [False, True])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("n", [0, 3], ids=["empty", "populated"])
def test_periodic_2d_is_unsupported(fill, exact, dtype, n):
    points = np.array([[0.1, 0.2], [0.3, 0.4], [0.7, 0.2]], dtype=dtype)[:n]
    # Neither the default 3D box nor an explicit 2D box enables a tiled fallback.
    with pytest.raises(NotImplementedError):
        fill(points, exact=exact)
    with pytest.raises(NotImplementedError):
        fill(points, exact, [0., 0.], [1., 1.])
