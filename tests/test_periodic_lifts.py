import itertools

import diode
import numpy as np
import pytest


def periodic_cloud(dim):
    # Seeded cloud with a pair whose shortest edge crosses the periodic boundary.
    n = 40 if dim == 2 else 250
    points = np.random.default_rng(700 + dim).random((n, dim))
    points[0] = [0.01, 0.5] if dim == 2 else [0.01, 0.5, 0.5]
    points[1] = [0.99, 0.5] if dim == 2 else [0.99, 0.5, 0.5]
    return points


def simplex_set(arrays):
    # Compare Delaunay outputs independently of their per-dimension ordering.
    return {tuple(int(v) for v in row) for array in arrays for row in np.sort(array, axis=1)}


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("exact", [False, True])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_periodic_delaunay_lifts_contract(dim, exact, dtype):
    points = periodic_cloud(dim).astype(dtype)
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
    points = periodic_cloud(2)
    vertices, offsets = diode.fill_periodic_delaunay_lifts_arrays(
        points, exact=exact, bbox_min=[0, 0], bbox_max=[1, 1]
    )
    row_index = np.flatnonzero(np.all(vertices[1] == [0, 1], axis=1))
    assert row_index.shape == (1,)
    lifted = points[vertices[1][row_index[0]]] + offsets[1][row_index[0]]
    np.testing.assert_allclose(lifted[1] - lifted[0], [-0.02, 0.0], atol=1e-12)


@pytest.mark.parametrize(
    "points,bbox_min,bbox_max,match",
    [
        (np.array([[0.1, 0.2], [np.nan, 0.3]]), [0, 0], [1, 1], "finite"),
        (np.array([[0.1, 0.2], [1.0, 0.3]]), [0, 0], [1, 1], "half-open"),
        (np.array([[0.1, 0.2], [0.3, 0.4]]), [1, 0], [0, 1], "empty or inverted"),
        (np.array([[0.1, 0.2], [0.3, 0.4]]), [-np.inf, 0], [1, 1], "finite"),
        (np.array([[0.1, 0.2], [0.3, 0.4]]), [-1e308, 0], [1e308, 1], "finite"),
    ],
)
def test_periodic_delaunay_lifts_validate_input(points, bbox_min, bbox_max, match):
    with pytest.raises(RuntimeError, match=match):
        diode.fill_periodic_delaunay_lifts_arrays(
            points, bbox_min=bbox_min, bbox_max=bbox_max
        )


def test_periodic_delaunay_lifts_function_exists():
    assert hasattr(diode, "fill_periodic_delaunay_lifts_arrays")
