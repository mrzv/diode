"""Additional fast-vs-slow coverage over non-uniform sample families.

The existing branch tests already compare the fast paths to the `_slow`
reference on uniform random clouds.  This file keeps the same oracle but covers
normal clouds, mixtures, and samples from simple shapes, with and without small
coordinate noise.  Noiseless lower-dimensional shape samples are intentionally
only checked with exact=True: exact=False may expose ordinary inexact-kernel
degeneracy rather than a diode wrapper regression.
"""
from itertools import combinations

import numpy as np
import pytest

import diode


def _value_dict(filtration):
    return {
        tuple(sorted(int(v) for v in verts)): float(alpha)
        for verts, alpha in filtration
    }


def _arrays_value_dict(verts_by_dim, vals_by_dim):
    out = {}
    for verts, vals in zip(verts_by_dim, vals_by_dim):
        verts = np.asarray(verts)
        vals = np.asarray(vals)
        assert verts.ndim == 2
        assert vals.ndim == 1
        assert verts.shape[0] == vals.shape[0]
        for row, value in zip(verts, vals):
            out[tuple(sorted(int(v) for v in row))] = float(value)
    return out


def _assert_close_dicts(actual, expected, *, exact):
    assert set(actual) == set(expected)
    rtol = 1e-12 if exact else 1e-6
    atol = 1e-10 if exact else 1e-8
    for simplex, actual_value in actual.items():
        expected_value = expected[simplex]
        assert np.isclose(actual_value, expected_value, rtol=rtol, atol=atol), (
            f"value mismatch at {simplex}: {actual_value} != {expected_value}"
        )


def _assert_filtered_complex(filtration):
    values = _value_dict(filtration)
    assert len(values) == len(filtration), "duplicate simplices emitted"
    for simplex, alpha in values.items():
        for size in range(1, len(simplex)):
            for face in combinations(simplex, size):
                assert face in values, f"missing face {face} of {simplex}"
                assert values[face] <= alpha + 1e-9, (
                    f"face {face} value {values[face]} > coface {simplex} value {alpha}"
                )


def _scale_to_unit_box(points):
    points = np.asarray(points, dtype=np.float64)
    lo = points.min(axis=0)
    span = points.max(axis=0) - lo
    span[span == 0.0] = 1.0
    return (points - lo) / span


def _noisy(points, noise, seed):
    points = np.asarray(points, dtype=np.float64)
    if noise == 0.0:
        return points
    rng = np.random.default_rng(seed)
    return points + rng.normal(scale=noise, size=points.shape)


def _normal_cloud(dim, seed):
    rng = np.random.default_rng(seed)
    return _scale_to_unit_box(rng.normal(loc=0.25, scale=1.7, size=(42, dim)))


def _uniform_cloud(dim, seed):
    rng = np.random.default_rng(seed)
    return rng.uniform(low=-2.0, high=3.0, size=(44, dim))


def _gaussian_uniform_mixture(dim, seed):
    rng = np.random.default_rng(seed)
    gaussian = rng.normal(loc=-0.4, scale=0.35, size=(28, dim))
    uniform = rng.uniform(low=0.15, high=1.15, size=(30, dim))
    return _scale_to_unit_box(np.vstack([gaussian, uniform]))


def _two_normal_mixture(dim, seed):
    rng = np.random.default_rng(seed)
    first = rng.normal(loc=-0.7, scale=0.20, size=(24, dim))
    second = rng.normal(loc=0.8, scale=0.45, size=(26, dim))
    return _scale_to_unit_box(np.vstack([first, second]))


def _weighted_data(kind, seed):
    rng = np.random.default_rng(seed)
    if kind == "normal_small_weights":
        points = _normal_cloud(3, seed)
        weights = rng.uniform(0.0, 0.01, size=(points.shape[0], 1))
    elif kind == "normal_mixed_weights":
        points = _normal_cloud(3, seed)
        weights = np.vstack([
            rng.uniform(0.0, 0.003, size=(points.shape[0] // 2, 1)),
            rng.uniform(0.02, 0.08, size=(points.shape[0] - points.shape[0] // 2, 1)),
        ])
    elif kind == "mixture_small_weights":
        points = _gaussian_uniform_mixture(3, seed)
        weights = rng.uniform(0.0, 0.015, size=(points.shape[0], 1))
    elif kind == "mixture_mixed_weights":
        points = _two_normal_mixture(3, seed)
        weights = rng.beta(2.0, 8.0, size=(points.shape[0], 1)) * 0.08
    else:
        raise AssertionError(kind)
    return np.hstack([points, weights])


def _rectangle(seed):
    rng = np.random.default_rng(seed)
    sides = []
    for t in np.linspace(0.0, 1.0, 12, endpoint=False):
        sides.extend([(t, 0.0), (1.0, t), (1.0 - t, 1.0), (0.0, 1.0 - t)])
    interior = rng.uniform(0.15, 0.85, size=(10, 2))
    return np.vstack([np.asarray(sides), interior])


def _triangle(seed):
    rng = np.random.default_rng(seed)
    a = np.array([0.0, 0.0])
    b = np.array([1.0, 0.0])
    c = np.array([0.25, 0.9])
    edges = []
    for t in np.linspace(0.0, 1.0, 16, endpoint=False):
        edges.extend([(1 - t) * a + t * b, (1 - t) * b + t * c, (1 - t) * c + t * a])
    interior = rng.dirichlet([2.0, 2.0, 2.0], size=10) @ np.vstack([a, b, c])
    return np.vstack([np.asarray(edges), interior])


def _circle(seed):
    rng = np.random.default_rng(seed)
    theta = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
    boundary = np.column_stack([np.cos(theta), np.sin(theta)])
    interior_theta = rng.uniform(0.0, 2.0 * np.pi, size=12)
    interior_r = np.sqrt(rng.uniform(0.0, 0.7, size=12))
    interior = np.column_stack([interior_r * np.cos(interior_theta),
                                interior_r * np.sin(interior_theta)])
    return _scale_to_unit_box(np.vstack([boundary, interior]))


def _box_surface(seed):
    rng = np.random.default_rng(seed)
    points = []
    grid = np.linspace(0.0, 1.0, 5)
    for x in grid:
        for y in grid:
            points.extend([(x, y, 0.0), (x, y, 1.0)])
    points.extend(rng.uniform(0.1, 0.9, size=(16, 3)))
    return np.asarray(points, dtype=np.float64)


def _sphere(seed):
    rng = np.random.default_rng(seed)
    u = np.linspace(-0.85, 0.85, 8)
    theta = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    points = []
    for z in u:
        r = np.sqrt(max(0.0, 1.0 - z * z))
        for t in theta:
            points.append((r * np.cos(t), r * np.sin(t), z))
    interior = rng.normal(size=(10, 3))
    interior /= np.linalg.norm(interior, axis=1)[:, None]
    interior *= rng.uniform(0.0, 0.65, size=(10, 1))
    return _scale_to_unit_box(np.vstack([np.asarray(points), interior]))


def _torus(seed):
    major = 1.0
    minor = 0.32
    points = []
    for u in np.linspace(0.0, 2.0 * np.pi, 9, endpoint=False):
        for v in np.linspace(0.0, 2.0 * np.pi, 7, endpoint=False):
            points.append(((major + minor * np.cos(v)) * np.cos(u),
                           (major + minor * np.cos(v)) * np.sin(u),
                           minor * np.sin(v)))
    points = np.asarray(points, dtype=np.float64)
    return _scale_to_unit_box(points)


DISTRIBUTIONS = {
    "normal": _normal_cloud,
    "uniform": _uniform_cloud,
    "gaussian_uniform_mixture": _gaussian_uniform_mixture,
    "two_normal_mixture": _two_normal_mixture,
}


WEIGHTED_DISTRIBUTIONS = [
    "normal_small_weights",
    "normal_mixed_weights",
    "mixture_small_weights",
    "mixture_mixed_weights",
]


SHAPES = {
    "rectangle": _rectangle,
    "triangle": _triangle,
    "circle": _circle,
    "box_surface": _box_surface,
    "sphere": _sphere,
    "torus": _torus,
}


@pytest.mark.parametrize("name", sorted(DISTRIBUTIONS))
@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("exact", [False, True])
def test_distribution_samples_match_slow(name, dim, exact):
    points = DISTRIBUTIONS[name](dim, seed=17 * len(name) + dim)
    fast = _value_dict(diode.fill_alpha_shapes(points, exact=exact))
    slow = _value_dict(diode.fill_alpha_shapes_slow(points, exact=exact))
    _assert_close_dicts(fast, slow, exact=exact)
    _assert_filtered_complex(diode.fill_alpha_shapes(points, exact=exact))


@pytest.mark.parametrize("name", sorted(DISTRIBUTIONS))
@pytest.mark.parametrize("dim", [2, 3])
def test_alpha_arrays_match_list_exporter(name, dim):
    points = DISTRIBUTIONS[name](dim, seed=100 + len(name) + dim)
    listed = _value_dict(diode.fill_alpha_shapes(points, exact=True))
    verts_by_dim, vals_by_dim = diode.fill_alpha_shapes_arrays(points, exact=True)
    arrays = _arrays_value_dict(verts_by_dim, vals_by_dim)
    _assert_close_dicts(arrays, listed, exact=True)


@pytest.mark.parametrize("name", WEIGHTED_DISTRIBUTIONS)
@pytest.mark.parametrize("exact", [False, True])
def test_weighted_distribution_samples_match_slow(name, exact):
    data = _weighted_data(name, seed=500 + len(name) + int(exact))
    fast = _value_dict(diode.fill_weighted_alpha_shapes(data, exact=exact))
    slow = _value_dict(diode.fill_weighted_alpha_shapes_slow(data, exact=exact))
    _assert_close_dicts(fast, slow, exact=exact)
    _assert_filtered_complex(diode.fill_weighted_alpha_shapes(data, exact=exact))


@pytest.mark.parametrize("name", sorted(SHAPES))
@pytest.mark.parametrize("noise", [0.0, 1e-5])
def test_shape_samples_match_slow_exact(name, noise):
    points = _noisy(SHAPES[name](seed=300 + len(name)), noise, seed=900 + len(name))
    fast = _value_dict(diode.fill_alpha_shapes(points, exact=True))
    slow = _value_dict(diode.fill_alpha_shapes_slow(points, exact=True))
    _assert_close_dicts(fast, slow, exact=True)
    _assert_filtered_complex(diode.fill_alpha_shapes(points, exact=True))


@pytest.mark.parametrize("name", sorted(SHAPES))
def test_shape_samples_with_noise_match_slow_inexact(name):
    points = _noisy(SHAPES[name](seed=700 + len(name)), 1e-5, seed=1100 + len(name))
    fast = _value_dict(diode.fill_alpha_shapes(points, exact=False))
    slow = _value_dict(diode.fill_alpha_shapes_slow(points, exact=False))
    _assert_close_dicts(fast, slow, exact=False)
    _assert_filtered_complex(diode.fill_alpha_shapes(points, exact=False))


NOISELESS_INEXACT_DEGENERATE = {"sphere", "torus", "triangle"}


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(
            name,
            marks=pytest.mark.xfail(
                strict=True,
                reason="noiseless shape sample is degenerate under exact=False",
            )
            if name in NOISELESS_INEXACT_DEGENERATE
            else (),
        )
        for name in sorted(SHAPES)
    ],
)
def test_shape_samples_without_noise_inexact_match_slow_when_non_degenerate(name):
    points = SHAPES[name](seed=1200 + len(name))
    fast = _value_dict(diode.fill_alpha_shapes(points, exact=False))
    slow = _value_dict(diode.fill_alpha_shapes_slow(points, exact=False))
    _assert_close_dicts(fast, slow, exact=False)
    _assert_filtered_complex(diode.fill_alpha_shapes(points, exact=False))
