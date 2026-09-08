import numpy as np
import pytest

import diode


def values(filtration):
    return {tuple(sorted(int(v) for v in simplex)): alpha
            for simplex, alpha in filtration}


def test_skinny_weighted_tetrahedron_retains_finite_power_radius():
    # Solving normal equations squares the condition number and can lose this
    # tetrahedron entirely. The original three power equations are nonsingular.
    epsilon = 1e-8
    weight = 1e-4
    points = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.],
                       [1., epsilon, 0., weight], [0., 0., 1., 0.]])
    result = values(diode.fill_weighted_alpha_shapes(points))
    center_y = (epsilon * epsilon - weight) / (2 * epsilon)
    assert result[(0, 1, 2, 3)] == pytest.approx(0.5 + center_y**2, rel=1e-12)


@pytest.mark.parametrize("periodic", [False, True])
def test_weight_shift_preserves_complex_and_shifts_every_alpha(periodic):
    rng = np.random.default_rng(20260903)
    points = np.column_stack((rng.random((300, 3)),
                              rng.integers(0, 1024, 300) / 1048576.0))
    fill = (diode.fill_weighted_periodic_alpha_shapes if periodic
            else diode.fill_weighted_alpha_shapes)
    baseline = values(fill(points))
    points[:, 3] -= 2.0
    shifted = values(fill(points))
    assert shifted.keys() == baseline.keys()
    for simplex, alpha in baseline.items():
        assert shifted[simplex] == pytest.approx(alpha + 2.0, rel=1e-7, abs=1e-10)


def test_periodic_large_common_weight_preserves_delaunay():
    coordinates = np.random.default_rng(732).random((300, 3))
    points = np.column_stack((coordinates, np.zeros(len(coordinates))))
    baseline = {tuple(sorted(s)) for s in diode.fill_weighted_periodic_delaunay(points)}
    points[:, 3] = -1e20
    shifted = {tuple(sorted(s)) for s in diode.fill_weighted_periodic_delaunay(points)}
    assert shifted == baseline


@pytest.mark.parametrize("exact", [False, True])
def test_periodic_hidden_sites_do_not_change_visible_complex(exact):
    rng = np.random.default_rng(1709)
    grid = np.indices((6, 6, 6)).reshape(3, -1).T
    coordinates = (grid + 0.25 + 0.5 * rng.random(grid.shape)) / 6
    points = np.column_stack((coordinates, rng.random(len(grid)) * 1e-4))
    baseline = values(diode.fill_weighted_periodic_alpha_shapes(points, exact))
    assert {s[0] for s in baseline if len(s) == 1} == set(range(len(points)))

    # Both a coincident lower-weight site and a distinct negative-weight site
    # have empty power cells. Prepending them also checks surviving input IDs.
    hidden = np.array([np.r_[points[0, :3], -4.0],
                       [0.125, 0.375, 0.625, -4.0]])
    augmented_points = np.vstack((hidden, points))
    augmented = values(diode.fill_weighted_periodic_alpha_shapes(
        augmented_points, exact))
    expected = {tuple(v + len(hidden) for v in simplex): alpha
                for simplex, alpha in baseline.items()}
    assert augmented.keys() == expected.keys()
    for simplex, alpha in expected.items():
        assert augmented[simplex] == pytest.approx(alpha, rel=1e-7, abs=1e-10)
    delaunay = diode.fill_weighted_periodic_delaunay_arrays(augmented_points, exact)
    assert {tuple(sorted(int(v) for v in simplex))
            for rows in delaunay for simplex in rows} == expected.keys()


def test_periodic_rectangular_domain_translation_and_scaling():
    rng = np.random.default_rng(762)
    extent = np.array([1., 1.25, 1.5])
    points = rng.random((800, 3)) * extent
    baseline = values(diode.fill_periodic_alpha_shapes(
        points, False, [0., 0., 0.], extent.tolist()))
    origin = np.array([-3., 5., 2.])
    scale = 2.0
    transformed = values(diode.fill_periodic_alpha_shapes(
        origin + scale * points, False, origin.tolist(),
        (origin + scale * extent).tolist()))
    assert transformed.keys() == baseline.keys()
    for simplex, alpha in baseline.items():
        assert transformed[simplex] == pytest.approx(scale**2 * alpha, rel=1e-7, abs=1e-10)
