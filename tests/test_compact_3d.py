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


def test_weight_shift_preserves_complex_and_shifts_every_alpha():
    rng = np.random.default_rng(20260903)
    points = np.column_stack((rng.random((300, 3)),
                              rng.integers(0, 1024, 300) / 1048576.0))
    baseline = values(diode.fill_weighted_alpha_shapes(points))
    points[:, 3] -= 2.0
    shifted = values(diode.fill_weighted_alpha_shapes(points))
    assert shifted.keys() == baseline.keys()
    for simplex, alpha in baseline.items():
        assert shifted[simplex] == pytest.approx(alpha + 2.0, rel=1e-7, abs=1e-10)


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
