"""Cross-check the fast Delaunay-direct alpha shapes against the slow reference.

diode.fill_alpha_shapes uses the fast paths (Delaunay_triangulation_3 + direct
squared-circumradius in 3D; Delaunay_triangulation_2 + face-info circumradius in
2D). diode.fill_alpha_shapes_slow keeps the original reference implementations
(CGAL::Alpha_shape_3 in 3D, the std::set-based path in 2D), unchanged, so we can
compare the two on every run.

We check, over many random clouds (2D/3D, exact True/False):
  * identical (simplex, alpha) sets,
  * for with_attachment, that the reported attacher tau is a valid Gabriel coface:
    squared_circumradius(tau) == alpha(sigma), for BOTH fast and slow.
"""
import numpy as np
import pytest
import diode


# ---------- helpers ----------
def to_value_dict(filtration):
    """list of (verts, alpha[, tau]) -> {sorted-vertex-tuple: alpha}."""
    d = {}
    for entry in filtration:
        verts, alpha = entry[0], entry[1]
        d[tuple(sorted(int(v) for v in verts))] = float(alpha)
    return d


def sq_circumradius_batch(P):
    """P: (m, k+1, D) simplex vertex coords -> (m,) squared circumradius.

    Circumcenter c = p0 + V^T t with V = (p_j - p0); solve (V V^T) t = 0.5 |V_j|^2,
    then R^2 = |c - p0|^2. Works for k in {1,2,3} (edges, triangles, tets).
    """
    m, kp1, D = P.shape
    k = kp1 - 1
    if k == 0:
        return np.zeros(m, dtype=float)
    if k == 1:                                 # edge: smallest sphere is diametric
        d = P[:, 1, :] - P[:, 0, :]
        return 0.25 * np.sum(d * d, axis=1)
    v0 = P[:, 0, :]
    V = P[:, 1:, :] - v0[:, None, :]          # (m, k, D)
    M = V @ np.transpose(V, (0, 2, 1))        # (m, k, k)
    b = 0.5 * np.sum(V * V, axis=2)           # (m, k)
    # numpy 2.x batched solve treats a 2-D b as a stack of matrices; pass an
    # explicit column vector (m, k, 1) and squeeze.
    t = np.linalg.solve(M, b[:, :, None])[:, :, 0]   # (m, k)
    c = np.einsum("mk,mkd->md", t, V)         # (m, D), = c - p0
    return np.sum(c * c, axis=1)


def check_attacher(filtration, points, rtol):
    """Assert squared_circumradius(tau) == alpha for every (sigma, alpha, tau)."""
    by_dim = {}
    for verts, alpha, tau in filtration:
        td = len(tau)
        by_dim.setdefault(td, ([], []))
        by_dim[td][0].append([int(x) for x in tau])
        by_dim[td][1].append(float(alpha))
    for td, (taus, alphas) in by_dim.items():
        taus = np.asarray(taus, dtype=np.int64)
        alphas = np.asarray(alphas, dtype=float)
        P = points[taus]                       # (m, td, D)
        r2 = sq_circumradius_batch(P)
        assert np.allclose(r2, alphas, rtol=rtol, atol=1e-9), (
            f"attacher circumradius mismatch at tau-dim {td-1}: "
            f"max |r2-alpha|={np.max(np.abs(r2-alphas)):.2e}")


SIZES = [10, 50, 200, 800]
DIMS = [2, 3]
EXACTS = [False, True]


@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("exact", EXACTS)
def test_fast_matches_slow_values(dim, n, exact):
    rng = np.random.default_rng(1000 * n + 10 * dim + int(exact))
    pts = rng.random((n, dim))
    fast = to_value_dict(diode.fill_alpha_shapes(pts, exact=exact))
    slow = to_value_dict(diode.fill_alpha_shapes_slow(pts, exact=exact))
    assert set(fast) == set(slow), (
        f"simplex sets differ (fast {len(fast)}, slow {len(slow)})")
    # same exact setting => same kernel => values match very tightly
    rtol = 1e-12 if exact else 1e-7
    for k, fv in fast.items():
        sv = slow[k]
        denom = max(abs(sv), 1e-12)
        assert abs(fv - sv) / denom < rtol or abs(fv - sv) < 1e-12, (
            f"value mismatch at {k}: fast {fv} slow {sv}")


@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("exact", EXACTS)
def test_fast_matches_slow_with_attachment(dim, n, exact):
    rng = np.random.default_rng(2000 * n + 10 * dim + int(exact))
    pts = rng.random((n, dim))
    fast = diode.fill_alpha_shapes(pts, exact=exact, with_attachment=True)
    slow = diode.fill_alpha_shapes_slow(pts, exact=exact, with_attachment=True)

    # same (simplex, alpha) set as each other and as the no-attachment path
    fd, sd = to_value_dict(fast), to_value_dict(slow)
    assert set(fd) == set(sd)
    plain = to_value_dict(diode.fill_alpha_shapes(pts, exact=exact))
    assert set(fd) == set(plain)
    for k in fd:
        assert abs(fd[k] - plain[k]) < (1e-12 if exact else 1e-6 * max(abs(plain[k]), 1e-9) + 1e-12)

    # the attacher contract (squared_circumradius(tau) == alpha) holds for both.
    # Tolerance is limited by our double-precision reference circumradius, not by
    # diode; a structurally wrong tau would be off by orders of magnitude.
    check_attacher(fast, pts, rtol=1e-5)
    check_attacher(slow, pts, rtol=1e-5)


def test_slow_function_exists():
    assert hasattr(diode, "fill_alpha_shapes_slow")
