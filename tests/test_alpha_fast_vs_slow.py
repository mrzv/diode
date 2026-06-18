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
    assert hasattr(diode, "fill_periodic_alpha_shapes_slow")


# ---- periodic 2D: fast (Delaunay-direct) vs slow (std::set) -----------------
# Periodic alpha picks an arbitrary periodic *offset representative* for a few
# near-degenerate edges, so a handful of alpha values can differ slightly --
# the old path is itself non-deterministic there. We therefore require the
# combinatorics (simplex set) to match exactly and the values to agree up to a
# small, bounded number of such edges. (We deliberately do NOT "fix" diode's
# existing periodic Gabriel test.)
@pytest.mark.parametrize("n", [20, 100, 500, 1500])
@pytest.mark.parametrize("exact", EXACTS)
def test_periodic_2d_fast_vs_slow_values(n, exact):
    rng = np.random.default_rng(5000 * n + int(exact))
    pts = rng.random((n, 2))
    fast = to_value_dict(diode.fill_periodic_alpha_shapes(pts, exact, [0., 0.], [1., 1.]))
    slow = to_value_dict(diode.fill_periodic_alpha_shapes_slow(pts, exact, [0., 0.], [1., 1.]))
    assert set(fast) == set(slow), "periodic simplex sets differ"
    diffs = np.array([abs(fast[k] - slow[k]) for k in fast])
    n_differ = int((diffs > 1e-7).sum())
    assert n_differ <= max(3, len(diffs) // 100), \
        f"{n_differ}/{len(diffs)} periodic values differ -- more than offset ambiguity explains"
    assert diffs.max(initial=0.0) < 1e-2, "periodic value difference too large"


def _gudhi_diagram(filtration, dim, gudhi):
    """Persistence diagram (per dimension) of an arbitrary (verts, value) filtration."""
    st = gudhi.SimplexTree()
    for verts, val in filtration:
        st.insert([int(v) for v in verts], float(val))
    st.make_filtration_non_decreasing()
    st.compute_persistence(persistence_dim_max=True)
    d = np.asarray(st.persistence_intervals_in_dimension(dim), dtype=float)
    return d[np.isfinite(d).all(axis=1)] if d.size else d   # finite bars only


# ---- periodic 3D: fast (Delaunay-direct + CGAL periodic is_Gabriel) vs slow
# (Alpha_shape_3) ------------------------------------------------------------
# Unlike the 2D direct path (which does a frame-mixed manual Gabriel test), the
# 3D direct path uses CGAL's Periodic_3_Delaunay_triangulation_3::is_Gabriel,
# which handles offsets internally -- the same predicate Alpha_shape_3 uses. So
# the values agree to ~1e-16, not just up to offset ambiguity. The periodic
# triangulation must be representable in 1 sheet (needs enough points in 3D);
# both paths raise the same error otherwise.
@pytest.mark.parametrize("n", [200, 800, 2000])
@pytest.mark.parametrize("exact", EXACTS)
def test_periodic_3d_fast_vs_slow_values(n, exact):
    rng = np.random.default_rng(7000 * n + int(exact))
    pts = rng.random((n, 3))
    try:
        slow = to_value_dict(diode.fill_periodic_alpha_shapes_slow(pts, exact, [0.]*3, [1.]*3))
    except RuntimeError:
        with pytest.raises(RuntimeError):
            diode.fill_periodic_alpha_shapes(pts, exact, [0.]*3, [1.]*3)
        pytest.skip("point cloud not representable in 1 sheet")
    fast = to_value_dict(diode.fill_periodic_alpha_shapes(pts, exact, [0.]*3, [1.]*3))
    assert set(fast) == set(slow), "periodic 3D simplex sets differ"
    diffs = np.array([abs(fast[k] - slow[k]) for k in fast])
    # CGAL is_Gabriel makes the direct values match Alpha_shape_3 to round-off.
    assert diffs.max(initial=0.0) < 1e-7, \
        f"periodic 3D value difference too large: {diffs.max(initial=0.0):.2e}"


@pytest.mark.parametrize("n", [100, 500, 1500])
def test_periodic_2d_fast_vs_slow_diagram(n):
    # Optional, stronger check: persistence diagrams agree under bottleneck
    # distance. Skipped if no PD library is installed (diode itself does not
    # depend on one).
    gudhi = pytest.importorskip("gudhi")
    rng = np.random.default_rng(9000 + n)
    pts = rng.random((n, 2))
    fast = diode.fill_periodic_alpha_shapes(pts, False, [0., 0.], [1., 1.])
    slow = diode.fill_periodic_alpha_shapes_slow(pts, False, [0., 0.], [1., 1.])
    for d in (0, 1):
        bd = gudhi.bottleneck_distance(_gudhi_diagram(fast, d, gudhi),
                                       _gudhi_diagram(slow, d, gudhi))
        assert bd < 1e-2, f"periodic dim-{d} bottleneck distance {bd} too large"
