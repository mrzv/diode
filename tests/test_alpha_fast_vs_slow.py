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
from collections import Counter
from itertools import combinations

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
    """Assert tau is a coface of sigma and squared_circumradius(tau) == alpha
    for every (sigma, alpha, tau)."""
    by_dim = {}
    for verts, alpha, tau in filtration:
        # tau must be a coface of sigma (sigma's vertices are a subset of tau's);
        # a structurally wrong tau with a coincidentally-right circumradius would
        # otherwise slip past the radius check below.
        assert set(int(v) for v in verts) <= set(int(t) for t in tau), (
            f"attacher tau {list(tau)} is not a coface of sigma {list(verts)}")
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
    assert hasattr(diode, "fill_weighted_alpha_shapes_slow")


# ---- regression: degenerate (lower-dimensional) inputs must not crash and must
# match the _slow reference (previously the fast paths segfaulted) --------------
DEGENERATE_2D = [
    np.array([[0., 0.]]),                                  # single point
    np.array([[0., 0.], [1., 0.]]),                        # 2 points (1-D)
    np.array([[0., 0.], [1., 0.], [2., 0.]]),              # collinear
]
DEGENERATE_3D = [
    np.array([[0., 0., 0.]]),                              # single point
    np.array([[0., 0., 0.], [1., 0., 0.]]),               # 2 points
    np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]),  # 3 points (2-D)
    np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.]]),  # coplanar
]


@pytest.mark.parametrize("pts", DEGENERATE_2D + DEGENERATE_3D)
@pytest.mark.parametrize("exact", EXACTS)
def test_degenerate_input_matches_slow(pts, exact):
    assert (to_value_dict(diode.fill_alpha_shapes(pts, exact=exact))
            == to_value_dict(diode.fill_alpha_shapes_slow(pts, exact=exact)))


@pytest.mark.parametrize("pts", DEGENERATE_2D + DEGENERATE_3D)
def test_degenerate_input_with_attachment_matches_slow(pts):
    assert (to_value_dict(diode.fill_alpha_shapes(pts, with_attachment=True))
            == to_value_dict(diode.fill_alpha_shapes_slow(pts, with_attachment=True)))


def test_degenerate_weighted_matches_slow():
    w = np.array([[0., 0., 0., .01], [1., 0., 0., .01]])  # 2 weighted points (1-D)
    assert (to_value_dict(diode.fill_weighted_alpha_shapes(w))
            == to_value_dict(diode.fill_weighted_alpha_shapes_slow(w)))


# ---- regression: duplicate coordinates must pick the same vertex representative
# as the _slow reference (last input index wins) --------------------------------
def _vertex_id_set(filtration):
    return {tuple(sorted(int(v) for v in entry[0])) for entry in filtration}


def _delaunay_vertex_id_set(simplices):
    return {tuple(sorted(int(v) for v in verts)) for verts in simplices}


def test_duplicate_coords_match_slow_vertex_ids():
    # index 0 and index 4 are the same point; _slow keeps index 4.
    p3 = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [0., 0., 0.]])
    assert _vertex_id_set(diode.fill_alpha_shapes(p3)) == _vertex_id_set(diode.fill_alpha_shapes_slow(p3))
    assert _delaunay_vertex_id_set(diode.fill_delaunay(p3)) == _vertex_id_set(diode.fill_alpha_shapes_slow(p3))

    p2 = np.array([[0., 0.], [1., 0.], [0., 1.], [0., 0.]])
    assert _vertex_id_set(diode.fill_alpha_shapes(p2)) == _vertex_id_set(diode.fill_alpha_shapes_slow(p2))

    w = np.array([[0., 0., 0., .01], [1., 0., 0., .01], [0., 1., 0., .01],
                  [0., 0., 1., .01], [0., 0., 0., .01]])
    assert _vertex_id_set(diode.fill_weighted_alpha_shapes(w)) == _vertex_id_set(diode.fill_weighted_alpha_shapes_slow(w))


# ---- regression: malformed from/to length must raise, not crash ---------------
def test_periodic_from_to_length_validated():
    rng = np.random.default_rng(0)
    with pytest.raises(RuntimeError):
        diode.fill_periodic_delaunay(rng.random((20, 2)), False, [0.0], [1.0])
    with pytest.raises(RuntimeError):
        diode.fill_periodic_alpha_shapes(rng.random((20, 3)), False, [0., 0.], [1., 1.])
    with pytest.raises(RuntimeError):
        diode.fill_weighted_periodic_delaunay(rng.random((20, 4)), False, [0., 0.], [1., 1.])
    # the length-3 default is accepted for 2D points (only the first 2 are used)
    diode.fill_periodic_alpha_shapes(rng.random((50, 2)))


def test_periodic_inverted_domain_raises():
    # An empty/inverted periodic box (from >= to on some axis) must raise, not
    # crash the interpreter deep inside CGAL.
    rng = np.random.default_rng(0)
    p3 = rng.random((50, 3))
    p2 = rng.random((50, 2))
    w = np.hstack([rng.random((50, 3)), rng.random((50, 1)) * 0.01])
    with pytest.raises(RuntimeError):
        diode.fill_periodic_alpha_shapes(p3, False, [1., 1., 1.], [0., 0., 0.])   # inverted 3D
    with pytest.raises(RuntimeError):
        diode.fill_periodic_alpha_shapes(p2, False, [1., 1.], [0., 0.])           # inverted 2D
    with pytest.raises(RuntimeError):
        diode.fill_periodic_delaunay(p3, False, [0., 0., 0.], [0., 1., 1.])       # zero extent on x
    with pytest.raises(RuntimeError):
        diode.fill_weighted_periodic_alpha_shapes(w, False, [1., 1., 1.], [0., 0., 0.])
    with pytest.raises(RuntimeError):
        diode.fill_weighted_periodic_delaunay(w, False, [1., 1., 1.], [0., 0., 0.])


# ---- third-party (gudhi) oracle for cases where the _slow reference is itself
# wrong: the fast paths are correct, so we pin them against gudhi, not _slow.
# See SLOW_REFERENCE_BUGS.md.
def test_2d_hull_vertices_match_gudhi():
    gudhi = pytest.importorskip("gudhi")
    # 7 distinct points on a hull config where the _slow 2D path drops a vertex
    # (uninitialized s[0]); the fast path keeps all 7 and matches gudhi.
    pts = np.array([[2., 1.], [0., 0.], [1., 0.], [2., 2.], [0., 2.], [1., 2.], [2., 0.]])
    fast = {frozenset(int(v) for v in s) for s, _ in diode.fill_alpha_shapes(pts, exact=True)}
    st = gudhi.AlphaComplex(points=pts.tolist()).create_simplex_tree()
    g = {frozenset(s) for s, _ in st.get_simplices()}
    assert fast == g, f"fast vs gudhi differ: {fast ^ g}"
    assert sum(1 for s in fast if len(s) == 1) == 7   # _slow drops one vertex here


def test_weighted_shared_coords_match_gudhi():
    gudhi = pytest.importorskip("gudhi")
    # two weighted sites at identical coordinates with different weights: the
    # regular triangulation keeps the larger-weight one, so the surviving vertex's
    # alpha is -max(weight) = -0.5. The fast path gets this right (matches gudhi);
    # the _slow reference mislabels it as -0.01.
    data = np.array([[0., 0., 0., 0.01],
                     [1., 0., 0., 0.01],
                     [0., 1., 0., 0.01],
                     [0., 0., 1., 0.01],
                     [0., 0., 0., 0.5]])
    fast = diode.fill_weighted_alpha_shapes(data)
    fv = sorted(float(a) for s, a in fast if len(s) == 1)
    st = gudhi.AlphaComplex(points=data[:, :3].tolist(),
                            weights=data[:, 3].tolist()).create_simplex_tree()
    gv = sorted(float(v) for s, v in st.get_simplices() if len(s) == 1)
    assert len(fv) == len(gv) and np.allclose(fv, gv, atol=1e-9), (fv, gv)
    assert min(fv) == pytest.approx(-0.5, abs=1e-9)   # the kept site has weight 0.5


# ---- weighted 3D: fast (Regular_triangulation_3 + Edelsbrunner) vs slow
# (Alpha_shape_3 on the regular triangulation) ------------------------------
# Input is a 4-column array (x, y, z, weight). The direct path uses CGAL's
# Regular_triangulation_3::is_Gabriel -- including is_Gabriel(vertex), since a
# weighted vertex (unlike an unweighted one) need not be Gabriel -- so the values
# match Alpha_shape_3 to round-off, not just up to an approximation.
@pytest.mark.parametrize("n", [50, 200, 800])
@pytest.mark.parametrize("exact", EXACTS)
@pytest.mark.parametrize("wscale", [0.01, 0.1])
def test_weighted_3d_fast_vs_slow_values(n, exact, wscale):
    rng = np.random.default_rng(8000 * n + 17 * int(exact) + int(1000 * wscale))
    data = np.hstack([rng.random((n, 3)), rng.random((n, 1)) * wscale])
    fast = to_value_dict(diode.fill_weighted_alpha_shapes(data, exact=exact))
    slow = to_value_dict(diode.fill_weighted_alpha_shapes_slow(data, exact=exact))
    assert set(fast) == set(slow), (
        f"weighted simplex sets differ (fast {len(fast)}, slow {len(slow)})")
    rtol = 1e-12 if exact else 1e-6
    for k, fv in fast.items():
        sv = slow[k]
        assert abs(fv - sv) <= rtol * max(abs(sv), 1.0) + 1e-9, (
            f"weighted value mismatch at {k}: fast {fv} slow {sv}")


# ---- weighted 3D PERIODIC: fast (Periodic_3_regular_triangulation_3 +
# Edelsbrunner) vs slow (Alpha_shape_3) -------------------------------------
# Weights must satisfy 0 <= w < 1/64 * domain^2 (CGAL's periodic-regular
# requirement). The triangulation must be 1-sheet representable, which needs
# enough points; the slow path raises otherwise. Like the unweighted periodic
# case, is_Gabriel handles offsets, so for non-degenerate clouds the values match
# Alpha_shape_3 to round-off (a few near-degenerate simplices may differ via the
# periodic Gabriel offset ambiguity).
@pytest.mark.parametrize("n", [1000, 2500])
@pytest.mark.parametrize("exact", EXACTS)
def test_weighted_periodic_3d_fast_vs_slow_values(n, exact):
    rng = np.random.default_rng(31 * n + int(exact))
    data = np.hstack([rng.random((n, 3)), rng.random((n, 1)) * 0.01])
    frm, to = [0.] * 3, [1.] * 3
    try:
        slow = to_value_dict(diode.fill_weighted_periodic_alpha_shapes_slow(data, exact, frm, to))
    except RuntimeError:
        with pytest.raises(RuntimeError):
            diode.fill_weighted_periodic_alpha_shapes(data, exact, frm, to)
        pytest.skip("point cloud not representable in 1 sheet")
    fast = to_value_dict(diode.fill_weighted_periodic_alpha_shapes(data, exact, frm, to))
    assert set(fast) == set(slow), "weighted periodic simplex sets differ"
    diffs = np.array([abs(fast[k] - slow[k]) for k in fast])
    n_big = int((diffs > 1e-7).sum())
    assert n_big <= max(3, len(diffs) // 100), \
        f"{n_big}/{len(diffs)} weighted-periodic values differ -- more than offset ambiguity explains"
    assert diffs.max(initial=0.0) < 1e-2


# The direct path always emits a valid filtered simplicial complex -- face-closed,
# values non-decreasing onto cofaces, no duplicate/garbage indices -- even for
# sparse clouds near the 1-sheet boundary where the slow Alpha_shape_3 path itself
# degenerates (emits out-of-order / garbage-indexed simplices).
@pytest.mark.parametrize("n", [200, 800, 1500])
@pytest.mark.parametrize("exact", EXACTS)
def test_weighted_periodic_3d_is_valid_complex(n, exact):
    rng = np.random.default_rng(99 * n + int(exact))
    data = np.hstack([rng.random((n, 3)), rng.random((n, 1)) * 0.01])
    try:
        f = diode.fill_weighted_periodic_alpha_shapes(data, exact, [0.] * 3, [1.] * 3)
    except RuntimeError:
        pytest.skip("point cloud not representable in 1 sheet")
    val = {tuple(sorted(int(x) for x in v)): a for v, a in f}
    assert len(val) == len(f), "duplicate index-simplices emitted"
    assert max(k[-1] for k in val) < n, "garbage (out-of-range) vertex index emitted"
    for verts, a in f:
        verts = sorted(int(x) for x in verts)
        for k in range(1, len(verts)):
            for face in combinations(verts, k):
                fa = val.get(tuple(face))
                assert fa is not None, f"missing face {face} of {verts}"
                assert fa <= a + 1e-9, f"face {face} value {fa} > coface {verts} value {a}"
    # NB the 3-torus Euler characteristic is 0 for non-degenerate clouds, but can
    # be off at the sparse 1-sheet boundary where CGAL drops degenerate vertices;
    # the set-match-vs-slow test above pins the topology at non-degenerate n.


# ---- periodic 2D: fast (Delaunay-direct) vs slow (tiling) -------------------
# 2D periodic alpha is resolved differently by the two paths (periodic offsets vs
# a finite tiling), and CGAL picks an arbitrary periodic *offset representative*,
# so for a few near-degenerate edges the alpha value differs by an unbounded amount
# AND varies run-to-run -- this is the documented, deliberately-unfixed 2D periodic
# ambiguity (CGAL's Periodic_2_Delaunay has no is_Gabriel to resolve it cleanly,
# unlike the 3D periodic case). The combinatorics are deterministic, so we pin the
# simplex set exactly and require all but a tiny fraction of values to agree; the
# persistence diagram (what actually matters) is checked under bottleneck distance
# in test_periodic_2d_fast_vs_slow_diagram.
@pytest.mark.parametrize("n", [20, 100, 500, 1500])
@pytest.mark.parametrize("exact", EXACTS)
def test_periodic_2d_fast_vs_slow_values(n, exact):
    rng = np.random.default_rng(5000 * n + int(exact))
    pts = rng.random((n, 2))
    fast = to_value_dict(diode.fill_periodic_alpha_shapes(pts, exact, [0., 0.], [1., 1.]))
    slow = to_value_dict(diode.fill_periodic_alpha_shapes_slow(pts, exact, [0., 0.], [1., 1.]))
    assert set(fast) == set(slow), "periodic simplex sets differ"
    diffs = np.array([abs(fast[k] - slow[k]) for k in fast])
    # The vast majority of edges agree exactly; only a few near-degenerate ones
    # differ (by a possibly large, run-dependent amount -- so we bound the COUNT,
    # not the magnitude, and rely on the diagram test for value correctness).
    assert np.median(diffs) < 1e-9, "most periodic values should match exactly"
    n_differ = int((diffs > 1e-7).sum())
    assert n_differ <= max(5, len(diffs) // 50), \
        f"{n_differ}/{len(diffs)} periodic values differ -- more than offset ambiguity explains"


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
