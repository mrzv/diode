"""Combinatorics-only Delaunay exporters: fill_delaunay / fill_delaunay_arrays.

These return the Delaunay simplices (the alpha-complex simplex set) WITHOUT any
alpha value, for consumers that recompute filtration values themselves. We check,
over random clouds (2D/3D, exact True/False):
  * the simplex SET equals that of fill_alpha_shapes (alpha complex == full
    Delaunay triangulation),
  * the list form and the arrays form agree and contain no duplicates,
  * the arrays are well-shaped (dimension d has d+1 vertex columns),
  * (optionally) the set matches GUDHI's alpha complex.
"""
import numpy as np
import pytest
import diode


def alpha_simplex_set(pts, exact):
    return {frozenset(int(v) for v in verts)
            for verts, _ in diode.fill_alpha_shapes(pts, exact=exact)}


def list_simplex_set(simplices):
    return {frozenset(int(v) for v in verts) for verts in simplices}


def arrays_simplex_set(verts_by_dim):
    s = set()
    for d, arr in enumerate(verts_by_dim):
        arr = np.asarray(arr)
        assert arr.dtype == np.int64
        assert arr.ndim == 2 and arr.shape[1] == d + 1, (d, arr.shape)
        for row in arr:
            s.add(frozenset(int(v) for v in row))
    return s


SIZES = [10, 50, 200, 800]
DIMS = [2, 3]
EXACTS = [False, True]


@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("n", SIZES)
@pytest.mark.parametrize("exact", EXACTS)
def test_delaunay_matches_alpha_simplex_set(dim, n, exact):
    rng = np.random.default_rng(3000 * n + 10 * dim + int(exact))
    pts = rng.random((n, dim))

    s_alpha = alpha_simplex_set(pts, exact)

    dl_list = diode.fill_delaunay(pts, exact=exact)
    s_list = list_simplex_set(dl_list)
    assert len(dl_list) == len(s_list), "fill_delaunay emitted duplicate simplices"
    assert s_list == s_alpha, (
        f"delaunay/alpha simplex sets differ (delaunay {len(s_list)}, alpha {len(s_alpha)})")

    dl_arr = diode.fill_delaunay_arrays(pts, exact=exact)
    s_arr = arrays_simplex_set(dl_arr)
    total = sum(np.asarray(a).shape[0] for a in dl_arr)
    assert total == len(s_arr), "fill_delaunay_arrays emitted duplicate simplices"
    assert s_arr == s_alpha, "arrays simplex set differs from alpha"


@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("n", [10, 200, 800])
@pytest.mark.parametrize("exact", EXACTS)
def test_delaunay_matches_gudhi(dim, n, exact):
    gudhi = pytest.importorskip("gudhi")
    rng = np.random.default_rng(4000 * n + 10 * dim + int(exact))
    pts = rng.random((n, dim))

    ac = gudhi.AlphaComplex(points=pts)
    st = ac.create_simplex_tree()
    s_gudhi = {frozenset(s) for s, _ in st.get_simplices()}

    s_arr = arrays_simplex_set(diode.fill_delaunay_arrays(pts, exact=exact))
    assert s_arr == s_gudhi, (
        f"delaunay/gudhi simplex sets differ (delaunay {len(s_arr)}, gudhi {len(s_gudhi)})")


def test_delaunay_arrays_empty_dims_and_dtype():
    # A handful of points still yields all dimensions present in 2D/3D.
    pts = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0]])
    verts_by_dim = diode.fill_delaunay_arrays(pts)
    assert len(verts_by_dim) == 4  # vertices, edges, triangles, tetrahedron
    assert [np.asarray(a).shape for a in verts_by_dim] == [(4, 1), (6, 2), (4, 3), (1, 4)]


def test_delaunay_functions_exist():
    assert hasattr(diode, "fill_delaunay")
    assert hasattr(diode, "fill_delaunay_arrays")


# ---- periodic: combinatorics vs the periodic alpha path ---------------------
# The periodic triangulation must be representable in 1 sheet; the alpha path
# raises "Cannot convert to 1-sheeted covering" otherwise (too few points for the
# domain). 3D needs noticeably more points than 2D. Both paths share this guard,
# so we compare only when the alpha path succeeds.
PERIODIC_SIZES = {2: [20, 100, 500], 3: [200, 800]}


@pytest.mark.parametrize("dim", DIMS)
@pytest.mark.parametrize("exact", EXACTS)
def test_periodic_delaunay_matches_alpha_simplex_set(dim, exact):
    frm, to = [0.0] * dim, [1.0] * dim
    for n in PERIODIC_SIZES[dim]:
        rng = np.random.default_rng(6000 * n + 10 * dim + int(exact))
        pts = rng.random((n, dim))
        try:
            s_alpha = {frozenset(int(v) for v in verts)
                       for verts, _ in diode.fill_periodic_alpha_shapes(pts, exact, frm, to)}
        except RuntimeError:
            # not representable in 1 sheet -- the delaunay path must also raise
            with pytest.raises(RuntimeError):
                diode.fill_periodic_delaunay(pts, exact, frm, to)
            continue

        dl_list = diode.fill_periodic_delaunay(pts, exact, frm, to)
        s_list = list_simplex_set(dl_list)
        assert len(dl_list) == len(s_list), "periodic fill_delaunay emitted duplicates"
        assert s_list == s_alpha, (
            f"periodic delaunay/alpha sets differ at dim={dim} n={n} "
            f"(delaunay {len(s_list)}, alpha {len(s_alpha)})")

        dl_arr = diode.fill_periodic_delaunay_arrays(pts, exact, frm, to)
        s_arr = arrays_simplex_set(dl_arr)
        assert sum(np.asarray(a).shape[0] for a in dl_arr) == len(s_arr), \
            "periodic fill_delaunay_arrays emitted duplicates"
        assert s_arr == s_alpha

        # every simplex is non-degenerate, and the alternating face count is the
        # Euler characteristic of the d-torus, which is 0.
        counts = [int(np.asarray(a).shape[0]) for a in dl_arr]
        euler = sum((-1) ** d * c for d, c in enumerate(counts))
        assert euler == 0, f"periodic Euler characteristic {euler} != 0 (d-torus) at dim={dim} n={n}"


def test_periodic_delaunay_functions_exist():
    assert hasattr(diode, "fill_periodic_delaunay")
    assert hasattr(diode, "fill_periodic_delaunay_arrays")
