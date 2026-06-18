"""Filtration-array exporters for the weighted / periodic / weighted-periodic
alpha shapes.

Each `fill_*_alpha_shapes_arrays` returns the same (simplex, alpha) filtration as
its `fill_*_alpha_shapes` list form, but as per-dimension NumPy arrays
(verts_by_dim[d] an (n_d, d+1) int64 array, vals_by_dim[d] an (n_d,) float64
array) with no Python object per simplex. Both forms run the identical fast
direct C++ path, so the arrays and the list must carry the same (sorted-vertex ->
alpha) map. The unweighted-non-periodic cell is already covered in
test_alpha_distribution_samples.py; this file covers the other three cells plus
the shared dtype / dedup / domain-validation behavior.
"""
import numpy as np
import pytest
import diode


EXACTS = [False, True]


def list_value_dict(filtration):
    """list of (verts, alpha) -> {sorted-vertex-tuple: alpha}."""
    return {tuple(sorted(int(v) for v in verts)): float(a) for verts, a in filtration}


def arrays_value_dict(res):
    """(verts_by_dim, vals_by_dim) -> {sorted-vertex-tuple: alpha}, checking shape,
    dtype and that no simplex is emitted twice."""
    verts_by_dim, vals_by_dim = res
    assert len(verts_by_dim) == len(vals_by_dim)
    out = {}
    n_rows = 0
    for d, (verts, vals) in enumerate(zip(verts_by_dim, vals_by_dim)):
        verts = np.asarray(verts)
        vals = np.asarray(vals)
        assert verts.dtype == np.int64, (d, verts.dtype)
        assert vals.dtype == np.float64, (d, vals.dtype)
        assert verts.ndim == 2 and verts.shape[1] == d + 1, (d, verts.shape)
        assert vals.ndim == 1 and vals.shape[0] == verts.shape[0], (d, vals.shape)
        n_rows += verts.shape[0]
        for row, value in zip(verts, vals):
            out[tuple(sorted(int(v) for v in row))] = float(value)
    assert n_rows == len(out), "arrays exporter emitted duplicate simplices"
    return out


def assert_maps_match(arrays, listed, *, exact, ambiguity=None):
    """The arrays map must equal the list map: both run the identical fast direct
    C++ path. The simplex set always matches exactly. Values match to round-off,
    EXCEPT for periodic clouds, where a few near-degenerate simplices pick a
    different periodic offset representative run-to-run (the documented periodic
    Gabriel offset ambiguity); for those we bound the COUNT of differing values,
    mirroring the fast-vs-slow tests. ambiguity is None / "2d" / "weighted_periodic".
    """
    assert set(arrays) == set(listed), (
        f"simplex sets differ (arrays {len(arrays)}, list {len(listed)})")
    diffs = np.array([abs(arrays[k] - listed[k]) for k in arrays]) if arrays else np.array([0.0])
    if ambiguity is None:
        scale = max(1.0, max((abs(v) for v in listed.values()), default=1.0))
        rtol = 1e-12 if exact else 1e-6
        assert diffs.max() <= rtol * scale + 1e-9, (
            f"arrays/list values differ by {diffs.max():.2e}")
    elif ambiguity == "2d":
        # 2D periodic: the offset representative differs by a possibly large,
        # run-dependent amount on a few edges -- bound the count, not the magnitude.
        assert np.median(diffs) < 1e-9, "most periodic values should match exactly"
        n_big = int((diffs > 1e-7).sum())
        assert n_big <= max(5, len(diffs) // 50), (
            f"{n_big}/{len(diffs)} values differ -- more than offset ambiguity explains")
    elif ambiguity == "weighted_periodic":
        assert np.median(diffs) < 1e-9, "most periodic values should match exactly"
        n_big = int((diffs > 1e-7).sum())
        assert n_big <= max(3, len(diffs) // 100), (
            f"{n_big}/{len(diffs)} values differ -- more than offset ambiguity explains")
        assert diffs.max() < 1e-2
    else:
        raise ValueError(ambiguity)


# ---- weighted (4-column x,y,z,weight), non-periodic -------------------------
@pytest.mark.parametrize("n", [50, 200, 800])
@pytest.mark.parametrize("exact", EXACTS)
def test_weighted_alpha_arrays_match_list(n, exact):
    rng = np.random.default_rng(8000 * n + int(exact))
    data = np.hstack([rng.random((n, 3)), rng.random((n, 1)) * 0.05])
    arrays = arrays_value_dict(diode.fill_weighted_alpha_shapes_arrays(data, exact=exact))
    listed = list_value_dict(diode.fill_weighted_alpha_shapes(data, exact=exact))
    assert_maps_match(arrays, listed, exact=exact)


def test_weighted_alpha_arrays_degenerate_matches_list():
    # Coplanar 4-column input is < 3D: the direct path defers to the reference
    # (CGAL::Alpha_shape_3, which needs a full-dimensional triangulation), so both
    # the list and the arrays form come back empty. Pins that the fallback flows
    # through the arrays sink too.
    data = np.array([[0., 0., 0., 0.01], [1., 0., 0., 0.02],
                     [0., 1., 0., 0.0], [1., 1., 0., 0.03], [0.5, 0.5, 0., 0.01]])
    listed = list_value_dict(diode.fill_weighted_alpha_shapes(data))
    arrays = arrays_value_dict(diode.fill_weighted_alpha_shapes_arrays(data))
    assert arrays == listed


# ---- unweighted periodic (2D and 3D) ---------------------------------------
# The cloud must be representable in one sheet of the periodic covering; both
# forms raise "Cannot convert to 1-sheeted covering" otherwise (3D needs more
# points than 2D). Compare only when the list form succeeds.
PERIODIC_SIZES = {2: [20, 100, 500], 3: [200, 800]}


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("exact", EXACTS)
def test_periodic_alpha_arrays_match_list(dim, exact):
    frm, to = [0.0] * dim, [1.0] * dim
    for n in PERIODIC_SIZES[dim]:
        rng = np.random.default_rng(9000 * n + 10 * dim + int(exact))
        pts = rng.random((n, dim))
        try:
            listed = list_value_dict(diode.fill_periodic_alpha_shapes(pts, exact, frm, to))
        except RuntimeError:
            with pytest.raises(RuntimeError):
                diode.fill_periodic_alpha_shapes_arrays(pts, exact, frm, to)
            continue
        arrays = arrays_value_dict(diode.fill_periodic_alpha_shapes_arrays(pts, exact, frm, to))
        # 2D periodic has the run-to-run offset ambiguity; 3D uses CGAL is_Gabriel
        # (offset-resolving), so the same direct path agrees to round-off.
        assert_maps_match(arrays, listed, exact=exact, ambiguity=("2d" if dim == 2 else None))


# ---- weighted periodic (4-column, 3D) --------------------------------------
# Weights must satisfy 0 <= w < 1/64 * domain^2 (CGAL's periodic-regular
# requirement); *0.01 keeps them small. The cloud must be 1-sheet representable.
@pytest.mark.parametrize("n", [1000, 2500])
@pytest.mark.parametrize("exact", EXACTS)
def test_weighted_periodic_alpha_arrays_match_list(n, exact):
    rng = np.random.default_rng(31 * n + int(exact))
    data = np.hstack([rng.random((n, 3)), rng.random((n, 1)) * 0.01])
    frm, to = [0.] * 3, [1.] * 3
    try:
        listed = list_value_dict(diode.fill_weighted_periodic_alpha_shapes(data, exact, frm, to))
    except RuntimeError:
        with pytest.raises(RuntimeError):
            diode.fill_weighted_periodic_alpha_shapes_arrays(data, exact, frm, to)
        pytest.skip("point cloud not representable in 1 sheet")
    arrays = arrays_value_dict(diode.fill_weighted_periodic_alpha_shapes_arrays(data, exact, frm, to))
    assert_maps_match(arrays, listed, exact=exact, ambiguity="weighted_periodic")


# ---- shared: dtype dispatch, domain validation, presence -------------------
def test_arrays_float32_input_gives_float64_values():
    rng = np.random.default_rng(123)
    data = np.hstack([rng.random((300, 3)), rng.random((300, 1)) * 0.05]).astype(np.float32)
    verts_by_dim, vals_by_dim = diode.fill_weighted_alpha_shapes_arrays(data)
    assert np.asarray(verts_by_dim[0]).dtype == np.int64
    assert np.asarray(vals_by_dim[0]).dtype == np.float64


def test_weighted_arrays_reject_non_4_column():
    pts = np.random.default_rng(0).random((50, 3))  # missing weight column
    with pytest.raises(RuntimeError):
        diode.fill_weighted_alpha_shapes_arrays(pts)


@pytest.mark.parametrize("frm,to", [([1., 1.], [0., 0.]),   # inverted
                                    ([0.], [1., 1.])])       # too few entries
def test_periodic_arrays_bad_domain_raises(frm, to):
    pts = np.random.default_rng(1).random((100, 2))
    with pytest.raises(RuntimeError):
        diode.fill_periodic_alpha_shapes_arrays(pts, False, frm, to)


@pytest.mark.parametrize("frm,to", [([1., 1., 1.], [0., 0., 0.]),  # inverted
                                    ([0., 0.], [1., 1., 1.])])      # too few entries
def test_weighted_periodic_arrays_bad_domain_raises(frm, to):
    rng = np.random.default_rng(2)
    data = np.hstack([rng.random((100, 3)), rng.random((100, 1)) * 0.01])
    with pytest.raises(RuntimeError):
        diode.fill_weighted_periodic_alpha_shapes_arrays(data, False, frm, to)


def test_alpha_arrays_functions_exist():
    for name in ("fill_weighted_alpha_shapes_arrays",
                 "fill_periodic_alpha_shapes_arrays",
                 "fill_weighted_periodic_alpha_shapes_arrays"):
        assert hasattr(diode, name)
