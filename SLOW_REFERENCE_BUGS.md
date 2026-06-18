# Known bugs in the `_slow` reference implementations

The `alpha-merged` branch replaced `CGAL::Alpha_shape_3` with direct
Delaunay/Regular-triangulation + Edelsbrunner code for the public
`fill_*` functions, and kept the original implementations as `*_slow`
(`fill_alpha_shapes_slow`, `fill_weighted_alpha_shapes_slow`,
`fill_periodic_alpha_shapes_slow`, `fill_weighted_periodic_alpha_shapes_slow`).
The test suite uses those `_slow` functions as the oracle.

Two independent reviews found that the `_slow` references contain
**pre-existing bugs** (present on `master`, not introduced by this branch). In
each case below the **fast path is correct** -- verified against GUDHI's
`AlphaComplex` as a neutral third party -- and the `_slow` reference is wrong.
This document records them so the oracle isn't mistaken for ground truth, and so
they can be fixed (or knowingly tolerated) in the reference code.

The test suite now pins these specific cases against GUDHI rather than `_slow`
(`tests/test_alpha_fast_vs_slow.py::test_2d_hull_vertices_match_gudhi`,
`::test_weighted_shared_coords_match_gudhi`).

Environment for the repros below:

```
PY=/path/to/python; SP=/path/to/site-packages   # with numpy, gudhi
PYTHONPATH=<diode build>/bindings/python:$SP $PY
```

---

## 1. 2D slow path drops a hull vertex

**Where:** `include/diode/diode.hpp`, `fill_alpha_shapes2d` (the unweighted 2D
`_slow` path, around line 1486). In its vertex-emit loop the simplex index
`s[0]` is assigned **only inside** a conditional that matches a face vertex to
the current vertex's point (`s[0] = point_map[cur->face()->vertex(i)]`, ~line
1623). For some hull vertices that match never fires, so `s[0]` keeps a stale/
default value, the 0-simplex collides with another in the dedup `std::set`, and
a vertex silently disappears. The same pattern appears in the 2D
with-attachment and 2D periodic `_slow` paths (~lines 1812, 2104).

**Repro (7 distinct points; exact kernel):**

```python
import numpy as np, diode, gudhi
pts = np.array([[2.,1.],[0.,0.],[1.,0.],[2.,2.],[0.,2.],[1.,2.],[2.,0.]])
nverts = lambda f: sorted(tuple(s) for s,*_ in f if len(s)==1)
slow = diode.fill_alpha_shapes_slow(pts, True)
fast = diode.fill_alpha_shapes(pts, True)
g = gudhi.AlphaComplex(points=pts.tolist()).create_simplex_tree()
print("slow  verts:", nverts(slow))                       # 6 -- vertex (1,) missing
print("fast  verts:", nverts(fast))                       # 7
print("gudhi verts:", sorted((s[0],) for s,_ in g.get_simplices() if len(s)==1))  # 7
```

`fast` and `gudhi` agree (7 vertices); `slow` drops vertex `(1,)`.

---

## 2. Weighted slow path mislabels coincident weighted points

**Where:** `include/diode/diode.hpp`, `fill_weighted_alpha_shapes` (the weighted
`_slow` path, line 688). When two input sites share coordinates but have
different weights, the regular triangulation keeps the larger-weight one, whose
vertex alpha should be `-max(weight)`. The `_slow` path reports the wrong
(smaller) weight. The fast path keys its duplicate handling on the full
`Weighted_point` and is correct.

**Repro:**

```python
import numpy as np, diode, gudhi
data = np.array([[0.,0.,0.,0.01],[1.,0.,0.,0.01],[0.,1.,0.,0.01],
                 [0.,0.,1.,0.01],[0.,0.,0.,0.5]])   # rows 0 and 4 coincide
vv = lambda f: sorted(round(float(a),4) for s,a in f if len(s)==1)
print("slow :", vv(diode.fill_weighted_alpha_shapes_slow(data)))   # [-0.01,-0.01,-0.01,-0.01]
print("fast :", vv(diode.fill_weighted_alpha_shapes(data)))        # [-0.5, -0.01,-0.01,-0.01]
st = gudhi.AlphaComplex(points=data[:,:3].tolist(), weights=data[:,3].tolist()).create_simplex_tree()
print("gudhi:", sorted(round(float(v),4) for s,v in st.get_simplices() if len(s)==1))  # [-0.5,...]
```

`fast` and `gudhi` agree (`-0.5` for the kept site); `slow` reports `-0.01`.

---

## 3. Weighted-periodic slow path emits index-duplicated simplices near the 1-sheet boundary

**Where:** `include/diode/diode.hpp`, `fill_weighted_periodic_alpha_shapes` (the
weighted-periodic `_slow` path, line 1072). For sparse clouds near the periodic
1-sheet boundary the slow path can emit the **same index simplex more than
once** (e.g. `(0,)` and `(0,164)` twice each) over a few near-degenerate input
indices -- a malformed, index-duplicated complex. The fast path deduplicates by
vertex-index set and always emits each canonical simplex once, so the two
disagree there. (Reported by review at seed 4, n=404, weights in `[0, 0.005]`;
exact reproduction is offset-dependent, so it is described rather than pinned by
a test. For non-degenerate clouds the two agree to round-off.)

This is the reason the `fill_weighted_periodic_alpha_shapes_direct` docstring was
softened from "same simplex set" to "same for non-degenerate clouds; may differ
at the sparse 1-sheet boundary, where the slow path may emit index-duplicated
simplices."

---

## Note: inverted periodic domain (fixed, was shared)

A separate issue -- an empty/inverted periodic box (`from >= to` on some axis) --
crashed the interpreter inside CGAL for **both** the fast and slow periodic
paths. This is now guarded at the binding layer (`check_periodic_domain` in
`bindings/python/diode.cpp`), which raises a `RuntimeError` before either path
runs, so it is no longer reachable.

---

## Recommendation

The fast paths are the correct implementations and are what users call. For the
`_slow` references you can either:

1. **Fix them** (1 and 2 are small, localized fixes; 1 is a one-line index fix,
   2 is in the weighted vertex value lookup), which restores them as a faithful
   oracle; or
2. **Leave them and treat them as known-imperfect**, relying on the GUDHI-anchored
   tests added for these specific cases.

Either way the bugs are confined to degenerate / coincident-point / sparse-
periodic inputs; they do not affect generic point clouds, where `_slow` and the
fast paths agree to round-off.
