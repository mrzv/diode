DioDe
=====

DioDe uses `Geogram <https://github.com/BrunoLevy/geogram>`_ to generate
ordinary, weighted, and periodic alpha-shape filtrations in a format that
Dionysus_ understands. Geogram is fetched and built automatically by CMake.

**Geometry backend:** Geogram 1.9.9 (BSD-3-Clause)

The plain 3D ordinary and weighted alpha exporters use compact cell/edge
storage and assign filtration values from tetrahedra down to vertices.
Ordinary triangulations use Geogram's serial ``Delaunay3d`` backend to avoid
PDEL dropping unfinished insertions between BRIO levels. Weighted triangulations
use BPOW, with hidden sites omitted. BRIO-Hilbert reordering remains enabled.
Weighted tetrahedral power centers are computed from the original linear
equations rather than a Gram system.

All periodic 3D exporters use Geogram's native periodic triangulation, including
weighted alpha shapes, Delaunay simplices, lattice offsets, and rectangular
domains. No explicit 27-copy triangulation is constructed by DioDe. Periodic
copies are identified by original vertex IDs only after checking that their
relative lattice offsets agree and that the resulting covering is closed.
Hidden weighted sites are omitted. The bundled Geogram build includes a
correction allowing periodic insertion to continue past hidden weighted sites;
common weight offsets are removed before triangulation to preserve precision.
Inputs that cannot be represented as a simplicial one-sheeted covering raise
an error. Periodic 2D still uses explicit tiling because Geogram 1.9.9 has no
native periodic 2D backend. Ordinary attachment exporters retain the general
implementation; periodic attachments remain unsupported.

These paths are shared by the list and NumPy-array alpha exporters; the
list exporter still sorts the filtration, while arrays remain unsorted.

Get, Build, Install
-------------------

The simplest way to install Diode as a Python package:

.. parsed-literal::

    pip install --verbose diode

or from this repository directly:

.. parsed-literal::

    pip install --verbose `git+https://github.com/mrzv/diode.git <https://github.com/mrzv/diode.git>`_

Alternatively, you can clone and build everything by hand.
To get Diode, either clone its `repository <https://github.com/mrzv/diode>`_:

.. parsed-literal::

    git clone `<https://github.com/mrzv/diode.git>`_

or download it as a `Zip archive <https://github.com/mrzv/diode/archive/master.zip>`_.

To build the project::

    mkdir build
    cd build
    cmake ..
    make

To use the Python bindings, either launch Python from ``.../build/bindings/python`` or add this directory to your ``PYTHONPATH`` variable, by adding::

    export PYTHONPATH=.../build/bindings/python:$PYTHONPATH

to your ``~/.bashrc`` or ``~/.zshrc``.


Usage
-----

See the `exactness note <#exactness>`_ below. Robust predicates are especially
important for degenerate point sets, such as repeated periodic domains.

See `examples/generate_alpha_shape.cpp <https://github.com/mrzv/diode/blob/master/examples/generate_alpha_shape.cpp>`_ and
`examples/generate_weighted_alpha_shape.cpp <https://github.com/mrzv/diode/blob/master/examples/generate_weighted_alpha_shape.cpp>`_ for C++ examples.

In Python, use ``diode.fill_alpha_shapes(...)`` and ``diode.fill_weighted_alpha_shapes(...)`` to fill a list of simplices, together with their alpha values::

    >>> import diode
    >>> import numpy as np

    >>> points = np.random.random((100,3))
    >>> simplices = diode.fill_alpha_shapes(points)

    >>> print(simplices)
     [([13L], 0.0),
      ([18L], 0.0),
      ([59L], 0.0),
      ([10L], 0.0),
      ([72L], 0.0),
      ...,
      ([91L, 4L, 16L, 49L], 546.991052812204),
      ([49L, 62L], 1933.2257381777533),
      ([62L, 34L, 49L], 1933.2257381777533),
      ([62L, 91L, 49L], 1933.2257381777533),
      ([62L, 91L, 34L, 49L], 1933.2257381777533)]

    >>> weighted_points = np.random.random((100,4))
    >>> simplices2 = diode.fill_weighted_alpha_shapes(weighted_points)
    >>> print(simplices2)
    [([24L], -0.987214836816236),
     ([35L], -0.968749877102265),
     ([50L], -0.9673151804059413),
     ([47L], -0.9640549893422644),
     ([71L], -0.9639978806827709),
     ([24L, 50L], -0.9540965704765515),
     ...
     ([54L, 10L], 29223.611044169364),
     ([10L, 54L, 43L], 29223.611044169364),
     ([13L, 10L, 54L], 29223.611044169364),
     ([13L, 10L, 54L, 43L], 29223.611044169364)]

The list can be passed to Dionysus_ to initialize a filtration::

    >>> import dionysus
    >>> f = dionysus.Filtration(simplices)
    >>> print(f)
    Filtration with 2287 simplices

DioDe also includes ``diode.fill_periodic_alpha_shapes(...)``, which generates
the alpha shape for a point set on a periodic cube, by default ``[0,0,0]
- [1,1,1]``. Canonical simplices are emitted once.

    >>> simplices_periodic = diode.fill_periodic_alpha_shapes(points)
    >>> f_periodic = dionysus.Filtration(simplices_periodic)
    >>> print(f_periodic)
    Filtration with 2912 simplices

    >>> for s in f_periodic: print(s)
    <0> 0
    <1> 0
    <2> 0
    <3> 0
    ...
    <77,94,97> 0.0704355
    <46,77,94,97> 0.0708062
    <30,77,94,97> 0.0708474
    <18,65,79> 0.0715833
    <18,64,65,79> 0.0715833
    <18,65,79,99> 0.0725366

.. _Dionysus:   http://mrzv.org/software/dionysus2

``diode.fill_weighted_periodic_alpha_shapes(...)`` generates the alpha shape
for a weighted point set on a periodic cube::

    >>> weighted_points[:,3] /= 64
    >>> simplices_weighted_periodic = diode.fill_weighted_periodic_alpha_shapes(weighted_points)

``diode.circumcenter(...)`` can be used to compute the circumcenter of a tetrahedron in 3D::

    >>> tet = np.random.random((4,3))
    >>> center = diode.circumcenter(tet)
    >>> print(center)
    [-0.00752673  0.14213101  1.0060982 ]


Delaunay combinatorics (no alpha values)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some consumers only need the simplicial complex (the Delaunay triangulation,
which for full-dimensional input is the same simplex set as the alpha complex)
and recompute their own filtration values -- for example a differentiable
Cech-Delaunay filtration that recomputes values as minimum-enclosing-ball radii.
For those,
``diode.fill_delaunay_arrays(...)`` returns just the combinatorics, avoiding
per-simplex Gabriel and orthosphere calculations::

    >>> verts_by_dim = diode.fill_delaunay_arrays(points)

The result is a list of per-dimension NumPy arrays, where ``verts_by_dim[d]`` is
an ``(n_d, d+1)`` int64 array of vertex ids (dimension 0 = vertices, 1 = edges,
and so on). ``diode.fill_delaunay(...)`` is the equivalent list-of-tuples form.
``diode.fill_periodic_delaunay_arrays(...)`` / ``diode.fill_periodic_delaunay(...)``
are the periodic counterparts (over the cube ``[from, to]``, default the unit
cube). All four take the same ``exact`` argument as the alpha-shape functions.

Consumers that need periodic geometry as well as combinatorics can use
``diode.fill_periodic_delaunay_lifts_arrays(...)``::

    >>> vertices, offsets = diode.fill_periodic_delaunay_lifts_arrays(
    ...     points, bbox_min=[0, 0, 0], bbox_max=[1, 1, 1])

``offsets[d]`` is aligned with ``vertices[d]`` and has shape
``(n_d, d+1, ambient_dim)``. A lifted coordinate is
``points[vertices[d]] + offsets[d] * (bbox_max - bbox_min)``. Vertex ids are
sorted within each simplex and the integer offsets are normalized by a common
lattice translation so that the first offset is zero. Points must be inside the
half-open domain ``[bbox_min, bbox_max)``. The exporter validates that each
vertex-id tuple has one coherent relative lattice lift.


Exactness
~~~~~~~~~

All functions retain the ``exact`` argument for API compatibility. Geogram uses
robust exact predicates for Delaunay and regular-triangulation decisions.
Circumcenters and alpha values are returned as doubles, not exact numbers.
The shared 3D sphere solver uses translated long-double intermediates and
adaptive expansion arithmetic for ill-conditioned facets and tetrahedra before
floating-point division. This improves conditioning without providing exact
alpha comparisons. ``exact=False`` and ``exact=True`` select the same backend
and numerical construction type.


License
-------

DioDe is distributed under the BSD 3-Clause License. See ``LICENSE`` for the
complete terms.

Geogram
~~~~~~~

Geogram is also distributed under the BSD 3-Clause License:

Copyright (c) 2000-2022 Inria. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
* Neither the name of Inria nor the names of its contributors may be used to
  endorse or promote products derived from this software without specific
  prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
