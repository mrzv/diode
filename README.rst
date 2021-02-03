DioDe
=====

DioDe uses CGAL to generate alpha shapes filtrations in a format that Dionysus_
understands. DioDe is not integrated into Dionysus_ because of licensing
restrictions (Dionysus is under BSD, DioDe is under GPL because of its
dependence on CGAL). It supports both ordinary and weighted alpha shapes.

**Dependencies:**

* `CGAL <http://www.cgal.org/>`_

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

NB: a remark below about `using exact computation <#exactness>`_. This issue is especially important when working with degenerate point sets
(e.g., repeated copies of a fundamental domain in a periodic point set).

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
- [1,1,1]``. (In the periodic case, it may happen that CGAL reports each
simplex multiple times. However, passing the result to
``dionysus.Filtration`` will take care of the duplicates.)::

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

When using CGAL version at least 4.11, DioDe includes
``diode.fill_weighted_periodic_alpha_shapes(...)``, which generates the alpha
shape for a weighted point set on a periodic cube::

    >>> weighted_points[:,3] /= 64
    >>> simplices_weighted_periodic = diode.fill_weighted_periodic_alpha_shapes(weighted_points)


Exactness
~~~~~~~~~

All functions take an argument ``exact``, set to ``False`` by default. The argument
determines a choice of the kernel in CGAL
(``Exact_predicates_inexact_constructions_kernel`` vs
``Exact_predicates_exact_constructions_kernel``). ``exact = True`` guarantees
correctness of the output; ``exact = False`` is faster, but can sometimes fail
(not even produce a simplicial complex). It's possible to run the two versions
adaptively by running the default ``exact = False`` version first, and if the
result is not a simplicial complex, then run ``exact = True``. This should be the
best of both worlds.
