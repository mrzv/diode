DioDe
=====

DioDe uses CGAL to generate alpha shapes filtrations in a format that Dionysus_
understands. DioDe is not integrated into Dionysus_ because of licensing
restrictions (Dionysus is under BSD, DioDe is under GPL because of its
dependence on CGAL).

**Dependencies:**

* `CGAL <http://www.cgal.org/>`_

Get, Build, Install
-------------------

The simplest way to install Dionysus as a Python package:

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

See `examples/generate_alpha_shape.cpp <https://github.com/mrzv/diode/blob/master/examples/generate_alpha_shape.cpp>`_ for a C++ example.

In Python, use ``diode.fill_alpha_shapes(...)`` to fill a list of simplices, together with their alpha values::

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

The list can be passed to Dionysus_ to initialize a filtration::

    >>> import dionysus
    >>> f = dionysus.Filtration(simplices)
    >>> print(f)
    Filtration with 2287 simplices

.. _Dionysus:   http://mrzv.org/software/dionysus

