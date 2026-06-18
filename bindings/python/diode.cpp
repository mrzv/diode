#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <diode/diode.h>

// Validate a periodic box before it reaches CGAL: it must have at least `dim`
// entries and a strictly positive extent on every axis. Otherwise CGAL's
// Iso_cuboid_3/Iso_rectangle is empty or inverted and the triangulation crashes
// the interpreter deep inside (no Python exception). Works for std::vector and
// std::array operands.
template<class Vec>
static void check_periodic_domain(const Vec& from_, const Vec& to_, std::size_t dim)
{
    if (from_.size() < dim || to_.size() < dim)
        throw std::runtime_error("from/to must have at least as many entries as the point dimension");
    for (std::size_t k = 0; k < dim; ++k)
        if (!(from_[k] < to_[k]))
            throw std::runtime_error("periodic domain is empty or inverted: require from[k] < to[k] on every axis");
}

// present a numpy array in a way that diode understands
template<class T>
struct ArrayWrapper
{
    using Real = T;

            ArrayWrapper(const py::array_t<T>& a_):
                a(a_)           {}

    size_t  size() const                                { return a.shape()[0]; }
    T       operator()(size_t i, size_t j) const        { return a.at(i,j); }

    const py::array_t<T>& a;
};

struct AddSimplex
{
    using Simplices = std::vector<std::tuple<std::vector<unsigned>, double>>;

        AddSimplex(Simplices* result_): result(result_)     {}

    template<unsigned long D>
    void operator()(const std::array<unsigned, D>& vertices, double a) const
    {
        std::vector<unsigned> verts(vertices.begin(), vertices.end());
        result->emplace_back(verts, a);
    }

    Simplices* result;
};

struct AddSimplexWithAttachment
{
    using Simplices = std::vector<std::tuple<std::vector<unsigned>, double, std::vector<unsigned>>>;

        AddSimplexWithAttachment(Simplices* result_): result(result_)     {}

    template<unsigned long D>
    void operator()(const std::array<unsigned, D>& vertices, double a, const std::vector<unsigned>& tau) const
    {
        std::vector<unsigned> verts(vertices.begin(), vertices.end());
        result->emplace_back(verts, a, tau);
    }

    Simplices* result;
};

template<class Simplices>
void sort_filtration(Simplices& filtration)
{
    using Simplex = typename Simplices::value_type;
    std::sort(filtration.begin(), filtration.end(), [](const Simplex& x, const Simplex& y)
    {
        auto xv = std::get<1>(x);
        auto yv = std::get<1>(y);

        if (xv < yv) return true;
        if (xv > yv) return false;

        auto& xvert = std::get<0>(x);
        auto& yvert = std::get<0>(y);

        if (xvert.size() < yvert.size()) return true;
        if (xvert.size() > yvert.size()) return false;

        return std::lexicographical_compare(xvert.begin(), xvert.end(), yvert.begin(), yvert.end());
    });
}

// ===========================================================================
// fill_alpha_shapes_arrays: hand the alpha filtration back as per-dimension NumPy
// arrays without materializing one Python object per simplex. Runs the fast
// Delaunay-direct traversal and buckets each simplex into flat, zero-copy buffers.
// ===========================================================================

// move a std::vector onto the heap and expose its buffer as a NumPy array with
// no copy: the array keeps the vector alive through a capsule base.
template<class T>
py::array_t<T> vector_to_numpy(std::vector<T>&& v, std::vector<py::ssize_t> shape)
{
    auto* held = new std::vector<T>(std::move(v));
    py::capsule owner(held, [](void* p) { delete reinterpret_cast<std::vector<T>*>(p); });
    return py::array_t<T>(std::move(shape), held->data(), owner);
}

// run the combinatorics-only Delaunay traversal (2D or 3D, dispatched on
// exact/dtype), feeding each simplex to a callback cb(std::array<unsigned,D>) with
// no alpha value. Mirrors run_alpha_direct_traversal but calls fill_delaunay /
// fill_delaunay2d, which skip all Gabriel/circumradius work.
template<class Cb>
void run_delaunay_traversal(py::array a, bool exact, const Cb& cb)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 3)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            if (exact) diode::AlphaShapes<true >::fill_delaunay(ArrayWrapper<float >(a), cb);
            else       diode::AlphaShapes<false>::fill_delaunay(ArrayWrapper<float >(a), cb);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            if (exact) diode::AlphaShapes<true >::fill_delaunay(ArrayWrapper<double>(a), cb);
            else       diode::AlphaShapes<false>::fill_delaunay(ArrayWrapper<double>(a), cb);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else if (a.shape()[1] == 2)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            if (exact) diode::fill_delaunay2d<true >(ArrayWrapper<float >(a), cb);
            else       diode::fill_delaunay2d<false>(ArrayWrapper<float >(a), cb);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            if (exact) diode::fill_delaunay2d<true >(ArrayWrapper<double>(a), cb);
            else       diode::fill_delaunay2d<false>(ArrayWrapper<double>(a), cb);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 2D or 3D Delaunay triangulations");
}

// run the fast Delaunay-direct traversal (2D or 3D, dispatched on exact/dtype),
// feeding each simplex to a generic callback cb(std::array<unsigned,D>, double).
template<class Cb>
void run_alpha_direct_traversal(py::array a, bool exact, const Cb& cb)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 3)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            if (exact) diode::AlphaShapes<true >::fill_alpha_shapes_direct(ArrayWrapper<float >(a), cb);
            else       diode::AlphaShapes<false>::fill_alpha_shapes_direct(ArrayWrapper<float >(a), cb);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            if (exact) diode::AlphaShapes<true >::fill_alpha_shapes_direct(ArrayWrapper<double>(a), cb);
            else       diode::AlphaShapes<false>::fill_alpha_shapes_direct(ArrayWrapper<double>(a), cb);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else if (a.shape()[1] == 2)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            if (exact) diode::fill_alpha_shapes2d_direct<true >(ArrayWrapper<float >(a), cb);
            else       diode::fill_alpha_shapes2d_direct<false>(ArrayWrapper<float >(a), cb);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            if (exact) diode::fill_alpha_shapes2d_direct<true >(ArrayWrapper<double>(a), cb);
            else       diode::fill_alpha_shapes2d_direct<false>(ArrayWrapper<double>(a), cb);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
}

// Bucket simplices by dimension into flat row-major vertex buffers (d+1 ids per
// simplex) + per-dimension value buffers. Grouping by dimension up front suits
// consumers that build a filtration dimension-by-dimension. Vertex ids are int64.
struct AddSimplexArrays
{
    using Verts = std::array<std::vector<long>, 4>;
    using Vals  = std::array<std::vector<double>, 4>;
    Verts* verts;
    Vals*  vals;

    template<unsigned long D>
    void operator()(const std::array<unsigned, D>& vertices, double a) const
    {
        auto& vb = (*verts)[D - 1];
        for (unsigned v : vertices)
            vb.push_back(static_cast<long>(v));
        (*vals)[D - 1].push_back(a);
    }
};

// returns (verts_by_dim, vals_by_dim): per-dimension NumPy arrays. verts_by_dim[d]
// is (n_d, d+1) int64, vals_by_dim[d] is (n_d,) float64. Unsorted within a
// dimension -- the consumer sorts by value.
py::object
fill_alpha_shapes_arrays(py::array a, bool exact)
{
    AddSimplexArrays::Verts verts;
    AddSimplexArrays::Vals  vals;
    run_alpha_direct_traversal(a, exact, AddSimplexArrays { &verts, &vals });

    int max_dim = -1;
    for (int d = 0; d < 4; ++d)
        if (!vals[d].empty())
            max_dim = d;

    py::list verts_by_dim, vals_by_dim;
    for (int d = 0; d <= max_dim; ++d)
    {
        py::ssize_t n = static_cast<py::ssize_t>(vals[d].size());
        py::ssize_t w = d + 1;
        verts_by_dim.append(vector_to_numpy(std::move(verts[d]), { n, w }));
        vals_by_dim.append (vector_to_numpy(std::move(vals[d]),  { n }));
    }
    return py::make_tuple(verts_by_dim, vals_by_dim);
}

// ===========================================================================
// Combinatorics-only Delaunay exporters (no alpha values). Same shape as the
// alpha exporters above but built on run_delaunay_traversal, which skips every
// Gabriel/circumradius evaluation. Intended for consumers that recompute the
// filtration values themselves (e.g. a differentiable Cech-Delaunay filtration).
// ===========================================================================

// per-dimension flat vertex buffers (d+1 ids per simplex), no value buffers.
struct AddSimplexArraysNoVal
{
    using Verts = std::array<std::vector<long>, 4>;
    Verts* verts;

    template<unsigned long D>
    void operator()(const std::array<unsigned, D>& vertices) const
    {
        auto& vb = (*verts)[D - 1];
        for (unsigned v : vertices)
            vb.push_back(static_cast<long>(v));
    }
};

// collect simplices as a flat list of vertex lists (one Python list per simplex).
struct AddSimplexNoVal
{
    using Simplices = std::vector<std::vector<unsigned>>;
    Simplices* result;

    template<unsigned long D>
    void operator()(const std::array<unsigned, D>& vertices) const
    {
        result->emplace_back(vertices.begin(), vertices.end());
    }
};

// returns verts_by_dim: per-dimension NumPy arrays, verts_by_dim[d] an (n_d, d+1)
// int64 array of vertex ids. No values. Unsorted within a dimension.
py::object
fill_delaunay_arrays(py::array a, bool exact)
{
    AddSimplexArraysNoVal::Verts verts;
    run_delaunay_traversal(a, exact, AddSimplexArraysNoVal { &verts });

    int max_dim = -1;
    for (int d = 0; d < 4; ++d)
        if (!verts[d].empty())
            max_dim = d;

    py::list verts_by_dim;
    for (int d = 0; d <= max_dim; ++d)
    {
        py::ssize_t w = d + 1;
        py::ssize_t n = static_cast<py::ssize_t>(verts[d].size()) / w;
        verts_by_dim.append(vector_to_numpy(std::move(verts[d]), { n, w }));
    }
    return verts_by_dim;
}

// returns a flat list of vertex lists (no values), one entry per Delaunay simplex.
py::object
fill_delaunay(py::array a, bool exact)
{
    AddSimplexNoVal::Simplices f;
    run_delaunay_traversal(a, exact, AddSimplexNoVal { &f });
    return py::cast(f);
}

// run the periodic combinatorics-only Delaunay traversal (2D or 3D), feeding each
// simplex to cb(std::array<unsigned,D>) with no alpha value.
template<class Cb>
void run_periodic_delaunay_traversal(py::array a, bool exact,
                                     std::vector<double> from_, std::vector<double> to_, const Cb& cb)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");
    auto cols = a.shape()[1];
    if (cols != 2 && cols != 3)
        throw std::runtime_error("Can only handle 2D or 3D Delaunay triangulations");
    check_periodic_domain(from_, to_, static_cast<std::size_t>(cols));
    bool is_float  = a.dtype().is(py::dtype::of<float>());
    bool is_double = a.dtype().is(py::dtype::of<double>());
    if (!is_float && !is_double)
        throw std::runtime_error("Unknown array dtype");

    if (cols == 3)
    {
        std::array<double,3> from { from_[0], from_[1], from_[2] }, to { to_[0], to_[1], to_[2] };
        auto run = [&](auto etag) {
            constexpr bool E = decltype(etag)::value;
            if (is_float) diode::AlphaShapes<E>::fill_periodic_delaunay(ArrayWrapper<float >(a), cb, from, to);
            else          diode::AlphaShapes<E>::fill_periodic_delaunay(ArrayWrapper<double>(a), cb, from, to);
        };
        if (exact) run(std::true_type{}); else run(std::false_type{});
    }
    else
    {
        std::array<double,2> from { from_[0], from_[1] }, to { to_[0], to_[1] };
        auto run = [&](auto etag) {
            constexpr bool E = decltype(etag)::value;
            if (is_float) diode::fill_periodic_delaunay2d<E>(ArrayWrapper<float >(a), cb, from, to);
            else          diode::fill_periodic_delaunay2d<E>(ArrayWrapper<double>(a), cb, from, to);
        };
        if (exact) run(std::true_type{}); else run(std::false_type{});
    }
}

// returns verts_by_dim (per-dim (n_d, d+1) int64 arrays) for the periodic Delaunay
// triangulation, no values.
py::object
fill_periodic_delaunay_arrays(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_)
{
    AddSimplexArraysNoVal::Verts verts;
    run_periodic_delaunay_traversal(a, exact, from_, to_, AddSimplexArraysNoVal { &verts });

    int max_dim = -1;
    for (int d = 0; d < 4; ++d)
        if (!verts[d].empty())
            max_dim = d;

    py::list verts_by_dim;
    for (int d = 0; d <= max_dim; ++d)
    {
        py::ssize_t w = d + 1;
        py::ssize_t n = static_cast<py::ssize_t>(verts[d].size()) / w;
        verts_by_dim.append(vector_to_numpy(std::move(verts[d]), { n, w }));
    }
    return verts_by_dim;
}

// returns a flat list of vertex lists for the periodic Delaunay triangulation.
py::object
fill_periodic_delaunay(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_)
{
    AddSimplexNoVal::Simplices f;
    run_periodic_delaunay_traversal(a, exact, from_, to_, AddSimplexNoVal { &f });
    return py::cast(f);
}

// ---- weighted combinatorics-only (regular triangulation, no alpha values) ----

template<class Cb>
void run_weighted_delaunay_traversal(py::array a, bool exact, const Cb& cb)
{
    if (a.ndim() != 2 || a.shape()[1] != 4)
        throw std::runtime_error("weighted input must be a 4-column array (x, y, z, weight)");
    if (a.dtype().is(py::dtype::of<float>())) {
        if (exact) diode::AlphaShapes<true >::fill_weighted_delaunay(ArrayWrapper<float >(a), cb);
        else       diode::AlphaShapes<false>::fill_weighted_delaunay(ArrayWrapper<float >(a), cb);
    } else if (a.dtype().is(py::dtype::of<double>())) {
        if (exact) diode::AlphaShapes<true >::fill_weighted_delaunay(ArrayWrapper<double>(a), cb);
        else       diode::AlphaShapes<false>::fill_weighted_delaunay(ArrayWrapper<double>(a), cb);
    } else
        throw std::runtime_error("Unknown array dtype");
}

py::object
fill_weighted_delaunay_arrays(py::array a, bool exact)
{
    AddSimplexArraysNoVal::Verts verts;
    run_weighted_delaunay_traversal(a, exact, AddSimplexArraysNoVal { &verts });
    int max_dim = -1;
    for (int d = 0; d < 4; ++d) if (!verts[d].empty()) max_dim = d;
    py::list verts_by_dim;
    for (int d = 0; d <= max_dim; ++d) {
        py::ssize_t w = d + 1;
        py::ssize_t n = static_cast<py::ssize_t>(verts[d].size()) / w;
        verts_by_dim.append(vector_to_numpy(std::move(verts[d]), { n, w }));
    }
    return verts_by_dim;
}

py::object
fill_weighted_delaunay(py::array a, bool exact)
{
    AddSimplexNoVal::Simplices f;
    run_weighted_delaunay_traversal(a, exact, AddSimplexNoVal { &f });
    return py::cast(f);
}

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
template<class Cb>
void run_weighted_periodic_delaunay_traversal(py::array a, bool exact,
                                              std::vector<double> from_, std::vector<double> to_, const Cb& cb)
{
    if (a.ndim() != 2 || a.shape()[1] != 4)
        throw std::runtime_error("weighted input must be a 4-column array (x, y, z, weight)");
    check_periodic_domain(from_, to_, 3);
    std::array<double,3> from { from_[0], from_[1], from_[2] }, to { to_[0], to_[1], to_[2] };
    auto run = [&](auto etag) {
        constexpr bool E = decltype(etag)::value;
        if (a.dtype().is(py::dtype::of<float>()))
            diode::AlphaShapes<E>::fill_weighted_periodic_delaunay(ArrayWrapper<float >(a), cb, from, to);
        else if (a.dtype().is(py::dtype::of<double>()))
            diode::AlphaShapes<E>::fill_weighted_periodic_delaunay(ArrayWrapper<double>(a), cb, from, to);
        else
            throw std::runtime_error("Unknown array dtype");
    };
    if (exact) run(std::true_type{}); else run(std::false_type{});
}

py::object
fill_weighted_periodic_delaunay_arrays(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_)
{
    AddSimplexArraysNoVal::Verts verts;
    run_weighted_periodic_delaunay_traversal(a, exact, from_, to_, AddSimplexArraysNoVal { &verts });
    int max_dim = -1;
    for (int d = 0; d < 4; ++d) if (!verts[d].empty()) max_dim = d;
    py::list verts_by_dim;
    for (int d = 0; d <= max_dim; ++d) {
        py::ssize_t w = d + 1;
        py::ssize_t n = static_cast<py::ssize_t>(verts[d].size()) / w;
        verts_by_dim.append(vector_to_numpy(std::move(verts[d]), { n, w }));
    }
    return verts_by_dim;
}

py::object
fill_weighted_periodic_delaunay(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_)
{
    AddSimplexNoVal::Simplices f;
    run_weighted_periodic_delaunay_traversal(a, exact, from_, to_, AddSimplexNoVal { &f });
    return py::cast(f);
}
#endif

// Slow == true selects the reference implementations kept for testing
// (CGAL::Alpha_shape_3 in 3D, the std::set-based path in 2D). Slow == false
// selects the fast Delaunay-direct paths. if constexpr keeps both compiled.
template<bool Slow, bool Exact, class T>
void run_alpha(const py::array& a, AddSimplex::Simplices& f, bool is2d)
{
    if (is2d) {
        if constexpr (Slow) diode::fill_alpha_shapes2d<Exact>(ArrayWrapper<T>(a), AddSimplex(&f));
        else                diode::fill_alpha_shapes2d_direct<Exact>(ArrayWrapper<T>(a), AddSimplex(&f));
    } else {
        if constexpr (Slow) diode::AlphaShapes<Exact>::fill_alpha_shapes(ArrayWrapper<T>(a), AddSimplex(&f));
        else                diode::AlphaShapes<Exact>::fill_alpha_shapes_direct(ArrayWrapper<T>(a), AddSimplex(&f));
    }
}

template<bool Slow, bool Exact, class T>
void run_alpha_attach(const py::array& a, AddSimplexWithAttachment::Simplices& f, bool is2d)
{
    if (is2d) {
        if constexpr (Slow) diode::fill_alpha_shapes2d_with_attachment<Exact>(ArrayWrapper<T>(a), AddSimplexWithAttachment(&f));
        else                diode::fill_alpha_shapes2d_direct_with_attachment<Exact>(ArrayWrapper<T>(a), AddSimplexWithAttachment(&f));
    } else {
        if constexpr (Slow) diode::AlphaShapes<Exact>::fill_alpha_shapes_with_attachment(ArrayWrapper<T>(a), AddSimplexWithAttachment(&f));
        else                diode::AlphaShapes<Exact>::fill_alpha_shapes_direct_with_attachment(ArrayWrapper<T>(a), AddSimplexWithAttachment(&f));
    }
}

template<bool Slow>
py::object
fill_alpha_shape_impl(py::array a, bool exact, bool with_attachment)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");
    auto cols = a.shape()[1];
    if (cols != 2 && cols != 3)
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
    bool is2d = (cols == 2);

    bool is_float  = a.dtype().is(py::dtype::of<float>());
    bool is_double = a.dtype().is(py::dtype::of<double>());
    if (!is_float && !is_double)
        throw std::runtime_error("Unknown array dtype");

    if (with_attachment)
    {
        AddSimplexWithAttachment::Simplices f;
        if (is_float) { if (exact) run_alpha_attach<Slow, true,  float >(a, f, is2d);
                        else       run_alpha_attach<Slow, false, float >(a, f, is2d); }
        else          { if (exact) run_alpha_attach<Slow, true,  double>(a, f, is2d);
                        else       run_alpha_attach<Slow, false, double>(a, f, is2d); }
        sort_filtration(f);
        return py::cast(f);
    }
    else
    {
        AddSimplex::Simplices f;
        if (is_float) { if (exact) run_alpha<Slow, true,  float >(a, f, is2d);
                        else       run_alpha<Slow, false, float >(a, f, is2d); }
        else          { if (exact) run_alpha<Slow, true,  double>(a, f, is2d);
                        else       run_alpha<Slow, false, double>(a, f, is2d); }
        sort_filtration(f);
        return py::cast(f);
    }
}

py::object
fill_alpha_shape(py::array a, bool exact, bool with_attachment)
{ return fill_alpha_shape_impl<false>(a, exact, with_attachment); }

py::object
fill_alpha_shape_slow(py::array a, bool exact, bool with_attachment)
{ return fill_alpha_shape_impl<true>(a, exact, with_attachment); }

// Slow == true selects the CGAL::Alpha_shape_3 reference (on the regular
// triangulation); Slow == false the fast weighted Delaunay-direct path.
template<bool Slow>
py::object
fill_weighted_alpha_shape_impl(py::array a, bool exact, bool with_attachment)
{
    if (with_attachment)
    {
        PyErr_SetString(PyExc_NotImplementedError,
                        "with_attachment is not yet supported for weighted alpha shapes");
        throw py::error_already_set();
    }

    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");
    if (a.shape()[1] != 4)
        throw std::runtime_error("Can only handle 3D alpha shapes (input must be a 4-column array: coordinates + weight)");

    bool is_float  = a.dtype().is(py::dtype::of<float>());
    bool is_double = a.dtype().is(py::dtype::of<double>());
    if (!is_float && !is_double)
        throw std::runtime_error("Unknown array dtype");

    AddSimplex::Simplices f;
    auto run = [&](auto etag) {
        constexpr bool E = decltype(etag)::value;
        if constexpr (Slow) {
            if (is_float) diode::AlphaShapes<E>::fill_weighted_alpha_shapes(ArrayWrapper<float >(a), AddSimplex(&f));
            else          diode::AlphaShapes<E>::fill_weighted_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&f));
        } else {
            if (is_float) diode::AlphaShapes<E>::fill_weighted_alpha_shapes_direct(ArrayWrapper<float >(a), AddSimplex(&f));
            else          diode::AlphaShapes<E>::fill_weighted_alpha_shapes_direct(ArrayWrapper<double>(a), AddSimplex(&f));
        }
    };
    if (exact) run(std::true_type{}); else run(std::false_type{});
    sort_filtration(f);
    return py::cast(f);
}

py::object
fill_weighted_alpha_shape(py::array a, bool exact, bool with_attachment)
{ return fill_weighted_alpha_shape_impl<false>(a, exact, with_attachment); }

py::object
fill_weighted_alpha_shape_slow(py::array a, bool exact, bool with_attachment)
{ return fill_weighted_alpha_shape_impl<true>(a, exact, with_attachment); }

// Slow == true selects the reference 2D periodic path (std::set); Slow == false
// the fast Delaunay-direct one. 3D periodic still uses CGAL::Alpha_shape_3 for
// both (no direct 3D periodic path yet).
template<bool Slow>
py::object
fill_periodic_alpha_shape_impl(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_, bool with_attachment)
{
    if (with_attachment)
    {
        PyErr_SetString(PyExc_NotImplementedError,
                        "with_attachment is not yet supported for periodic alpha shapes");
        throw py::error_already_set();
    }
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");
    auto cols = a.shape()[1];
    if (cols != 2 && cols != 3)
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
    check_periodic_domain(from_, to_, static_cast<std::size_t>(cols));
    bool is_float  = a.dtype().is(py::dtype::of<float>());
    bool is_double = a.dtype().is(py::dtype::of<double>());
    if (!is_float && !is_double)
        throw std::runtime_error("Unknown array dtype");

    AddSimplex::Simplices f;
    if (cols == 3)
    {
        std::array<double,3> from { from_[0], from_[1], from_[2] }, to { to_[0], to_[1], to_[2] };
        auto run = [&](auto etag) {
            constexpr bool E = decltype(etag)::value;
            if constexpr (Slow) {
                if (is_float) diode::AlphaShapes<E>::fill_periodic_alpha_shapes(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
                else          diode::AlphaShapes<E>::fill_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
            } else {
                if (is_float) diode::AlphaShapes<E>::fill_periodic_alpha_shapes_direct(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
                else          diode::AlphaShapes<E>::fill_periodic_alpha_shapes_direct(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
            }
        };
        if (exact) run(std::true_type{}); else run(std::false_type{});
    }
    else
    {
        std::array<double,2> from { from_[0], from_[1] }, to { to_[0], to_[1] };
        auto run = [&](auto etag) {
            constexpr bool E = decltype(etag)::value;
            if constexpr (Slow) {
                if (is_float) diode::fill_periodic_alpha_shapes2d<E>(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
                else          diode::fill_periodic_alpha_shapes2d<E>(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
            } else {
                if (is_float) diode::fill_periodic_alpha_shapes2d_direct<E>(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
                else          diode::fill_periodic_alpha_shapes2d_direct<E>(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
            }
        };
        if (exact) run(std::true_type{}); else run(std::false_type{});
    }
    sort_filtration(f);
    return py::cast(f);
}

py::object
fill_periodic_alpha_shape(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_, bool with_attachment)
{ return fill_periodic_alpha_shape_impl<false>(a, exact, from_, to_, with_attachment); }

py::object
fill_periodic_alpha_shape_slow(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_, bool with_attachment)
{ return fill_periodic_alpha_shape_impl<true>(a, exact, from_, to_, with_attachment); }

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
// Slow == true selects CGAL::Alpha_shape_3 on the periodic regular triangulation;
// Slow == false the fast weighted periodic Delaunay-direct path.
template<bool Slow>
py::object
fill_weighted_periodic_alpha_shape_impl(py::array a, bool exact, std::array<double,3> from, std::array<double,3> to, bool with_attachment)
{
    if (with_attachment)
    {
        PyErr_SetString(PyExc_NotImplementedError,
                        "with_attachment is not yet supported for weighted periodic alpha shapes");
        throw py::error_already_set();
    }

    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");
    if (a.shape()[1] != 4)
        throw std::runtime_error("Can only handle 3D alpha shapes (input must be a 4-column array: coordinates + weight)");
    check_periodic_domain(from, to, 3);

    bool is_float  = a.dtype().is(py::dtype::of<float>());
    bool is_double = a.dtype().is(py::dtype::of<double>());
    if (!is_float && !is_double)
        throw std::runtime_error("Unknown array dtype");

    AddSimplex::Simplices f;
    auto run = [&](auto etag) {
        constexpr bool E = decltype(etag)::value;
        if constexpr (Slow) {
            if (is_float) diode::AlphaShapes<E>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
            else          diode::AlphaShapes<E>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
        } else {
            if (is_float) diode::AlphaShapes<E>::fill_weighted_periodic_alpha_shapes_direct(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
            else          diode::AlphaShapes<E>::fill_weighted_periodic_alpha_shapes_direct(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
        }
    };
    if (exact) run(std::true_type{}); else run(std::false_type{});
    sort_filtration(f);
    return py::cast(f);
}

py::object
fill_weighted_periodic_alpha_shape(py::array a, bool exact, std::array<double,3> from, std::array<double,3> to, bool with_attachment)
{ return fill_weighted_periodic_alpha_shape_impl<false>(a, exact, from, to, with_attachment); }

py::object
fill_weighted_periodic_alpha_shape_slow(py::array a, bool exact, std::array<double,3> from, std::array<double,3> to, bool with_attachment)
{ return fill_weighted_periodic_alpha_shape_impl<true>(a, exact, from, to, with_attachment); }
#endif

py::array
circumcenter(py::array a, bool exact)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (!(a.shape()[0] == 3 || a.shape()[0] == 4) || a.shape()[1] != 3)
        throw std::runtime_error("Expected 3 or 4 points in 3D");

    if (a.dtype().is(py::dtype::of<float>()))
    {
        std::array<float,3> center;
        if (exact)
            center = diode::AlphaShapes<true>::circumcenter(ArrayWrapper<float>(a));
        else
            center = diode::AlphaShapes<false>::circumcenter(ArrayWrapper<float>(a));
        return py::array_t<float>(3, center.data());
    }
    else if (a.dtype().is(py::dtype::of<double>()))
    {
        std::array<double,3> center;
        if (exact)
            center = diode::AlphaShapes<true>::circumcenter(ArrayWrapper<double>(a));
        else
            center = diode::AlphaShapes<false>::circumcenter(ArrayWrapper<double>(a));
        return py::array_t<double>(3, center.data());
    }
    else
        throw std::runtime_error("Unknown array dtype");
}


PYBIND11_MODULE(diode, m)
{
    m.doc() = "DioDe python bindings";

    using namespace pybind11::literals;
    m.def("fill_alpha_shapes",  &fill_alpha_shape,
          "data"_a, "exact"_a = false, "with_attachment"_a = false,
          "Returns the (sorted) alpha shape filtration of the input points.\n"
          "\n"
          "If with_attachment=False (default), each entry is a pair\n"
          "    (sigma_vertices, alpha)\n"
          "where alpha is the filtration value of simplex sigma.\n"
          "\n"
          "If with_attachment=True, each entry is a triple\n"
          "    (sigma_vertices, alpha, tau_vertices)\n"
          "where tau_vertices are the vertices of a Gabriel coface tau of sigma\n"
          "such that the squared circumradius of tau equals alpha. For Gabriel\n"
          "sigma, tau == sigma. The triple form is intended for differentiable\n"
          "alpha-filtration pipelines where alpha values must be expressed as\n"
          "smooth functions of point coordinates.");
    m.def("fill_alpha_shapes_slow",  &fill_alpha_shape_slow,
          "data"_a, "exact"_a = false, "with_attachment"_a = false,
          "Reference implementation of fill_alpha_shapes kept for testing: 3D uses\n"
          "CGAL::Alpha_shape_3, 2D the std::set-based path. Same result as\n"
          "fill_alpha_shapes but much slower; used to cross-check the fast path.");
    m.def("fill_alpha_shapes_arrays", &fill_alpha_shapes_arrays,
          "data"_a, "exact"_a = false,
          "Alpha shape filtration as per-dimension NumPy arrays: returns\n"
          "(verts_by_dim, vals_by_dim) where verts_by_dim[d] is an (n_d, d+1) int64\n"
          "array of vertex ids and vals_by_dim[d] an (n_d,) float64 array of alpha\n"
          "values. Uses the fast Delaunay-direct traversal; avoids one Python object\n"
          "per simplex. The consumer sorts within each dimension.");
    m.def("fill_delaunay_arrays", &fill_delaunay_arrays,
          "data"_a, "exact"_a = false,
          "Delaunay simplices (== the alpha-complex simplex set for full-dimensional\n"
          "input) as per-dimension NumPy arrays WITHOUT alpha values: returns\n"
          "verts_by_dim where verts_by_dim[d] is an (n_d, d+1) int64 array of vertex\n"
          "ids. Skips all Gabriel/circumradius work, so it is faster than\n"
          "fill_alpha_shapes_arrays. Intended for consumers that recompute filtration\n"
          "values themselves (e.g. a differentiable Cech-Delaunay filtration). For\n"
          "degenerate (collinear/coplanar) input this returns the lower-dimensional\n"
          "Delaunay complex (fill_alpha_shapes returns nothing there). Unsorted within\n"
          "each dimension.");
    m.def("fill_delaunay", &fill_delaunay,
          "data"_a, "exact"_a = false,
          "Delaunay simplices as a flat list of vertex lists, WITHOUT alpha values.\n"
          "List form of fill_delaunay_arrays (see it for the degenerate-input note).");
    m.def("fill_periodic_delaunay_arrays", &fill_periodic_delaunay_arrays,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "Periodic Delaunay simplices as per-dimension NumPy arrays WITHOUT alpha\n"
          "values (periodic counterpart of fill_delaunay_arrays). Each canonical\n"
          "simplex is emitted once. verts_by_dim[d] is an (n_d, d+1) int64 array.");
    m.def("fill_periodic_delaunay", &fill_periodic_delaunay,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "Periodic Delaunay simplices as a flat list of vertex lists, WITHOUT alpha\n"
          "values. List form of fill_periodic_delaunay_arrays.");
    m.def("fill_weighted_delaunay_arrays", &fill_weighted_delaunay_arrays,
          "data"_a, "exact"_a = false,
          "Regular-triangulation (weighted Delaunay) simplices as per-dimension NumPy\n"
          "arrays WITHOUT alpha values, for a 4-column input (x, y, z, weight). Same\n"
          "shape as fill_delaunay_arrays. Hidden (redundant) weighted points are absent.");
    m.def("fill_weighted_delaunay", &fill_weighted_delaunay,
          "data"_a, "exact"_a = false,
          "Weighted Delaunay (regular triangulation) simplices as a flat list of vertex\n"
          "lists, WITHOUT alpha values. List form of fill_weighted_delaunay_arrays.");
#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
    m.def("fill_weighted_periodic_delaunay_arrays", &fill_weighted_periodic_delaunay_arrays,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "Periodic weighted Delaunay (regular triangulation) simplices as per-dimension\n"
          "NumPy arrays WITHOUT alpha values. Each canonical simplex appears once.");
    m.def("fill_weighted_periodic_delaunay", &fill_weighted_periodic_delaunay,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "Periodic weighted Delaunay simplices as a flat list of vertex lists, WITHOUT\n"
          "alpha values. List form of fill_weighted_periodic_delaunay_arrays.");
#endif
    m.def("fill_weighted_alpha_shapes",  &fill_weighted_alpha_shape,
          "data"_a, "exact"_a = false, "with_attachment"_a = false,
          "returns (sorted) alpha shape filtration of the weighted input points "
          "(4-column array x,y,z,weight). Uses the fast weighted Delaunay-direct "
          "path (Regular_triangulation_3 + Edelsbrunner). with_attachment=True is "
          "not yet supported.");
    m.def("fill_weighted_alpha_shapes_slow",  &fill_weighted_alpha_shape_slow,
          "data"_a, "exact"_a = false, "with_attachment"_a = false,
          "Reference implementation of fill_weighted_alpha_shapes kept for testing: "
          "uses CGAL::Alpha_shape_3 on the regular triangulation. Same result as "
          "fill_weighted_alpha_shapes but slower.");
    m.def("fill_periodic_alpha_shapes",  &fill_periodic_alpha_shape,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "with_attachment"_a = false,
          "returns (sorted) alpha shape filtration of the input points on a periodic domain (with_attachment=True is not yet supported)");
    m.def("fill_periodic_alpha_shapes_slow",  &fill_periodic_alpha_shape_slow,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "with_attachment"_a = false,
          "Reference periodic alpha filtration kept for testing: 2D uses the\n"
          "std::set-based path (3D uses CGAL::Alpha_shape_3, same as the fast one).");
#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
    m.def("fill_weighted_periodic_alpha_shapes",  &fill_weighted_periodic_alpha_shape,
          "data"_a, "exact"_a = false,
          "from"_a = std::array<double,3> {0.,0.,0.},
          "to"_a   = std::array<double,3> {1.,1.,1.},
          "with_attachment"_a = false,
          "returns (sorted) alpha shape filtration of the weighted input points on a "
          "periodic domain. Uses the fast weighted periodic Delaunay-direct path "
          "(Periodic_3_regular_triangulation_3 + Edelsbrunner). with_attachment=True "
          "is not yet supported.");
    m.def("fill_weighted_periodic_alpha_shapes_slow",  &fill_weighted_periodic_alpha_shape_slow,
          "data"_a, "exact"_a = false,
          "from"_a = std::array<double,3> {0.,0.,0.},
          "to"_a   = std::array<double,3> {1.,1.,1.},
          "with_attachment"_a = false,
          "Reference implementation of fill_weighted_periodic_alpha_shapes kept for "
          "testing: uses CGAL::Alpha_shape_3 on the periodic regular triangulation.");
#endif
    m.def("circumcenter", &circumcenter, "points"_a, "exact"_a = false, "returns circumcenter of the intput points");
}

