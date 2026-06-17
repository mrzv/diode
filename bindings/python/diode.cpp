#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <diode/diode.h>

// Optional capsule path: build Oineus CellWithValue<Simplex> objects directly.
// Enabled when the build is pointed at Oineus's headers (-DOINEUS_INCLUDE_DIR),
// which defines DIODE_HAVE_OINEUS. The NumPy array path needs none of this.
#ifdef DIODE_HAVE_OINEUS
#include <oineus/simplex.h>
#include <oineus/cell_with_value.h>
#endif

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
// Exporters that hand an alpha filtration to Oineus without materializing one
// Python object per simplex. Both run the fast Delaunay-direct traversal and
// differ only in the per-simplex sink:
//   * fill_alpha_shapes_arrays -- per-dimension NumPy arrays. Zero-copy,
//       framework/library agnostic, needs no Oineus headers. Recommended.
//   * fill_alpha_shapes_cells  -- a PyCapsule owning a std::vector of Oineus
//       CellWithValue<Simplex>; only built when DIODE_HAVE_OINEUS is defined
//       and requires a matching Int/Real/allocator ABI.
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
// simplex) + per-dimension value buffers. Oineus sorts by dimension first, so
// this is exactly the grouping it wants. Vertex ids are int64 (Oineus's Int).
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
// dimension -- Oineus re-sorts by value.
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

#ifdef DIODE_HAVE_OINEUS
// MUST match Oineus's bindings: OINEUS_PYTHON_INT (long), OINEUS_PYTHON_REAL
// (double). Oineus copies the cells out of the capsule, so allocator backends
// need not match (capsule frees its own with the system allocator).
using OinInt     = long;
using OinReal    = double;
using OinSimplex = oineus::Simplex<OinInt>;
using OinCell    = oineus::CellWithValue<OinSimplex, OinReal>;
using OinCellVec = std::vector<OinCell, oineus::JeAllocator<OinCell>>;

static const char* k_oineus_capsule_name = "oineus_cellvector_simplex_long_double_v1";

struct AddSimplexCells
{
    std::array<OinCellVec, 4>* buckets;

    template<unsigned long D>
    void operator()(const std::array<unsigned, D>& vertices, double a) const
    {
        typename OinSimplex::IdxVector vs;
        vs.reserve(D);
        for (unsigned v : vertices)
            vs.push_back(static_cast<OinInt>(v));
        (*buckets)[D - 1].emplace_back(OinSimplex(vs), static_cast<OinReal>(a));
    }
};

// returns a PyCapsule owning a dimension-grouped std::vector<CellWithValue<...>>.
py::object
fill_alpha_shapes_cells(py::array a, bool exact)
{
    auto buckets = std::make_unique<std::array<OinCellVec, 4>>();
    run_alpha_direct_traversal(a, exact, AddSimplexCells { buckets.get() });

    auto* cells = new OinCellVec();
    size_t total = 0;
    for (auto& b : *buckets) total += b.size();
    cells->reserve(total);
    for (auto& b : *buckets)
        for (auto& c : b)
            cells->emplace_back(std::move(c));

    return py::capsule(cells, k_oineus_capsule_name, [](PyObject* o) {
        void* p = PyCapsule_GetPointer(o, "oineus_cellvector_simplex_long_double_v1");
        delete reinterpret_cast<OinCellVec*>(p);
    });
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

py::object
fill_weighted_alpha_shape(py::array a, bool exact, bool with_attachment)
{
    if (with_attachment)
    {
        PyErr_SetString(PyExc_NotImplementedError,
                        "with_attachment is not yet supported for weighted alpha shapes");
        throw py::error_already_set();
    }

    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 4)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_weighted_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration));
            else
                diode::AlphaShapes<false>::fill_weighted_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration));
            return py::cast(filtration);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_weighted_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration));
            else
                diode::AlphaShapes<false>::fill_weighted_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration));
            return py::cast(filtration);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 3D alpha shapes (input must be a 4-column array: coordinates + weight)");
}

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
            if (is_float) diode::AlphaShapes<E>::fill_periodic_alpha_shapes(ArrayWrapper<float >(a), AddSimplex(&f), from, to);
            else          diode::AlphaShapes<E>::fill_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&f), from, to);
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
py::object
fill_weighted_periodic_alpha_shape(py::array a, bool exact, std::array<double,3> from, std::array<double,3> to, bool with_attachment)
{
    if (with_attachment)
    {
        PyErr_SetString(PyExc_NotImplementedError,
                        "with_attachment is not yet supported for weighted periodic alpha shapes");
        throw py::error_already_set();
    }

    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 4)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration), from, to);
            else
                diode::AlphaShapes<false>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration), from, to);
            return py::cast(filtration);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            else
                diode::AlphaShapes<false>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            return py::cast(filtration);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 3D alpha shapes");
}
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
#ifdef DIODE_HAVE_OINEUS
    m.def("fill_alpha_shapes_cells", &fill_alpha_shapes_cells,
          "data"_a, "exact"_a = false,
          "Alpha shape filtration as a PyCapsule owning a\n"
          "std::vector<oineus::CellWithValue<Simplex<long>, double>> (dimension-\n"
          "grouped). Requires a diode built against Oineus headers.");
    m.attr("has_oineus_cells") = true;
#else
    m.attr("has_oineus_cells") = false;
#endif
    m.def("fill_weighted_alpha_shapes",  &fill_weighted_alpha_shape,
          "data"_a, "exact"_a = false, "with_attachment"_a = false,
          "returns (sorted) alpha shape filtration of the weighted input points (with_attachment=True is not yet supported)");
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
          "returns (sorted) alpha shape filtration of the weighted input points on a periodic domain (with_attachment=True is not yet supported)");
#endif
    m.def("circumcenter", &circumcenter, "points"_a, "exact"_a = false, "returns circumcenter of the intput points");
}

