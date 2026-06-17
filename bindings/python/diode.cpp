#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <diode/diode.h>

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

py::object
fill_alpha_shape(py::array a, bool exact, bool with_attachment)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 3)
    {
        if (with_attachment)
        {
            AddSimplexWithAttachment::Simplices filtration;
            if (a.dtype().is(py::dtype::of<float>()))
            {
                if (exact)
                    diode::AlphaShapes<true>::fill_alpha_shapes_with_attachment(ArrayWrapper<float>(a), AddSimplexWithAttachment(&filtration));
                else
                    diode::AlphaShapes<false>::fill_alpha_shapes_with_attachment(ArrayWrapper<float>(a), AddSimplexWithAttachment(&filtration));
            }
            else if (a.dtype().is(py::dtype::of<double>()))
            {
                if (exact)
                    diode::AlphaShapes<true>::fill_alpha_shapes_with_attachment(ArrayWrapper<double>(a), AddSimplexWithAttachment(&filtration));
                else
                    diode::AlphaShapes<false>::fill_alpha_shapes_with_attachment(ArrayWrapper<double>(a), AddSimplexWithAttachment(&filtration));
            }
            else
                throw std::runtime_error("Unknown array dtype");
            return py::cast(filtration);
        }
        else
        {
            // Faster Delaunay-direct path (plain Delaunay_triangulation_3 + direct
            // squared-circumradius), equivalent to the Alpha_shape_3 path but
            // emits simplices unsorted, so sort to keep the filtration order.
            AddSimplex::Simplices filtration;
            if (a.dtype().is(py::dtype::of<float>()))
            {
                if (exact)
                    diode::AlphaShapes<true>::fill_alpha_shapes_direct(ArrayWrapper<float>(a), AddSimplex(&filtration));
                else
                    diode::AlphaShapes<false>::fill_alpha_shapes_direct(ArrayWrapper<float>(a), AddSimplex(&filtration));
            }
            else if (a.dtype().is(py::dtype::of<double>()))
            {
                if (exact)
                    diode::AlphaShapes<true>::fill_alpha_shapes_direct(ArrayWrapper<double>(a), AddSimplex(&filtration));
                else
                    diode::AlphaShapes<false>::fill_alpha_shapes_direct(ArrayWrapper<double>(a), AddSimplex(&filtration));
            }
            else
                throw std::runtime_error("Unknown array dtype");
            sort_filtration(filtration);
            return py::cast(filtration);
        }
    } else if (a.shape()[1] == 2)
    {
        if (with_attachment)
        {
            AddSimplexWithAttachment::Simplices filtration;
            if (a.dtype().is(py::dtype::of<float>()))
                diode::fill_alpha_shapes2d_with_attachment(ArrayWrapper<float>(a), AddSimplexWithAttachment(&filtration));
            else if (a.dtype().is(py::dtype::of<double>()))
                diode::fill_alpha_shapes2d_with_attachment(ArrayWrapper<double>(a), AddSimplexWithAttachment(&filtration));
            else
                throw std::runtime_error("Unknown array dtype");
            sort_filtration(filtration);
            return py::cast(filtration);
        }
        else
        {
            AddSimplex::Simplices filtration;
            if (a.dtype().is(py::dtype::of<float>()))
                diode::fill_alpha_shapes2d(ArrayWrapper<float>(a), AddSimplex(&filtration));
            else if (a.dtype().is(py::dtype::of<double>()))
                diode::fill_alpha_shapes2d(ArrayWrapper<double>(a), AddSimplex(&filtration));
            else
                throw std::runtime_error("Unknown array dtype");
            sort_filtration(filtration);
            return py::cast(filtration);
        }
    }
    else
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
}

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

py::object
fill_periodic_alpha_shape(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_, bool with_attachment)
{
    if (with_attachment)
    {
        PyErr_SetString(PyExc_NotImplementedError,
                        "with_attachment is not yet supported for periodic alpha shapes");
        throw py::error_already_set();
    }

    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 3)
    {
        std::array<double,3> from { from_[0], from_[1], from_[2] },
                             to   { to_[0],   to_[1],   to_[2]   };

        if (a.dtype().is(py::dtype::of<float>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_periodic_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration), from, to);
            else
                diode::AlphaShapes<false>::fill_periodic_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration), from, to);
            return py::cast(filtration);
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            else
                diode::AlphaShapes<false>::fill_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            return py::cast(filtration);
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else if (a.shape()[1] == 2)
    {
        std::array<double,2> from { from_[0], from_[1] },
                             to   { to_[0],   to_[1]   };

        AddSimplex::Simplices filtration;
        if (a.dtype().is(py::dtype::of<float>()))
            diode::fill_periodic_alpha_shapes2d(ArrayWrapper<float>(a), AddSimplex(&filtration),from,to);
        else if (a.dtype().is(py::dtype::of<double>()))
            diode::fill_periodic_alpha_shapes2d(ArrayWrapper<double>(a), AddSimplex(&filtration),from,to);
        else
            throw std::runtime_error("Unknown array dtype");

        sort_filtration(filtration);
        return py::cast(filtration);
    }
    else
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
}

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
    m.def("fill_weighted_alpha_shapes",  &fill_weighted_alpha_shape,
          "data"_a, "exact"_a = false, "with_attachment"_a = false,
          "returns (sorted) alpha shape filtration of the weighted input points (with_attachment=True is not yet supported)");
    m.def("fill_periodic_alpha_shapes",  &fill_periodic_alpha_shape,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "with_attachment"_a = false,
          "returns (sorted) alpha shape filtration of the input points on a periodic domain (with_attachment=True is not yet supported)");
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

