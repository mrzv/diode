#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <diode/diode.h>

// present a numpy array in a way that diode understands
template<class T>
struct ArrayWrapper
{
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

void sort_filtration(AddSimplex::Simplices& filtration)
{
    using Simplex = AddSimplex::Simplices::value_type;
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

AddSimplex::Simplices
fill_alpha_shape(py::array a, bool exact)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 3)
    {
        if (a.dtype().is(py::dtype::of<float>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration));
            else
                diode::AlphaShapes<false>::fill_alpha_shapes(ArrayWrapper<float>(a), AddSimplex(&filtration));
            return filtration;
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration));
            else
                diode::AlphaShapes<false>::fill_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration));
            return filtration;
        }
        else
            throw std::runtime_error("Unknown array dtype");
    } else if (a.shape()[1] == 2)
    {
        AddSimplex::Simplices filtration;
        if (a.dtype().is(py::dtype::of<float>()))
            diode::fill_alpha_shapes2d(ArrayWrapper<float>(a), AddSimplex(&filtration));
        else if (a.dtype().is(py::dtype::of<double>()))
            diode::fill_alpha_shapes2d(ArrayWrapper<double>(a), AddSimplex(&filtration));
        else
            throw std::runtime_error("Unknown array dtype");

        sort_filtration(filtration);
        return filtration;
    }
    else
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
}

AddSimplex::Simplices
fill_weighted_alpha_shape(py::array a, bool exact)
{
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
            return filtration;
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_weighted_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration));
            else
                diode::AlphaShapes<false>::fill_weighted_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration));
            return filtration;
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 3D alpha shapes (input must be a 4-column array: coordinates + weight)");
}

AddSimplex::Simplices
fill_periodic_alpha_shape(py::array a, bool exact, std::vector<double> from_, std::vector<double> to_)
{
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
            return filtration;
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            else
                diode::AlphaShapes<false>::fill_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            return filtration;
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
        return filtration;
    }
    else
        throw std::runtime_error("Can only handle 2D or 3D alpha shapes");
}

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
AddSimplex::Simplices
fill_weighted_periodic_alpha_shape(py::array a, bool exact, std::array<double,3> from, std::array<double,3> to)
{
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
            return filtration;
        }
        else if (a.dtype().is(py::dtype::of<double>()))
        {
            AddSimplex::Simplices filtration;
            if (exact)
                diode::AlphaShapes<true>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            else
                diode::AlphaShapes<false>::fill_weighted_periodic_alpha_shapes(ArrayWrapper<double>(a), AddSimplex(&filtration), from, to);
            return filtration;
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 3D alpha shapes");
}
#endif


PYBIND11_MODULE(diode, m)
{
    m.doc() = "DioDe pythonn bindings";

    using namespace pybind11::literals;
    m.def("fill_alpha_shapes",  &fill_alpha_shape,
          "data"_a, "exact"_a = false,
          "returns (sorted) alpha shape filtration of the input points");
    m.def("fill_weighted_alpha_shapes",  &fill_weighted_alpha_shape,
          "data"_a, "exact"_a = false,
          "returns (sorted) alpha shape filtration of the weighted input points");
    m.def("fill_periodic_alpha_shapes",  &fill_periodic_alpha_shape,
          "data"_a, "exact"_a = false,
          "from"_a = std::vector<double> {0.,0.,0.},
          "to"_a   = std::vector<double> {1.,1.,1.},
          "returns (sorted) alpha shape filtration of the input points on a periodic domain");
#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
    m.def("fill_weighted_periodic_alpha_shapes",  &fill_weighted_periodic_alpha_shape,
          "data"_a, "exact"_a = false,
          "from"_a = std::array<double,3> {0.,0.,0.},
          "to"_a   = std::array<double,3> {1.,1.,1.},
          "returns (sorted) alpha shape filtration of the weighted input points on a periodic domain");
#endif
}

