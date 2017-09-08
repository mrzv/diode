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

template<class T>
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

typename AddSimplex<void>::Simplices
fill_alpha_shape(py::array a)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 3)
    {
        if (a.dtype() == py::dtype::of<float>())
        {
            AddSimplex<float>::Simplices filtration;
            diode::fill_alpha_shapes(ArrayWrapper<float>(a), AddSimplex<float>(&filtration));
            return filtration;
        }
        else if (a.dtype() == py::dtype::of<double>())
        {
            AddSimplex<double>::Simplices filtration;
            diode::fill_alpha_shapes(ArrayWrapper<double>(a), AddSimplex<double>(&filtration));
            return filtration;
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 3D alpha shapes");
}

typename AddSimplex<void>::Simplices
fill_weighted_alpha_shape(py::array a)
{
    if (a.ndim() != 2)
        throw std::runtime_error("Unknown input dimension: can only process 2D arrays");

    if (a.shape()[1] == 4)
    {
        if (a.dtype() == py::dtype::of<float>())
        {
            AddSimplex<float>::Simplices filtration;
            diode::fill_weighted_alpha_shapes(ArrayWrapper<float>(a), AddSimplex<float>(&filtration));
            return filtration;
        }
        else if (a.dtype() == py::dtype::of<double>())
        {
            AddSimplex<double>::Simplices filtration;
            diode::fill_weighted_alpha_shapes(ArrayWrapper<double>(a), AddSimplex<double>(&filtration));
            return filtration;
        }
        else
            throw std::runtime_error("Unknown array dtype");
    }
    else
        throw std::runtime_error("Can only handle 3D alpha shapes (input must be a 4-column array: coordinates + weight)");
}

PYBIND11_PLUGIN(diode)
{
    py::module m("diode", "DioDe pythonn bindings");

    using namespace pybind11::literals;
    m.def("fill_alpha_shapes",  &fill_alpha_shape,
          "data"_a,
          "returns (sorted) alpha shape filtration of the input points");
    m.def("fill_weighted_alpha_shapes",  &fill_weighted_alpha_shape,
          "data"_a,
          "returns (sorted) alpha shape filtration of the weighted input points");

    return m.ptr();
}

