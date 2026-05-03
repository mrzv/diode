#pragma once

#include <array>
#include <tuple>
#include <vector>

#include <boost/range/adaptor/map.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>

#include <CGAL/version_macros.h>

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
#include <CGAL/Periodic_3_regular_triangulation_traits_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#endif


namespace diode
{

// SimplexCallback contract:
//   - For functions WITHOUT _with_attachment, the callback is invoked as
//         add_simplex(sigma_vertices, alpha)
//     where sigma_vertices is a std::array<unsigned, D> for D in {1,2,3,4}
//     and alpha is a double (the filtration value).
//   - For functions WITH _with_attachment, the callback is invoked as
//         add_simplex(sigma_vertices, alpha, tau_vertices)
//     where tau_vertices is a std::array<unsigned, D'> for D' in {1,2,3,4}
//     listing the vertices of a simplex tau whose own squared circumradius
//     (smallest enclosing sphere through tau's own vertices) equals alpha.
//     For Gabriel sigma, tau == sigma. For non-Gabriel sigma, tau is a Gabriel
//     coface of sigma.
template<bool exact = false>
struct AlphaShapes
{
    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_with_attachment(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_weighted_alpha_shapes(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);


#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
    template<class Points, class SimplexCallback>
    static void fill_weighted_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                                    std::array<double, 3> from, std::array<double, 3> to);
#endif

    template<class Points>
    static std::array<typename Points::Real, 3> circumcenter(const Points& points);
};

template<class Points, class SimplexCallback>
void fill_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex);

template<class Points, class SimplexCallback>
void fill_alpha_shapes2d_with_attachment(const Points& points, const SimplexCallback& add_simplex);

template<class Points, class SimplexCallback>
void fill_periodic_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex,
                                std::array<double, 2> from, std::array<double, 2> to);



}

#include "diode.hpp"
