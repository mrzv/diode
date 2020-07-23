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

template<bool exact = false>
struct AlphaShapes
{
    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes(const Points& points, const SimplexCallback& add_simplex);

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
};

template<class Points, class SimplexCallback>
void fill_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex);

template<class Points, class SimplexCallback>
void fill_periodic_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex,
                                std::array<double, 2> from, std::array<double, 2> to);



}

#include "diode.hpp"
