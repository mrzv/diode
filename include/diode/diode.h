#pragma once

#include <array>
#include <tuple>
#include <vector>

#include <boost/range/adaptor/map.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

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

    // Faster equivalent of fill_alpha_shapes for 3D unweighted alpha shapes: builds
    // a plain Delaunay_triangulation_3 (vertex index stored in vertex info, O(1)
    // lookup) and assigns alpha = squared circumradius directly via Edelsbrunner's
    // algorithm (is_Gabriel + min-over-cofaces), instead of constructing the full
    // CGAL::Alpha_shape_3 spectrum. Produces the same (simplex, alpha) set. Simplices
    // are emitted unsorted; the caller sorts if a filtration order is needed.
    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_direct(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_with_attachment(const Points& points, const SimplexCallback& add_simplex);

    // Faster equivalent of fill_alpha_shapes_with_attachment (3D unweighted), built
    // on the Delaunay-direct path. Emits the same attacher tau: a Gabriel coface
    // whose own squared circumradius equals alpha. Computed as a byproduct of the
    // Edelsbrunner propagation (Gabriel sigma -> tau = sigma; non-Gabriel -> the
    // determining coface). Simplices emitted unsorted.
    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_direct_with_attachment(const Points& points, const SimplexCallback& add_simplex);

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

// 2D alpha shapes select the kernel from `exact`, like the 3D AlphaShapes<exact>
// methods: exact=true uses CGAL's exact-construction kernel (EPECK), exact=false
// the inexact-construction kernel (EPICK; exact predicates, fast double values).
template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex);

// Faster equivalent of fill_alpha_shapes2d: Delaunay_triangulation_2 with the
// input index in vertex info (O(1) lookup) and the face circumradius cached in
// face info, instead of a std::set<Simplex2D> with per-edge recomputation.
// Produces the same (simplex, alpha) set. Simplices are emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d_direct(const Points& points, const SimplexCallback& add_simplex);

template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d_with_attachment(const Points& points, const SimplexCallback& add_simplex);

// Faster equivalent of fill_alpha_shapes2d_with_attachment, built on the 2D
// Delaunay-direct path. Same attacher contract. Simplices emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d_direct_with_attachment(const Points& points, const SimplexCallback& add_simplex);

template<bool exact, class Points, class SimplexCallback>
void fill_periodic_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex,
                                std::array<double, 2> from, std::array<double, 2> to);



}

#include "diode.hpp"
