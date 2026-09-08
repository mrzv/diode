#pragma once

#include <array>
#include <tuple>
#include <vector>

#include <geogram/basic/common.h>


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

    // Geogram-backed implementation. The `exact` template parameter is retained
    // for source compatibility; Geogram always uses exact Delaunay predicates
    // and double-precision constructions.
    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_direct(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_with_attachment(const Points& points, const SimplexCallback& add_simplex);

    // Attachment output records a Gabriel coface whose orthosphere determines
    // the simplex filtration value.
    template<class Points, class SimplexCallback>
    static void fill_alpha_shapes_direct_with_attachment(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_weighted_alpha_shapes(const Points& points, const SimplexCallback& add_simplex);

    // Weighted regular triangulation; input columns are x, y, z, weight.
    template<class Points, class SimplexCallback>
    static void fill_weighted_alpha_shapes_direct(const Points& points, const SimplexCallback& add_simplex);

    template<class Points, class SimplexCallback>
    static void fill_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);

    // Periodic 3D alpha complex over the supplied rectangular domain.
    template<class Points, class SimplexCallback>
    static void fill_periodic_alpha_shapes_direct(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);

    // Combinatorics-only 3D export. Emits every simplex by input vertex index
    // without computing alpha values. Full-dimensional input has the same
    // simplex set as fill_alpha_shapes.
    template<class Points, class SimplexCallback>
    static void fill_delaunay(const Points& points, const SimplexCallback& add_simplex);

    // Combinatorics-only export (3D, weighted): the regular-triangulation simplices
    // (== the weighted alpha-complex simplex set) by vertex index, without alpha
    // values. Input is a 4-column array (x, y, z, weight). Redundant (hidden)
    // weighted points are absent. Callback: add_simplex(vertices).
    template<class Points, class SimplexCallback>
    static void fill_weighted_delaunay(const Points& points, const SimplexCallback& add_simplex);

    // Combinatorics-only export (3D, weighted, periodic): like
    // fill_weighted_delaunay, using Geogram's native periodic regular
    // triangulation over [from, to]. Each canonical simplex is emitted once.
    // Hidden sites are omitted; non-one-sheeted coverings raise an error.
    // Callback: add_simplex(vertices).
    template<class Points, class SimplexCallback>
    static void fill_weighted_periodic_delaunay(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);

    // Combinatorics-only export (3D, unweighted, periodic): like fill_delaunay,
    // using Geogram's native periodic Delaunay triangulation over [from, to].
    // Each canonical simplex is emitted once, matching the periodic alpha paths.
    // No alpha values. Callback: add_simplex(vertices).
    template<class Points, class SimplexCallback>
    static void fill_periodic_delaunay(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);

    // Offset-aware periodic Delaunay export (3D, unweighted). The callback is
    // add_simplex(vertices, offsets); sorted vertex ids and aligned offsets are
    // normalized so the first sorted vertex has offset zero.
    template<class Points, class SimplexCallback>
    static void fill_periodic_delaunay_lifts(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);


    template<class Points, class SimplexCallback>
    static void fill_weighted_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                                    std::array<double, 3> from, std::array<double, 3> to);

    template<class Points, class SimplexCallback>
    static void fill_weighted_periodic_alpha_shapes_direct(const Points& points, const SimplexCallback& add_simplex,
                                                    std::array<double, 3> from, std::array<double, 3> to);

    template<class Points>
    static std::array<typename Points::Real, 3> circumcenter(const Points& points);
};

// `exact` is retained for API compatibility. Geogram's triangulations use
// robust exact predicates and double-precision constructions.
// Only ordinary 2D triangulations are supported; Geogram's native periodic
// triangulations are 3D-only.
template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex);

// Geogram-backed ordinary 2D alpha complex. fill_alpha_shapes2d delegates to
// this implementation. Simplices are emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d_direct(const Points& points, const SimplexCallback& add_simplex);

template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d_with_attachment(const Points& points, const SimplexCallback& add_simplex);

// Ordinary 2D alpha complex with the same attacher contract as the 3D path.
// fill_alpha_shapes2d_with_attachment delegates here. Simplices emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_alpha_shapes2d_direct_with_attachment(const Points& points, const SimplexCallback& add_simplex);

// Combinatorics-only export (2D, unweighted): builds the same ordinary Geogram
// Delaunay triangulation as fill_alpha_shapes2d_direct and emits every finite
// simplex (faces, edges, vertices) by vertex index, WITHOUT computing alpha
// values. Same simplex set as fill_alpha_shapes2d.
// Callback: add_simplex(vertices). Simplices emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_delaunay2d(const Points& points, const SimplexCallback& add_simplex);



}

#include "diode.hpp"
