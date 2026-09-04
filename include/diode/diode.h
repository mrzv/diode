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
    // fill_weighted_delaunay, using a 3x3x3 tiled regular triangulation over
    // [from, to]. Each canonical simplex is emitted once. Repeated real vertices
    // within a tiled cell are omitted. No alpha values.
    // Callback: add_simplex(vertices).
    template<class Points, class SimplexCallback>
    static void fill_weighted_periodic_delaunay(const Points& points, const SimplexCallback& add_simplex,
                                    std::array<double, 3> from, std::array<double, 3> to);

    // Combinatorics-only export (3D, unweighted, periodic): like fill_delaunay,
    // using a 3x3x3 tiled Delaunay triangulation over the cuboid [from, to].
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

// Faster equivalent of fill_periodic_alpha_shapes2d: same periodic geometry and
// Gabriel test, but caches face circumradii in a hash map keyed by vertex set
// (first-wins, matching the std::set dedup) instead of a std::set with per-edge
// recomputation. Produces the same (simplex, alpha) set. Emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_periodic_alpha_shapes2d_direct(const Points& points, const SimplexCallback& add_simplex,
                                std::array<double, 2> from, std::array<double, 2> to);

// Combinatorics-only export (2D, unweighted): builds the same
// Delaunay_triangulation_2 as fill_alpha_shapes2d_direct (vertex index in vertex
// info) and emits every finite simplex (faces, edges, vertices) by vertex index,
// WITHOUT computing any alpha value. Same simplex set as fill_alpha_shapes2d.
// Callback: add_simplex(vertices). Simplices emitted unsorted.
template<bool exact, class Points, class SimplexCallback>
void fill_delaunay2d(const Points& points, const SimplexCallback& add_simplex);

// Combinatorics-only export (2D, unweighted, periodic): like fill_delaunay2d but
// on a Periodic_2_Delaunay_triangulation_2 over [from, to]. Each canonical simplex
// is emitted once (deduplicated by vertex-index set). No alpha values.
// Callback: add_simplex(vertices).
template<bool exact, class Points, class SimplexCallback>
void fill_periodic_delaunay2d(const Points& points, const SimplexCallback& add_simplex,
                                std::array<double, 2> from, std::array<double, 2> to);

// Offset-aware counterpart of fill_periodic_delaunay2d. Vertex ids are sorted
// and aligned offsets are normalized by a common lattice translation.
template<bool exact, class Points, class SimplexCallback>
void fill_periodic_delaunay2d_lifts(const Points& points, const SimplexCallback& add_simplex,
                                std::array<double, 2> from, std::array<double, 2> to);



}

#include "diode.hpp"
