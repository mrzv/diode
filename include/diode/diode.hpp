#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <array>
#include <CGAL/Alpha_shape_vertex_base_3.h>

namespace diode
{
namespace detail
{

template<class NT>
double to_floating_point(const NT& x, typename std::enable_if< !std::is_floating_point<NT>::value >::type* = 0)
{ return CGAL::to_double(x.exact()); }

// do nothing, if we already have a floating point type
template<class T>
T to_floating_point(T x, typename std::enable_if< std::is_floating_point<T>::value >::type* = 0)
{ return x; }

template<bool exact>
using Kernel = typename std::conditional<exact,
                                         CGAL::Exact_predicates_exact_constructions_kernel,
                                         CGAL::Exact_predicates_inexact_constructions_kernel>::type;

template<class Delaunay_, class Point_ = typename Delaunay_::Point, class Vertex_ = unsigned>
struct AlphaShapeWrapper
{
    using Vertex     = Vertex_;
    using Delaunay   = Delaunay_;
    using AlphaShape = CGAL::Alpha_shape_3<Delaunay>;
    using Point      = Point_;
    using FT         = typename AlphaShape::FT;
    using PointsMap  = std::map<Point, Vertex>;

    // Helpers: build the vertex list of a CGAL simplex (cell/facet/edge/vertex).
    static std::vector<Vertex> verts_of_cell(typename AlphaShape::Cell_handle c, const PointsMap& points)
    {
        std::vector<Vertex> vs(4);
        for (size_t i = 0; i < 4; ++i)
            vs[i] = points.find(c->vertex(i)->point())->second;
        return vs;
    }

    static std::vector<Vertex> verts_of_facet(const typename AlphaShape::Facet& f, const PointsMap& points)
    {
        std::vector<Vertex> vs;
        vs.reserve(3);
        typename AlphaShape::Cell_handle c = f.first;
        for (size_t i = 0; i < 4; ++i)
            if (static_cast<int>(i) != f.second)
                vs.push_back(points.find(c->vertex(i)->point())->second);
        return vs;
    }

    static std::vector<Vertex> verts_of_edge(const typename AlphaShape::Edge& e, const PointsMap& points)
    {
        typename AlphaShape::Cell_handle c = e.first;
        std::vector<Vertex> vs(2);
        vs[0] = points.find(c->vertex(e.second)->point())->second;
        vs[1] = points.find(c->vertex(e.third)->point())->second;
        return vs;
    }

    static std::vector<Vertex> verts_of_vertex(typename AlphaShape::Vertex_handle v, const PointsMap& points)
    {
        return std::vector<Vertex> { points.find(v->point())->second };
    }

    // Find a Gabriel coface tau of facet f such that squared_circumradius(tau) == alpha.
    // For non-Gabriel facets, tau is one of the at-most-two adjacent finite cells.
    static std::vector<Vertex> find_attacher_for_facet(
        const typename AlphaShape::Facet& f,
        const AlphaShape& as,
        const FT& alpha,
        const PointsMap& points)
    {
        auto fas = as.get_alpha_status(f);
        if (fas.is_Gabriel())
            return verts_of_facet(f, points);

        typename AlphaShape::Cell_handle c0 = f.first;
        typename AlphaShape::Cell_handle c1 = c0->neighbor(f.second);
        if (!as.is_infinite(c0) && c0->get_alpha() == alpha)
            return verts_of_cell(c0, points);
        if (!as.is_infinite(c1) && c1->get_alpha() == alpha)
            return verts_of_cell(c1, points);
        throw std::runtime_error("find_attacher_for_facet: no matching coface");
    }

    // Find a Gabriel coface tau of edge e such that squared_circumradius(tau) == alpha.
    // Search order: Gabriel facets (lower-dim) first, then incident cells.
    static std::vector<Vertex> find_attacher_for_edge(
        const typename AlphaShape::Edge& e,
        const AlphaShape& as,
        const FT& alpha,
        const PointsMap& points)
    {
        auto eas = as.get_alpha_status(e);
        if (eas.is_Gabriel())
            return verts_of_edge(e, points);

        typename AlphaShape::Cell_handle c = e.first;
        int i = e.second;
        int j = e.third;

        // Gabriel facets first (they correspond to Gabriel-facet-contributions in compute_edge_status).
        typename AlphaShape::Facet_circulator fcirc = as.incident_facets(c, i, j);
        typename AlphaShape::Facet_circulator fdone = fcirc;
        do {
            if (!as.is_infinite(*fcirc)) {
                auto fas_it = (*fcirc).first->get_facet_status((*fcirc).second);
                if (fas_it->is_Gabriel() && fas_it->alpha_min() == alpha)
                    return verts_of_facet(*fcirc, points);
            }
        } while (++fcirc != fdone);

        // Then incident cells.
        typename AlphaShape::Cell_circulator ccirc = as.incident_cells(c, i, j);
        typename AlphaShape::Cell_circulator cdone = ccirc;
        do {
            if (!as.is_infinite(ccirc)) {
                if (ccirc->get_alpha() == alpha)
                    return verts_of_cell(ccirc, points);
            }
        } while (++ccirc != cdone);

        throw std::runtime_error("find_attacher_for_edge: no matching coface");
    }

    // Find a Gabriel coface tau of vertex v such that squared_circumradius(tau) == alpha.
    // For unweighted alpha shapes (Tag_false), every vertex is itself Gabriel and
    // the search reduces to tau = v.
    static std::vector<Vertex> find_attacher_for_vertex(
        typename AlphaShape::Vertex_handle v,
        const AlphaShape& as,
        const FT& alpha,
        const PointsMap& points)
    {
        // Fast path for the unweighted case: vertices are always Gabriel.
        // Detected via is_Gabriel() to keep the helper safe in weighted/extended contexts.
        auto* vas = v->get_alpha_status();
        if (vas->is_Gabriel())
            return verts_of_vertex(v, points);

        // Weighted / general case: search incident Gabriel edges, then Gabriel facets,
        // then incident cells. Lower-dim attacher preferred.
        std::vector<typename AlphaShape::Edge> incident_es;
        as.incident_edges(v, std::back_inserter(incident_es));
        for (const auto& e : incident_es)
        {
            if (as.is_infinite(e)) continue;
            auto eas = as.get_alpha_status(e);
            if (eas.is_Gabriel() && eas.alpha_min() == alpha)
                return verts_of_edge(e, points);
        }

        std::vector<typename AlphaShape::Facet> incident_fs;
        as.incident_facets(v, std::back_inserter(incident_fs));
        for (const auto& f : incident_fs)
        {
            if (as.is_infinite(f)) continue;
            auto fas_it = f.first->get_facet_status(f.second);
            if (fas_it->is_Gabriel() && fas_it->alpha_min() == alpha)
                return verts_of_facet(f, points);
        }

        std::vector<typename AlphaShape::Cell_handle> incident_cs;
        as.incident_cells(v, std::back_inserter(incident_cs));
        for (auto c : incident_cs)
        {
            if (as.is_infinite(c)) continue;
            if (c->get_alpha() == alpha)
                return verts_of_cell(c, points);
        }

        throw std::runtime_error("find_attacher_for_vertex: no matching coface");
    }

    // CGAL's peculiar design for its output of the alpha shape filtration requires
    // that once dereferenced this output iterator can be assigned both a
    // CGAL::Object and K::FT. The two are assigned in sequence. This requires
    // keeping state, namely, the last assigned object.
    template<class AddSimplex>
    struct ASOutputIterator3
    {
                          ASOutputIterator3(const AddSimplex& add_simplex_, const PointsMap& points_):
                              add_simplex(add_simplex_), points(points_)       {}

        ASOutputIterator3& operator*()                       { return *this; }
        ASOutputIterator3& operator=(const CGAL::Object& o_) { o = o_; return *this; }       // store the object for later
        ASOutputIterator3& operator=(const FT& a)
        {
            using V = typename AlphaShape::Vertex_handle;
            using E = typename AlphaShape::Edge;
            using F = typename AlphaShape::Facet;
            using C = typename AlphaShape::Cell_handle;

            if (const V* v = CGAL::object_cast<V>(&o))
            {
                auto u = points.find((*v)->point())->second;
                std::array<Vertex, 1> vertices { u };
                add_simplex(vertices, to_floating_point(a));
            }
            else if (const E* e = CGAL::object_cast<E>(&o))
            {
                C c = e->first;
                auto u = points.find(c->vertex(e->second)->point())->second;
                auto v = points.find(c->vertex(e->third)->point())->second;
                std::array<Vertex, 2> vertices { u, v };
                add_simplex(vertices, to_floating_point(a));
            }
            else if (const F* f = CGAL::object_cast<F>(&o))
            {
                std::array<Vertex, 3> vertices;
                size_t j = 0;
                C c = f->first;
                for (size_t i = 0; i < 4; ++i)
                    if (i != f->second)
                        vertices[j++] = points.find(c->vertex(i)->point())->second;
                add_simplex(vertices, to_floating_point(a));
            } else if (const C* c = CGAL::object_cast<C>(&o))
            {
                std::array<Vertex, 4> vertices;
                for (size_t i = 0; i < 4; ++i)
                    vertices[i] = points.find((*c)->vertex(i)->point())->second;
                add_simplex(vertices, to_floating_point(a));
            } else
                throw std::runtime_error("Unknown object type in ASOutputIterator3");

            return *this;
        }

        ASOutputIterator3& operator++()                      { return *this; }
        ASOutputIterator3& operator++(int)                   { return *this; }

        const AddSimplex&   add_simplex;
        const PointsMap&    points;
        CGAL::Object        o;
    };

    // Output iterator that, in addition to the simplex and its alpha value,
    // computes a Gabriel coface tau (the "attacher") and emits its vertex tuple.
    // Callback signature: add_simplex(sigma_vertices, alpha, tau_vertices).
    template<class AddSimplex>
    struct ASOutputIterator3WithAttachment
    {
                          ASOutputIterator3WithAttachment(const AddSimplex& add_simplex_,
                                                          const PointsMap& points_,
                                                          const AlphaShape& as_):
                              add_simplex(add_simplex_), points(points_), as(&as_) {}

        ASOutputIterator3WithAttachment& operator*()                       { return *this; }
        ASOutputIterator3WithAttachment& operator=(const CGAL::Object& o_) { o = o_; return *this; }
        ASOutputIterator3WithAttachment& operator=(const FT& a)
        {
            using V = typename AlphaShape::Vertex_handle;
            using E = typename AlphaShape::Edge;
            using F = typename AlphaShape::Facet;
            using C = typename AlphaShape::Cell_handle;

            if (const V* v = CGAL::object_cast<V>(&o))
            {
                auto u = points.find((*v)->point())->second;
                std::array<Vertex, 1> vertices { u };
                auto tau = find_attacher_for_vertex(*v, *as, a, points);
                add_simplex(vertices, to_floating_point(a), tau);
            }
            else if (const E* e = CGAL::object_cast<E>(&o))
            {
                C c = e->first;
                auto u = points.find(c->vertex(e->second)->point())->second;
                auto vv = points.find(c->vertex(e->third)->point())->second;
                std::array<Vertex, 2> vertices { u, vv };
                auto tau = find_attacher_for_edge(*e, *as, a, points);
                add_simplex(vertices, to_floating_point(a), tau);
            }
            else if (const F* f = CGAL::object_cast<F>(&o))
            {
                std::array<Vertex, 3> vertices;
                size_t j = 0;
                C c = f->first;
                for (size_t i = 0; i < 4; ++i)
                    if (i != f->second)
                        vertices[j++] = points.find(c->vertex(i)->point())->second;
                auto tau = find_attacher_for_facet(*f, *as, a, points);
                add_simplex(vertices, to_floating_point(a), tau);
            } else if (const C* c = CGAL::object_cast<C>(&o))
            {
                std::array<Vertex, 4> vertices;
                for (size_t i = 0; i < 4; ++i)
                    vertices[i] = points.find((*c)->vertex(i)->point())->second;
                auto tau = verts_of_cell(*c, points);   // cells are always Gabriel
                add_simplex(vertices, to_floating_point(a), tau);
            } else
                throw std::runtime_error("Unknown object type in ASOutputIterator3WithAttachment");

            return *this;
        }

        ASOutputIterator3WithAttachment& operator++()                      { return *this; }
        ASOutputIterator3WithAttachment& operator++(int)                   { return *this; }

        const AddSimplex&   add_simplex;
        const PointsMap&    points;
        const AlphaShape*   as;
        CGAL::Object        o;
    };

    template<class AddSimplex>
    static void fill_filtration(Delaunay& dt, const PointsMap& points_map, const AddSimplex& add_simplex)
    {
        AlphaShape as(dt, std::numeric_limits<FT>::infinity(), AlphaShape::GENERAL);
        as.filtration_with_alpha_values(ASOutputIterator3<AddSimplex>(add_simplex, points_map));
    }

    template<class AddSimplex>
    static void fill_filtration_with_attachment(Delaunay& dt, const PointsMap& points_map, const AddSimplex& add_simplex)
    {
        AlphaShape as(dt, std::numeric_limits<FT>::infinity(), AlphaShape::GENERAL);
        as.filtration_with_alpha_values(ASOutputIterator3WithAttachment<AddSimplex>(add_simplex, points_map, as));
    }
};



}   // detail
}   // diode

template<bool exact>
template<class Points, class SimplexCallback>
void
diode::AlphaShapes<exact>::
fill_alpha_shapes(const Points& points, const SimplexCallback& add_simplex)
{
    using K          = detail::Kernel<exact>;
    using Vb         = CGAL::Alpha_shape_vertex_base_3<K>;
    using Fb         = CGAL::Alpha_shape_cell_base_3<K>;
    using TDS        = CGAL::Triangulation_data_structure_3<Vb,Fb>;
    using Delaunay   = CGAL::Delaunay_triangulation_3<K,TDS,CGAL::Fast_location>;

    using ASWrapper  = detail::AlphaShapeWrapper<Delaunay>;
    using PointsMap  = typename ASWrapper::PointsMap;
    using AlphaShape = typename ASWrapper::AlphaShape;
    using Vertex     = typename ASWrapper::Vertex;
    using Point      = typename ASWrapper::Point;

    PointsMap points_map;
    for (Vertex i = 0; i < points.size(); ++i)
    {
        Point p(points(i,0), points(i,1), points(i,2));
        points_map[p] = i;
    }

    auto points_range = points_map | boost::adaptors::map_keys;
    Delaunay dt(std::begin(points_range), std::end(points_range));
    ASWrapper::fill_filtration(dt, points_map, add_simplex);
}

template<bool exact>
template<class Points, class SimplexCallback>
void
diode::AlphaShapes<exact>::
fill_alpha_shapes_with_attachment(const Points& points, const SimplexCallback& add_simplex)
{
    using K          = detail::Kernel<exact>;
    using Vb         = CGAL::Alpha_shape_vertex_base_3<K>;
    using Fb         = CGAL::Alpha_shape_cell_base_3<K>;
    using TDS        = CGAL::Triangulation_data_structure_3<Vb,Fb>;
    using Delaunay   = CGAL::Delaunay_triangulation_3<K,TDS,CGAL::Fast_location>;

    using ASWrapper  = detail::AlphaShapeWrapper<Delaunay>;
    using PointsMap  = typename ASWrapper::PointsMap;
    using AlphaShape = typename ASWrapper::AlphaShape;
    using Vertex     = typename ASWrapper::Vertex;
    using Point      = typename ASWrapper::Point;

    PointsMap points_map;
    for (Vertex i = 0; i < points.size(); ++i)
    {
        Point p(points(i,0), points(i,1), points(i,2));
        points_map[p] = i;
    }

    auto points_range = points_map | boost::adaptors::map_keys;
    Delaunay dt(std::begin(points_range), std::end(points_range));
    ASWrapper::fill_filtration_with_attachment(dt, points_map, add_simplex);
}

template<bool exact>
template<class Points, class SimplexCallback>
void
diode::AlphaShapes<exact>::
fill_weighted_alpha_shapes(const Points& points, const SimplexCallback& add_simplex)
{
    using K         = detail::Kernel<exact>;

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR <= 9) || (CGAL_VERSION_MAJOR < 4)
    using Gt        = CGAL::Regular_triangulation_euclidean_traits_3<K>;
    using Vb        = CGAL::Alpha_shape_vertex_base_3<Gt>;
    using Fb        = CGAL::Alpha_shape_cell_base_3<Gt>;
    using TDS       = CGAL::Triangulation_data_structure_3<Vb,Fb>;
    using Delaunay  = CGAL::Regular_triangulation_3<Gt,TDS>;
#else
    using Rvb       = CGAL::Regular_triangulation_vertex_base_3<K>;
    using Vb        = CGAL::Alpha_shape_vertex_base_3<K, Rvb>;
    using Rcb       = CGAL::Regular_triangulation_cell_base_3<K>;
    using Cb        = CGAL::Alpha_shape_cell_base_3<K, Rcb>;
    using TDS       = CGAL::Triangulation_data_structure_3<Vb,Cb>;
    using Delaunay  = CGAL::Regular_triangulation_3<K,TDS>;
#endif

    using ASWrapper     = detail::AlphaShapeWrapper<Delaunay, typename Delaunay::Weighted_point>;
    using PointsMap     = typename ASWrapper::PointsMap;
    using AlphaShape    = typename ASWrapper::AlphaShape;
    using Vertex        = typename ASWrapper::Vertex;
    using Point         = typename ASWrapper::Point;

    PointsMap points_map;
    for (Vertex i = 0; i < points.size(); ++i)
    {
        Point p({points(i,0), points(i,1), points(i,2)}, points(i,3));
        points_map[p] = i;
    }

    auto points_range = points_map | boost::adaptors::map_keys;
    Delaunay dt(std::begin(points_range), std::end(points_range));
    ASWrapper::fill_filtration(dt, points_map, add_simplex);
}

template<bool exact>
template<class Points, class SimplexCallback>
void
diode::AlphaShapes<exact>::
fill_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                           std::array<double, 3> from, std::array<double, 3> to)
{
    using K          = detail::Kernel<exact>;
    using PK         = CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>;

    using DsVb       = CGAL::Periodic_3_triangulation_ds_vertex_base_3<>;
    using Vb         = CGAL::Triangulation_vertex_base_3<PK,DsVb>;
    using AsVb       = CGAL::Alpha_shape_vertex_base_3<PK,Vb>;
    using DsCb       = CGAL::Periodic_3_triangulation_ds_cell_base_3<>;
    using Cb         = CGAL::Triangulation_cell_base_3<PK,DsCb>;
    using AsCb       = CGAL::Alpha_shape_cell_base_3<PK,Cb>;
    using TDS        = CGAL::Triangulation_data_structure_3<AsVb,AsCb>;
    using Delaunay   = CGAL::Periodic_3_Delaunay_triangulation_3<PK,TDS>;

    using ASWrapper  = detail::AlphaShapeWrapper<Delaunay>;
    using PointsMap  = typename ASWrapper::PointsMap;
    using AlphaShape = typename ASWrapper::AlphaShape;
    using Vertex     = typename ASWrapper::Vertex;
    using Point      = typename ASWrapper::Point;

    PointsMap points_map;
    for (Vertex i = 0; i < points.size(); ++i)
    {
        Point p(points(i,0), points(i,1), points(i,2));
        points_map[p] = i;
    }

    Delaunay pdt(typename PK::Iso_cuboid_3(from[0], from[1], from[2], to[0], to[1], to[2]));
    auto points_range = points_map | boost::adaptors::map_keys;
    pdt.insert(std::begin(points_range), std::end(points_range), true);
    if (pdt.is_triangulation_in_1_sheet())
        pdt.convert_to_1_sheeted_covering();
    else
        throw std::runtime_error("Cannot convert to 1-sheeted covering");
    ASWrapper::fill_filtration(pdt, points_map, add_simplex);
}

template<bool exact>
template<class Points>
std::array<typename Points::Real, 3>
diode::AlphaShapes<exact>::
circumcenter(const Points& points)
{
    using K          = detail::Kernel<exact>;
    using Point      = CGAL::Point_3<K>;
    using Real       = typename Points::Real;

    Point p_0(points(0,0), points(0,1), points(0,2));
    Point p_1(points(1,0), points(1,1), points(1,2));
    Point p_2(points(2,0), points(2,1), points(2,2));

    Point center;

    if (points.size() == 3)
        center = CGAL::circumcenter(p_0,p_1,p_2);
    else     // points.size() == 4
    {
        Point p_3(points(3,0), points(3,1), points(3,2));
        center = CGAL::circumcenter(p_0,p_1,p_2,p_3);
    }

    return std::array<Real, 3> { Real(CGAL::to_double(center[0])),
                                 Real(CGAL::to_double(center[1])),
                                 Real(CGAL::to_double(center[2])) };
}

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
template<bool exact>
template<class Points, class SimplexCallback>
void
diode::AlphaShapes<exact>::
fill_weighted_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                           std::array<double, 3> from, std::array<double, 3> to)
{
    using K          = detail::Kernel<exact>;
    using PK         = CGAL::Periodic_3_regular_triangulation_traits_3<K>;

    using DsVb       = CGAL::Periodic_3_triangulation_ds_vertex_base_3<>;
    using Vb         = CGAL::Regular_triangulation_vertex_base_3<PK,DsVb>;
    using AsVb       = CGAL::Alpha_shape_vertex_base_3<PK,Vb>;
    using DsCb       = CGAL::Periodic_3_triangulation_ds_cell_base_3<>;
    using Cb         = CGAL::Regular_triangulation_cell_base_3<PK,DsCb>;
    using AsCb       = CGAL::Alpha_shape_cell_base_3<PK,Cb>;
    using TDS        = CGAL::Triangulation_data_structure_3<AsVb,AsCb>;
    using Delaunay   = CGAL::Periodic_3_regular_triangulation_3<PK,TDS>;

    using ASWrapper  = detail::AlphaShapeWrapper<Delaunay>;
    using PointsMap  = typename ASWrapper::PointsMap;
    using AlphaShape = typename ASWrapper::AlphaShape;
    using Vertex     = typename ASWrapper::Vertex;
    using Point      = typename ASWrapper::Point;
    using FT         = typename ASWrapper::FT;

    auto domain_size = to[0] - from[0];
    auto upper_bound = FT(0.015625) * domain_size * domain_size;

    PointsMap points_map;
    for (Vertex i = 0; i < points.size(); ++i)
    {
        auto x = points(i,0);
        auto y = points(i,1);
        auto z = points(i,2);
        auto w = points(i,3);
        if (w < 0 || w >= upper_bound)
        {
            std::ostringstream oss;
            oss << "Point weight w must satisfy: 0 <= w < 1/64 * domain_size * domain_size; but got point"
                << " (" << x << ", " << y << ", " << z << ") weight = " << w;
            throw std::runtime_error(oss.str());
        }

        Point p({x,y,z}, w);
        points_map[p] = i;
    }

    Delaunay pdt(typename PK::Iso_cuboid_3(from[0], from[1], from[2], to[0], to[1], to[2]));
    auto points_range = points_map | boost::adaptors::map_keys;
    pdt.insert(std::begin(points_range), std::end(points_range), true);
    if (pdt.is_triangulation_in_1_sheet())
        pdt.convert_to_1_sheeted_covering();
    else
        throw std::runtime_error("Cannot convert to 1-sheeted covering");
    ASWrapper::fill_filtration(pdt, points_map, add_simplex);
}
#endif


template<class Points, class SimplexCallback>
void
diode::
fill_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex)
{
    using K             = CGAL::Exact_predicates_exact_constructions_kernel;
    using Delaunay2D    = CGAL::Delaunay_triangulation_2<K>;
    using Vertex_handle = Delaunay2D::Vertex_handle;
    using Point         = Delaunay2D::Point;
    using Face_handle   = Delaunay2D::Face_handle;

    using ASPointMap    = std::unordered_map<Vertex_handle, unsigned>;

    Delaunay2D  Dt;
    ASPointMap  point_map;
    for (unsigned i = 0; i < points.size(); ++i)
    {
        auto x = points(i,0);
        auto y = points(i,1);
        point_map[Dt.insert(Point(x,y))] = i;
    }

    // fill simplex set
    struct Simplex2D: public std::array<unsigned, 3>
    {
        double value;
        unsigned dimension;

        void sort() { std::sort(begin(), begin() + dimension + 1); }

        bool operator<(const Simplex2D& s) const
        {
            if (dimension < s.dimension) return true;
            if (dimension == s.dimension)
                return std::lexicographical_compare(begin(), begin() + dimension + 1, s.begin(), s.begin() + dimension + 1);

            return false;
        }


        bool operator==(const Simplex2D& s) const
        {
            if (dimension != s.dimension)
                return false;
            for (unsigned i = 0; i < dimension + 1; ++i)
                if ((*this)[i] != s[i])
                    return false;
            return true;
        }
    };

    auto simplex_from_face = [&](const Delaunay2D::Face& f)
    {
        Simplex2D s;

        s.dimension = 2;
        for (int i = 0; i < 3; ++i) s[i] = point_map[f.vertex(i)];
        auto p1    = f.vertex(0)->point();
        auto p2    = f.vertex(1)->point();
        auto p3    = f.vertex(2)->point();
        auto alpha = CGAL::squared_radius(p1, p2, p3);
        s.value = detail::to_floating_point(alpha);

        s.sort();

        return s;
    };

    std::set<Simplex2D> simplices;

    // faces
    for(auto cur = Dt.finite_faces_begin(); cur != Dt.finite_faces_end(); ++cur)
    {
        Simplex2D s = simplex_from_face(*cur);
        simplices.emplace(s);
    }

    // edges
    for(auto cur = Dt.finite_edges_begin(); cur != Dt.finite_edges_end(); ++cur)
    {
        auto e = *cur;
        Simplex2D s; s.dimension = 1;

        std::array<Point,2> points;
        unsigned j = 0;
        Face_handle f = e.first;
        for (int i = 0; i < 3; ++i)
            if (i != e.second)
            {
                points[j] = f->vertex(i)->point();
                s[j++] = point_map[f->vertex(i)];
            }
        auto& p1 = points[0];
        auto& p2 = points[1];

        Face_handle o = f->neighbor(e.second);
        if (o == Face_handle())
        {
            s.value = detail::to_floating_point(CGAL::squared_radius(p1, p2));
        } else
        {
            int oi = o->index(f);

            bool attached = false;
            if (!Dt.is_infinite(f->vertex(e.second)) &&
                CGAL::side_of_bounded_circle(p1, p2,
                                             f->vertex(e.second)->point()) == CGAL::ON_BOUNDED_SIDE)
                attached = true;
            else if (!Dt.is_infinite(o->vertex(oi)) &&
                     CGAL::side_of_bounded_circle(p1, p2,
                                                  o->vertex(oi)->point()) == CGAL::ON_BOUNDED_SIDE)
                attached = true;
            else
                s.value = detail::to_floating_point(CGAL::squared_radius(p1, p2));

            if (attached)
            {
                if (Dt.is_infinite(f))
                    s.value = simplices.find(simplex_from_face(*o))->value;
                else if (Dt.is_infinite(o))
                    s.value = simplices.find(simplex_from_face(*f))->value;
                else
                    s.value = std::min(simplices.find(simplex_from_face(*f))->value,
                                       simplices.find(simplex_from_face(*o))->value);
            }
        }

        s.sort();
        simplices.emplace(s);
    }

    // vertices
    for(auto cur = Dt.finite_vertices_begin(); cur != Dt.finite_vertices_end(); ++cur)
    {
        Simplex2D s;

        s.dimension = 0;
        s.value = 0;
        for (int i = 0; i < 3; ++i)
            if (cur->face()->vertex(i) != Vertex_handle() && cur->face()->vertex(i)->point() == cur->point())
                s[0] = point_map[cur->face()->vertex(i)];

        simplices.emplace(s);
    }

    // invoke callback with the simplices
    for (auto& s : simplices)
    {
        if (s.dimension == 0)
        {
            std::array<unsigned, 1> vertices { s[0] };
            add_simplex(vertices, s.value);
        } else if (s.dimension == 1)
        {
            std::array<unsigned, 2> vertices { s[0], s[1] };
            add_simplex(vertices, s.value);
        } else if (s.dimension == 2)
        {
            std::array<unsigned, 3> vertices { s[0], s[1], s[2] };
            add_simplex(vertices, s.value);
        }
    }
}


// Variant that emits attacher info for each simplex.
// Callback signature: add_simplex(sigma_verts, alpha, tau_verts).
// For Gabriel sigma, tau == sigma. For non-Gabriel sigma, tau is a Gabriel
// coface whose own squared circumradius equals alpha.
template<class Points, class SimplexCallback>
void
diode::
fill_alpha_shapes2d_with_attachment(const Points& points, const SimplexCallback& add_simplex)
{
    using K             = CGAL::Exact_predicates_exact_constructions_kernel;
    using Delaunay2D    = CGAL::Delaunay_triangulation_2<K>;
    using Vertex_handle = Delaunay2D::Vertex_handle;
    using Point         = Delaunay2D::Point;
    using Face_handle   = Delaunay2D::Face_handle;

    using ASPointMap    = std::unordered_map<Vertex_handle, unsigned>;

    Delaunay2D  Dt;
    ASPointMap  point_map;
    for (unsigned i = 0; i < points.size(); ++i)
    {
        auto x = points(i,0);
        auto y = points(i,1);
        point_map[Dt.insert(Point(x,y))] = i;
    }

    // Same Simplex2D as the no-attachment version, plus attacher fields:
    //   attacher_dim:     0 (vertex), 1 (edge), 2 (face)
    //   attacher_vertices: meaningful prefix of length attacher_dim + 1
    struct Simplex2D: public std::array<unsigned, 3>
    {
        double value;
        unsigned dimension;
        unsigned attacher_dim;
        std::array<unsigned, 3> attacher_vertices;

        void sort() { std::sort(begin(), begin() + dimension + 1); }

        bool operator<(const Simplex2D& s) const
        {
            if (dimension < s.dimension) return true;
            if (dimension == s.dimension)
                return std::lexicographical_compare(begin(), begin() + dimension + 1, s.begin(), s.begin() + dimension + 1);

            return false;
        }


        bool operator==(const Simplex2D& s) const
        {
            if (dimension != s.dimension)
                return false;
            for (unsigned i = 0; i < dimension + 1; ++i)
                if ((*this)[i] != s[i])
                    return false;
            return true;
        }
    };

    auto simplex_from_face = [&](const Delaunay2D::Face& f)
    {
        Simplex2D s;

        s.dimension = 2;
        for (int i = 0; i < 3; ++i) s[i] = point_map[f.vertex(i)];
        auto p1    = f.vertex(0)->point();
        auto p2    = f.vertex(1)->point();
        auto p3    = f.vertex(2)->point();
        auto alpha = CGAL::squared_radius(p1, p2, p3);
        s.value = detail::to_floating_point(alpha);

        // Faces are always self-attaching (their own squared circumradius
        // is what defines their alpha value).
        s.attacher_dim = 2;
        s.attacher_vertices = { s[0], s[1], s[2] };

        s.sort();

        return s;
    };

    std::set<Simplex2D> simplices;

    // faces
    for(auto cur = Dt.finite_faces_begin(); cur != Dt.finite_faces_end(); ++cur)
    {
        Simplex2D s = simplex_from_face(*cur);
        simplices.emplace(s);
    }

    // edges
    for(auto cur = Dt.finite_edges_begin(); cur != Dt.finite_edges_end(); ++cur)
    {
        auto e = *cur;
        Simplex2D s; s.dimension = 1;

        std::array<Point,2> points;
        unsigned j = 0;
        Face_handle f = e.first;
        for (int i = 0; i < 3; ++i)
            if (i != e.second)
            {
                points[j] = f->vertex(i)->point();
                s[j++] = point_map[f->vertex(i)];
            }
        auto& p1 = points[0];
        auto& p2 = points[1];

        // Default: edge is its own attacher (Gabriel case).
        s.attacher_dim = 1;
        s.attacher_vertices = { s[0], s[1], 0 };

        Face_handle o = f->neighbor(e.second);
        if (o == Face_handle())
        {
            s.value = detail::to_floating_point(CGAL::squared_radius(p1, p2));
        } else
        {
            int oi = o->index(f);

            bool attached = false;
            if (!Dt.is_infinite(f->vertex(e.second)) &&
                CGAL::side_of_bounded_circle(p1, p2,
                                             f->vertex(e.second)->point()) == CGAL::ON_BOUNDED_SIDE)
                attached = true;
            else if (!Dt.is_infinite(o->vertex(oi)) &&
                     CGAL::side_of_bounded_circle(p1, p2,
                                                  o->vertex(oi)->point()) == CGAL::ON_BOUNDED_SIDE)
                attached = true;
            else
                s.value = detail::to_floating_point(CGAL::squared_radius(p1, p2));

            if (attached)
            {
                Simplex2D winning_face;
                if (Dt.is_infinite(f))
                    winning_face = *simplices.find(simplex_from_face(*o));
                else if (Dt.is_infinite(o))
                    winning_face = *simplices.find(simplex_from_face(*f));
                else
                {
                    Simplex2D sf = *simplices.find(simplex_from_face(*f));
                    Simplex2D so = *simplices.find(simplex_from_face(*o));
                    winning_face = (sf.value <= so.value) ? sf : so;
                }
                s.value = winning_face.value;
                s.attacher_dim = 2;
                s.attacher_vertices = winning_face.attacher_vertices;
            }
        }

        s.sort();
        simplices.emplace(s);
    }

    // vertices: always Gabriel for unweighted alpha shapes (alpha = 0); tau = vertex itself.
    for(auto cur = Dt.finite_vertices_begin(); cur != Dt.finite_vertices_end(); ++cur)
    {
        Simplex2D s;

        s.dimension = 0;
        s.value = 0;
        for (int i = 0; i < 3; ++i)
            if (cur->face()->vertex(i) != Vertex_handle() && cur->face()->vertex(i)->point() == cur->point())
                s[0] = point_map[cur->face()->vertex(i)];

        s.attacher_dim = 0;
        s.attacher_vertices = { s[0], 0, 0 };

        simplices.emplace(s);
    }

    // invoke callback with the simplices
    for (auto& s : simplices)
    {
        std::vector<unsigned> tau_verts(s.attacher_vertices.begin(),
                                        s.attacher_vertices.begin() + s.attacher_dim + 1);
        if (s.dimension == 0)
        {
            std::array<unsigned, 1> vertices { s[0] };
            add_simplex(vertices, s.value, tau_verts);
        } else if (s.dimension == 1)
        {
            std::array<unsigned, 2> vertices { s[0], s[1] };
            add_simplex(vertices, s.value, tau_verts);
        } else if (s.dimension == 2)
        {
            std::array<unsigned, 3> vertices { s[0], s[1], s[2] };
            add_simplex(vertices, s.value, tau_verts);
        }
    }
}


template<class Points, class SimplexCallback>
void
diode::
fill_periodic_alpha_shapes2d(const Points& points, const SimplexCallback& add_simplex,std::array<double, 2> from, std::array<double, 2> to)
{
    using K             = CGAL::Exact_predicates_exact_constructions_kernel;
    using GT            = CGAL::Periodic_2_Delaunay_triangulation_traits_2<K>;
    using PDelaunay2D   = CGAL::Periodic_2_Delaunay_triangulation_2<GT>;
    using Vertex_handle = PDelaunay2D::Vertex_handle;
    using Point         = PDelaunay2D::Point;
    using Face_handle   = PDelaunay2D::Face_handle;
    using Iso_rectangle = PDelaunay2D::Iso_rectangle;

    using ASPointMap    = std::unordered_map<Vertex_handle, unsigned>;

    Iso_rectangle domain(from[0],from[1],to[0],to[1]);
    PDelaunay2D  pdt(domain);


    ASPointMap  point_map;

    for (unsigned i = 0; i < points.size(); ++i)
    {
        auto x = points(i,0);
        auto y = points(i,1);
        point_map[pdt.insert(Point(x,y))] = i;
    }
    
    if (pdt.is_triangulation_in_1_sheet())
        pdt.convert_to_1_sheeted_covering();
    else
        throw std::runtime_error("Cannot convert to 1-sheeted covering");

    // fill simplex set
    struct Simplex2D: public std::array<unsigned, 3>
    {
        double value;
        unsigned dimension;

        void sort() { std::sort(begin(), begin() + dimension + 1); }

        bool operator<(const Simplex2D& s) const
        {
            if (dimension < s.dimension) return true;
            if (dimension == s.dimension)
                return std::lexicographical_compare(begin(), begin() + dimension + 1, s.begin(), s.begin() + dimension + 1);

            return false;
        }


        bool operator==(const Simplex2D& s) const
        {
            if (dimension != s.dimension)
                return false;
            for (unsigned i = 0; i < dimension + 1; ++i)
                if ((*this)[i] != s[i])
                    return false;
            return true;
        }
    };

     auto simplex_from_face = [&](const PDelaunay2D::Face_handle f)
     {
        Simplex2D s;

        s.dimension = 2;
        for (int i = 0; i < 3; ++i) s[i] = point_map[f->vertex(i)];
        
        auto T =pdt.triangle(pdt.periodic_triangle(f));

        auto p1    = T.vertex(0);
        auto p2    = T.vertex(1);
        auto p3    = T.vertex(2);

        auto alpha = CGAL::squared_radius(p1, p2, p3);
        s.value = detail::to_floating_point(alpha);

        s.sort();

        return s;
    };
    
    std::set<Simplex2D> simplices;


    // faces
    for(auto cur = pdt.finite_faces_begin(); cur != pdt.finite_faces_end(); ++cur)
    {
        Simplex2D s = simplex_from_face(cur);
        simplices.emplace(s);
    }

    // edges
    for(auto cur = pdt.finite_edges_begin(); cur != pdt.finite_edges_end(); ++cur)
    {
        auto e = *cur;
        Simplex2D s; s.dimension = 1;

        std::array<Point,2> points;
        unsigned j = 0;
        Face_handle f = e.first;
        auto t = pdt.triangle(f);
        for (int i = 0; i < 3; ++i)
            if (i != e.second)
            {
                points[j] = t.vertex(i);
                s[j++] = point_map[f->vertex(i)];
            }
        auto& p1 = points[0];
        auto& p2 = points[1];

        Face_handle o = f->neighbor(e.second);
        if (o == Face_handle())
        {
            s.value = detail::to_floating_point(CGAL::squared_radius(p1, p2));
        } else
        {
            int oi = o->index(f);

            bool attached = false;
            if (!pdt.is_infinite(f->vertex(e.second)) &&
                CGAL::side_of_bounded_circle(p1, p2,
                                             f->vertex(e.second)->point()) == CGAL::ON_BOUNDED_SIDE)
                attached = true;
            else if (!pdt.is_infinite(o->vertex(oi)) &&
                     CGAL::side_of_bounded_circle(p1, p2,
                                                  o->vertex(oi)->point()) == CGAL::ON_BOUNDED_SIDE)
                attached = true;
            else
                s.value = detail::to_floating_point(CGAL::squared_radius(p1, p2));

            if (attached)
            {
                if (pdt.is_infinite(f)){

                    s.value = simplices.find(simplex_from_face(o))->value;

                    }
                else if (pdt.is_infinite(o)){
                    s.value = simplices.find(simplex_from_face(f))->value;

                }
                else{
                    s.value = std::min(simplices.find(simplex_from_face(f))->value,
                                       simplices.find(simplex_from_face(o))->value);
                }
                    
            }
        }

        s.sort();
        simplices.emplace(s);
    }

    // vertices
    for(auto cur = pdt.finite_vertices_begin(); cur != pdt.finite_vertices_end(); ++cur)
    {
        Simplex2D s;

        s.dimension = 0;
        s.value = 0;
        for (int i = 0; i < 3; ++i)
            if (cur->face()->vertex(i) != Vertex_handle() && cur->face()->vertex(i)->point() == cur->point())
                s[0] = point_map[cur->face()->vertex(i)];

        simplices.emplace(s);
    }

    // invoke callback with the simplices
    for (auto& s : simplices)
    {
        if (s.dimension == 0)
        {
            std::array<unsigned, 1> vertices { s[0] };
            add_simplex(vertices, s.value);
        } else if (s.dimension == 1)
        {
            std::array<unsigned, 2> vertices { s[0], s[1] };
            add_simplex(vertices, s.value);
        } else if (s.dimension == 2)
        {
            std::array<unsigned, 3> vertices { s[0], s[1], s[2] };
            add_simplex(vertices, s.value);
        }
    }
}
