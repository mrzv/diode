#include <sstream>
#include <unordered_map>
#include <algorithm>
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

    template<class AddSimplex>
    static void fill_filtration(Delaunay& dt, const PointsMap& points_map, const AddSimplex& add_simplex)
    {
        AlphaShape as(dt, std::numeric_limits<FT>::infinity(), AlphaShape::GENERAL);
        as.filtration_with_alpha_values(ASOutputIterator3<AddSimplex>(add_simplex, points_map));
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
