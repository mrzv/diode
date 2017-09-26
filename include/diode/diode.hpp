#include <sstream>

namespace diode
{
namespace detail
{

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
                add_simplex(vertices, a);
            }
            else if (const E* e = CGAL::object_cast<E>(&o))
            {
                C c = e->first;
                auto u = points.find(c->vertex(e->second)->point())->second;
                auto v = points.find(c->vertex(e->third)->point())->second;
                std::array<Vertex, 2> vertices { u, v };
                add_simplex(vertices, a);
            }
            else if (const F* f = CGAL::object_cast<F>(&o))
            {
                std::array<Vertex, 3> vertices;
                size_t j = 0;
                C c = f->first;
                for (size_t i = 0; i < 4; ++i)
                    if (i != f->second)
                        vertices[j++] = points.find(c->vertex(i)->point())->second;
                add_simplex(vertices, a);
            } else if (const C* c = CGAL::object_cast<C>(&o))
            {
                std::array<Vertex, 4> vertices;
                for (size_t i = 0; i < 4; ++i)
                    vertices[i] = points.find((*c)->vertex(i)->point())->second;
                add_simplex(vertices, a);
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

template<class Points, class SimplexCallback>
void
diode::fill_alpha_shapes(const Points& points, const SimplexCallback& add_simplex)
{
    using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
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

template<class Points, class SimplexCallback>
void
diode::fill_weighted_alpha_shapes(const Points& points, const SimplexCallback& add_simplex)
{
    using K         = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Rvb       = CGAL::Regular_triangulation_vertex_base_3<K>;
    using Vb        = CGAL::Alpha_shape_vertex_base_3<K, Rvb>;
    using Rcb       = CGAL::Regular_triangulation_cell_base_3<K>;
    using Cb        = CGAL::Alpha_shape_cell_base_3<K, Rcb>;
    using TDS       = CGAL::Triangulation_data_structure_3<Vb,Cb>;
    using Delaunay  = CGAL::Regular_triangulation_3<K,TDS>;

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

template<class Points, class SimplexCallback>
void
diode::fill_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                  std::array<double, 3> from, std::array<double, 3> to)
{
    using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
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

    Delaunay pdt(PK::Iso_cuboid_3(from[0], from[1], from[2], to[0], to[1], to[2]));
    auto points_range = points_map | boost::adaptors::map_keys;
    pdt.insert(std::begin(points_range), std::end(points_range), true);
    if (pdt.is_triangulation_in_1_sheet())
        pdt.convert_to_1_sheeted_covering();
    else
        throw std::runtime_error("Cannot convert to 1-sheeted covering");
    ASWrapper::fill_filtration(pdt, points_map, add_simplex);
}

#if (CGAL_VERSION_MAJOR == 4 && CGAL_VERSION_MINOR >= 11) || (CGAL_VERSION_MAJOR > 4)
template<class Points, class SimplexCallback>
void
diode::fill_weighted_periodic_alpha_shapes(const Points& points, const SimplexCallback& add_simplex,
                                           std::array<double, 3> from, std::array<double, 3> to)
{
    using K          = CGAL::Exact_predicates_inexact_constructions_kernel;
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

    Delaunay pdt(PK::Iso_cuboid_3(from[0], from[1], from[2], to[0], to[1], to[2]));
    auto points_range = points_map | boost::adaptors::map_keys;
    pdt.insert(std::begin(points_range), std::end(points_range), true);
    if (pdt.is_triangulation_in_1_sheet())
        pdt.convert_to_1_sheeted_covering();
    else
        throw std::runtime_error("Cannot convert to 1-sheeted covering");
    ASWrapper::fill_filtration(pdt, points_map, add_simplex);
}
#endif

