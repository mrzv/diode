namespace diode
{
namespace detail
{

template<class Delaunay_, class Vertex_ = unsigned>
struct AlphaShapeWrapper
{
    using Vertex     = Vertex_;
    using Delaunay   = Delaunay_;
    using AlphaShape = CGAL::Alpha_shape_3<Delaunay>;
    using Point      = typename Delaunay::Point;
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
    static void fill_filtration(const PointsMap& points_map, const AddSimplex& add_simplex)
    {
        auto points_range = points_map | boost::adaptors::map_keys;
        AlphaShape as(std::begin(points_range), std::end(points_range), std::numeric_limits<FT>::infinity(), AlphaShape::GENERAL);
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

    ASWrapper::fill_filtration(points_map, add_simplex);
}
