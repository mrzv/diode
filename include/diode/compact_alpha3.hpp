// Included inside diode::detail by diode.hpp.

inline bool compact_alpha3_full_dimension(const std::vector<Vertex<3>>& vertices) {
    if(vertices.size()<4) return false;
    const auto& a=vertices[0].point;
    const auto& b=vertices[1].point;
    std::size_t third=2;
    for(;third<vertices.size();++third) {
        const auto& c=vertices[third].point;
        bool collinear=true;
        for(unsigned omitted=0;omitted<3;++omitted) {
            const unsigned x=(omitted+1)%3, y=(omitted+2)%3;
            const double pa[2]={a[x],a[y]}, pb[2]={b[x],b[y]}, pc[2]={c[x],c[y]};
            if(GEO::PCK::orient_2d(pa,pb,pc)!=GEO::ZERO) {
                collinear=false;
                break;
            }
        }
        if(!collinear) break;
    }
    if(third==vertices.size()) return false;
    for(std::size_t i=third+1;i<vertices.size();++i)
        if(GEO::PCK::orient_3d(a.data(),b.data(),vertices[third].point.data(),
                             vertices[i].point.data())!=GEO::ZERO) return true;
    return false;
}

inline std::uint64_t compact_alpha3_edge_key(unsigned a,unsigned b) {
    if(a>b) std::swap(a,b);
    return (static_cast<std::uint64_t>(a)<<32)|b;
}

// Structure-of-arrays open addressing avoids a node and allocator call per edge.
class CompactAlpha3Edges {
public:
    explicit CompactAlpha3Edges(std::size_t expected) {
        std::size_t capacity=16;
        while(capacity-capacity/4<expected) {
            if(capacity>std::numeric_limits<std::size_t>::max()/2)
                throw std::runtime_error("alpha edge table capacity overflow");
            capacity*=2;
        }
        reset(capacity);
    }

    void insert(std::uint64_t key,double alpha,bool witness) {
        if(size_+1>keys_.size()-keys_.size()/4) {
            if(keys_.size()>std::numeric_limits<std::size_t>::max()/2)
                throw std::runtime_error("alpha edge table capacity overflow");
            grow(keys_.size()*2);
        }
        insert_without_growth(key,alpha,witness);
    }

    std::size_t size() const { return size_; }

    template<class Callback> void for_each(const Callback& callback) const {
        for(std::size_t i=0;i<keys_.size();++i)
            if(keys_[i]!=empty_key)
                callback(static_cast<unsigned>(keys_[i]>>32),
                         static_cast<unsigned>(keys_[i]),alpha_[i],gabriel_[i]!=0);
    }

private:
    static constexpr std::uint64_t empty_key=std::numeric_limits<std::uint64_t>::max();

    void reset(std::size_t capacity) {
        keys_.assign(capacity,empty_key);
        alpha_.resize(capacity);
        gabriel_.resize(capacity);
        size_=0;
    }

    void insert_without_growth(std::uint64_t key,double alpha,bool witness) {
        std::uint64_t hash=key;
        hash=(hash^(hash>>30))*UINT64_C(0xbf58476d1ce4e5b9);
        hash=(hash^(hash>>27))*UINT64_C(0x94d049bb133111eb);
        hash^=hash>>31;
        const std::size_t mask=keys_.size()-1;
        std::size_t slot=static_cast<std::size_t>(hash)&mask;
        while(keys_[slot]!=empty_key && keys_[slot]!=key) slot=(slot+1)&mask;
        if(keys_[slot]==empty_key) {
            keys_[slot]=key;
            alpha_[slot]=alpha;
            gabriel_[slot]=!witness;
            ++size_;
        } else {
            alpha_[slot]=std::min(alpha_[slot],alpha);
            if(witness) gabriel_[slot]=0;
        }
    }

    void grow(std::size_t capacity) {
        auto keys=std::move(keys_);
        auto alpha=std::move(alpha_);
        auto gabriel=std::move(gabriel_);
        reset(capacity);
        for(std::size_t i=0;i<keys.size();++i)
            if(keys[i]!=empty_key) insert_without_growth(keys[i],alpha[i],!gabriel[i]);
    }

    std::vector<std::uint64_t> keys_;
    std::vector<double> alpha_;
    std::vector<unsigned char> gabriel_;
    std::size_t size_=0;
};

struct CompactAlpha3Sphere {
    Point<3> center{};
    double alpha=0.0;
};

inline double compact_alpha3_power(const Point<3>& center,const Vertex<3>& vertex) {
    const double x=center[0]-vertex.point[0];
    const double y=center[1]-vertex.point[1];
    const double z=center[2]-vertex.point[2];
    return x*x+y*y+z*z-vertex.weight;
}

template<std::size_t N>
inline double compact_alpha3_sphere(
    std::array<const Vertex<3>*,N> vertices,Point<3>* center=nullptr
) {
    // Match generic records' canonical order, including their translation
    // origin. Ill-conditioned simplices must not depend on cell-local order.
    std::sort(vertices.begin(),vertices.end(),[](const auto* a,const auto* b) {
        return a->id<b->id;
    });
    std::array<Point<3>,N> points;
    std::array<double,N> weights;
    for(std::size_t i=0;i<N;++i) {
        points[i]=vertices[i]->point;
        weights[i]=vertices[i]->weight;
    }
    return sphere(points,weights,center);
}

template<std::size_t N,class CB>
void compact_alpha3_emit(const std::vector<Vertex<3>>& vertices,
                         std::array<unsigned,N> indices,double alpha,const CB& cb) {
    if(!std::isfinite(alpha)) throw std::runtime_error("non-finite alpha radius");
    for(auto& index:indices) index=vertices[index].id;
    std::sort(indices.begin(),indices.end());
    cb(indices,alpha);
}

template<class Points,class CB>
bool try_compact_alpha3(const Points& points,const CB& cb,bool weighted) {
    if(points.size()>std::numeric_limits<unsigned>::max())
        throw std::runtime_error("too many alpha points");
    initialize_geogram();
    const auto vertices=unique_points<Points,3>(points,weighted);
    // False is reserved for geometric rank deficiency, before any callback.
    if(!compact_alpha3_full_dimension(vertices)) return false;
    const std::size_t count=vertices.size();
    if(count>=std::numeric_limits<GEO::index_t>::max())
        throw std::runtime_error("too many Geogram alpha vertices");
    const std::size_t stride=weighted?4:3;
    std::vector<double> coordinates;
    coordinates.reserve(stride*count);
    double max_weight=vertices.front().weight;
    if(weighted) for(const auto& vertex:vertices) max_weight=std::max(max_weight,vertex.weight);
    for(const auto& vertex:vertices) {
        for(double value:vertex.point) coordinates.push_back(value);
        if(weighted) {
            // Subtract in wider precision so finite, opposite-sign weights need
            // not overflow before taking a representable square root.
            const long double difference=static_cast<long double>(max_weight)-vertex.weight;
            coordinates.push_back(static_cast<double>(std::sqrt(difference)));
        }
    }
    GEO::Delaunay_var triangulation;
    {
        std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
        GEO::Numeric::random_reset();
        triangulation=weighted
            ? GEO::Delaunay::create(4,"BPOW")
            : static_cast<GEO::Delaunay*>(new GEO::Delaunay3d(3));
        if(triangulation.is_null()) throw std::runtime_error("Geogram alpha backend unavailable");
        triangulation->set_reorder(true);
        triangulation->set_vertices(static_cast<GEO::index_t>(count),coordinates.data());
    }
    const GEO::index_t cells=triangulation->nb_cells();
    if(triangulation->cell_size()!=4 || triangulation->nb_vertices()!=count || cells==0)
        throw std::runtime_error("Geogram did not produce a full 3D alpha triangulation");
    std::vector<double> cell_alpha(cells);
    std::vector<unsigned char> visible(count,0);
    std::size_t facets=0;
    for(GEO::index_t cell=0;cell<cells;++cell) {
        std::array<unsigned,4> ids{};
        for(unsigned local=0;local<4;++local) {
            const GEO::index_t vertex=triangulation->cell_vertex(cell,local);
            if(vertex>=count) throw std::runtime_error("invalid Geogram alpha vertex");
            ids[local]=static_cast<unsigned>(vertex);
            visible[vertex]=1;
            const GEO::index_t adjacent=triangulation->cell_adjacent(cell,local);
            if(adjacent!=GEO::NO_INDEX && adjacent>=cells)
                throw std::runtime_error("invalid Geogram alpha adjacency");
            if(adjacent==GEO::NO_INDEX || adjacent>cell) ++facets;
        }
        cell_alpha[cell]=compact_alpha3_sphere<4>({
            &vertices[ids[0]],&vertices[ids[1]],&vertices[ids[2]],&vertices[ids[3]]});
        if(!std::isfinite(cell_alpha[cell])) throw std::runtime_error("non-finite alpha cell radius");
    }
    const std::size_t visible_count=static_cast<std::size_t>(
        std::count(visible.begin(),visible.end(),static_cast<unsigned char>(1)));
    if(!weighted && visible_count!=count)
        throw std::runtime_error("Geogram omitted an ordinary alpha vertex");
    if(facets>std::numeric_limits<std::size_t>::max()-visible_count || visible_count+facets<=cells)
        throw std::runtime_error("invalid alpha edge count");
    const std::size_t expected_edges=visible_count+facets-cells-1;
    CompactAlpha3Edges edges(expected_edges);
    for(GEO::index_t cell=0;cell<cells;++cell) {
        std::array<unsigned,4> ids{};
        for(unsigned local=0;local<4;++local)
            ids[local]=static_cast<unsigned>(triangulation->cell_vertex(cell,local));
        compact_alpha3_emit(vertices,ids,cell_alpha[cell],cb);
        for(unsigned opposite=0;opposite<4;++opposite) {
            const GEO::index_t adjacent=triangulation->cell_adjacent(cell,opposite);
            if(adjacent!=GEO::NO_INDEX && adjacent<cell) continue;
            std::array<unsigned,3> face{};
            unsigned next=0;
            for(unsigned local=0;local<4;++local) if(local!=opposite) face[next++]=ids[local];
            CompactAlpha3Sphere face_sphere;
            if(weighted) face_sphere.alpha=compact_alpha3_sphere<3>({
                &vertices[face[0]],&vertices[face[1]],&vertices[face[2]]},&face_sphere.center);
            const auto interior=[&](unsigned witness) {
                if(weighted) return compact_alpha3_power(face_sphere.center,vertices[witness])<face_sphere.alpha;
                const double* a=vertices[face[0]].point.data();
                const double* b=vertices[face[1]].point.data();
                const double* c=vertices[face[2]].point.data();
                return GEO::PCK::side3_SOS(a,b,c,vertices[witness].point.data(),a,b,c,3)==GEO::NEGATIVE;
            };
            bool gabriel=!interior(ids[opposite]);
            double alpha=cell_alpha[cell];
            if(adjacent!=GEO::NO_INDEX) {
                unsigned local=0;
                while(local<4 && triangulation->cell_adjacent(adjacent,local)!=cell) ++local;
                if(local==4) throw std::runtime_error("inconsistent Geogram alpha adjacency");
                gabriel=gabriel && !interior(static_cast<unsigned>(triangulation->cell_vertex(adjacent,local)));
                alpha=std::min(alpha,cell_alpha[adjacent]);
            }
            // Clamp own radii to cofaces to preserve filtration monotonicity
            // even when independent floating-point constructions differ by ulps.
            if(gabriel) alpha=std::min(alpha,weighted?face_sphere.alpha:
                compact_alpha3_sphere<3>({&vertices[face[0]],&vertices[face[1]],&vertices[face[2]]}));
            compact_alpha3_emit(vertices,face,alpha,cb);
            for(unsigned omit=0;omit<3;++omit) {
                const unsigned a=face[(omit+1)%3],b=face[(omit+2)%3],witness=face[omit];
                bool inside;
                if(weighted) {
                    CompactAlpha3Sphere edge;
                    edge.alpha=compact_alpha3_sphere<2>({&vertices[a],&vertices[b]},&edge.center);
                    inside=compact_alpha3_power(edge.center,vertices[witness])<edge.alpha;
                } else {
                    const double* pa=vertices[a].point.data();
                    const double* pb=vertices[b].point.data();
                    inside=GEO::PCK::side2_SOS(pa,pb,vertices[witness].point.data(),pa,pb,3)==GEO::NEGATIVE;
                }
                edges.insert(compact_alpha3_edge_key(a,b),alpha,inside);
            }
        }
    }
    if(edges.size()!=expected_edges) throw std::runtime_error("alpha Delaunay Euler characteristic check failed");
    std::vector<double> vertex_alpha(weighted?count:0,std::numeric_limits<double>::infinity());
    std::vector<unsigned char> vertex_gabriel(weighted?count:0,1);
    edges.for_each([&](unsigned a,unsigned b,double alpha,bool gabriel) {
        if(gabriel) {
            const double own=compact_alpha3_sphere<2>({&vertices[a],&vertices[b]});
            if(!std::isfinite(own)) throw std::runtime_error("non-finite alpha edge radius");
            alpha=std::min(alpha,own);
        }
        compact_alpha3_emit(vertices,std::array<unsigned,2>{a,b},alpha,cb);
        if(weighted) {
            vertex_alpha[a]=std::min(vertex_alpha[a],alpha);
            vertex_alpha[b]=std::min(vertex_alpha[b],alpha);
            if(compact_alpha3_power(vertices[a].point,vertices[b])<-vertices[a].weight) vertex_gabriel[a]=0;
            if(compact_alpha3_power(vertices[b].point,vertices[a])<-vertices[b].weight) vertex_gabriel[b]=0;
        }
    });
    for(unsigned i=0;i<count;++i) if(visible[i]) {
        double alpha=0.0;
        if(weighted) alpha=vertex_gabriel[i]?std::min(-vertices[i].weight,vertex_alpha[i]):vertex_alpha[i];
        compact_alpha3_emit(vertices,std::array<unsigned,1>{i},alpha,cb);
    }
    return true;
}
