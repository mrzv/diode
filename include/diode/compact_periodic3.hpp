// Included inside diode::detail, after the generic geometry helpers.

[[noreturn]] inline void periodic3_covering_error() {
    throw std::runtime_error("Cannot convert periodic triangulation to a one-sheeted covering");
}

inline Point<3> periodic3_extent(const Point<3>& from,const Point<3>& to) {
    Point<3> period{};
    for(std::size_t d=0;d<3;++d) {
        period[d]=to[d]-from[d];
        if(!std::isfinite(from[d]) || !std::isfinite(to[d]) || !std::isfinite(period[d]))
            throw std::runtime_error("periodic domain bounds and extents must be finite");
        if(!(period[d]>0.0))
            throw std::runtime_error("periodic domain is empty or inverted: require from[k] < to[k] on every axis");
    }
    return period;
}

struct Periodic3Alpha {
    Point<3> center{};
    double own=std::numeric_limits<double>::infinity();
    double alpha=std::numeric_limits<double>::infinity();
    bool gabriel=true;
};
struct Periodic3NoAlpha {};
template<bool Alpha> using Periodic3Value=std::conditional_t<Alpha,Periodic3Alpha,Periodic3NoAlpha>;
template<bool Alpha> struct Periodic3Edge: Periodic3Value<Alpha> {
    std::array<int,3> offset{};
};
template<bool Alpha> struct Periodic3Facet: Periodic3Value<Alpha> {
    unsigned cofaces=0;
};
struct Periodic3CellAlpha { double alpha=std::numeric_limits<double>::infinity(); };

template<bool Alpha> struct Periodic3Complex {
    std::vector<Vertex<3>> vertices;
    Point<3> period{};
    std::vector<unsigned char> visible;
    std::unordered_map<std::array<unsigned,2>,Periodic3Edge<Alpha>,VertexKeyHash<2>> edges;
    std::unordered_map<std::array<unsigned,3>,Periodic3Facet<Alpha>,VertexKeyHash<3>> facets;
    std::unordered_map<std::array<unsigned,4>,std::conditional_t<Alpha,Periodic3CellAlpha,Periodic3NoAlpha>,VertexKeyHash<4>> cells;

    template<std::size_t N>
    std::array<Vertex<3>,N> simplex(const std::array<unsigned,N>& key) const {
        std::array<Vertex<3>,N> result{};
        for(std::size_t i=0;i<N;++i) {
            result[i]=vertices[key[i]];
            if(i!=0) result[i].offset=edges.at({key[0],key[i]}).offset;
            else result[i].offset={};
            for(std::size_t d=0;d<3;++d)
                result[i].point[d]+=result[i].offset[d]*period[d];
        }
        return result;
    }

    template<std::size_t N>
    std::array<std::array<int,3>,N> offsets(const std::array<unsigned,N>& key) const {
        std::array<std::array<int,3>,N> result{};
        for(std::size_t i=1;i<N;++i) {
            const auto& native=edges.at({key[0],key[i]}).offset;
            for(std::size_t d=0;d<3;++d) {
                const long long value=static_cast<long long>(native[d])+
                    vertices[key[i]].offset[d]-vertices[key[0]].offset[d];
                if(value<std::numeric_limits<int>::min() || value>std::numeric_limits<int>::max())
                    periodic3_covering_error();
                result[i][d]=static_cast<int>(value);
            }
        }
        return result;
    }
};

template<std::size_t N>
std::array<unsigned,N-1> periodic3_face(const std::array<unsigned,N>& key,std::size_t omit) {
    std::array<unsigned,N-1> face{};
    std::size_t next=0;
    for(std::size_t i=0;i<N;++i) if(i!=omit) face[next++]=key[i];
    return face;
}

template<std::size_t N>
void periodic3_sphere(Periodic3Alpha& value,const std::array<Vertex<3>,N>& simplex) {
    std::array<Point<3>,N> points{};
    std::array<double,N> weights{};
    for(std::size_t i=0;i<N;++i) { points[i]=simplex[i].point; weights[i]=simplex[i].weight; }
    value.own=sphere(points,weights,&value.center);
    if(!std::isfinite(value.own)) throw std::runtime_error("non-finite alpha radius");
}

inline void periodic3_witness(Periodic3Alpha& value,const Vertex<3>& witness,
                              const Point<3>& shift) {
    if(!value.gabriel) return;
    long double power=-static_cast<long double>(witness.weight);
    for(std::size_t d=0;d<3;++d) {
        const long double delta=static_cast<long double>(value.center[d])-witness.point[d]-shift[d];
        power+=delta*delta;
    }
    const long double tolerance=128.0L*std::numeric_limits<double>::epsilon()*
        std::max(std::abs(static_cast<long double>(value.own)),std::abs(power));
    if(power<value.own-tolerance) value.gabriel=false;
}

template<bool Alpha>
Periodic3Complex<Alpha> acquire_periodic3(std::vector<Vertex<3>> vertices,bool weighted,
                                         const Point<3>& period) {
    Periodic3Complex<Alpha> out;
    out.period=period;
    // Preserve the established empty result for fewer than four unique sites.
    if(vertices.size()<4) return out;
    bool wrapped=false;
    for(auto& vertex:vertices) for(std::size_t d=0;d<3;++d) {
        if(vertex.point[d]>=0.0 && vertex.point[d]<period[d]) continue;
        wrapped=true;
        const long double coordinate=vertex.point[d];
        const long double quotient=std::floor(coordinate/period[d]);
        if(!std::isfinite(coordinate) || -quotient<std::numeric_limits<int>::min() ||
           -quotient>std::numeric_limits<int>::max()) periodic3_covering_error();
        vertex.offset[d]=static_cast<int>(-quotient);
        vertex.point[d]=static_cast<double>(coordinate-quotient*period[d]);
        // Keep a rounded upper endpoint inside the half-open native domain
        // without identifying a nearby negative site with the exact origin.
        if(vertex.point[d]>=period[d]) vertex.point[d]=std::nextafter(period[d],0.0);
    }
    if(wrapped) {
        // Period-equivalent sites must obey the same strongest-weight/latest-ID
        // selection rule as coincident input coordinates.
        std::map<Point<3>,Vertex<3>> selected;
        for(const auto& vertex:vertices) {
            auto inserted=selected.emplace(vertex.point,vertex);
            if(!inserted.second) {
                auto& previous=inserted.first->second;
                if(vertex.weight>previous.weight ||
                   (vertex.weight==previous.weight && vertex.id>previous.id)) previous=vertex;
            }
        }
        vertices.clear();
        for(const auto& item:selected) vertices.push_back(item.second);
        if(vertices.size()<4) return out;
    }
    // Native initialization requires a tetrahedron of real sites. A flat cloud
    // cannot define the unique simplicial one-sheeted covering exported here.
    if(!compact_alpha3_full_dimension(vertices)) periodic3_covering_error();
    if(vertices.size()>std::numeric_limits<GEO::index_t>::max()/27)
        throw std::runtime_error("too many Geogram periodic vertices");
    // Internal keys have the same ordering as the public input IDs.
    std::sort(vertices.begin(),vertices.end(),[](const auto& a,const auto& b) { return a.id<b.id; });
    out.vertices=std::move(vertices);
    const std::size_t count=out.vertices.size();
    out.visible.assign(count,0);
    std::vector<double> coordinates;
    std::vector<double> weights;
    coordinates.reserve(3*count);
    if(weighted) weights.reserve(count);
    double max_weight=0.0;
    if(weighted) {
        max_weight=out.vertices.front().weight;
        for(const auto& vertex:out.vertices) max_weight=std::max(max_weight,vertex.weight);
    }
    for(const auto& vertex:out.vertices) {
        for(double coordinate:vertex.point) coordinates.push_back(coordinate);
        // A common weight offset cannot change the regular triangulation.
        // Remove it before Geogram forms rounded lifted coordinates.
        if(weighted) weights.push_back(vertex.weight-max_weight);
    }
    {
        std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
        GEO::Numeric::random_reset();
        GEO::PeriodicDelaunay3d triangulation(GEO::vec3(period[0],period[1],period[2]));
        triangulation.set_vertices(static_cast<GEO::index_t>(count),coordinates.data());
        if(weighted) triangulation.set_weights(weights.data());
        triangulation.compute();
        // An aborted computation is not a triangulation with hidden sites.
        if(triangulation.has_empty_cells() || triangulation.cell_size()!=4 || triangulation.nb_cells()==0)
            throw std::runtime_error("Geogram did not produce a complete periodic triangulation");
        out.cells.reserve(triangulation.nb_cells());
        out.edges.reserve(triangulation.nb_cells());
        for(GEO::index_t cell=0;cell<triangulation.nb_cells();++cell) {
            std::array<std::pair<unsigned,std::array<int,3>>,4> sites{};
            for(GEO::index_t local=0;local<4;++local) {
                const auto pv=triangulation.cell_vertex(cell,local);
                if(pv>=27*count) throw std::runtime_error("invalid Geogram periodic vertex");
                sites[local].first=triangulation.periodic_vertex_real(pv);
                auto& offset=sites[local].second;
                triangulation.periodic_vertex_get_T(pv,offset[0],offset[1],offset[2]);
            }
            std::sort(sites.begin(),sites.end(),[](const auto& a,const auto& b) { return a.first<b.first; });
            std::array<unsigned,4> key{};
            for(std::size_t i=0;i<4;++i) {
                key[i]=sites[i].first;
                if(i!=0 && key[i]==key[i-1]) periodic3_covering_error();
                out.visible[key[i]]=1;
                for(std::size_t j=0;j<i;++j) {
                    std::array<int,3> offset{};
                    for(std::size_t d=0;d<3;++d) offset[d]=sites[i].second[d]-sites[j].second[d];
                    auto inserted=out.edges.try_emplace({key[j],key[i]});
                    if(inserted.second) {
                        inserted.first->second.offset=offset;
                        // Validate input-relative lift range before any emitter.
                        (void)out.offsets(std::array<unsigned,2>{key[j],key[i]});
                    }
                    else if(inserted.first->second.offset!=offset) periodic3_covering_error();
                }
            }
            out.cells.try_emplace(key);
        }
    }
    // Edge lifts determine all simplex lifts. Validate the entire covering
    // before any callback, including all faces of compressed boundary copies.
    out.facets.reserve(2*out.cells.size());
    for(const auto& cell:out.cells) for(std::size_t omit=0;omit<4;++omit)
        ++out.facets[periodic3_face(cell.first,omit)].cofaces;
    for(const auto& face:out.facets) if(face.second.cofaces!=2) periodic3_covering_error();
    const auto visible=static_cast<std::size_t>(std::count(out.visible.begin(),out.visible.end(),1));
    if(!weighted && visible!=count) throw std::runtime_error("Geogram omitted an ordinary periodic vertex");
    if(visible+out.facets.size()!=out.edges.size()+out.cells.size()) periodic3_covering_error();
    return out;
}

template<class Points,bool Alpha>
Periodic3Complex<Alpha> make_periodic3(const Points& points,bool weighted,
                                      const Point<3>& from,const Point<3>& to) {
    initialize_geogram();
    const auto period=periodic3_extent(from,to);
    if(points.size()>std::numeric_limits<unsigned>::max()) throw std::runtime_error("too many periodic points");
    return acquire_periodic3<Alpha>(unique_points<Points,3>(points,weighted,&from),weighted,period);
}

inline std::vector<Periodic3Alpha> periodic3_alpha(Periodic3Complex<true>& out) {
    std::vector<Periodic3Alpha> vertices(out.vertices.size());
    for(std::size_t i=0;i<vertices.size();++i) if(out.visible[i])
        periodic3_sphere(vertices[i],std::array<Vertex<3>,1>{out.vertices[i]});
    for(auto& edge:out.edges) periodic3_sphere(edge.second,out.simplex(edge.first));
    for(auto& face:out.facets) periodic3_sphere(face.second,out.simplex(face.first));
    for(auto& cell:out.cells) {
        const auto simplex=out.simplex(cell.first);
        Periodic3Alpha value;
        periodic3_sphere(value,simplex);
        cell.second.alpha=value.own;
        for(std::size_t omit=0;omit<4;++omit) {
            const auto key=periodic3_face(cell.first,omit);
            auto& face=out.facets.at(key);
            Point<3> shift{};
            const std::size_t first=omit==0 ? 1 : 0;
            for(std::size_t d=0;d<3;++d) shift[d]=out.vertices[key[0]].point[d]-simplex[first].point[d];
            periodic3_witness(face,simplex[omit],shift);
            face.alpha=std::min(face.alpha,value.own);
        }
    }
    for(auto& face:out.facets) {
        if(face.second.gabriel) face.second.alpha=std::min(face.second.alpha,face.second.own);
        const auto simplex=out.simplex(face.first);
        for(std::size_t omit=0;omit<3;++omit) {
            const auto key=periodic3_face(face.first,omit);
            auto& edge=out.edges.at(key);
            Point<3> shift{};
            const std::size_t first=omit==0 ? 1 : 0;
            for(std::size_t d=0;d<3;++d) shift[d]=out.vertices[key[0]].point[d]-simplex[first].point[d];
            periodic3_witness(edge,simplex[omit],shift);
            edge.alpha=std::min(edge.alpha,face.second.alpha);
        }
    }
    for(auto& edge:out.edges) {
        if(edge.second.gabriel) edge.second.alpha=std::min(edge.second.alpha,edge.second.own);
        const auto simplex=out.simplex(edge.first);
        for(std::size_t i=0;i<2;++i) {
            auto& vertex=vertices[edge.first[i]];
            Point<3> shift{};
            for(std::size_t d=0;d<3;++d) shift[d]=out.vertices[edge.first[i]].point[d]-simplex[i].point[d];
            periodic3_witness(vertex,simplex[1-i],shift);
            vertex.alpha=std::min(vertex.alpha,edge.second.alpha);
        }
    }
    for(std::size_t i=0;i<vertices.size();++i) if(out.visible[i]) {
        if(vertices[i].gabriel) vertices[i].alpha=std::min(vertices[i].alpha,vertices[i].own);
        if(!std::isfinite(vertices[i].alpha)) throw std::runtime_error("non-finite alpha radius");
    }
    return vertices;
}

template<std::size_t N,bool Alpha>
std::array<unsigned,N> periodic3_ids(const Periodic3Complex<Alpha>& out,std::array<unsigned,N> key) {
    for(auto& index:key) index=out.vertices[index].id;
    return key;
}

template<class Points,class CB>
void fill_compact_periodic3_alpha(const Points& points,const CB& cb,bool weighted,
                                  const Point<3>& from,const Point<3>& to) {
    auto out=make_periodic3<Points,true>(points,weighted,from,to);
    const auto vertices=periodic3_alpha(out);
    for(std::size_t i=0;i<vertices.size();++i) if(out.visible[i])
        cb(std::array<unsigned,1>{out.vertices[i].id},vertices[i].alpha);
    for(const auto& edge:out.edges) cb(periodic3_ids(out,edge.first),edge.second.alpha);
    for(const auto& face:out.facets) cb(periodic3_ids(out,face.first),face.second.alpha);
    for(const auto& cell:out.cells) cb(periodic3_ids(out,cell.first),cell.second.alpha);
}

template<bool Lifts,class Points,class CB>
void fill_compact_periodic3_delaunay(const Points& points,const CB& cb,bool weighted,
                                     const Point<3>& from,const Point<3>& to) {
    const auto out=make_periodic3<Points,false>(points,weighted,from,to);
    const auto emit=[&](const auto& key) {
        const auto ids=periodic3_ids(out,key);
        if constexpr(Lifts) {
            cb(ids,out.offsets(key));
        } else cb(ids);
    };
    for(std::size_t i=0;i<out.vertices.size();++i) if(out.visible[i]) emit(std::array<unsigned,1>{static_cast<unsigned>(i)});
    for(const auto& edge:out.edges) emit(edge.first);
    for(const auto& face:out.facets) emit(face.first);
    for(const auto& cell:out.cells) emit(cell.first);
}

inline Complex<3> triangulate3_periodic(std::vector<Vertex<3>> vertices,bool weighted,
                                       const Point<3>& period,bool collect_alpha) {
    const auto native=acquire_periodic3<false>(std::move(vertices),weighted,period);
    Complex<3> out(collect_alpha);
    out.reserve(native.vertices.size(),native.cells.size());
    for(const auto& cell:native.cells) {
        auto simplex=native.simplex(cell.first);
        const auto offsets=native.offsets(cell.first);
        for(std::size_t i=0;i<4;++i) simplex[i].offset=offsets[i];
        add_tetrahedron(out,simplex);
    }
    return out;
}
