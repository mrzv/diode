// Included inside diode::detail, after the generic geometry helpers.

struct compact_periodic_sphere {
    Point<3> center{};
    double alpha=std::numeric_limits<double>::infinity();
};

inline Point<3> compact_periodic_subtract(const Point<3>& a,const Point<3>& b) {
    return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};
}

inline double compact_periodic_dot(const Point<3>& a,const Point<3>& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline Point<3> compact_periodic_point(
    const GEO::PeriodicDelaunay3d& triangulation,GEO::index_t vertex
) {
    const GEO::vec3 p=triangulation.vertex(vertex);
    return {p.x,p.y,p.z};
}

inline compact_periodic_sphere compact_periodic_edge_sphere(
    const Point<3>& p0,const Point<3>& p1
) {
    const auto delta=compact_periodic_subtract(p1,p0);
    return {{p0[0]+0.5*delta[0],p0[1]+0.5*delta[1],p0[2]+0.5*delta[2]},
            0.25*compact_periodic_dot(delta,delta)};
}

inline compact_periodic_sphere compact_periodic_triangle_sphere(
    const Point<3>& p0,const Point<3>& p1,const Point<3>& p2
) {
    const auto d1=compact_periodic_subtract(p1,p0);
    const auto d2=compact_periodic_subtract(p2,p0);
    const double g11=compact_periodic_dot(d1,d1);
    const double g12=compact_periodic_dot(d1,d2);
    const double g22=compact_periodic_dot(d2,d2);
    const Point<3> normal{d1[1]*d2[2]-d1[2]*d2[1],
                          d1[2]*d2[0]-d1[0]*d2[2],
                          d1[0]*d2[1]-d1[1]*d2[0]};
    const double determinant=compact_periodic_dot(normal,normal);
    if(!(determinant>0.0)) return {};
    const double a=0.5*g22*(g11-g12)/determinant;
    const double b=0.5*g11*(g22-g12)/determinant;
    const Point<3> delta{a*d1[0]+b*d2[0],a*d1[1]+b*d2[1],a*d1[2]+b*d2[2]};
    return {{p0[0]+delta[0],p0[1]+delta[1],p0[2]+delta[2]},
            compact_periodic_dot(delta,delta)};
}

inline compact_periodic_sphere compact_periodic_tetra_sphere(
    const std::array<Point<3>,4>& points
) {
    // Solve in a vertex-relative frame, avoiding differences of large norms.
    std::array<std::array<double,4>,3> matrix{};
    for(std::size_t row=0;row<3;++row) {
        const auto delta=compact_periodic_subtract(points[row+1],points[0]);
        for(std::size_t d=0;d<3;++d) matrix[row][d]=delta[d];
        matrix[row][3]=0.5*compact_periodic_dot(delta,delta);
    }
    for(std::size_t column=0;column<3;++column) {
        std::size_t pivot=column;
        for(std::size_t row=column+1;row<3;++row)
            if(std::abs(matrix[row][column])>std::abs(matrix[pivot][column])) pivot=row;
        if(matrix[pivot][column]==0.0) return {};
        std::swap(matrix[pivot],matrix[column]);
        const double divisor=matrix[column][column];
        for(std::size_t j=column;j<4;++j) matrix[column][j]/=divisor;
        for(std::size_t row=0;row<3;++row) {
            if(row==column) continue;
            const double factor=matrix[row][column];
            for(std::size_t j=column;j<4;++j) matrix[row][j]-=factor*matrix[column][j];
        }
    }
    const Point<3> delta{matrix[0][3],matrix[1][3],matrix[2][3]};
    return {{points[0][0]+delta[0],points[0][1]+delta[1],points[0][2]+delta[2]},
            compact_periodic_dot(delta,delta)};
}

inline bool compact_periodic_inside(
    const compact_periodic_sphere& sphere,const Point<3>& witness
) {
    const auto delta=compact_periodic_subtract(sphere.center,witness);
    return compact_periodic_dot(delta,delta)<sphere.alpha;
}

template<std::size_t N>
std::array<unsigned,N> compact_periodic_key(std::array<unsigned,N> vertices) {
    std::sort(vertices.begin(),vertices.end());
    return vertices;
}

template<std::size_t N>
using compact_periodic_values=
    std::unordered_map<std::array<unsigned,N>,double,VertexKeyHash<N>>;

template<std::size_t N>
void compact_periodic_relax(
    compact_periodic_values<N>& values,const std::array<unsigned,N>& key,double alpha
) {
    auto inserted=values.emplace(key,alpha);
    if(!inserted.second) inserted.first->second=std::min(inserted.first->second,alpha);
}

struct compact_periodic_edge {
    std::array<int,3> offset{};
    double own_alpha=std::numeric_limits<double>::infinity();
    double minimum_facet_alpha=std::numeric_limits<double>::infinity();
    bool gabriel=true;
};

template<class Points,class CB>
bool try_compact_periodic3(
    const Points& p,const CB& cb,std::array<double,3> from,std::array<double,3> to
) {
    initialize_geogram();
    const auto vertices=unique_points<Points,3>(p,false,&from);
    Point<3> period{};
    for(std::size_t d=0;d<3;++d) {
        period[d]=to[d]-from[d];
        if(!std::isfinite(from[d]) || !std::isfinite(to[d]) ||
           !std::isfinite(period[d]) || !(period[d]>0.0)) return false;
        // Leave out-of-domain input to the existing API implementation rather
        // than changing its contract by wrapping or discarding sites.
        for(const auto& vertex:vertices)
            if(!(vertex.point[d]>=0.0 && vertex.point[d]<period[d])) return false;
    }
    // The native backend starts with an ordinary, full-dimensional tetrahedron.
    if(!compact_alpha3_full_dimension(vertices)) return false;
    if(vertices.size()>std::numeric_limits<GEO::index_t>::max()/27) return false;
    std::vector<double> coordinates;
    coordinates.reserve(3*vertices.size());
    for(const auto& vertex:vertices)
        for(double coordinate:vertex.point) coordinates.push_back(coordinate);

    compact_periodic_values<4> cells;
    compact_periodic_values<3> facets;
    std::unordered_map<std::array<unsigned,2>,compact_periodic_edge,VertexKeyHash<2>> edges;
    {
        std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
        GEO::Numeric::random_reset();
        GEO::PeriodicDelaunay3d triangulation(GEO::vec3(period[0],period[1],period[2]));
        triangulation.set_vertices(static_cast<GEO::index_t>(vertices.size()),coordinates.data());
        triangulation.compute();
        if(triangulation.has_empty_cells() || triangulation.cell_size()!=4 ||
           triangulation.nb_cells()==0) return false;

        const std::size_t occurrence_count=triangulation.nb_cells();
        std::vector<double> occurrence_alpha(occurrence_count);
        cells.reserve(occurrence_count);
        edges.reserve(occurrence_count);
        for(GEO::index_t cell=0;cell<triangulation.nb_cells();++cell) {
            std::array<unsigned,4> ids{};
            std::array<Point<3>,4> points{};
            std::array<std::array<int,3>,4> offsets{};
            for(GEO::index_t local=0;local<4;++local) {
                const auto pv=triangulation.cell_vertex(cell,local);
                ids[local]=vertices[triangulation.periodic_vertex_real(pv)].id;
                points[local]=compact_periodic_point(triangulation,pv);
                triangulation.periodic_vertex_get_T(
                    pv,offsets[local][0],offsets[local][1],offsets[local][2]
                );
            }
            const auto key=compact_periodic_key(ids);
            if(std::adjacent_find(key.begin(),key.end())!=key.end()) return false;
            const double alpha=compact_periodic_tetra_sphere(points).alpha;
            if(!std::isfinite(alpha)) return false;
            occurrence_alpha[cell]=alpha;
            compact_periodic_relax(cells,key,alpha);

            for(std::size_t i=0;i<4;++i) for(std::size_t j=i+1;j<4;++j) {
                const std::size_t a=ids[i]<ids[j] ? i : j;
                const std::size_t b=ids[i]<ids[j] ? j : i;
                const std::array<unsigned,2> edge_key{ids[a],ids[b]};
                std::array<int,3> offset{};
                for(std::size_t d=0;d<3;++d) offset[d]=offsets[b][d]-offsets[a][d];
                auto inserted=edges.try_emplace(edge_key);
                auto& edge=inserted.first->second;
                if(inserted.second) edge.offset=offset;
                else if(edge.offset!=offset) return false;
                // Equality of every edge lift implies equality of the lift of
                // every higher simplex. Never min-reduce distinct coverings.
                const auto sphere=compact_periodic_edge_sphere(points[a],points[b]);
                if(!std::isfinite(sphere.alpha)) return false;
                edge.own_alpha=std::min(edge.own_alpha,sphere.alpha);
                if(edge.gabriel)
                    for(std::size_t k=0;k<4;++k)
                        if(k!=a && k!=b && compact_periodic_inside(sphere,points[k]))
                            edge.gabriel=false;
            }
        }

        facets.reserve(2*cells.size());
        for(GEO::index_t cell=0;cell<triangulation.nb_cells();++cell) {
            for(GEO::index_t opposite=0;opposite<4;++opposite) {
                const auto adjacent=triangulation.cell_adjacent(cell,opposite);
                // Compressed boundary copies may lack an opposite cell. An
                // equivalent complete copy supplies both Gabriel witnesses.
                if(adjacent==GEO::NO_INDEX || adjacent<cell) continue;
                std::array<unsigned,3> ids{};
                std::array<Point<3>,3> points{};
                std::size_t next=0;
                for(GEO::index_t local=0;local<4;++local) if(local!=opposite) {
                    const auto pv=triangulation.cell_vertex(cell,local);
                    ids[next]=vertices[triangulation.periodic_vertex_real(pv)].id;
                    points[next++]=compact_periodic_point(triangulation,pv);
                }
                const auto sphere=compact_periodic_triangle_sphere(points[0],points[1],points[2]);
                if(!std::isfinite(sphere.alpha)) return false;
                GEO::index_t adjacent_opposite=0;
                while(adjacent_opposite<4 &&
                      triangulation.cell_adjacent(adjacent,adjacent_opposite)!=cell)
                    ++adjacent_opposite;
                if(adjacent_opposite==4) return false;
                const bool gabriel=!compact_periodic_inside(
                    sphere,compact_periodic_point(triangulation,triangulation.cell_vertex(cell,opposite))
                ) && !compact_periodic_inside(
                    sphere,compact_periodic_point(
                        triangulation,triangulation.cell_vertex(adjacent,adjacent_opposite)
                    )
                );
                const double coface_alpha=std::min(occurrence_alpha[cell],occurrence_alpha[adjacent]);
                compact_periodic_relax(
                    facets,compact_periodic_key(ids),gabriel ? std::min(sphere.alpha,coface_alpha) : coface_alpha
                );
            }
        }
        // Clamp after canonical copy reduction: translated sphere solves can
        // differ by ulps, but the face/coface inequality must hold exactly.
        for(const auto& cell:cells) for(std::size_t removed=0;removed<4;++removed) {
            std::array<unsigned,3> face{};
            std::size_t next=0;
            for(std::size_t i=0;i<4;++i) if(i!=removed) face[next++]=cell.first[i];
            auto found=facets.find(face);
            if(found==facets.end()) return false;
            found->second=std::min(found->second,cell.second);
        }
        for(const auto& face:facets) for(std::size_t removed=0;removed<3;++removed) {
            std::array<unsigned,2> edge_key{};
            std::size_t next=0;
            for(std::size_t i=0;i<3;++i) if(i!=removed) edge_key[next++]=face.first[i];
            auto& edge=edges.at(edge_key);
            edge.minimum_facet_alpha=std::min(edge.minimum_facet_alpha,face.second);
        }
        if(vertices.size()+facets.size()!=edges.size()+cells.size()) return false;
    }
    // No user callback runs until dimensional and covering eligibility, and
    // all assignments, have succeeded. The global backend lock is released.
    for(const auto& vertex:vertices) cb(std::array<unsigned,1>{vertex.id},0.0);
    for(const auto& item:edges) {
        const auto& edge=item.second;
        cb(item.first,edge.gabriel ? std::min(edge.own_alpha,edge.minimum_facet_alpha)
                                  : edge.minimum_facet_alpha);
    }
    for(const auto& item:facets) cb(item.first,item.second);
    for(const auto& item:cells) cb(item.first,item.second);
    return true;
}
