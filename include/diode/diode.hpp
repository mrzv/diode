#include <geogram/basic/numeric.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/delaunay/delaunay_3d.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace diode {
namespace detail {

inline void initialize_geogram() {
    static std::once_flag once;
    std::call_once(once, [] { GEO::initialize(); });
}
inline std::mutex& geogram_triangulation_mutex() {
    static std::mutex mutex;
    return mutex;
}


template<std::size_t D> using Point = std::array<double, D>;

template<std::size_t D> struct Vertex {
    unsigned id;
    Point<D> point;
    std::array<int, D> offset{};
    double weight = 0.0;
};

template<std::size_t N, std::size_t D> struct LiftKey {
    std::array<unsigned, N> vertices{};
    std::array<std::array<int, D>, N> offsets{};
    bool operator<(const LiftKey& rhs) const {
        return vertices < rhs.vertices || (!(rhs.vertices < vertices) && offsets < rhs.offsets);
    }
};

template<std::size_t N, std::size_t D> struct Record {
    LiftKey<N,D> key;
    std::array<Point<D>,N> points{};
    std::array<double,N> weights{};
    std::vector<std::pair<Point<D>,double>> witnesses;
    double alpha = std::numeric_limits<double>::infinity();
    std::vector<unsigned> tau;
};

template<std::size_t N, std::size_t D>
LiftKey<N,D> canonicalize(std::array<Vertex<D>,N>& simplex) {
    std::sort(simplex.begin(), simplex.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
    LiftKey<N,D> key;
    const auto base = simplex[0].offset;
    for(std::size_t i=0; i<N; ++i) {
        key.vertices[i] = simplex[i].id;
        for(std::size_t d=0; d<D; ++d) key.offsets[i][d] = simplex[i].offset[d] - base[d];
    }
    return key;
}

template<std::size_t N, std::size_t D>
Record<N,D>& add_record(std::map<LiftKey<N,D>,Record<N,D>>& records,
                        std::array<Vertex<D>,N> simplex,
                        Point<D>* frame_shift=nullptr) {
    auto key = canonicalize(simplex);
    auto [it, inserted] = records.emplace(key, Record<N,D>{});
    if(inserted) {
        it->second.key = key;
        for(std::size_t i=0; i<N; ++i) {
            it->second.points[i] = simplex[i].point;
            it->second.weights[i] = simplex[i].weight;
        }
    }
    if(frame_shift)
        for(std::size_t d=0; d<D; ++d)
            (*frame_shift)[d]=it->second.points[0][d]-simplex[0].point[d];
    return it->second;
}

template<std::size_t N, std::size_t D>
Record<N,D>& add_record_with_witness(
    std::map<LiftKey<N,D>,Record<N,D>>& records,
    std::array<Vertex<D>,N> simplex,
    const Vertex<D>& witness
) {
    Point<D> shift{};
    auto& record=add_record(records,std::move(simplex),&shift);
    Point<D> point=witness.point;
    for(std::size_t d=0; d<D; ++d) point[d]+=shift[d];
    record.witnesses.emplace_back(point,witness.weight);
    return record;
}

template<std::size_t M>
bool solve(std::array<std::array<double,M>,M> a, std::array<double,M> b,
           std::array<double,M>& x) {
    double matrix_scale=0.0;
    for(const auto& row:a) for(double value:row)
        matrix_scale=std::max(matrix_scale,std::abs(value));
    if(matrix_scale==0.0) return false;
    for(std::size_t k=0; k<M; ++k) {
        std::size_t pivot=k;
        for(std::size_t i=k+1; i<M; ++i) if(std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot=i;
        if(std::abs(a[pivot][k]) <=
           64.0*std::numeric_limits<double>::epsilon()*matrix_scale) return false;
        std::swap(a[k],a[pivot]); std::swap(b[k],b[pivot]);
        for(std::size_t i=k+1; i<M; ++i) {
            const double q=a[i][k]/a[k][k];
            for(std::size_t j=k; j<M; ++j) a[i][j]-=q*a[k][j];
            b[i]-=q*b[k];
        }
    }
    for(std::size_t ii=M; ii-- > 0;) {
        double v=b[ii];
        for(std::size_t j=ii+1; j<M; ++j) v-=a[ii][j]*x[j];
        x[ii]=v/a[ii][ii];
    }
    return true;
}

template<std::size_t N, std::size_t D>
double sphere(const std::array<Point<D>,N>& p, const std::array<double,N>& w,
              Point<D>* center_out=nullptr) {
    Point<D> center=p[0];
    if constexpr(N > 1) {
        constexpr std::size_t M=N-1;
        std::array<Point<D>,M> v{};
        std::array<std::array<double,M>,M> gram{};
        std::array<double,M> rhs{}, coeff{};
        for(std::size_t i=0; i<M; ++i) {
            for(std::size_t d=0; d<D; ++d) v[i][d]=p[i+1][d]-p[0][d];
            double n2=0; for(double z:v[i]) n2+=z*z;
            rhs[i]=0.5*(n2-w[i+1]+w[0]);
        }
        for(std::size_t i=0; i<M; ++i) for(std::size_t j=0; j<M; ++j)
            for(std::size_t d=0; d<D; ++d) gram[i][j]+=v[i][d]*v[j][d];
        if(!solve<M>(gram,rhs,coeff)) return std::numeric_limits<double>::infinity();
        for(std::size_t i=0; i<M; ++i) for(std::size_t d=0; d<D; ++d) center[d]+=coeff[i]*v[i][d];
    }
    if(center_out) *center_out=center;
    double r=-w[0];
    for(std::size_t d=0; d<D; ++d) { double q=center[d]-p[0][d]; r+=q*q; }
    return r;
}

template<std::size_t N, std::size_t D>
void own_alpha(Record<N,D>& r) {
    Point<D> center{};
    const double radius=sphere(r.points,r.weights,&center);
    bool gabriel=std::isfinite(radius);
    for(const auto& witness:r.witnesses) {
        double power=-witness.second;
        for(std::size_t d=0; d<D; ++d) { double q=center[d]-witness.first[d]; power+=q*q; }
        const double tolerance=128.0*std::numeric_limits<double>::epsilon()*
            std::max(std::abs(radius),std::abs(power));
        if(power < radius-tolerance) { gabriel=false; break; }
    }
    if(gabriel) { r.alpha=radius; r.tau.assign(r.key.vertices.begin(),r.key.vertices.end()); }
}

template<std::size_t N, std::size_t D>
void take_coface(Record<N,D>& r, double alpha, const std::vector<unsigned>& tau) {
    if(alpha < r.alpha) { r.alpha=alpha; r.tau=tau; }
}

template<class Points, std::size_t D>
std::vector<Vertex<D>> unique_points(const Points& points, bool weighted=false,
                                     const Point<D>* from=nullptr) {
    struct Selected {
        unsigned id;
        double weight;
    };
    std::map<Point<D>,Selected> selected;
    for(unsigned i=0; i<points.size(); ++i) {
        Point<D> p{};
        for(std::size_t d=0; d<D; ++d) {
            p[d]=static_cast<double>(points(i,d));
            if(!std::isfinite(p[d])) throw std::runtime_error("points must be finite");
        }
        const double weight=weighted ? static_cast<double>(points(i,D)) : 0.0;
        if(!std::isfinite(weight)) throw std::runtime_error("weights must be finite");
        auto it=selected.find(p);
        if(it==selected.end())
            selected.emplace(p,Selected{i,weight});
        else if(!weighted || weight>=it->second.weight)
            it->second={i,weight};
    }
    std::vector<Vertex<D>> result;
    result.reserve(selected.size());
    for(const auto& item:selected) {
        Vertex<D> vertex;
        vertex.id=item.second.id;
        vertex.point=item.first;
        vertex.weight=item.second.weight;
        if(from) for(std::size_t d=0; d<D; ++d) vertex.point[d]-=(*from)[d];
        result.push_back(vertex);
    }
    return result;
}

template<std::size_t D> struct Complex;
template<> struct Complex<2> {
    std::map<LiftKey<1,2>,Record<1,2>> v;
    std::map<LiftKey<2,2>,Record<2,2>> e;
    std::map<LiftKey<3,2>,Record<3,2>> f;
    std::vector<std::pair<LiftKey<1,2>,LiftKey<2,2>>> ve;
    std::vector<std::pair<LiftKey<2,2>,LiftKey<3,2>>> ef;
};
template<> struct Complex<3> {
    std::map<LiftKey<1,3>,Record<1,3>> v;
    std::map<LiftKey<2,3>,Record<2,3>> e;
    std::map<LiftKey<3,3>,Record<3,3>> f;
    std::map<LiftKey<4,3>,Record<4,3>> c;
    std::vector<std::pair<LiftKey<1,3>,LiftKey<2,3>>> ve;
    std::vector<std::pair<LiftKey<2,3>,LiftKey<3,3>>> ef;
    std::vector<std::pair<LiftKey<3,3>,LiftKey<4,3>>> fc;
};

inline Record<2,2>& add_edge(Complex<2>& out, const std::array<Vertex<2>,2>& edge) {
    auto& record=add_record(out.e,edge);
    for(int i=0;i<2;++i) {
        std::array<Vertex<2>,1> one{edge[i]};
        auto& vertex=add_record_with_witness(out.v,one,edge[1-i]);
        out.ve.emplace_back(vertex.key,record.key);
    }
    return record;
}

inline void add_triangle(Complex<2>& out, const std::array<Vertex<2>,3>& cell) {
    auto& face=add_record(out.f,cell);
    for(int skip=0; skip<3; ++skip) {
        std::array<Vertex<2>,2> edge{cell[(skip+1)%3],cell[(skip+2)%3]};
        auto& er=add_edge(out,edge);
        Point<2> shift{};
        add_record(out.e,edge,&shift);
        Point<2> witness=cell[skip].point;
        for(std::size_t d=0;d<2;++d) witness[d]+=shift[d];
        er.witnesses.emplace_back(witness,cell[skip].weight);
        out.ef.emplace_back(er.key,face.key);
    }
}

inline void add_triangle(Complex<3>& out, const std::array<Vertex<3>,3>& cell) {
    auto& face=add_record(out.f,cell);
    for(int skip=0; skip<3; ++skip) {
        std::array<Vertex<3>,2> edge{cell[(skip+1)%3],cell[(skip+2)%3]};
        auto& er=add_record(out.e,edge);
        out.ef.emplace_back(er.key,face.key);
        for(int i=0;i<2;++i) {
            std::array<Vertex<3>,1> one{edge[i]};
            auto& vr=add_record(out.v,one);
            out.ve.emplace_back(vr.key,er.key);
        }
    }
}

inline Complex<3> lower_dimensional_complex(
    std::vector<Vertex<3>> vertices, bool weighted
) {
    Complex<3> out;
    if(vertices.empty()) return out;
    const auto origin=vertices.front().point;
    Point<3> axis{};
    double axis_norm=0.0;
    for(const auto& vertex:vertices) {
        Point<3> delta{};
        double norm=0.0;
        for(int d=0;d<3;++d) {
            delta[d]=vertex.point[d]-origin[d];
            norm+=delta[d]*delta[d];
        }
        if(norm>axis_norm) {
            axis_norm=norm;
            axis=delta;
        }
    }
    if(axis_norm==0.0) {
        add_record(out.v,std::array<Vertex<3>,1>{vertices.front()});
        return out;
    }
    for(double& value:axis) value/=std::sqrt(axis_norm);

    Point<3> second{};
    double second_norm=0.0;
    for(const auto& vertex:vertices) {
        Point<3> delta{};
        double projection=0.0;
        for(int d=0;d<3;++d) {
            delta[d]=vertex.point[d]-origin[d];
            projection+=delta[d]*axis[d];
        }
        double norm=0.0;
        for(int d=0;d<3;++d) {
            delta[d]-=projection*axis[d];
            norm+=delta[d]*delta[d];
        }
        if(norm>second_norm) {
            second_norm=norm;
            second=delta;
        }
    }
    const double relative_tolerance=
        1024.0*std::numeric_limits<double>::epsilon();
    const double tolerance=relative_tolerance*relative_tolerance*axis_norm;
    if(second_norm<=tolerance) {
        struct LineSite {
            double coordinate;
            Vertex<3> vertex;
        };
        std::vector<LineSite> sites;
        sites.reserve(vertices.size());
        for(const auto& vertex:vertices) {
            double coordinate=0.0;
            for(int d=0;d<3;++d)
                coordinate+=(vertex.point[d]-origin[d])*axis[d];
            sites.push_back({coordinate,vertex});
        }
        std::sort(sites.begin(),sites.end(),[](const auto& lhs,const auto& rhs) {
            return lhs.coordinate<rhs.coordinate;
        });
        std::vector<LineSite> hull;
        std::vector<double> starts;
        for(const auto& site:sites) {
            double start=-std::numeric_limits<double>::infinity();
            while(!hull.empty()) {
                const auto& previous=hull.back();
                start=((site.coordinate*site.coordinate-site.vertex.weight)-
                       (previous.coordinate*previous.coordinate-previous.vertex.weight))/
                      (2.0*(site.coordinate-previous.coordinate));
                if(hull.size()==1 || start>starts.back()) break;
                hull.pop_back();
                starts.pop_back();
            }
            hull.push_back(site);
            starts.push_back(start);
        }
        for(const auto& site:hull)
            add_record(out.v,std::array<Vertex<3>,1>{site.vertex});
        for(std::size_t i=1;i<hull.size();++i) {
            auto& edge=add_record(
                out.e,std::array<Vertex<3>,2>{hull[i-1].vertex,hull[i].vertex}
            );
            for(int endpoint=0;endpoint<2;++endpoint) {
                const auto& vertex=endpoint==0 ? hull[i-1].vertex : hull[i].vertex;
                auto& record=add_record(out.v,std::array<Vertex<3>,1>{vertex});
                out.ve.emplace_back(record.key,edge.key);
            }
        }
        return out;
    }

    for(double& value:second) value/=std::sqrt(second_norm);
    std::vector<double> coordinates;
    coordinates.reserve(vertices.size()*(weighted ? 3 : 2));
    double max_weight=0.0;
    if(weighted) {
        max_weight=vertices.front().weight;
        for(const auto& vertex:vertices)
            max_weight=std::max(max_weight,vertex.weight);
    }
    for(const auto& vertex:vertices) {
        double x=0.0;
        double y=0.0;
        for(int d=0;d<3;++d) {
            const double delta=vertex.point[d]-origin[d];
            x+=delta*axis[d];
            y+=delta*second[d];
        }
        coordinates.push_back(x);
        coordinates.push_back(y);
        if(weighted)
            coordinates.push_back(std::sqrt(std::max(0.0,max_weight-vertex.weight)));
    }
    std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
    GEO::Numeric::random_reset();
    GEO::SmartPointer<GEO::Delaunay> triangulation=weighted
        ? static_cast<GEO::Delaunay*>(new GEO::RegularWeightedDelaunay2d())
        : static_cast<GEO::Delaunay*>(new GEO::Delaunay2d());
    triangulation->set_reorder(false);
    triangulation->set_vertices(vertices.size(),coordinates.data());
    for(GEO::index_t cell_index=0;cell_index<triangulation->nb_cells();++cell_index) {
        std::array<Vertex<3>,3> cell;
        for(int i=0;i<3;++i)
            cell[i]=vertices[triangulation->cell_vertex(cell_index,i)];
        add_triangle(out,cell);
    }
    return out;
}

inline bool spans_three_dimensions(const std::vector<Vertex<3>>& vertices) {
    if(vertices.size()<4) return false;
    const auto origin=vertices.front().point;
    Point<3> first{};
    double first_norm=0.0;
    for(const auto& vertex:vertices) {
        Point<3> q{}; double n=0.0;
        for(int d=0;d<3;++d) { q[d]=vertex.point[d]-origin[d]; n+=q[d]*q[d]; }
        if(n>first_norm) { first_norm=n; first=q; }
    }
    Point<3> normal{};
    double normal_norm=0.0;
    for(const auto& vertex:vertices) {
        Point<3> q{};
        for(int d=0;d<3;++d) q[d]=vertex.point[d]-origin[d];
        Point<3> n{
            first[1]*q[2]-first[2]*q[1],
            first[2]*q[0]-first[0]*q[2],
            first[0]*q[1]-first[1]*q[0]
        };
        double n2=0.0; for(double value:n) n2+=value*value;
        if(n2>normal_norm) { normal_norm=n2; normal=n; }
    }
    if(normal_norm==0.0) return false;
    double height=0.0;
    for(const auto& vertex:vertices) {
        double value=0.0;
        for(int d=0;d<3;++d) value+=(vertex.point[d]-origin[d])*normal[d];
        height=std::max(height,std::abs(value));
    }
    const double scale=std::sqrt(first_norm*normal_norm);
    return height>1024.0*std::numeric_limits<double>::epsilon()*scale;
}

inline void add_tetrahedron(Complex<3>& out, const std::array<Vertex<3>,4>& cell) {
    auto& cr=add_record(out.c,cell);
    for(int omit=0; omit<4; ++omit) {
        std::array<Vertex<3>,3> face{}; int k=0;
        for(int i=0;i<4;++i) if(i!=omit) face[k++]=cell[i];
        auto& fr=add_record_with_witness(out.f,face,cell[omit]);
        out.fc.emplace_back(fr.key,cr.key);
        for(int a=0;a<3;++a) for(int b=a+1;b<3;++b) {
            std::array<Vertex<3>,2> edge{face[a],face[b]};
            auto& er=add_record_with_witness(out.e,edge,face[3-a-b]);
            out.ef.emplace_back(er.key,fr.key);
            for(int i=0;i<2;++i) {
                std::array<Vertex<3>,1> one{edge[i]};
                auto& vr=add_record_with_witness(out.v,one,edge[1-i]);
                out.ve.emplace_back(vr.key,er.key);
            }
        }
    }
}

inline void compute_alpha(Complex<2>& x) {
    for(auto& p:x.f) own_alpha(p.second);
    for(auto& p:x.e) own_alpha(p.second);
    for(const auto& r:x.ef) take_coface(x.e.at(r.first),x.f.at(r.second).alpha,x.f.at(r.second).tau);
    for(auto& p:x.v) own_alpha(p.second);
    for(const auto& r:x.ve) take_coface(x.v.at(r.first),x.e.at(r.second).alpha,x.e.at(r.second).tau);
}
inline void compute_alpha(Complex<3>& x) {
    for(auto& p:x.c) own_alpha(p.second);
    for(auto& p:x.f) own_alpha(p.second);
    for(const auto& r:x.fc) take_coface(x.f.at(r.first),x.c.at(r.second).alpha,x.c.at(r.second).tau);
    for(auto& p:x.e) own_alpha(p.second);
    for(const auto& r:x.ef) take_coface(x.e.at(r.first),x.f.at(r.second).alpha,x.f.at(r.second).tau);
    for(auto& p:x.v) own_alpha(p.second);
    for(const auto& r:x.ve) take_coface(x.v.at(r.first),x.e.at(r.second).alpha,x.e.at(r.second).tau);
}

inline Complex<2> lower_dimensional_complex(std::vector<Vertex<2>> vertices) {
    Complex<2> out;
    if(vertices.empty()) return out;
    if(vertices.size()==1) {
        add_record(out.v,std::array<Vertex<2>,1>{vertices.front()});
        return out;
    }
    const auto origin=vertices.front().point;
    Point<2> axis{};
    double axis_norm=0.0;
    for(const auto& vertex:vertices) {
        Point<2> delta{};
        double norm=0.0;
        for(int d=0;d<2;++d) {
            delta[d]=vertex.point[d]-origin[d];
            norm+=delta[d]*delta[d];
        }
        if(norm>axis_norm) {
            axis_norm=norm;
            axis=delta;
        }
    }
    for(double& value:axis) value/=std::sqrt(axis_norm);
    std::sort(vertices.begin(),vertices.end(),[&](const auto& lhs,const auto& rhs) {
        double left=0.0;
        double right=0.0;
        for(int d=0;d<2;++d) {
            left+=(lhs.point[d]-origin[d])*axis[d];
            right+=(rhs.point[d]-origin[d])*axis[d];
        }
        return left<right;
    });
    for(const auto& vertex:vertices)
        add_record(out.v,std::array<Vertex<2>,1>{vertex});
    for(std::size_t i=1;i<vertices.size();++i)
        add_edge(out,std::array<Vertex<2>,2>{vertices[i-1],vertices[i]});
    return out;
}

template<class Points>
Complex<2> triangulate2(const Points& points, bool periodic=false,
                        Point<2> from={0,0}, Point<2> to={1,1}) {
    initialize_geogram();
    auto base=unique_points<Points,2>(points,false,periodic?&from:nullptr);
    std::vector<Vertex<2>> vertices;
    if(periodic) {
        vertices.reserve(base.size()*9);
        for(int ox=-1;ox<=1;++ox) for(int oy=-1;oy<=1;++oy) for(auto v:base) {
            v.offset={ox,oy}; v.point[0]+=ox*(to[0]-from[0]); v.point[1]+=oy*(to[1]-from[1]); vertices.push_back(v);
        }
    } else vertices=base;
    Complex<2> out;
    if(base.size()<3)
        return periodic ? out : lower_dimensional_complex(std::move(base));
    std::vector<double> coords; coords.reserve(vertices.size()*2);
    for(const auto& v:vertices) { coords.push_back(v.point[0]); coords.push_back(v.point[1]); }
    std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
    GEO::Numeric::random_reset();
    GEO::SmartPointer<GEO::Delaunay2d> dt=new GEO::Delaunay2d();
    dt->set_reorder(false);
    dt->set_vertices(vertices.size(),coords.data());
    for(GEO::index_t c=0;c<dt->nb_cells();++c) {
        std::array<Vertex<2>,3> cell;
        for(int i=0;i<3;++i) cell[i]=vertices[dt->cell_vertex(c,i)];
        if(periodic) {
            std::array<Point<2>,3> p{cell[0].point,cell[1].point,cell[2].point};
            std::array<double,3> w{}; Point<2> center{}; sphere(p,w,&center);
            const double tolerance_x=1024.0*std::numeric_limits<double>::epsilon()*
                (to[0]-from[0]);
            const double tolerance_y=1024.0*std::numeric_limits<double>::epsilon()*
                (to[1]-from[1]);
            if(center[0]<-tolerance_x || center[0]>to[0]-from[0]+tolerance_x ||
               center[1]<-tolerance_y || center[1]>to[1]-from[1]+tolerance_y) continue;
        }
        add_triangle(out,cell);
    }
    if(!periodic && out.f.empty()) return lower_dimensional_complex(std::move(base));
    return out;
}


inline Complex<3> triangulate3_regular(
    const std::vector<Vertex<3>>& vertices, bool weighted,
    const Point<3>* periodic_extent=nullptr
) {
    std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
    GEO::Numeric::random_reset();
    Complex<3> out;
    std::vector<double> coordinates;
    coordinates.reserve(vertices.size()*(weighted ? 4 : 3));
    double max_weight=0.0;
    if(weighted) {
        max_weight=vertices.front().weight;
        for(const auto& vertex:vertices) max_weight=std::max(max_weight,vertex.weight);
    }
    for(const auto& vertex:vertices) {
        for(double x:vertex.point) coordinates.push_back(x);
        if(weighted) coordinates.push_back(std::sqrt(std::max(0.0,max_weight-vertex.weight)));
    }
    GEO::SmartPointer<GEO::Delaunay> dt = weighted
        ? static_cast<GEO::Delaunay*>(new GEO::RegularWeightedDelaunay3d())
        : static_cast<GEO::Delaunay*>(new GEO::Delaunay3d());
    dt->set_reorder(false);
    dt->set_vertices(vertices.size(),coordinates.data());
    for(GEO::index_t c=0;c<dt->nb_cells();++c) {
        std::array<Vertex<3>,4> cell;
        bool repeated=false;
        for(int i=0;i<4;++i) {
            cell[i]=vertices[dt->cell_vertex(c,i)];
            for(int j=0;j<i;++j) if(cell[j].id==cell[i].id) repeated=true;
        }
        if(repeated) continue;
        if(periodic_extent) {
            std::array<Point<3>,4> p{};
            std::array<double,4> w{};
            Point<3> center{};
            for(int i=0;i<4;++i) { p[i]=cell[i].point; w[i]=cell[i].weight; }
            sphere(p,w,&center);
            bool central=true;
            for(int d=0;d<3;++d) {
                const double tolerance=1024.0*std::numeric_limits<double>::epsilon()*
                    (*periodic_extent)[d];
                central=central && center[d]>=-tolerance &&
                    center[d]<=(*periodic_extent)[d]+tolerance;
            }
            if(!central) continue;
        }
        add_tetrahedron(out,cell);
    }
    return out;
}

inline Complex<3> triangulate3_tiled(
    const std::vector<Vertex<3>>& base, bool weighted, const Point<3>& period
) {
    std::vector<Vertex<3>> vertices;
    vertices.reserve(base.size()*27);
    for(int ox=-1;ox<=1;++ox) for(int oy=-1;oy<=1;++oy) for(int oz=-1;oz<=1;++oz)
        for(auto vertex:base) {
            vertex.offset={ox,oy,oz};
            vertex.point[0]+=ox*period[0];
            vertex.point[1]+=oy*period[1];
            vertex.point[2]+=oz*period[2];
            vertices.push_back(vertex);
        }
    return triangulate3_regular(vertices,weighted,&period);
}

template<class Points>
Complex<3> triangulate3(const Points& points, bool weighted=false, bool periodic=false,
                        Point<3> from={0,0,0}, Point<3> to={1,1,1},
                        bool preserve_lower_dimension=false) {
    initialize_geogram();
    auto vertices=unique_points<Points,3>(points,weighted,periodic?&from:nullptr);
    if(!periodic && !spans_three_dimensions(vertices))
        return preserve_lower_dimension
            ? lower_dimensional_complex(std::move(vertices),weighted)
            : Complex<3>{};
    if(vertices.size()<4) return Complex<3>{};
    if(!periodic) return triangulate3_regular(vertices,weighted);
    Point<3> period{to[0]-from[0],to[1]-from[1],to[2]-from[2]};
    return triangulate3_tiled(vertices,weighted,period);
}

template<std::size_t N, std::size_t D, class Callback>
void emit_values_plain(const std::map<LiftKey<N,D>,Record<N,D>>& records, const Callback& cb) {
    std::map<std::array<unsigned,N>,const Record<N,D>*> unique;
    for(const auto& item:records) {
        auto [it,inserted]=unique.emplace(item.first.vertices,&item.second);
        if(!inserted)
            throw std::runtime_error("Cannot convert periodic triangulation to a one-sheeted covering");
    }
    for(const auto& item:unique) if(std::isfinite(item.second->alpha)) cb(item.first,item.second->alpha);
}

template<std::size_t N, std::size_t D, class Callback>
void emit_values_attachment(const std::map<LiftKey<N,D>,Record<N,D>>& records, const Callback& cb) {
    std::map<std::array<unsigned,N>,const Record<N,D>*> unique;
    for(const auto& item:records) {
        auto [it,inserted]=unique.emplace(item.first.vertices,&item.second);
        if(!inserted)
            throw std::runtime_error("Cannot convert periodic triangulation to a one-sheeted covering");
    }
    for(const auto& item:unique) if(std::isfinite(item.second->alpha)) cb(item.first,item.second->alpha,item.second->tau);
}

template<std::size_t N, std::size_t D, class Callback>
void emit_combinatorics(const std::map<LiftKey<N,D>,Record<N,D>>& records, const Callback& cb) {
    std::set<std::array<unsigned,N>> seen;
    for(const auto& item:records) {
        if(!seen.insert(item.first.vertices).second)
            throw std::runtime_error("Cannot convert periodic triangulation to a one-sheeted covering");
        cb(item.first.vertices);
    }
}

template<std::size_t N, std::size_t D, class Callback>
void emit_lifts(const std::map<LiftKey<N,D>,Record<N,D>>& records, const Callback& cb) {
    std::map<std::array<unsigned,N>,std::array<std::array<int,D>,N>> seen;
    for(const auto& item:records) {
        auto [it,inserted]=seen.emplace(item.first.vertices,item.first.offsets);
        if(!inserted && it->second!=item.first.offsets)
            throw std::runtime_error("Cannot convert periodic triangulation to a one-sheeted covering");
    }
    for(const auto& item:seen) cb(item.first,item.second);
}

template<class C, class CB> void emit_alpha(C& x,const CB& cb) {
    compute_alpha(x);
    emit_values_plain(x.v,cb); emit_values_plain(x.e,cb); emit_values_plain(x.f,cb);
}
template<class CB> void emit_alpha(Complex<3>& x,const CB& cb) {
    compute_alpha(x);
    emit_values_plain(x.v,cb); emit_values_plain(x.e,cb); emit_values_plain(x.f,cb); emit_values_plain(x.c,cb);
}
template<class C, class CB> void emit_attachment(C& x,const CB& cb) {
    compute_alpha(x);
    emit_values_attachment(x.v,cb); emit_values_attachment(x.e,cb); emit_values_attachment(x.f,cb);
}
template<class CB> void emit_attachment(Complex<3>& x,const CB& cb) {
    compute_alpha(x);
    emit_values_attachment(x.v,cb); emit_values_attachment(x.e,cb); emit_values_attachment(x.f,cb); emit_values_attachment(x.c,cb);
}
template<class C, class CB> void emit_delaunay(const C& x,const CB& cb) {
    emit_combinatorics(x.v,cb); emit_combinatorics(x.e,cb); emit_combinatorics(x.f,cb);
}
template<class CB> void emit_delaunay(const Complex<3>& x,const CB& cb) {
    emit_combinatorics(x.v,cb); emit_combinatorics(x.e,cb); emit_combinatorics(x.f,cb); emit_combinatorics(x.c,cb);
}
template<class C, class CB> void emit_periodic_lifts(const C& x,const CB& cb) {
    emit_lifts(x.v,cb); emit_lifts(x.e,cb); emit_lifts(x.f,cb);
}
template<class CB> void emit_periodic_lifts(const Complex<3>& x,const CB& cb) {
    emit_lifts(x.v,cb); emit_lifts(x.e,cb); emit_lifts(x.f,cb); emit_lifts(x.c,cb);
}

} // namespace detail

template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes(const Points& p,const SimplexCallback& cb) { fill_alpha_shapes_direct(p,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes_direct(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p); detail::emit_alpha(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes_with_attachment(const Points& p,const SimplexCallback& cb) { fill_alpha_shapes_direct_with_attachment(p,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes_direct_with_attachment(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p); detail::emit_attachment(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_alpha_shapes(const Points& p,const SimplexCallback& cb) { fill_weighted_alpha_shapes_direct(p,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_alpha_shapes_direct(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p,true); detail::emit_alpha(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_alpha_shapes(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { fill_periodic_alpha_shapes_direct(p,cb,a,b); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_alpha_shapes_direct(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,false,true,a,b); detail::emit_alpha(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_delaunay(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p,false,false,{0,0,0},{1,1,1},true); detail::emit_delaunay(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_delaunay(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p,true,false,{0,0,0},{1,1,1},true); detail::emit_delaunay(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_delaunay(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,false,true,a,b); detail::emit_delaunay(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_delaunay_lifts(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,false,true,a,b); detail::emit_periodic_lifts(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_periodic_alpha_shapes(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { fill_weighted_periodic_alpha_shapes_direct(p,cb,a,b); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_periodic_alpha_shapes_direct(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,true,true,a,b); detail::emit_alpha(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_periodic_delaunay(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,true,true,a,b); detail::emit_delaunay(x,cb); }

template<bool exact> template<class Points>
std::array<typename Points::Real,3> AlphaShapes<exact>::circumcenter(const Points& p) {
    if(p.size()!=3 && p.size()!=4) throw std::runtime_error("circumcenter expects three or four 3D points");
    if(p.size()==3) {
        std::array<detail::Point<3>,3> q{}; std::array<double,3>w{}; detail::Point<3> c{};
        for(int i=0;i<3;++i) for(int d=0;d<3;++d) q[i][d]=p(i,d); detail::sphere(q,w,&c);
        return {typename Points::Real(c[0]),typename Points::Real(c[1]),typename Points::Real(c[2])};
    }
    std::array<detail::Point<3>,4> q{}; std::array<double,4>w{}; detail::Point<3> c{};
    for(int i=0;i<4;++i) for(int d=0;d<3;++d) q[i][d]=p(i,d); detail::sphere(q,w,&c);
    return {typename Points::Real(c[0]),typename Points::Real(c[1]),typename Points::Real(c[2])};
}

template<bool exact,class Points,class SimplexCallback>
void fill_alpha_shapes2d(const Points& p,const SimplexCallback& cb) { fill_alpha_shapes2d_direct<exact>(p,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_alpha_shapes2d_direct(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate2(p); detail::emit_alpha(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_alpha_shapes2d_with_attachment(const Points& p,const SimplexCallback& cb) { fill_alpha_shapes2d_direct_with_attachment<exact>(p,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_alpha_shapes2d_direct_with_attachment(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate2(p); detail::emit_attachment(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_periodic_alpha_shapes2d(const Points& p,const SimplexCallback& cb,std::array<double,2> a,std::array<double,2> b) { fill_periodic_alpha_shapes2d_direct<exact>(p,cb,a,b); }
template<bool exact,class Points,class SimplexCallback>
void fill_periodic_alpha_shapes2d_direct(const Points& p,const SimplexCallback& cb,std::array<double,2> a,std::array<double,2> b) { auto x=detail::triangulate2(p,true,a,b); detail::emit_alpha(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_delaunay2d(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate2(p); detail::emit_delaunay(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_periodic_delaunay2d(const Points& p,const SimplexCallback& cb,std::array<double,2> a,std::array<double,2> b) { auto x=detail::triangulate2(p,true,a,b); detail::emit_delaunay(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_periodic_delaunay2d_lifts(const Points& p,const SimplexCallback& cb,std::array<double,2> a,std::array<double,2> b) { auto x=detail::triangulate2(p,true,a,b); detail::emit_periodic_lifts(x,cb); }

} // namespace diode
