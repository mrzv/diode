#include <geogram/basic/numeric.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/delaunay/delaunay_3d.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>
#include <geogram/numerics/expansion_nt.h>
#include <geogram/delaunay/periodic_delaunay_3d.h>
#include <cstdint>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
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
    bool operator==(const LiftKey& rhs) const {
        return vertices==rhs.vertices && offsets==rhs.offsets;
    }
};

inline std::size_t mix_hash(std::size_t seed, std::size_t value) {
    value^=value>>30;
    value*=static_cast<std::size_t>(0xbf58476d1ce4e5b9ULL);
    value^=value>>27;
    value*=static_cast<std::size_t>(0x94d049bb133111ebULL);
    value^=value>>31;
    return seed^(value+static_cast<std::size_t>(0x9e3779b97f4a7c15ULL)+(seed<<6)+(seed>>2));
}

template<std::size_t N, std::size_t D> struct LiftKeyHash {
    std::size_t operator()(const LiftKey<N,D>& key) const {
        std::size_t hash=0;
        for(unsigned vertex:key.vertices) hash=mix_hash(hash,vertex);
        for(const auto& offset:key.offsets)
            for(int value:offset)
                hash=mix_hash(hash,static_cast<std::size_t>(static_cast<unsigned>(value)));
        return hash;
    }
};

template<std::size_t N> struct VertexKeyHash {
    std::size_t operator()(const std::array<unsigned,N>& key) const {
        std::size_t hash=0;
        for(unsigned vertex:key) hash=mix_hash(hash,vertex);
        return hash;
    }
};

template<std::size_t N, std::size_t D>
double sphere(const std::array<Point<D>,N>& p, const std::array<double,N>& w,
              Point<D>* center_out);

template<std::size_t N, std::size_t D> struct Record {
    std::array<Point<D>,N> points{};
    std::array<double,N> weights{};
    Point<D> center{};
    double own_radius=std::numeric_limits<double>::infinity();
    double alpha=std::numeric_limits<double>::infinity();
    std::array<unsigned,D+1> tau{};
    unsigned char tau_size=0;
    bool gabriel=false;
};

template<std::size_t N, std::size_t D>
using RecordMap=std::unordered_map<LiftKey<N,D>,Record<N,D>,LiftKeyHash<N,D>>;

template<std::size_t N, std::size_t D> struct RecordInsertion {
    Record<N,D>* record;
    Point<D> frame_shift{};
    bool inserted;
};

template<std::size_t N, std::size_t D>
LiftKey<N,D> canonicalize(std::array<Vertex<D>,N>& simplex) {
    std::sort(simplex.begin(),simplex.end(),[](const auto& a,const auto& b) {
        return a.id<b.id;
    });
    LiftKey<N,D> key;
    const auto base=simplex[0].offset;
    for(std::size_t i=0;i<N;++i) {
        key.vertices[i]=simplex[i].id;
        for(std::size_t d=0;d<D;++d)
            key.offsets[i][d]=simplex[i].offset[d]-base[d];
    }
    return key;
}

template<std::size_t N, std::size_t D>
RecordInsertion<N,D> add_record(
    RecordMap<N,D>& records, std::array<Vertex<D>,N> simplex, bool collect_alpha
) {
    auto key=canonicalize(simplex);
    auto [it,inserted]=records.try_emplace(key);
    auto& record=it->second;
    if(inserted) {
        for(std::size_t i=0;i<N;++i) {
            record.points[i]=simplex[i].point;
            record.weights[i]=simplex[i].weight;
            record.tau[i]=key.vertices[i];
        }
        if(collect_alpha) {
            record.own_radius=sphere(record.points,record.weights,&record.center);
            record.gabriel=std::isfinite(record.own_radius);
        }
    }
    RecordInsertion<N,D> result{&record,{},inserted};
    for(std::size_t d=0;d<D;++d)
        result.frame_shift[d]=record.points[0][d]-simplex[0].point[d];
    return result;
}

template<std::size_t N, std::size_t D>
void add_witness(const RecordInsertion<N,D>& insertion,const Vertex<D>& witness) {
    auto& record=*insertion.record;
    if(!record.gabriel) return;
    double power=-witness.weight;
    for(std::size_t d=0;d<D;++d) {
        const double coordinate=witness.point[d]+insertion.frame_shift[d];
        const double delta=record.center[d]-coordinate;
        power+=delta*delta;
    }
    const double tolerance=128.0*std::numeric_limits<double>::epsilon()*
        std::max(std::abs(record.own_radius),std::abs(power));
    if(power<record.own_radius-tolerance) record.gabriel=false;
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

template<std::size_t N>
bool sphere3_expansion_offset(const std::array<Point<3>,N>& p,
                              const std::array<double,N>& w,
                              std::array<long double,3>& offset) {
    // Only ill-conditioned constructions take this allocating path. Exact
    // numerators and denominators retain full-rank simplices even when their
    // determinant is smaller than floating-point elimination can resolve.
    initialize_geogram();
    using Exact=GEO::expansion_nt;
    std::array<std::array<Exact,3>,N-1> v;
    std::array<Exact,N-1> rhs;
    for(std::size_t i=0;i<N-1;++i) {
        Exact norm(0.0);
        for(std::size_t d=0;d<3;++d) {
            v[i][d]=Exact(p[i+1][d])-Exact(p[0][d]);
            norm+=GEO::expansion_nt_square(v[i][d]);
        }
        rhs[i]=(norm+w[0]-w[i+1])*0.5;
    }
    if constexpr(N==3) {
        std::array<Exact,3> normal,q;
        Exact determinant(0.0);
        for(std::size_t d=0;d<3;++d) {
            const std::size_t j=(d+1)%3,k=(d+2)%3;
            normal[d]=v[0][j]*v[1][k]-v[0][k]*v[1][j];
            determinant+=GEO::expansion_nt_square(normal[d]);
            q[d]=rhs[0]*v[1][d]-rhs[1]*v[0][d];
        }
        if(determinant.sign()==GEO::ZERO) return false;
        const long double denominator=determinant.estimate();
        for(std::size_t d=0;d<3;++d) {
            const std::size_t j=(d+1)%3,k=(d+2)%3;
            offset[d]=static_cast<long double>(
                (q[j]*normal[k]-q[k]*normal[j]).estimate())/denominator;
        }
    } else {
        static_assert(N==4,"exact sphere construction expects a facet or tetrahedron");
        const Exact determinant=GEO::det3x3(
            v[0][0],v[0][1],v[0][2],
            v[1][0],v[1][1],v[1][2],
            v[2][0],v[2][1],v[2][2]);
        if(determinant.sign()==GEO::ZERO) return false;
        const long double denominator=determinant.estimate();
        for(std::size_t d=0;d<3;++d) {
            const Exact numerator=GEO::det3x3(
                d==0 ? rhs[0] : v[0][0],d==1 ? rhs[0] : v[0][1],d==2 ? rhs[0] : v[0][2],
                d==0 ? rhs[1] : v[1][0],d==1 ? rhs[1] : v[1][1],d==2 ? rhs[1] : v[1][2],
                d==0 ? rhs[2] : v[2][0],d==1 ? rhs[2] : v[2][1],d==2 ? rhs[2] : v[2][2]);
            offset[d]=static_cast<long double>(numerator.estimate())/denominator;
        }
    }
    return true;
}

template<std::size_t N, std::size_t D>
double sphere(const std::array<Point<D>,N>& p,const std::array<double,N>& w,
              Point<D>* center_out) {
    if constexpr(D==3 && N>1) {
        // Work in the vertex-relative frame and round only the constructed
        // center/radius. Normal equations square the condition number of skinny
        // tetrahedra; extended intermediates also retain small weight differences.
        std::array<std::array<long double,3>,N-1> v{};
        std::array<long double,N-1> rhs{},norms{};
        for(std::size_t i=0;i<N-1;++i) {
            long double norm=0.0L;
            for(std::size_t d=0;d<3;++d) {
                v[i][d]=static_cast<long double>(p[i+1][d])-p[0][d];
                norm+=v[i][d]*v[i][d];
            }
            norms[i]=norm;
            rhs[i]=0.5L*(norm+(static_cast<long double>(w[0])-w[i+1]));
        }
        std::array<long double,3> offset{};
        if constexpr(N==2) {
            if(norms[0]==0.0L) return std::numeric_limits<double>::infinity();
            const long double scale=rhs[0]/norms[0];
            for(std::size_t d=0;d<3;++d) offset[d]=scale*v[0][d];
        } else if constexpr(N==3) {
            // The cross-product norm does not subtract nearly equal Gram
            // products when the triangle is almost collinear.
            std::array<long double,3> normal{},q{};
            long double determinant=0.0L;
            for(std::size_t d=0;d<3;++d) {
                const std::size_t j=(d+1)%3,k=(d+2)%3;
                normal[d]=v[0][j]*v[1][k]-v[0][k]*v[1][j];
                determinant+=normal[d]*normal[d];
                q[d]=rhs[0]*v[1][d]-rhs[1]*v[0][d];
            }
            if(determinant<=1e-12L*norms[0]*norms[1]) {
                if(!sphere3_expansion_offset(p,w,offset))
                    return std::numeric_limits<double>::infinity();
            } else {
                for(std::size_t d=0;d<3;++d) {
                    const std::size_t j=(d+1)%3,k=(d+2)%3;
                    offset[d]=(q[j]*normal[k]-q[k]*normal[j])/determinant;
                }
            }
        } else {
            static_assert(N==4,"a 3D sphere has at most four defining points");
            // Solve the original power equations instead of a Gram system.
            long double matrix_scale=0.0L;
            for(const auto& row:v) for(long double value:row)
                matrix_scale=std::max(matrix_scale,std::abs(value));
            bool use_expansion=false;
            for(std::size_t column=0;column<3;++column) {
                std::size_t pivot=column;
                for(std::size_t row=column+1;row<3;++row)
                    if(std::abs(v[row][column])>std::abs(v[pivot][column])) pivot=row;
                // This is an accuracy filter, not a degeneracy tolerance.
                if(std::abs(v[pivot][column])<=1e-6L*matrix_scale) {
                    use_expansion=true;
                    break;
                }
                std::swap(v[column],v[pivot]);
                std::swap(rhs[column],rhs[pivot]);
                for(std::size_t row=column+1;row<3;++row) {
                    const long double factor=v[row][column]/v[column][column];
                    for(std::size_t d=column+1;d<3;++d)
                        v[row][d]-=factor*v[column][d];
                    rhs[row]-=factor*rhs[column];
                }
            }
            if(use_expansion) {
                if(!sphere3_expansion_offset(p,w,offset))
                    return std::numeric_limits<double>::infinity();
            } else {
                for(std::size_t row=3;row-- >0;) {
                    long double value=rhs[row];
                    for(std::size_t d=row+1;d<3;++d) value-=v[row][d]*offset[d];
                    offset[row]=value/v[row][row];
                }
            }
        }
        long double radius=-static_cast<long double>(w[0]);
        for(std::size_t d=0;d<3;++d) {
            if(center_out) (*center_out)[d]=static_cast<double>(p[0][d]+offset[d]);
            radius+=offset[d]*offset[d];
        }
        return static_cast<double>(radius);
    }
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
void own_alpha(Record<N,D>& record) {
    if(record.gabriel) {
        record.alpha=record.own_radius;
        record.tau_size=static_cast<unsigned char>(N);
    }
}

template<std::size_t N, std::size_t M, std::size_t D>
void take_coface(Record<N,D>& record,const Record<M,D>& coface) {
    if(coface.alpha<record.alpha) {
        record.alpha=coface.alpha;
        record.tau=coface.tau;
        record.tau_size=coface.tau_size;
    }
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
    explicit Complex(bool collect=true): collect_alpha(collect) {
        v.max_load_factor(0.8f); e.max_load_factor(0.8f); f.max_load_factor(0.8f);
    }
    Complex(const Complex&)=delete;
    Complex& operator=(const Complex&)=delete;
    Complex(Complex&&)=default;
    Complex& operator=(Complex&&)=default;
    void reserve(std::size_t vertices,std::size_t cells) {
        v.reserve(vertices); e.reserve(2*cells); f.reserve(cells);
        ve.reserve(4*cells); ef.reserve(3*cells);
    }
    bool collect_alpha;
    RecordMap<1,2> v;
    RecordMap<2,2> e;
    RecordMap<3,2> f;
    std::vector<std::pair<Record<1,2>*,Record<2,2>*>> ve;
    std::vector<std::pair<Record<2,2>*,Record<3,2>*>> ef;
};
template<> struct Complex<3> {
    explicit Complex(bool collect=true): collect_alpha(collect) {
        v.max_load_factor(0.8f); e.max_load_factor(0.8f);
        f.max_load_factor(0.8f); c.max_load_factor(0.8f);
    }
    Complex(const Complex&)=delete;
    Complex& operator=(const Complex&)=delete;
    Complex(Complex&&)=default;
    Complex& operator=(Complex&&)=default;
    void reserve(std::size_t vertices,std::size_t cells) {
        v.reserve(vertices); e.reserve(2*cells); f.reserve(2*cells); c.reserve(cells);
        ve.reserve(4*cells); ef.reserve(6*cells); fc.reserve(4*cells);
    }
    bool collect_alpha;
    RecordMap<1,3> v;
    RecordMap<2,3> e;
    RecordMap<3,3> f;
    RecordMap<4,3> c;
    std::vector<std::pair<Record<1,3>*,Record<2,3>*>> ve;
    std::vector<std::pair<Record<2,3>*,Record<3,3>*>> ef;
    std::vector<std::pair<Record<3,3>*,Record<4,3>*>> fc;
};

inline RecordInsertion<2,2> add_edge(
    Complex<2>& out,const std::array<Vertex<2>,2>& edge
) {
    auto result=add_record(out.e,edge,out.collect_alpha);
    if(!result.inserted) return result;
    for(int i=0;i<2;++i) {
        auto vertex=add_record(
            out.v,std::array<Vertex<2>,1>{edge[i]},out.collect_alpha
        );
        add_witness(vertex,edge[1-i]);
        out.ve.emplace_back(vertex.record,result.record);
    }
    return result;
}

inline void add_triangle(Complex<2>& out,const std::array<Vertex<2>,3>& cell) {
    auto face=add_record(out.f,cell,out.collect_alpha);
    if(!face.inserted) return;
    for(int skip=0;skip<3;++skip) {
        std::array<Vertex<2>,2> edge{cell[(skip+1)%3],cell[(skip+2)%3]};
        auto edge_record=add_edge(out,edge);
        add_witness(edge_record,cell[skip]);
        out.ef.emplace_back(edge_record.record,face.record);
    }
}

inline void add_triangle(Complex<3>& out,const std::array<Vertex<3>,3>& cell) {
    auto face=add_record(out.f,cell,out.collect_alpha);
    if(!face.inserted) return;
    for(int skip=0;skip<3;++skip) {
        std::array<Vertex<3>,2> edge{cell[(skip+1)%3],cell[(skip+2)%3]};
        auto edge_record=add_record(out.e,edge,out.collect_alpha);
        add_witness(edge_record,cell[skip]);
        out.ef.emplace_back(edge_record.record,face.record);
        if(!edge_record.inserted) continue;
        for(int i=0;i<2;++i) {
            auto vertex=add_record(
                out.v,std::array<Vertex<3>,1>{edge[i]},out.collect_alpha
            );
            add_witness(vertex,edge[1-i]);
            out.ve.emplace_back(vertex.record,edge_record.record);
        }
    }
}

inline Complex<3> lower_dimensional_complex(
    std::vector<Vertex<3>> vertices,bool weighted,bool collect_alpha=true
) {
    Complex<3> out(collect_alpha);
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
        add_record(
            out.v,std::array<Vertex<3>,1>{vertices.front()},out.collect_alpha
        );
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
        out.reserve(hull.size(),hull.size());
        for(const auto& site:hull)
            add_record(
                out.v,std::array<Vertex<3>,1>{site.vertex},out.collect_alpha
            );
        for(std::size_t i=1;i<hull.size();++i) {
            auto edge=add_record(
                out.e,std::array<Vertex<3>,2>{hull[i-1].vertex,hull[i].vertex},
                out.collect_alpha
            );
            for(int endpoint=0;endpoint<2;++endpoint) {
                const auto& vertex=endpoint==0 ? hull[i-1].vertex : hull[i].vertex;
                const auto& witness=endpoint==0 ? hull[i].vertex : hull[i-1].vertex;
                auto record=add_record(
                    out.v,std::array<Vertex<3>,1>{vertex},out.collect_alpha
                );
                add_witness(record,witness);
                out.ve.emplace_back(record.record,edge.record);
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
    triangulation->set_reorder(true);
    triangulation->set_vertices(vertices.size(),coordinates.data());
    out.reserve(vertices.size(),triangulation->nb_cells());
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

inline void add_tetrahedron(Complex<3>& out,const std::array<Vertex<3>,4>& cell) {
    auto tetrahedron=add_record(out.c,cell,out.collect_alpha);
    if(!tetrahedron.inserted) return;
    for(int omit=0;omit<4;++omit) {
        std::array<Vertex<3>,3> face{};
        int k=0;
        for(int i=0;i<4;++i) if(i!=omit) face[k++]=cell[i];
        auto face_record=add_record(out.f,face,out.collect_alpha);
        add_witness(face_record,cell[omit]);
        out.fc.emplace_back(face_record.record,tetrahedron.record);
        if(!face_record.inserted) continue;
        for(int skip=0;skip<3;++skip) {
            std::array<Vertex<3>,2> edge{
                face[(skip+1)%3],face[(skip+2)%3]
            };
            auto edge_record=add_record(out.e,edge,out.collect_alpha);
            add_witness(edge_record,face[skip]);
            out.ef.emplace_back(edge_record.record,face_record.record);
            if(!edge_record.inserted) continue;
            for(int i=0;i<2;++i) {
                auto vertex=add_record(
                    out.v,std::array<Vertex<3>,1>{edge[i]},out.collect_alpha
                );
                add_witness(vertex,edge[1-i]);
                out.ve.emplace_back(vertex.record,edge_record.record);
            }
        }
    }
}

inline void compute_alpha(Complex<2>& complex) {
    for(auto& item:complex.f) own_alpha(item.second);
    for(auto& item:complex.e) own_alpha(item.second);
    for(const auto& relation:complex.ef)
        take_coface(*relation.first,*relation.second);
    for(auto& item:complex.v) own_alpha(item.second);
    for(const auto& relation:complex.ve)
        take_coface(*relation.first,*relation.second);
}
inline void compute_alpha(Complex<3>& complex) {
    for(auto& item:complex.c) own_alpha(item.second);
    for(auto& item:complex.f) own_alpha(item.second);
    for(const auto& relation:complex.fc)
        take_coface(*relation.first,*relation.second);
    for(auto& item:complex.e) own_alpha(item.second);
    for(const auto& relation:complex.ef)
        take_coface(*relation.first,*relation.second);
    for(auto& item:complex.v) own_alpha(item.second);
    for(const auto& relation:complex.ve)
        take_coface(*relation.first,*relation.second);
}

inline Complex<2> lower_dimensional_complex(
    std::vector<Vertex<2>> vertices,bool collect_alpha=true
) {
    Complex<2> out(collect_alpha);
    if(vertices.empty()) return out;
    if(vertices.size()==1) {
        add_record(
            out.v,std::array<Vertex<2>,1>{vertices.front()},out.collect_alpha
        );
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
    out.reserve(vertices.size(),vertices.size());
    for(const auto& vertex:vertices)
        add_record(
            out.v,std::array<Vertex<2>,1>{vertex},out.collect_alpha
        );
    for(std::size_t i=1;i<vertices.size();++i)
        add_edge(out,std::array<Vertex<2>,2>{vertices[i-1],vertices[i]});
    return out;
}

template<class Points>
Complex<2> triangulate2(
    const Points& points,bool periodic=false,Point<2> from={0,0},
    Point<2> to={1,1},bool collect_alpha=true
) {
    initialize_geogram();
    auto base=unique_points<Points,2>(points,false,periodic?&from:nullptr);
    std::vector<Vertex<2>> vertices;
    if(periodic) {
        vertices.reserve(base.size()*9);
        for(int ox=-1;ox<=1;++ox) for(int oy=-1;oy<=1;++oy) for(auto v:base) {
            v.offset={ox,oy}; v.point[0]+=ox*(to[0]-from[0]); v.point[1]+=oy*(to[1]-from[1]); vertices.push_back(v);
        }
    } else vertices=base;
    Complex<2> out(collect_alpha);
    if(base.size()<3)
        return periodic
            ? std::move(out)
            : lower_dimensional_complex(std::move(base),collect_alpha);
    std::vector<double> coords; coords.reserve(vertices.size()*2);
    for(const auto& v:vertices) { coords.push_back(v.point[0]); coords.push_back(v.point[1]); }
    std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
    GEO::Numeric::random_reset();
    GEO::SmartPointer<GEO::Delaunay2d> dt=new GEO::Delaunay2d();
    dt->set_reorder(true);
    dt->set_vertices(vertices.size(),coords.data());
    const std::size_t expected_cells=periodic
        ? static_cast<std::size_t>(dt->nb_cells())/9+8
        : static_cast<std::size_t>(dt->nb_cells());
    out.reserve(base.size(),expected_cells);
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
    if(!periodic && out.f.empty())
        return lower_dimensional_complex(std::move(base),collect_alpha);
    return out;
}


inline Complex<3> triangulate3_regular(
    const std::vector<Vertex<3>>& vertices,bool weighted,
    const Point<3>* periodic_extent=nullptr,bool collect_alpha=true
) {
    std::lock_guard<std::mutex> lock(geogram_triangulation_mutex());
    GEO::Numeric::random_reset();
    Complex<3> out(collect_alpha);
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
    dt->set_reorder(true);
    dt->set_vertices(vertices.size(),coordinates.data());
    const std::size_t tile_count=periodic_extent ? 27 : 1;
    const std::size_t expected_cells=
        static_cast<std::size_t>(dt->nb_cells())/tile_count+8;
    out.reserve(vertices.size()/tile_count,expected_cells);
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
    const std::vector<Vertex<3>>& base,bool weighted,const Point<3>& period,
    bool collect_alpha=true
) {
    std::vector<Vertex<3>> vertices;
    vertices.reserve(base.size()*27);
    for(int ox=-1;ox<=1;++ox)
        for(int oy=-1;oy<=1;++oy)
            for(int oz=-1;oz<=1;++oz)
                for(auto vertex:base) {
                    vertex.offset={ox,oy,oz};
                    vertex.point[0]+=ox*period[0];
                    vertex.point[1]+=oy*period[1];
                    vertex.point[2]+=oz*period[2];
                    vertices.push_back(vertex);
                }
    return triangulate3_regular(vertices,weighted,&period,collect_alpha);
}

template<class Points>
Complex<3> triangulate3(
    const Points& points,bool weighted=false,bool periodic=false,
    Point<3> from={0,0,0},Point<3> to={1,1,1},
    bool preserve_lower_dimension=false,bool collect_alpha=true
) {
    initialize_geogram();
    auto vertices=unique_points<Points,3>(
        points,weighted,periodic ? &from : nullptr
    );
    if(!periodic && !spans_three_dimensions(vertices))
        return preserve_lower_dimension
            ? lower_dimensional_complex(
                std::move(vertices),weighted,collect_alpha
            )
            : Complex<3>(collect_alpha);
    if(vertices.size()<4) return Complex<3>(collect_alpha);
    if(!periodic)
        return triangulate3_regular(
            vertices,weighted,nullptr,collect_alpha
        );
    Point<3> period{to[0]-from[0],to[1]-from[1],to[2]-from[2]};
    return triangulate3_tiled(vertices,weighted,period,collect_alpha);
}

template<std::size_t N,std::size_t D>
bool has_nonzero_offsets(const RecordMap<N,D>& records) {
    for(const auto& item:records)
        for(const auto& offset:item.first.offsets)
            for(int value:offset)
                if(value!=0) return true;
    return false;
}

template<std::size_t N,std::size_t D>
void validate_one_sheeted(const RecordMap<N,D>& records) {
    if(!has_nonzero_offsets(records)) return;
    std::unordered_set<std::array<unsigned,N>,VertexKeyHash<N>> seen;
    seen.reserve(records.size());
    for(const auto& item:records)
        if(!seen.insert(item.first.vertices).second)
            throw std::runtime_error(
                "Cannot convert periodic triangulation to a one-sheeted covering"
            );
}

template<std::size_t N,std::size_t D,class Callback>
void emit_values_plain(const RecordMap<N,D>& records,const Callback& callback) {
    validate_one_sheeted(records);
    for(const auto& item:records)
        if(std::isfinite(item.second.alpha))
            callback(item.first.vertices,item.second.alpha);
}

template<std::size_t N,std::size_t D,class Callback>
void emit_values_attachment(
    const RecordMap<N,D>& records,const Callback& callback
) {
    validate_one_sheeted(records);
    for(const auto& item:records) {
        if(std::isfinite(item.second.alpha)) {
            const auto& record=item.second;
            std::vector<unsigned> tau(
                record.tau.begin(),record.tau.begin()+record.tau_size
            );
            callback(item.first.vertices,record.alpha,tau);
        }
    }
}

template<std::size_t N,std::size_t D,class Callback>
void emit_combinatorics(
    const RecordMap<N,D>& records,const Callback& callback
) {
    validate_one_sheeted(records);
    for(const auto& item:records) callback(item.first.vertices);
}

template<std::size_t N,std::size_t D,class Callback>
void emit_lifts(const RecordMap<N,D>& records,const Callback& callback) {
    using Offsets=std::array<std::array<int,D>,N>;
    std::unordered_map<std::array<unsigned,N>,Offsets,VertexKeyHash<N>> seen;
    seen.reserve(records.size());
    for(const auto& item:records) {
        auto [iterator,inserted]=seen.emplace(
            item.first.vertices,item.first.offsets
        );
        if(!inserted && iterator->second!=item.first.offsets)
            throw std::runtime_error(
                "Cannot convert periodic triangulation to a one-sheeted covering"
            );
    }
    for(const auto& item:seen) callback(item.first,item.second);
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

#include "compact_alpha3.hpp"
#include "compact_periodic3.hpp"

} // namespace detail

template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes(const Points& p,const SimplexCallback& cb) { fill_alpha_shapes_direct(p,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes_direct(const Points& p,const SimplexCallback& cb) {
    if(detail::try_compact_alpha3(p,cb,false)) return;
    auto x=detail::triangulate3(p); detail::emit_alpha(x,cb);
}
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes_with_attachment(const Points& p,const SimplexCallback& cb) { fill_alpha_shapes_direct_with_attachment(p,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_alpha_shapes_direct_with_attachment(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p); detail::emit_attachment(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_alpha_shapes(const Points& p,const SimplexCallback& cb) { fill_weighted_alpha_shapes_direct(p,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_alpha_shapes_direct(const Points& p,const SimplexCallback& cb) {
    if(detail::try_compact_alpha3(p,cb,true)) return;
    auto x=detail::triangulate3(p,true); detail::emit_alpha(x,cb);
}
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_alpha_shapes(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { fill_periodic_alpha_shapes_direct(p,cb,a,b); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_alpha_shapes_direct(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) {
    if(detail::try_compact_periodic3(p,cb,a,b)) return;
    auto x=detail::triangulate3(p,false,true,a,b); detail::emit_alpha(x,cb);
}
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_delaunay(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p,false,false,{0,0,0},{1,1,1},true,false); detail::emit_delaunay(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_delaunay(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate3(p,true,false,{0,0,0},{1,1,1},true,false); detail::emit_delaunay(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_delaunay(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,false,true,a,b,false,false); detail::emit_delaunay(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_periodic_delaunay_lifts(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,false,true,a,b,false,false); detail::emit_periodic_lifts(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_periodic_alpha_shapes(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { fill_weighted_periodic_alpha_shapes_direct(p,cb,a,b); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_periodic_alpha_shapes_direct(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,true,true,a,b); detail::emit_alpha(x,cb); }
template<bool exact> template<class Points,class SimplexCallback>
void AlphaShapes<exact>::fill_weighted_periodic_delaunay(const Points& p,const SimplexCallback& cb,std::array<double,3> a,std::array<double,3> b) { auto x=detail::triangulate3(p,true,true,a,b,false,false); detail::emit_delaunay(x,cb); }

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
void fill_delaunay2d(const Points& p,const SimplexCallback& cb) { auto x=detail::triangulate2(p,false,{0,0},{1,1},false); detail::emit_delaunay(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_periodic_delaunay2d(const Points& p,const SimplexCallback& cb,std::array<double,2> a,std::array<double,2> b) { auto x=detail::triangulate2(p,true,a,b,false); detail::emit_delaunay(x,cb); }
template<bool exact,class Points,class SimplexCallback>
void fill_periodic_delaunay2d_lifts(const Points& p,const SimplexCallback& cb,std::array<double,2> a,std::array<double,2> b) { auto x=detail::triangulate2(p,true,a,b,false); detail::emit_periodic_lifts(x,cb); }

} // namespace diode
