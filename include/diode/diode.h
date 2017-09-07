#pragma once

#include <array>
#include <tuple>
#include <vector>

#include <boost/range/adaptor/map.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Alpha_shape_3.h>

namespace diode
{

template<class Points, class SimplexCallback>
void fill_alpha_shapes(const Points& points, const SimplexCallback& add_simplex);

}

#include "diode.hpp"
