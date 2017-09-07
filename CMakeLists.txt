project                     (DioDe)
cmake_minimum_required      (VERSION 2.8)

# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories         (include)

find_package                (CGAL)
include                     (${CGAL_USE_FILE})
add_definitions             (${CGAL_CXX_FLAGS_INIT})
include_directories         (${CGAL_INCLUDE_DIRS})

add_subdirectory            (examples)
