add_executable          (generate_alpha_shape   generate_alpha_shape.cpp)
target_link_libraries   (generate_alpha_shape   ${CGAL_LIBRARY} ${GMP_LIBRARIES} ${MPFR_LIBRARIES})
