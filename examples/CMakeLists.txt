add_executable          (generate_alpha_shape   generate_alpha_shape.cpp)
target_link_libraries   (generate_alpha_shape PRIVATE diode_backend)

add_executable          (generate_weighted_alpha_shape   generate_weighted_alpha_shape.cpp)
target_link_libraries   (generate_weighted_alpha_shape PRIVATE diode_backend)

add_executable          (generate_periodic_alpha_shape   generate_periodic_alpha_shape.cpp)
target_link_libraries   (generate_periodic_alpha_shape PRIVATE diode_backend)
