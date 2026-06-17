set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module         (diode diode.cpp)
target_link_libraries       (diode PUBLIC CGAL::CGAL)

# Optional: enable the CellWithValue<Simplex> capsule exporter
# (diode.fill_alpha_shapes_cells) by pointing at Oineus's include directory.
# The NumPy exporter (fill_alpha_shapes_arrays) is always built and needs nothing.
set(OINEUS_INCLUDE_DIR "" CACHE PATH "Path to Oineus include/ to enable the CellWithValue capsule export")
if(OINEUS_INCLUDE_DIR)
    target_include_directories(diode PRIVATE "${OINEUS_INCLUDE_DIR}")
    target_compile_definitions(diode PRIVATE DIODE_HAVE_OINEUS)
    message(STATUS "DioDe: Oineus cell capsule export ENABLED (OINEUS_INCLUDE_DIR=${OINEUS_INCLUDE_DIR})")
else()
    message(STATUS "DioDe: Oineus cell capsule export disabled (set OINEUS_INCLUDE_DIR to enable)")
endif()

install(TARGETS diode DESTINATION .)
