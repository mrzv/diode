set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module         (diode diode.cpp)
target_link_libraries       (diode PRIVATE diode_backend)

install(TARGETS diode DESTINATION .)
