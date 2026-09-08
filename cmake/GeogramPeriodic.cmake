# Geogram 1.9.9 treats hidden weighted sites as fatal, although its boundary
# classification already skips sites without incident tetrahedra. Continue
# insertion for these sites just as the nonperiodic regular triangulation does.
# Keep the patch pinned and fail configuration if upstream context changes.
set(_periodic_source "${geogram_SOURCE_DIR}/src/lib/geogram/delaunay/periodic_delaunay_3d.cpp")
file(READ "${_periodic_source}" _periodic_code)
set(_hidden_before [=[                has_empty_cells_ = true;
                return true;]=])
set(_hidden_after [=[                // DioDe: a hidden weighted site is a successful insertion.
                return true;]=])
set(_boundary_before [=[	// Test for empty cells
	// TODO: I tested, it really occurs that there is still an empty
	// cell here, why is it not detected before ? To be understood.
        update_v_to_cell();
        for(index_t v=0; v<nb_vertices_non_periodic_; ++v) {
            if(v_to_cell_[v] == NO_INDEX) {
                has_empty_cells_ = true;
                return;
            }
        }]=])
set(_boundary_after [=[        // DioDe: hidden sites have no incident cells and generate no copies.
        update_v_to_cell();]=])
foreach(_part hidden boundary)
    string(FIND "${_periodic_code}" "${_${_part}_before}" _before_position)
    string(FIND "${_periodic_code}" "${_${_part}_after}" _after_position)
    if(NOT _before_position EQUAL -1)
        string(REPLACE "${_${_part}_before}" "${_${_part}_after}" _periodic_code "${_periodic_code}")
    elseif(_after_position EQUAL -1)
        message(FATAL_ERROR "Unexpected Geogram periodic source: cannot apply ${_part} correction")
    endif()
endforeach()
file(READ "${_periodic_source}" _periodic_original)
if(NOT _periodic_code STREQUAL _periodic_original)
    file(WRITE "${_periodic_source}" "${_periodic_code}")
endif()
