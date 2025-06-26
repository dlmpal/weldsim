#=============================================================================
target_sources(weldsim_${spacedim}d PRIVATE 
               ${CMAKE_CURRENT_SOURCE_DIR}/params.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/material_properties.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/heat_source.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/bc.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/advance_implicit.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/advance_explicit.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/single_level.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/thermal_solver.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/level_build.cpp
               ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)