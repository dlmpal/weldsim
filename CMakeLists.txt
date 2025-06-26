cmake_minimum_required(VERSION 3.18)
project(WELDSIM VERSION "0.0.1")

# Build type
if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

# Physical dimension
set(spacedim 3)

# Sources
add_executable(weldsim_${spacedim}d)
add_subdirectory(src)

# AMReX
find_package(AMReX REQUIRED ${spacedim}D AMRLEVEL LSOLVERS)
target_link_libraries(weldsim_${spacedim}d PRIVATE AMReX::amrex_${spacedim}d)

# Release flags
if(CMAKE_BUILD_TYPE MATCHES RELEASE)
    list(APPEND cxx_release_flags -O3)
endif()

# Debug flags
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    list(APPEND  cxx_debug_flags -g;-Wall;-Wextra;-Wpedantic;-Wconversion)
endif()

# Set flags (depending on build type)
target_compile_options(weldsim_${spacedim}d PRIVATE
                       $<$<CONFIG:RELEASE>:${cxx_release_flags}>
                       $<$<CONFIG:DEBUG>:${cxx_debug_flags}>)

# Installation
install(FILES
       ${CMAKE_CURRENT_BINARY_DIR}/weldsim_${spacedim}d
       DESTINATION ${PROJECT_SOURCE_DIR}/bin/)