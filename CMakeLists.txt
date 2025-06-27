cmake_minimum_required(VERSION 3.18)
project(WELDSIM VERSION "0.0.1")

# Build type
if (CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

# Sources
add_executable(weldsim_${SPACEDIM}d)
add_subdirectory(src)

# AMReX
find_package(AMReX REQUIRED ${SPACEDIM}D AMRLEVEL LSOLVERS)
target_link_libraries(weldsim_${SPACEDIM}d PRIVATE AMReX::amrex_${SPACEDIM}d)

# Release flags
if(CMAKE_BUILD_TYPE MATCHES RELEASE)
    list(APPEND cxx_release_flags -O3)
endif()

# Debug flags
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    list(APPEND  cxx_debug_flags -g;-Wall;-Wextra;-Wpedantic;-Wconversion)
endif()

# Set flags (depending on build type)
target_compile_options(weldsim_${SPACEDIM}d PRIVATE
                       $<$<CONFIG:RELEASE>:${cxx_release_flags}>
                       $<$<CONFIG:DEBUG>:${cxx_debug_flags}>)

# Installation
install(FILES
       ${CMAKE_CURRENT_BINARY_DIR}/weldsim_${SPACEDIM}d
       DESTINATION ${PROJECT_SOURCE_DIR}/bin/
       PERMISSIONS OWNER_READ OWNER_EXECUTE)