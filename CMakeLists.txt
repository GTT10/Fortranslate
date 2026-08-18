cmake_minimum_required(VERSION 3.23)

project(
  PeleF
  VERSION 0.9.0
  DESCRIPTION "Modern Fortran reimplementation of selected PeleC algorithms"
  LANGUAGES Fortran
)

option(PELEF_ENABLE_TESTS "Build PeleF tests" ON)
option(PELEF_ENABLE_CANTERA_REFERENCE "Enable live Cantera parity tests" OFF)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(
  pelef_core
  src/core/precision_mod.F90
  src/core/constants_mod.F90
  src/core/state_indices_mod.F90
  src/core/multispecies_state_mod.F90
  src/core/mesh_mod.F90
  src/core/mesh_2d_mod.F90
  src/physics/eos_ideal_mod.F90
  src/physics/nasa7_thermo_mod.F90
  src/physics/thermo_database_mod.F90
  src/physics/mixture_thermo_mod.F90
  src/chemistry/elementary_kinetics_mod.F90
  src/generated/h2o2_elementary_mechanism_mod.F90
  src/chemistry/constant_volume_reactor_mod.F90
  src/chemistry/isomerization_reactor_mod.F90
  src/hydro/state_conversion_mod.F90
  src/hydro/boundary_conditions_mod.F90
  src/hydro/boundary_conditions_multispecies_mod.F90
  src/hydro/slope_limiter_mod.F90
  src/hydro/reconstruction_plm_mod.F90
  src/hydro/reconstruction_pelec_plm_mod.F90
  src/hydro/riemann_rusanov_mod.F90
  src/hydro/riemann_pelec_mod.F90
  src/hydro/riemann_flux_mod.F90
  src/hydro/directional_flux_mod.F90
  src/hydro/multispecies_flux_mod.F90
  src/hydro/ctu_2d_mod.F90
  src/hydro/ctu_multispecies_2d_mod.F90
  src/hydro/reconstruction_multispecies_mod.F90
  src/hydro/finite_volume_multispecies_mod.F90
  src/hydro/finite_volume_mod.F90
  src/driver/time_integrator_mod.F90
  src/driver/time_integrator_multispecies_mod.F90
  src/driver/diagnostics_mod.F90
  src/driver/diagnostics_2d_mod.F90
  src/driver/simulation_config_mod.F90
  src/driver/simulation_config_2d_mod.F90
  src/driver/simulation_config_multispecies_mod.F90
  src/driver/simulation_config_reactor_mod.F90
  src/driver/simulation_config_h2o2_reactor_mod.F90
  src/problems/sod_problem_mod.F90
  src/problems/shu_osher_problem_mod.F90
  src/problems/sedov_problem_mod.F90
  src/problems/isentropic_vortex_problem_mod.F90
  src/problems/multispec_sod_problem_mod.F90
  src/io/csv_io_mod.F90
  src/io/csv_io_2d_mod.F90
  src/io/csv_io_multispecies_mod.F90
)

set_target_properties(
  pelef_core
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

target_include_directories(
  pelef_core
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>"
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(
    pelef_core
    PRIVATE
      -Wall
      -Wextra
      -Wpedantic
      -Wimplicit-interface
      -Wconversion-extra
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(
      pelef_core
      PRIVATE
        -fcheck=all
        -fbacktrace
        -ffpe-trap=invalid,zero,overflow
    )
  endif()
endif()

add_executable(pelef app/pelef.F90)
target_link_libraries(pelef PRIVATE pelef_core)

add_executable(pelef2d app/pelef2d.F90)
target_link_libraries(pelef2d PRIVATE pelef_core)

add_executable(pelef_ms app/pelef_ms.F90)
target_link_libraries(pelef_ms PRIVATE pelef_core)

add_executable(pelef0d app/pelef0d.F90)
target_link_libraries(pelef0d PRIVATE pelef_core)

add_executable(pelef0d_h2o2 app/pelef0d_h2o2.F90)
target_link_libraries(pelef0d_h2o2 PRIVATE pelef_core)

install(
  TARGETS pelef pelef2d pelef_ms pelef0d pelef0d_h2o2
  RUNTIME DESTINATION bin
)

if(PELEF_ENABLE_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()
