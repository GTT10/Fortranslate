cmake_minimum_required(VERSION 3.23)

project(
  PeleF
  VERSION 0.155.0
  DESCRIPTION "Modern Fortran reimplementation of selected PeleC algorithms"
  LANGUAGES Fortran
)

option(PELEF_ENABLE_TESTS "Build PeleF tests" ON)
option(PELEF_ENABLE_CANTERA_REFERENCE "Enable live Cantera parity tests" OFF)
option(PELEF_ENABLE_MPI "Build MPI verification targets" OFF)

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
  src/eb/eb_geometry_2d_mod.F90
  src/physics/eos_ideal_mod.F90
  src/physics/nasa7_thermo_mod.F90
  src/physics/thermo_database_mod.F90
  src/physics/h2o2_full_thermo_mod.F90
  src/physics/mixture_thermo_mod.F90
  src/transport/transport_database_mod.F90
  src/transport/mixture_transport_mod.F90
  src/driver/simulation_config_reactive_1d_mod.F90
  src/driver/simulation_config_reactive_2d_mod.F90
  src/driver/simulation_config_reactive_eb_2d_mod.F90
  src/driver/simulation_config_reactive_eb_amr_2d_mod.F90
  src/hydro/reconstruction_weno_mod.F90
  src/reactive/reactive_1d_mod.F90
  src/eb/reactive_eb_cfl_2d_mod.F90
  src/eb/eb_reactive_wall_flux_2d_mod.F90
  src/eb/eb_reactive_redistribution_2d_mod.F90
  src/amr/amr_reactive_1d_mod.F90
  src/amr/amr_eb_hierarchy_2d_mod.F90
  src/amr/amr_eb_patch_tree_2d_mod.F90
  src/amr/amr_eb_multilevel_2d_mod.F90
  src/amr/amr_eb_flux_register_2d_mod.F90
  src/amr/amr_eb_reactive_2d_mod.F90
  src/amr/amr_eb_patch_tree_reactive_2d_mod.F90
  src/amr/amr_eb_multilevel_reactive_2d_mod.F90
  src/amr/amr_eb_regrid_2d_mod.F90
  src/boundary/reactive_boundary_2d_mod.F90
  src/transport/reactive_transport_2d_mod.F90
  src/reactive/reactive_2d_mod.F90
  src/eb/eb_reactive_reconstruction_2d_mod.F90
  src/eb/eb_reactive_hydro_2d_mod.F90
  src/eb/eb_reactive_transport_2d_mod.F90
  src/amr/amr_eb_transport_2d_mod.F90
  src/amr/amr_eb_multilevel_transport_2d_mod.F90
  src/amr/amr_eb_multipatch_transport_2d_mod.F90
  src/driver/reactive_eb_2d_driver_mod.F90
  src/driver/reactive_eb_amr_2d_driver_mod.F90
  src/chemistry/elementary_kinetics_mod.F90
  src/generated/h2o2_elementary_mechanism_mod.F90
  src/generated/h2o2_full_mechanism_mod.F90
  src/chemistry/constant_volume_reactor_mod.F90
  src/chemistry/isomerization_reactor_mod.F90
  src/hydro/state_conversion_mod.F90
  src/hydro/boundary_conditions_mod.F90
  src/hydro/boundary_conditions_multispecies_mod.F90
  src/hydro/slope_limiter_mod.F90
  src/amr/amr_hierarchy_1d_mod.F90
  src/amr/amr_multipatch_1d_mod.F90
  src/amr/amr_patch_tree_1d_mod.F90
  src/amr/amr_patch_tree_reactive_1d_mod.F90
  src/amr/amr_regrid_1d_mod.F90
  src/amr/amr_multilevel_reactive_1d_mod.F90
  src/amr/amr_multipatch_reactive_1d_mod.F90
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
  src/driver/simulation_config_h2o2_full_mod.F90
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

add_executable(pelef0d_h2o2_full app/pelef0d_h2o2_full.F90)
target_link_libraries(pelef0d_h2o2_full PRIVATE pelef_core)

add_executable(pelef_reactive_1d app/pelef_reactive_1d.F90)
target_link_libraries(pelef_reactive_1d PRIVATE pelef_core)

add_executable(pelef_amr_reactive_1d app/pelef_amr_reactive_1d.F90)
target_link_libraries(pelef_amr_reactive_1d PRIVATE pelef_core)

add_executable(pelef_reactive_2d app/pelef_reactive_2d.F90)
target_link_libraries(pelef_reactive_2d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_2d app/pelef_reactive_eb_2d.F90)
target_link_libraries(pelef_reactive_eb_2d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_amr_2d app/pelef_reactive_eb_amr_2d.F90)
target_link_libraries(pelef_reactive_eb_amr_2d PRIVATE pelef_core)

add_executable(pelef_transport_probe app/pelef_transport_probe.F90)
target_link_libraries(pelef_transport_probe PRIVATE pelef_core)

install(
  TARGETS pelef pelef2d pelef_ms pelef0d pelef0d_h2o2 pelef_reactive_1d
    pelef_amr_reactive_1d
    pelef0d_h2o2_full pelef_reactive_2d pelef_reactive_eb_2d
    pelef_reactive_eb_amr_2d
    pelef_transport_probe
  RUNTIME DESTINATION bin
)

if(PELEF_ENABLE_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()


if(PELEF_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  add_library(
    pelef_mpi_support
    src/parallel/mpi_domain_1d_mod.F90
    src/parallel/mpi_amr_patch_1d_mod.F90
    src/parallel/mpi_amr_sparse_patch_1d_mod.F90
    src/parallel/mpi_amr_eb_patch_2d_mod.F90
    src/parallel/mpi_amr_eb_io_2d_mod.F90
    src/parallel/mpi_reactive_transport_1d_mod.F90
    src/parallel/mpi_reactive_1d_mod.F90
  )
  target_link_libraries(pelef_mpi_support PUBLIC pelef_core MPI::MPI_Fortran)
  set_target_properties(pelef_mpi_support PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mpi_modules")
  target_include_directories(pelef_mpi_support PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mpi_modules>")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(pelef_mpi_support PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(pelef_mpi_support PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
  endif()
  add_executable(pelef_mpi_1d app/pelef_mpi_1d.F90)
  target_link_libraries(pelef_mpi_1d PRIVATE pelef_mpi_support)
  add_executable(pelef_mpi_multispecies_1d app/pelef_mpi_multispecies_1d.F90)
  target_link_libraries(pelef_mpi_multispecies_1d PRIVATE pelef_mpi_support)
  add_executable(pelef_mpi_h2o2_batch app/pelef_mpi_h2o2_batch.F90)
  target_link_libraries(pelef_mpi_h2o2_batch PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_reactive_transport_1d
    app/pelef_mpi_reactive_transport_1d.F90
  )
  target_link_libraries(
    pelef_mpi_reactive_transport_1d PRIVATE pelef_mpi_support
  )
  add_executable(pelef_mpi_reactive_1d app/pelef_mpi_reactive_1d.F90)
  target_link_libraries(pelef_mpi_reactive_1d PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_amr_patch_1d
    app/pelef_mpi_amr_patch_1d.F90
  )
  target_link_libraries(pelef_mpi_amr_patch_1d PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_eb_amr_patch_2d
    app/pelef_mpi_eb_amr_patch_2d.F90
  )
  target_link_libraries(
    pelef_mpi_eb_amr_patch_2d PRIVATE pelef_mpi_support
  )
  add_executable(
    pelef_mpi_amr_reactive_1d
    app/pelef_mpi_amr_reactive_1d.F90
  )
  target_link_libraries(pelef_mpi_amr_reactive_1d PRIVATE pelef_mpi_support)
  install(
    TARGETS
      pelef_mpi_1d
      pelef_mpi_multispecies_1d
      pelef_mpi_h2o2_batch
      pelef_mpi_reactive_transport_1d
      pelef_mpi_reactive_1d
      pelef_mpi_amr_patch_1d
      pelef_mpi_eb_amr_patch_2d
      pelef_mpi_amr_reactive_1d
    RUNTIME DESTINATION bin
  )
endif()
