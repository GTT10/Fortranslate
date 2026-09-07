cmake_minimum_required(VERSION 3.23)

project(
  PeleF
  VERSION 0.245.0
  DESCRIPTION "Modern Fortran reimplementation of selected PeleC algorithms"
  LANGUAGES Fortran
)

option(PELEF_ENABLE_TESTS "Build PeleF tests" ON)
option(PELEF_ENABLE_CANTERA_REFERENCE "Enable live Cantera parity tests" OFF)
option(PELEF_ENABLE_MPI "Build MPI verification targets" OFF)
option(
  PELEF_ENABLE_SUNDIALS
  "Build the optional official-Fortran CVODE selected-reactor backend"
  OFF
)
set(
  PELEF_MECHANISM_BUNDLE "" CACHE FILEPATH
  "Optional normalized JSON mechanism bundle to compile and probe"
)
set(
  PELEF_MECHANISM_SOURCE "" CACHE FILEPATH
  "Optional source YAML whose SHA-256 must match the selected bundle"
)
if(PELEF_MECHANISM_SOURCE AND NOT PELEF_MECHANISM_BUNDLE)
  message(FATAL_ERROR "PELEF_MECHANISM_SOURCE requires PELEF_MECHANISM_BUNDLE")
endif()

if(PELEF_ENABLE_TESTS OR PELEF_MECHANISM_BUNDLE)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

include(cmake/RequireMpiF08.cmake)
include(cmake/MechanismBundle.cmake)
include(cmake/RequireSundialsDoublePrecision.cmake)
if(PELEF_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  pelef_require_mpi_f08()
endif()
if(PELEF_ENABLE_SUNDIALS)
  find_package(
    SUNDIALS 7.2.0 EXACT CONFIG REQUIRED
    COMPONENTS
      fcvode_mod
      fnvecserial_mod
      fsunmatrixdense_mod
      fsunlinsoldense_mod
  )
  foreach(
      sundials_target
      SUNDIALS::fcvode_mod_static
      SUNDIALS::fnvecserial_mod_static
      SUNDIALS::fsunmatrixdense_mod_static
      SUNDIALS::fsunlinsoldense_mod_static
    )
    if(NOT TARGET "${sundials_target}")
      message(
        FATAL_ERROR
        "SUNDIALS official Fortran target is unavailable: ${sundials_target}"
      )
    endif()
  endforeach()
  get_target_property(
    sundials_include_directories
    SUNDIALS::fcvode_mod_static
    INTERFACE_INCLUDE_DIRECTORIES
  )
  set(sundials_config_header "")
  foreach(sundials_include_directory IN LISTS sundials_include_directories)
    if(
        EXISTS
          "${sundials_include_directory}/sundials/sundials_config.h"
      )
      set(
        sundials_config_header
        "${sundials_include_directory}/sundials/sundials_config.h"
      )
      break()
    endif()
  endforeach()
  if(NOT sundials_config_header)
    message(
      FATAL_ERROR
      "SUNDIALS official precision header is unavailable from "
      "SUNDIALS::fcvode_mod_static"
    )
  endif()
  pelef_require_sundials_double_precision("${sundials_config_header}")
endif()
if(PELEF_ENABLE_TESTS)
  include(CTest)
  enable_testing()
endif()

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
  src/core/mesh_3d_mod.F90
  src/eb/eb_geometry_2d_mod.F90
  src/eb/eb_geometry_3d_mod.F90
  src/physics/eos_ideal_mod.F90
  src/physics/nasa7_thermo_mod.F90
  src/physics/thermo_database_mod.F90
  src/physics/h2o2_full_thermo_mod.F90
  src/physics/mixture_thermo_mod.F90
  src/transport/gas_transport_mod.F90
  src/transport/transport_database_mod.F90
  src/transport/mixture_transport_mod.F90
  src/driver/selected_composition_mod.F90
  src/driver/simulation_config_reactive_1d_mod.F90
  src/driver/simulation_config_reactive_2d_mod.F90
  src/driver/simulation_config_reactive_3d_mod.F90
  src/driver/simulation_config_amr_reactive_3d_mod.F90
  src/driver/simulation_config_reactive_eb_3d_mod.F90
  src/driver/simulation_config_reactive_eb_2d_mod.F90
  src/driver/simulation_config_reactive_eb_amr_2d_mod.F90
  src/hydro/reconstruction_weno_mod.F90
  src/reactive/reactive_1d_mod.F90
  src/reactive/reactive_directional_flux_3d_mod.F90
  src/eb/reactive_eb_cfl_2d_mod.F90
  src/eb/eb_reactive_wall_flux_2d_mod.F90
  src/eb/eb_reactive_wall_flux_3d_mod.F90
  src/eb/reactive_eb_cfl_3d_mod.F90
  src/eb/eb_reactive_redistribution_2d_mod.F90
  src/eb/eb_reactive_redistribution_3d_mod.F90
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
  src/transport/reactive_transport_3d_mod.F90
  src/eb/eb_reactive_transport_3d_mod.F90
  src/reactive/reactive_2d_mod.F90
  src/driver/reactive_2d_application_mod.F90
  src/reactive/reactive_3d_mod.F90
  src/eb/eb_reactive_hydro_3d_mod.F90
  src/amr/amr_hierarchy_3d_mod.F90
  src/amr/amr_reactive_3d_mod.F90
  src/amr/amr_reactive_transport_3d_mod.F90
  src/io/amr_reactive_3d_checkpoint_mod.F90
  src/driver/amr_reactive_3d_application_mod.F90
  src/io/reactive_eb_3d_checkpoint_mod.F90
  src/eb/eb_reactive_reconstruction_2d_mod.F90
  src/eb/eb_reactive_hydro_2d_mod.F90
  src/eb/eb_reactive_transport_2d_mod.F90
  src/amr/amr_eb_transport_2d_mod.F90
  src/amr/amr_eb_multilevel_transport_2d_mod.F90
  src/amr/amr_eb_multipatch_transport_2d_mod.F90
  src/driver/reactive_eb_2d_driver_mod.F90
  src/driver/reactive_eb_2d_application_mod.F90
  src/driver/reactive_eb_amr_2d_driver_mod.F90
  src/driver/reactive_eb_amr_2d_application_mod.F90
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
  src/driver/amr_reactive_1d_application_mod.F90
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
  src/hydro/finite_volume_3d_mod.F90
  src/driver/time_integrator_mod.F90
  src/driver/time_integrator_multispecies_mod.F90
  src/driver/diagnostics_mod.F90
  src/driver/diagnostics_2d_mod.F90
  src/driver/diagnostics_3d_mod.F90
  src/driver/reactive_eb_3d_driver_mod.F90
  src/driver/reactive_eb_3d_application_mod.F90
  src/driver/simulation_config_mod.F90
  src/driver/simulation_config_2d_mod.F90
  src/driver/simulation_config_3d_mod.F90
  src/driver/simulation_config_multispecies_mod.F90
  src/driver/simulation_config_reactor_mod.F90
  src/driver/simulation_config_h2o2_reactor_mod.F90
  src/driver/simulation_config_h2o2_full_mod.F90
  src/problems/sod_problem_mod.F90
  src/problems/shu_osher_problem_mod.F90
  src/problems/sedov_problem_mod.F90
  src/problems/isentropic_vortex_problem_mod.F90
  src/problems/entropy_wave_3d_problem_mod.F90
  src/problems/reactive_entropy_wave_3d_problem_mod.F90
  src/problems/multispec_sod_problem_mod.F90
  src/io/csv_io_mod.F90
  src/io/csv_io_2d_mod.F90
  src/io/csv_io_3d_mod.F90
  src/io/reactive_csv_io_3d_mod.F90
  src/driver/reactive_3d_application_mod.F90
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

if(PELEF_ENABLE_TESTS OR PELEF_MECHANISM_BUNDLE)
  add_library(
    pelef_selected_runtime_core
    STATIC
    src/core/precision_mod.F90
    src/core/constants_mod.F90
    src/physics/nasa7_thermo_mod.F90
    src/physics/mixture_thermo_mod.F90
    src/chemistry/elementary_kinetics_mod.F90
    src/chemistry/constant_volume_reactor_mod.F90
    src/driver/selected_composition_mod.F90
    src/transport/gas_transport_mod.F90
    src/driver/selected_mechanism_runtime_mod.F90
    src/driver/simulation_config_selected_reactor_mod.F90
  )
  set_target_properties(
    pelef_selected_runtime_core
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_runtime_modules"
  )
  target_include_directories(
    pelef_selected_runtime_core
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_runtime_modules>"
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_runtime_core
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_runtime_core
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  add_library(
    pelef_selected_reactive_1d_runtime
    STATIC
    src/core/state_indices_mod.F90
    src/transport/mixture_transport_mod.F90
    src/driver/simulation_config_reactive_1d_mod.F90
    src/hydro/slope_limiter_mod.F90
    src/hydro/reconstruction_weno_mod.F90
    src/reactive/reactive_1d_mod.F90
  )
  set_target_properties(
    pelef_selected_reactive_1d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_runtime_modules"
  )
  target_include_directories(
    pelef_selected_reactive_1d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_runtime_modules>"
  )
  target_link_libraries(
    pelef_selected_reactive_1d_runtime PUBLIC pelef_selected_runtime_core
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_reactive_1d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_reactive_1d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  if(PELEF_ENABLE_MPI)
    add_library(
      pelef_selected_mpi_reactive_1d_runtime
      STATIC
      src/parallel/mpi_domain_1d_mod.F90
      src/parallel/mpi_reactive_transport_1d_mod.F90
      src/parallel/mpi_reactive_1d_mod.F90
      src/driver/mpi_reactive_1d_application_mod.F90
    )
    set_target_properties(
      pelef_selected_mpi_reactive_1d_runtime
      PROPERTIES
        Fortran_MODULE_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_1d_modules"
    )
    target_include_directories(
      pelef_selected_mpi_reactive_1d_runtime
      PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_1d_modules>"
    )
    target_link_libraries(
      pelef_selected_mpi_reactive_1d_runtime
      PUBLIC pelef_selected_reactive_1d_runtime MPI::MPI_Fortran
    )
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(
        pelef_selected_mpi_reactive_1d_runtime
        PRIVATE
          -Wall
          -Wextra
          -Wpedantic
          -Wimplicit-interface
          -Wconversion-extra
      )
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(
          pelef_selected_mpi_reactive_1d_runtime
          PRIVATE
            -fcheck=all
            -fbacktrace
            -ffpe-trap=invalid,zero,overflow
        )
      endif()
    endif()
  endif()

  add_library(
    pelef_selected_amr_reactive_1d_runtime
    STATIC
    src/amr/amr_hierarchy_1d_mod.F90
    src/amr/amr_multipatch_1d_mod.F90
    src/amr/amr_patch_tree_1d_mod.F90
    src/amr/amr_regrid_1d_mod.F90
    src/amr/amr_reactive_1d_mod.F90
    src/amr/amr_patch_tree_reactive_1d_mod.F90
    src/amr/amr_multilevel_reactive_1d_mod.F90
    src/amr/amr_multipatch_reactive_1d_mod.F90
    src/driver/amr_reactive_1d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_amr_reactive_1d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_runtime_modules"
  )
  target_include_directories(
    pelef_selected_amr_reactive_1d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_runtime_modules>"
  )
  target_link_libraries(
    pelef_selected_amr_reactive_1d_runtime
    PUBLIC pelef_selected_reactive_1d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_amr_reactive_1d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_amr_reactive_1d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  if(PELEF_ENABLE_MPI)
    add_library(
      pelef_selected_mpi_amr_reactive_1d_runtime
      STATIC
      src/parallel/mpi_amr_patch_1d_mod.F90
      src/parallel/mpi_amr_sparse_patch_1d_mod.F90
      src/driver/mpi_amr_reactive_1d_application_mod.F90
    )
    set_target_properties(
      pelef_selected_mpi_amr_reactive_1d_runtime
      PROPERTIES
        Fortran_MODULE_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_1d_modules"
    )
    target_include_directories(
      pelef_selected_mpi_amr_reactive_1d_runtime
      PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_1d_modules>"
    )
    target_link_libraries(
      pelef_selected_mpi_amr_reactive_1d_runtime
      PUBLIC pelef_selected_amr_reactive_1d_runtime MPI::MPI_Fortran
    )
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(
        pelef_selected_mpi_amr_reactive_1d_runtime
        PRIVATE
          -Wall
          -Wextra
          -Wpedantic
          -Wimplicit-interface
          -Wconversion-extra
      )
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(
          pelef_selected_mpi_amr_reactive_1d_runtime
          PRIVATE
            -fcheck=all
            -fbacktrace
            -ffpe-trap=invalid,zero,overflow
        )
      endif()
    endif()
  endif()

  add_library(
    pelef_selected_reactive_2d_runtime
    STATIC
    src/core/mesh_mod.F90
    src/driver/simulation_config_reactive_2d_mod.F90
    src/boundary/reactive_boundary_2d_mod.F90
    src/transport/reactive_transport_2d_mod.F90
    src/reactive/reactive_2d_mod.F90
    src/driver/reactive_2d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_reactive_2d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_runtime_modules"
  )
  target_include_directories(
    pelef_selected_reactive_2d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_runtime_modules>"
  )
  target_link_libraries(
    pelef_selected_reactive_2d_runtime
    PUBLIC pelef_selected_reactive_1d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_reactive_2d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_reactive_2d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  add_library(
    pelef_selected_reactive_eb_2d_runtime
    STATIC
    src/eb/eb_geometry_2d_mod.F90
    src/driver/simulation_config_reactive_eb_2d_mod.F90
    src/eb/reactive_eb_cfl_2d_mod.F90
    src/eb/eb_reactive_wall_flux_2d_mod.F90
    src/eb/eb_reactive_redistribution_2d_mod.F90
    src/eb/eb_reactive_reconstruction_2d_mod.F90
    src/eb/eb_reactive_hydro_2d_mod.F90
    src/eb/eb_reactive_transport_2d_mod.F90
    src/driver/reactive_eb_2d_driver_mod.F90
    src/driver/reactive_eb_2d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_reactive_eb_2d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_runtime_modules"
  )
  target_include_directories(
    pelef_selected_reactive_eb_2d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_runtime_modules>"
  )
  target_link_libraries(
    pelef_selected_reactive_eb_2d_runtime
    PUBLIC pelef_selected_reactive_2d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_reactive_eb_2d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_reactive_eb_2d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  add_library(
    pelef_selected_reactive_eb_amr_2d_runtime
    STATIC
    src/driver/simulation_config_reactive_eb_amr_2d_mod.F90
    src/amr/amr_eb_hierarchy_2d_mod.F90
    src/amr/amr_eb_patch_tree_2d_mod.F90
    src/amr/amr_eb_multilevel_2d_mod.F90
    src/amr/amr_eb_flux_register_2d_mod.F90
    src/amr/amr_eb_reactive_2d_mod.F90
    src/amr/amr_eb_patch_tree_reactive_2d_mod.F90
    src/amr/amr_eb_multilevel_reactive_2d_mod.F90
    src/amr/amr_eb_regrid_2d_mod.F90
    src/amr/amr_eb_transport_2d_mod.F90
    src/amr/amr_eb_multilevel_transport_2d_mod.F90
    src/amr/amr_eb_multipatch_transport_2d_mod.F90
    src/driver/reactive_eb_amr_2d_driver_mod.F90
    src/driver/reactive_eb_amr_2d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_reactive_eb_amr_2d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_runtime_modules"
  )
  target_include_directories(
    pelef_selected_reactive_eb_amr_2d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_runtime_modules>"
  )
  target_link_libraries(
    pelef_selected_reactive_eb_amr_2d_runtime
    PUBLIC pelef_selected_reactive_eb_2d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_reactive_eb_amr_2d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_reactive_eb_amr_2d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  if(PELEF_ENABLE_MPI)
    add_library(
      pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
      STATIC
      src/parallel/mpi_amr_eb_patch_tree_2d_mod.F90
      src/parallel/mpi_amr_eb_patch_tree_io_2d_mod.F90
      src/driver/mpi_reactive_eb_patch_tree_2d_application_mod.F90
    )
    set_target_properties(
      pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
      PROPERTIES
        Fortran_MODULE_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_eb_patch_tree_2d_modules"
    )
    target_include_directories(
      pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
      PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_eb_patch_tree_2d_modules>"
    )
    target_link_libraries(
      pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
      PUBLIC pelef_selected_reactive_eb_amr_2d_runtime MPI::MPI_Fortran
    )
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(
        pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
        PRIVATE
          -Wall
          -Wextra
          -Wpedantic
          -Wimplicit-interface
          -Wconversion-extra
      )
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(
          pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
          PRIVATE
            -fcheck=all
            -fbacktrace
            -ffpe-trap=invalid,zero,overflow
        )
      endif()
    endif()
  endif()

  add_library(
    pelef_selected_reactive_3d_runtime
    STATIC
    src/core/mesh_3d_mod.F90
    src/driver/simulation_config_reactive_3d_mod.F90
    src/reactive/reactive_directional_flux_3d_mod.F90
    src/transport/reactive_transport_3d_mod.F90
    src/reactive/reactive_3d_mod.F90
    src/problems/reactive_entropy_wave_3d_problem_mod.F90
    src/io/reactive_csv_io_3d_mod.F90
    src/driver/reactive_3d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_reactive_3d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_3d_private_modules"
  )
  target_include_directories(
    pelef_selected_reactive_3d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_3d_private_modules>"
  )
  target_link_libraries(
    pelef_selected_reactive_3d_runtime
    PUBLIC pelef_selected_reactive_2d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_reactive_3d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_reactive_3d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  add_library(
    pelef_selected_amr_reactive_3d_runtime
    STATIC
    src/driver/simulation_config_amr_reactive_3d_mod.F90
    src/amr/amr_hierarchy_3d_mod.F90
    src/amr/amr_reactive_3d_mod.F90
    src/amr/amr_reactive_transport_3d_mod.F90
    src/io/amr_reactive_3d_checkpoint_mod.F90
    src/driver/amr_reactive_3d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_amr_reactive_3d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_modules"
  )
  target_include_directories(
    pelef_selected_amr_reactive_3d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_modules>"
  )
  target_link_libraries(
    pelef_selected_amr_reactive_3d_runtime
    PUBLIC pelef_selected_reactive_3d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_amr_reactive_3d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_amr_reactive_3d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  if(PELEF_ENABLE_MPI)
    add_library(
      pelef_selected_mpi_amr_reactive_3d_runtime
      STATIC
      src/parallel/mpi_amr_sparse_reactive_3d_mod.F90
      src/driver/mpi_amr_reactive_3d_application_mod.F90
    )
    set_target_properties(
      pelef_selected_mpi_amr_reactive_3d_runtime
      PROPERTIES
        Fortran_MODULE_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_modules"
    )
    target_include_directories(
      pelef_selected_mpi_amr_reactive_3d_runtime
      PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_modules>"
    )
    target_link_libraries(
      pelef_selected_mpi_amr_reactive_3d_runtime
      PUBLIC pelef_selected_amr_reactive_3d_runtime MPI::MPI_Fortran
    )
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(
        pelef_selected_mpi_amr_reactive_3d_runtime
        PRIVATE
          -Wall
          -Wextra
          -Wpedantic
          -Wimplicit-interface
          -Wconversion-extra
      )
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(
          pelef_selected_mpi_amr_reactive_3d_runtime
          PRIVATE
            -fcheck=all
            -fbacktrace
            -ffpe-trap=invalid,zero,overflow
        )
      endif()
    endif()
  endif()

  add_library(
    pelef_selected_reactive_eb_3d_runtime
    STATIC
    src/eb/eb_geometry_3d_mod.F90
    src/driver/simulation_config_reactive_eb_3d_mod.F90
    src/eb/reactive_eb_cfl_3d_mod.F90
    src/eb/eb_reactive_wall_flux_3d_mod.F90
    src/eb/eb_reactive_redistribution_3d_mod.F90
    src/eb/eb_reactive_hydro_3d_mod.F90
    src/eb/eb_reactive_transport_3d_mod.F90
    src/io/reactive_eb_3d_checkpoint_mod.F90
    src/driver/reactive_eb_3d_driver_mod.F90
    src/driver/reactive_eb_3d_application_mod.F90
  )
  set_target_properties(
    pelef_selected_reactive_eb_3d_runtime
    PROPERTIES
      Fortran_MODULE_DIRECTORY
        "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_runtime_modules"
  )
  target_include_directories(
    pelef_selected_reactive_eb_3d_runtime
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_runtime_modules>"
  )
  target_link_libraries(
    pelef_selected_reactive_eb_3d_runtime
    PUBLIC pelef_selected_reactive_3d_runtime
  )
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
      pelef_selected_reactive_eb_3d_runtime
      PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
    )
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(
        pelef_selected_reactive_eb_3d_runtime
        PRIVATE
          -fcheck=all
          -fbacktrace
          -ffpe-trap=invalid,zero,overflow
      )
    endif()
  endif()

  if(PELEF_ENABLE_SUNDIALS)
    add_library(
      pelef_selected_cvode_runtime
      STATIC src/chemistry/sundials_constant_volume_reactor_mod.F90
    )
    set_target_properties(
      pelef_selected_cvode_runtime
      PROPERTIES
        Fortran_MODULE_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_runtime_modules"
    )
    target_include_directories(
      pelef_selected_cvode_runtime
      PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_runtime_modules>"
    )
    target_link_libraries(
      pelef_selected_cvode_runtime
      PUBLIC
        pelef_selected_runtime_core
        SUNDIALS::fcvode_mod_static
        SUNDIALS::fnvecserial_mod_static
        SUNDIALS::fsunmatrixdense_mod_static
        SUNDIALS::fsunlinsoldense_mod_static
    )
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(
        pelef_selected_cvode_runtime
        PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra
      )
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(
          pelef_selected_cvode_runtime
          PRIVATE
            -fcheck=all
            -fbacktrace
            -ffpe-trap=invalid,zero,overflow
        )
      endif()
    endif()
  endif()
endif()

if(PELEF_MECHANISM_BUNDLE)
  set(pelef_selected_test_option)
  set(pelef_selected_mpi_option)
  if(PELEF_ENABLE_TESTS)
    list(APPEND pelef_selected_test_option ADD_TEST)
  endif()
  if(PELEF_ENABLE_MPI)
    list(
      APPEND pelef_selected_mpi_option
      MPI_REACTIVE_1D_NAME pelef_mpi_reactive_1d_selected
      MPI_AMR_REACTIVE_1D_NAME pelef_mpi_amr_reactive_1d_selected
      MPI_AMR_REACTIVE_3D_NAME pelef_mpi_amr_reactive_3d_selected
      MPI_REACTIVE_EB_PATCH_TREE_2D_NAME
        pelef_mpi_reactive_eb_patch_tree_2d_selected
    )
  endif()
  if(PELEF_MECHANISM_SOURCE)
    pelef_add_mechanism_bundle(
      NAME pelef_selected
      BUNDLE "${PELEF_MECHANISM_BUNDLE}"
      SOURCE "${PELEF_MECHANISM_SOURCE}"
      PROBE_NAME pelef_mechanism_probe
      REACTOR_NAME pelef0d_selected
      REACTIVE_1D_NAME pelef_reactive_1d_selected
      AMR_REACTIVE_1D_NAME pelef_amr_reactive_1d_selected
      REACTIVE_2D_NAME pelef_reactive_2d_selected
      REACTIVE_EB_2D_NAME pelef_reactive_eb_2d_selected
      REACTIVE_EB_AMR_2D_NAME pelef_reactive_eb_amr_2d_selected
      REACTIVE_3D_NAME pelef_reactive_3d_selected
      AMR_REACTIVE_3D_NAME pelef_amr_reactive_3d_selected
      REACTIVE_EB_3D_NAME pelef_reactive_eb_3d_selected
      ${pelef_selected_mpi_option}
      INSTALL
      ${pelef_selected_test_option}
    )
  else()
    pelef_add_mechanism_bundle(
      NAME pelef_selected
      BUNDLE "${PELEF_MECHANISM_BUNDLE}"
      PROBE_NAME pelef_mechanism_probe
      REACTOR_NAME pelef0d_selected
      REACTIVE_1D_NAME pelef_reactive_1d_selected
      AMR_REACTIVE_1D_NAME pelef_amr_reactive_1d_selected
      REACTIVE_2D_NAME pelef_reactive_2d_selected
      REACTIVE_EB_2D_NAME pelef_reactive_eb_2d_selected
      REACTIVE_EB_AMR_2D_NAME pelef_reactive_eb_amr_2d_selected
      REACTIVE_3D_NAME pelef_reactive_3d_selected
      AMR_REACTIVE_3D_NAME pelef_amr_reactive_3d_selected
      REACTIVE_EB_3D_NAME pelef_reactive_eb_3d_selected
      ${pelef_selected_mpi_option}
      INSTALL
      ${pelef_selected_test_option}
    )
  endif()
endif()

add_executable(pelef app/pelef.F90)
target_link_libraries(pelef PRIVATE pelef_core)

add_executable(pelef2d app/pelef2d.F90)
target_link_libraries(pelef2d PRIVATE pelef_core)

add_executable(pelef3d app/pelef3d.F90)
target_link_libraries(pelef3d PRIVATE pelef_core)

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

add_executable(pelef_reactive_3d app/pelef_reactive_3d.F90)
target_link_libraries(pelef_reactive_3d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_3d app/pelef_reactive_eb_3d.F90)
target_link_libraries(pelef_reactive_eb_3d PRIVATE pelef_core)

add_executable(pelef_amr_reactive_3d app/pelef_amr_reactive_3d.F90)
target_link_libraries(pelef_amr_reactive_3d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_2d app/pelef_reactive_eb_2d.F90)
target_link_libraries(pelef_reactive_eb_2d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_amr_2d app/pelef_reactive_eb_amr_2d.F90)
target_link_libraries(pelef_reactive_eb_amr_2d PRIVATE pelef_core)

add_executable(
  pelef_reactive_eb_patch_tree_2d
  app/pelef_reactive_eb_patch_tree_2d.F90
)
target_link_libraries(pelef_reactive_eb_patch_tree_2d PRIVATE pelef_core)

add_executable(pelef_transport_probe app/pelef_transport_probe.F90)
target_link_libraries(pelef_transport_probe PRIVATE pelef_core)

install(
  TARGETS pelef pelef2d pelef3d pelef_ms pelef0d pelef0d_h2o2
    pelef_reactive_1d
    pelef_amr_reactive_1d
    pelef0d_h2o2_full pelef_reactive_2d pelef_reactive_3d
    pelef_reactive_eb_3d
    pelef_amr_reactive_3d
    pelef_reactive_eb_2d
    pelef_reactive_eb_amr_2d
    pelef_reactive_eb_patch_tree_2d
    pelef_transport_probe
  RUNTIME DESTINATION bin
)

if(PELEF_ENABLE_TESTS)
  add_subdirectory(tests)
endif()


if(PELEF_ENABLE_MPI)
  add_library(
    pelef_mpi_support
    src/parallel/mpi_domain_1d_mod.F90
    src/parallel/mpi_amr_reactive_3d_mod.F90
    src/parallel/mpi_amr_sparse_reactive_3d_mod.F90
    src/parallel/mpi_amr_patch_1d_mod.F90
    src/parallel/mpi_amr_sparse_patch_1d_mod.F90
    src/parallel/mpi_amr_eb_patch_tree_2d_mod.F90
    src/parallel/mpi_amr_eb_patch_tree_io_2d_mod.F90
    src/driver/mpi_reactive_eb_patch_tree_2d_application_mod.F90
    src/parallel/mpi_amr_eb_patch_2d_mod.F90
    src/parallel/mpi_amr_eb_io_2d_mod.F90
    src/parallel/mpi_reactive_transport_1d_mod.F90
    src/parallel/mpi_reactive_1d_mod.F90
    src/driver/mpi_reactive_1d_application_mod.F90
    src/driver/mpi_amr_reactive_1d_application_mod.F90
    src/driver/mpi_amr_reactive_3d_application_mod.F90
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
  pelef_require_coherent_mpi_runtime(pelef_mpi_reactive_1d)
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
    pelef_mpi_amr_eb_patch_tree_2d
    app/pelef_mpi_amr_eb_patch_tree_2d.F90
  )
  target_link_libraries(
    pelef_mpi_amr_eb_patch_tree_2d PRIVATE pelef_mpi_support
  )
  add_executable(
    pelef_mpi_reactive_eb_patch_tree_2d
    app/pelef_mpi_reactive_eb_patch_tree_2d.F90
  )
  target_link_libraries(
    pelef_mpi_reactive_eb_patch_tree_2d PRIVATE pelef_mpi_support
  )
  add_executable(
    pelef_mpi_amr_reactive_1d
    app/pelef_mpi_amr_reactive_1d.F90
  )
  target_link_libraries(pelef_mpi_amr_reactive_1d PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_amr_reactive_3d
    app/pelef_mpi_amr_reactive_3d.F90
  )
  target_link_libraries(pelef_mpi_amr_reactive_3d PRIVATE pelef_mpi_support)
  if(PELEF_ENABLE_TESTS)
    add_executable(
      test_mpi_reactive_integrator_policy
      tests/mpi/test_mpi_reactive_integrator_policy.F90
    )
    target_link_libraries(
      test_mpi_reactive_integrator_policy PRIVATE pelef_mpi_support
    )
    foreach(mpi_integrator_policy_ranks IN ITEMS 1 2 4)
      add_test(
        NAME
          unit_mpi_reactive_integrator_policy_np${mpi_integrator_policy_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_integrator_policy_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_reactive_integrator_policy>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_reactive_integrator_policy_np${mpi_integrator_policy_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_integrator_policy_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_reactive_integrator_policy: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_sparse_integrator_policy
      tests/mpi/test_mpi_amr_sparse_integrator_policy.F90
    )
    target_link_libraries(
      test_mpi_amr_sparse_integrator_policy PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_integrator_policy_ranks IN ITEMS 1 2 4)
      add_test(
        NAME
          unit_mpi_amr_sparse_integrator_policy_np${mpi_amr_integrator_policy_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_integrator_policy_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_sparse_integrator_policy>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_amr_sparse_integrator_policy_np${mpi_amr_integrator_policy_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_integrator_policy_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_sparse_integrator_policy: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_eb_sparse_integrator_policy
      tests/mpi/test_mpi_amr_eb_sparse_integrator_policy.F90
    )
    target_link_libraries(
      test_mpi_amr_eb_sparse_integrator_policy PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_eb_integrator_policy_ranks IN ITEMS 1 2 4)
      add_test(
        NAME
          unit_mpi_amr_eb_sparse_integrator_policy_np${mpi_amr_eb_integrator_policy_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_eb_integrator_policy_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_eb_sparse_integrator_policy>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_amr_eb_sparse_integrator_policy_np${mpi_amr_eb_integrator_policy_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_eb_integrator_policy_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_eb_sparse_integrator_policy: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_selected_context_mismatch
      tests/mpi/test_mpi_amr_selected_context_mismatch.F90
    )
    target_link_libraries(
      test_mpi_amr_selected_context_mismatch PRIVATE pelef_mpi_support
    )
    set(
      MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_amr_selected_context_mismatch"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_selected_context_mismatch
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
        --launcher "${MPIEXEC_EXECUTABLE}"
        "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --ranks 2
        --executable $<TARGET_FILE:test_mpi_amr_selected_context_mismatch>
        --expected "Selected sparse MPI AMR rank context mismatch"
        --forbidden-output
          "${MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}/mpi_amr_context_mismatch.csv"
    )
    set_tests_properties(
      regression_mpi_amr_selected_context_mismatch
      PROPERTIES
        WORKING_DIRECTORY
          "${MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        TIMEOUT 60
    )

    add_executable(
      test_mpi_amr_selected_context_3d
      tests/mpi/test_mpi_amr_selected_context_3d.F90
    )
    target_link_libraries(
      test_mpi_amr_selected_context_3d PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_selected_context_3d_ranks IN ITEMS 2 4)
      add_test(
        NAME
          unit_mpi_amr_selected_context_3d_np${mpi_amr_selected_context_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_selected_context_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_selected_context_3d>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_amr_selected_context_3d_np${mpi_amr_selected_context_3d_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_selected_context_3d_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_selected_context_3d: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_eb_selected_context_mismatch
      tests/mpi/test_mpi_amr_eb_selected_context_mismatch.F90
    )
    target_link_libraries(
      test_mpi_amr_eb_selected_context_mismatch PRIVATE pelef_mpi_support
    )
    set(
      MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_amr_eb_selected_context_mismatch"
    )
    file(
      MAKE_DIRECTORY
      "${MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}"
    )
    add_test(
      NAME regression_mpi_amr_eb_selected_context_mismatch
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
        --launcher "${MPIEXEC_EXECUTABLE}"
        "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --ranks 2
        --executable $<TARGET_FILE:test_mpi_amr_eb_selected_context_mismatch>
        --expected "Selected sparse MPI EB rank context mismatch"
        --forbidden-output
          "${MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}/mpi_amr_eb_context_mismatch.csv"
    )
    set_tests_properties(
      regression_mpi_amr_eb_selected_context_mismatch
      PROPERTIES
        WORKING_DIRECTORY
          "${MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        TIMEOUT 60
    )

    add_executable(
      test_mpi_amr_eb_selected_checkpoint_presence
      tests/mpi/test_mpi_amr_eb_selected_checkpoint_presence.F90
    )
    target_link_libraries(
      test_mpi_amr_eb_selected_checkpoint_presence PRIVATE pelef_mpi_support
    )
    add_test(
      NAME unit_mpi_amr_eb_selected_checkpoint_presence_np2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_eb_selected_checkpoint_presence>
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      unit_mpi_amr_eb_selected_checkpoint_presence_np2
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION
          "test_mpi_amr_eb_selected_checkpoint_presence: PASS"
        TIMEOUT 60
    )

    set(
      SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/tests/selected_mpi_reactive_1d_full"
    )
    file(MAKE_DIRECTORY "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_reactive_1d_fixed_np1_reference
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_reactive_1d>
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/fixed_np1.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_mpi_reactive_1d_fixed_np1_reference
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 1
        PASS_REGULAR_EXPRESSION "MPI ranks: 1"
        TIMEOUT 300
    )
    add_test(
      NAME regression_selected_mpi_reactive_1d_full_fixed_parity
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_reactive_1d.py"
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/fixed_np1.csv"
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np1.csv"
    )
    set_tests_properties(
      regression_selected_mpi_reactive_1d_full_fixed_parity
      PROPERTIES
        DEPENDS
          "regression_mpi_reactive_1d_fixed_np1_reference;regression_selected_mpi_reactive_1d_full_np1"
    )

    add_executable(
      test_mpi_amr_reactive_3d
      tests/mpi/test_mpi_amr_reactive_3d.F90
    )
    target_link_libraries(test_mpi_amr_reactive_3d PRIVATE pelef_mpi_support)
    foreach(mpi_amr_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_reactive_3d> ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_3d_ranks}
          PASS_REGULAR_EXPRESSION "test_mpi_amr_reactive_3d: PASS"
          TIMEOUT 300
      )
    endforeach()

    add_executable(
      test_mpi_amr_reactive_chemistry_3d
      tests/mpi/test_mpi_amr_reactive_chemistry_3d.F90
    )
    target_link_libraries(
      test_mpi_amr_reactive_chemistry_3d PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_chemistry_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_chemistry_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_reactive_chemistry_3d>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_3d_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_chemistry_3d_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_reactive_chemistry_3d: PASS"
          TIMEOUT 300
      )
    endforeach()

    set(
      MPI_AMR_REACTIVE_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "static two-level 3D reactive AMR"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_3D_RUN_TESTS "")
    foreach(mpi_amr_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave.nml"
          "mpi_amr_np${mpi_amr_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_3d_ranks}
          PASS_REGULAR_EXPRESSION "distributed-slab static two-level 3D"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_3D_RUN_TESTS
        regression_mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np1_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np1_fine.csv"
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --time 1.0e-6
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_3d_np1
    )
    add_test(
      NAME regression_mpi_amr_reactive_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np8"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_3d_serial_reference;${MPI_AMR_REACTIVE_3D_RUN_TESTS}"
    )

    set(
      MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_transport_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_transport_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "Completed transport fine substeps: 16"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_TRANSPORT_3D_RUN_TESTS "")
    foreach(mpi_amr_transport_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_transport_3d_np${mpi_amr_transport_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_transport_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport.nml"
          "mpi_amr_transport_np${mpi_amr_transport_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_transport_3d_np${mpi_amr_transport_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_transport_3d_ranks}
          PASS_REGULAR_EXPRESSION "Completed transport fine substeps: 16"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_TRANSPORT_3D_RUN_TESTS
        regression_mpi_amr_reactive_transport_3d_np${mpi_amr_transport_3d_ranks}
      )
    endforeach()
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_np8
      PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
    )
    add_test(
      NAME regression_mpi_amr_reactive_transport_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np8_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np8_fine.csv"
        --species H2 H O O2 OH H2O HO2 H2O2 AR N2
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --final-time 1.0e-9
        --activity-species H2
        --initial-mass-fraction 0.028502723329253223
        --minimum-change 1.0e-12
    )
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_transport_3d_np8
    )
    add_test(
      NAME regression_mpi_amr_reactive_transport_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/amr_hotspot_transport_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np8"
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --expected-coarse-sha256
          $<IF:$<CONFIG:Debug>,5687cb97f834b759b3932bbd3f0b6e30fc5b48e6029e3f47cadfb3ec6a8d74fb,f0713f3a68f48b6537afaafd57209ea8ed4767f9c7d093fe2961bd3a1fa37fae>
        --expected-fine-sha256
          eaf132bbf89974d739c26a35f6b76137f9774438cfe9eebfc12436aa364cff58
    )
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_transport_3d_serial_reference;${MPI_AMR_REACTIVE_TRANSPORT_3D_RUN_TESTS}"
    )

    set(
      MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_chemistry_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_chemistry_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_chemistry.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_chemistry_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "Chemistry:  T"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_CHEMISTRY_3D_RUN_TESTS "")
    foreach(mpi_amr_chemistry_app_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_app_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_chemistry_app_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_chemistry.nml"
          "mpi_amr_chemistry_np${mpi_amr_chemistry_app_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_app_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_chemistry_app_ranks}
          PASS_REGULAR_EXPRESSION "Chemistry:  T"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_CHEMISTRY_3D_RUN_TESTS
        regression_mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_app_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_chemistry_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np8_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np8_fine.csv"
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --time 1.0e-9
    )
    set_tests_properties(
      regression_mpi_amr_reactive_chemistry_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_chemistry_3d_np8
    )
    add_test(
      NAME regression_mpi_amr_reactive_chemistry_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/amr_hotspot_chemistry_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np8"
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --expected-coarse-sha256
          3313a9a6eff5d6175625e4bbe3613a782cdd794314fd9fbdd53ddaa309fc9516
        --expected-fine-sha256
          4ff0785c3462a49648fd694d090eb23bffd5af05a43a33360cd187417bd8a31d
    )
    set_tests_properties(
      regression_mpi_amr_reactive_chemistry_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_chemistry_3d_serial_reference;${MPI_AMR_REACTIVE_CHEMISTRY_3D_RUN_TESTS}"
    )

    add_test(
      NAME regression_mpi_amr_reactive_3d_checkpoint_np2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml"
        "mpi_amr_checkpoint_np2"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_3D_RESTART_TESTS "")
    foreach(mpi_amr_3d_ranks IN ITEMS 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_3d_restart_np${mpi_amr_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_restart.nml"
          "mpi_amr_restart_np${mpi_amr_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_3d_restart_np${mpi_amr_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_3d_ranks}
          PASS_REGULAR_EXPRESSION "Restored coarse steps: 4"
          DEPENDS regression_mpi_amr_reactive_3d_checkpoint_np2
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_3D_RESTART_TESTS
        regression_mpi_amr_reactive_3d_restart_np${mpi_amr_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_3d_restart_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_restart_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_restart_np8"
        --checkpoint
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d.chk"
        --final-time 1.0e-6
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_restart_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_3d_serial_reference;${MPI_AMR_REACTIVE_3D_RESTART_TESTS}"
    )

    set(
      MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_plm_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_PLM_3D_RUN_TESTS "")
    foreach(mpi_amr_plm_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_plm_3d_np${mpi_amr_plm_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_plm_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm.nml"
          "mpi_amr_plm_np${mpi_amr_plm_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_plm_3d_np${mpi_amr_plm_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_plm_3d_ranks}
          PASS_REGULAR_EXPRESSION "rank-local sparse distributed-slab"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_PLM_3D_RUN_TESTS
        regression_mpi_amr_reactive_plm_3d_np${mpi_amr_plm_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np1_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np1_fine.csv"
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --time 1.0e-6
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_plm_3d_np1
    )
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np8"
        --expected-coarse-sha256
          864e9cabe8f5c2bee19f97e9279bdada7647305171955c6c767cc83e4db15c21
        --expected-fine-sha256
          8f764e38491e9cf62f8c0f17eac2470c4705f7d80bfe7e686db4ba54de2f58d7
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_plm_3d_serial_reference;${MPI_AMR_REACTIVE_PLM_3D_RUN_TESTS}"
    )
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_checkpoint_np2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm_checkpoint.nml"
        "mpi_amr_plm_checkpoint_np2"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_PLM_3D_RESTART_TESTS "")
    foreach(mpi_amr_plm_3d_ranks IN ITEMS 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_plm_3d_restart_np${mpi_amr_plm_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_plm_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm_restart.nml"
          "mpi_amr_plm_restart_np${mpi_amr_plm_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_plm_3d_restart_np${mpi_amr_plm_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_plm_3d_ranks}
          PASS_REGULAR_EXPRESSION "Restored coarse steps: 4"
          DEPENDS regression_mpi_amr_reactive_plm_3d_checkpoint_np2
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_PLM_3D_RESTART_TESTS
        regression_mpi_amr_reactive_plm_3d_restart_np${mpi_amr_plm_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_restart_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_restart_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_restart_np8"
        --checkpoint
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d.chk"
        --final-time 1.0e-6
        --reconstruction characteristic_plm
        --limiter mc
        --expected-checkpoint-sha256
          8f51be7f192dd61ab0e1f4eccf665870c2103205d13ea562fb52d32e981dbf9d
        --expected-coarse-sha256
          864e9cabe8f5c2bee19f97e9279bdada7647305171955c6c767cc83e4db15c21
        --expected-fine-sha256
          8f764e38491e9cf62f8c0f17eac2470c4705f7d80bfe7e686db4ba54de2f58d7
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_restart_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_plm_3d_serial_reference;${MPI_AMR_REACTIVE_PLM_3D_RESTART_TESTS}"
    )
  endif()
  install(
    TARGETS
      pelef_mpi_1d
      pelef_mpi_multispecies_1d
      pelef_mpi_h2o2_batch
      pelef_mpi_reactive_transport_1d
      pelef_mpi_reactive_1d
      pelef_mpi_amr_patch_1d
      pelef_mpi_eb_amr_patch_2d
      pelef_mpi_amr_eb_patch_tree_2d
      pelef_mpi_reactive_eb_patch_tree_2d
      pelef_mpi_amr_reactive_1d
      pelef_mpi_amr_reactive_3d
    RUNTIME DESTINATION bin
  )
endif()

if(PELEF_ENABLE_TESTS)
  # Exercise the installed 3D AMR transport/restart contract in a private
  # prefix.  The ordinary 0.245 tests above run build-tree executables; these
  # tests deliberately resolve the binaries through the install tree so that
  # install rules and restart schema selection are covered as well.
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT
    "${CMAKE_CURRENT_BINARY_DIR}/installed_amr_reactive_3d_transport_restart_0245"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/prefix"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/fixed"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/selected"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/fixed-mpi"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/selected-mpi"
  )
  file(
    MAKE_DIRECTORY
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}"
  )
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256
      "6a1b1ba7bece94e4f59e0da95fb6dc8ddc5bc5d356f3d88f56d37e184e0b573d"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256
      "0b0955e45d7d5e6296bdd69270cd78ed02eba4f814a3789d7218d861fd0eac6f"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256
      "45d51c690e8f3b1d20bb224deae24cbbb5ec4dc0c648a9b1367d6963acaeb37b"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256
      "6fa40282c6a12598bf7147efe2253540e49a66abb5e907aac8b5226443fdf342"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256
      "34d781bd9faef30a0fe5f859300fac459380403d67b00795f1c3da95e294cb39"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256
      "e4c3f6fe27ed9c8e62d33912c42adc6bdef4ee1c3e8cd0b6aebfdf34c42956ed"
    )
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256
      "15c2e5a3a087a38fe34ece1fa575bbc7857ec7c651ae9ee84f3fd8d52615afc5"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256
      "c81eb517ab4111f66eafa09b70d1f9543251e206d62143cfe02e50ca74e4f441"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256
      "45d51c690e8f3b1d20bb224deae24cbbb5ec4dc0c648a9b1367d6963acaeb37b"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256
      "52526217eb88e869167b8ef8b958f3cf6176b5346746a5d9ff1d1757ca526df3"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256
      "21470b7ac704274637266990f8bd77bbeacdba5c2808c2943fdc17fb106ee35a"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256
      "e4c3f6fe27ed9c8e62d33912c42adc6bdef4ee1c3e8cd0b6aebfdf34c42956ed"
    )
  else()
    message(
      FATAL_ERROR
      "installed 0.245 transport restart tests require "
      "CMAKE_BUILD_TYPE=Debug or Release"
    )
  endif()

  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST
    regression_installed_amr_reactive_3d_transport_restart_0245_install
  )
  add_test(
    NAME ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
    COMMAND
      "${CMAKE_COMMAND}"
      "-DTEST_BUILD_DIR=${CMAKE_BINARY_DIR}"
      "-DTEST_INSTALL_PREFIX=${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}"
      "-DTEST_CONFIG=$<CONFIG>"
      -P "${PROJECT_SOURCE_DIR}/tests/cmake/install_clean_prefix.cmake"
  )
  set_tests_properties(
    ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
    PROPERTIES
      WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}"
      TIMEOUT 600
  )

  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/reference.log"
      --
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d"
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_reference.nml"
  )
  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/checkpoint.log"
      --
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d"
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
  )
  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/restart.log"
      --
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d"
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run
    PROPERTIES
      WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}"
      TIMEOUT 300
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    PROPERTIES DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run
    PROPERTIES
      DEPENDS
        regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
      PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
  )
  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_exact
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_reference_0245"
      --candidate-prefix
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_restarted_0245"
      --checkpoint
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
      --expected-checkpoint-sha256
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256}"
      --expected-coarse-sha256
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256}"
      --expected-fine-sha256
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256}"
      --expected-schema 4
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --reference-log
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/reference.log"
      --candidate-log
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/restart.log"
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_exact
    PROPERTIES DEPENDS
      "regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference;regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint;regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run"
  )

  set(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED OFF)
  if(PELEF_MECHANISM_BUNDLE)
    file(SHA256 "${PELEF_MECHANISM_BUNDLE}"
      pelef_installed_transport_0245_bundle_sha256)
    if(pelef_installed_transport_0245_bundle_sha256 STREQUAL
        "f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3")
      set(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED ON)
    endif()
  endif()

  if(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED)
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/reference.log"
        --
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_reference.nml"
    )
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/checkpoint.log"
        --
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
    )
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_run
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/restart.log"
        --
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_run
      PROPERTIES
        WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}"
        TIMEOUT 600
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      PROPERTIES DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_run
      PROPERTIES
        DEPENDS
          regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
        PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    )
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_reference_0245"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_restarted_0245"
        --checkpoint
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_0245.chk"
        --expected-checkpoint-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256}"
        --expected-coarse-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256}"
        --expected-fine-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256}"
        --expected-schema 5
        --expected-bundle-sha256
          f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
        --expected-integrator implicit
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --reference-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/reference.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/restart.log"
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference;regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint;regression_installed_amr_reactive_3d_transport_restart_0245_selected_run"
    )
  endif()

  if(PELEF_ENABLE_MPI)
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/checkpoint-np2.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d"
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
          installed_fixed_transport_checkpoint_np2
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
        TIMEOUT 300
    )
    set(PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_TESTS)
    foreach(installed_fixed_mpi_rank IN ITEMS 1 2 4 8)
      set(
        installed_fixed_mpi_restart_test
        regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_np${installed_fixed_mpi_rank}
      )
      add_test(
        NAME ${installed_fixed_mpi_restart_test}
        COMMAND
          "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
          --output
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np${installed_fixed_mpi_rank}.log"
          -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          ${installed_fixed_mpi_rank} ${MPIEXEC_PREFLAGS}
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d"
            "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
            "installed_fixed_transport_restart_mpi_np${installed_fixed_mpi_rank}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        ${installed_fixed_mpi_restart_test}
        PROPERTIES
          WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${installed_fixed_mpi_rank}
          PASS_REGULAR_EXPRESSION "MPI ranks: ${installed_fixed_mpi_rank}"
          DEPENDS
            regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2
          TIMEOUT 300
      )
      list(APPEND PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_TESTS
        ${installed_fixed_mpi_restart_test})
    endforeach()
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_np8
      PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_reference_0245"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np1"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np2"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np4"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np8"
        --checkpoint
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
        --expected-checkpoint-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256}"
        --expected-coarse-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256}"
        --expected-fine-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256}"
        --expected-schema 4
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --reference-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/reference.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np1.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np2.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np4.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np8.log"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2;${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_TESTS}"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_exact
      COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2"
    )
  endif()

  if(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED AND PELEF_ENABLE_MPI)
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/checkpoint-np2.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
          installed_selected_transport_checkpoint_np2
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
        TIMEOUT 600
    )
    set(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_TESTS)
    foreach(installed_selected_mpi_rank IN ITEMS 1 2 4 8)
      set(
        installed_selected_mpi_restart_test
        regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_np${installed_selected_mpi_rank}
      )
      add_test(
        NAME ${installed_selected_mpi_restart_test}
        COMMAND
          "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
          --output
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np${installed_selected_mpi_rank}.log"
          -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          ${installed_selected_mpi_rank} ${MPIEXEC_PREFLAGS}
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d_selected"
            "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
            "installed_selected_transport_restart_mpi_np${installed_selected_mpi_rank}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        ${installed_selected_mpi_restart_test}
        PROPERTIES
          WORKING_DIRECTORY
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${installed_selected_mpi_rank}
          PASS_REGULAR_EXPRESSION "MPI ranks: ${installed_selected_mpi_rank}"
          DEPENDS
            regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2
          TIMEOUT 600
      )
      list(APPEND PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_TESTS
        ${installed_selected_mpi_restart_test})
    endforeach()
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_np8
      PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_reference_0245"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np1"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np2"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np4"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np8"
        --checkpoint
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/selected_h2o2_full_transport_restart_0245.chk"
        --expected-checkpoint-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256}"
        --expected-coarse-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256}"
        --expected-fine-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256}"
        --expected-schema 5
        --expected-bundle-sha256
          f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
        --expected-integrator implicit
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --reference-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/reference.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np1.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np2.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np4.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np8.log"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2;${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_TESTS}"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_exact
      COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_0245.chk"
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/selected_h2o2_full_transport_restart_0245.chk"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2"
    )
  endif()

  add_test(
    NAME quality_project_contract
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_project_contract.py"
      --source "${PROJECT_SOURCE_DIR}"
  )
  add_test(
    NAME quality_mpi_f08_contract_accept
    COMMAND "${CMAKE_COMMAND}" -DTEST_MPI_F08=ON -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_mpi_f08_contract.cmake"
  )
  add_test(
    NAME quality_mpi_f08_contract_reject
    COMMAND "${CMAKE_COMMAND}" -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_mpi_f08_rejection.cmake"
  )
  add_test(
    NAME quality_mpi_launcher_contract_reject
    COMMAND "${CMAKE_COMMAND}" -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_mpi_launcher_rejection.cmake"
  )
  add_test(
    NAME quality_sundials_double_precision_accept
    COMMAND "${CMAKE_COMMAND}"
      -DTEST_PRECISION=double
      -DTEST_BINARY_DIR=${PROJECT_BINARY_DIR}/tests/cmake
      -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_sundials_precision_contract.cmake"
  )
  add_test(
    NAME quality_sundials_non_double_precision_reject
    COMMAND "${CMAKE_COMMAND}"
      -DTEST_BINARY_DIR=${PROJECT_BINARY_DIR}/tests/cmake
      -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_sundials_precision_rejection.cmake"
  )

  set(
    PELEF_STACK_CHECK_TARGETS
    $<TARGET_FILE:pelef>
    $<TARGET_FILE:pelef2d>
    $<TARGET_FILE:pelef3d>
    $<TARGET_FILE:pelef_ms>
    $<TARGET_FILE:pelef0d>
    $<TARGET_FILE:pelef0d_h2o2>
    $<TARGET_FILE:pelef0d_h2o2_full>
    $<TARGET_FILE:pelef_reactive_1d>
    $<TARGET_FILE:pelef_amr_reactive_1d>
    $<TARGET_FILE:pelef_reactive_2d>
    $<TARGET_FILE:pelef_reactive_3d>
    $<TARGET_FILE:pelef_reactive_eb_3d>
    $<TARGET_FILE:pelef_amr_reactive_3d>
    $<TARGET_FILE:pelef_reactive_eb_2d>
    $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    $<TARGET_FILE:pelef_transport_probe>
    $<TARGET_FILE:test_amr_eb_multilevel_2d>
    $<TARGET_FILE:test_selected_reactor_config>
    $<TARGET_FILE:test_mechanism_bundle>
    $<TARGET_FILE:test_selected_fixture_reactor>
    $<TARGET_FILE:test_selected_fixture_reactive_1d>
    $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    $<TARGET_FILE:test_selected_fixture_reactive_2d>
    $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    $<TARGET_FILE:test_selected_fixture_reactive_3d>
    $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    $<TARGET_FILE:test_selected_full_bundle>
    $<TARGET_FILE:test_selected_full_reactor>
    $<TARGET_FILE:test_selected_full_reactive_1d>
    $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    $<TARGET_FILE:test_selected_full_reactive_2d>
    $<TARGET_FILE:test_selected_full_reactive_eb_2d>
    $<TARGET_FILE:test_selected_full_reactive_3d>
    $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    $<TARGET_FILE:test_selected_full_reactive_eb_3d>
    $<TARGET_FILE:test_selected_elementary_bundle>
    $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
    $<TARGET_FILE:test_selected_unsupported_element_bundle>
    $<TARGET_FILE:test_selected_unsupported_element_reactive_eb_3d>
  )
  if(PELEF_MECHANISM_BUNDLE)
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:pelef_mechanism_probe>
      $<TARGET_FILE:pelef0d_selected>
      $<TARGET_FILE:pelef_reactive_1d_selected>
      $<TARGET_FILE:pelef_amr_reactive_1d_selected>
      $<TARGET_FILE:pelef_reactive_2d_selected>
      $<TARGET_FILE:pelef_reactive_eb_2d_selected>
      $<TARGET_FILE:pelef_reactive_eb_amr_2d_selected>
      $<TARGET_FILE:pelef_reactive_3d_selected>
      $<TARGET_FILE:pelef_amr_reactive_3d_selected>
      $<TARGET_FILE:pelef_reactive_eb_3d_selected>
    )
    if(PELEF_ENABLE_MPI)
      list(
        APPEND PELEF_STACK_CHECK_TARGETS
        $<TARGET_FILE:pelef_mpi_reactive_1d_selected>
        $<TARGET_FILE:pelef_mpi_amr_reactive_1d_selected>
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d_selected>
        $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d_selected>
      )
    endif()
  endif()
  if(PELEF_ENABLE_SUNDIALS)
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:test_sundials_constant_volume_reactor>
      $<TARGET_FILE:test_sundials_constant_volume_multi_context>
    )
  endif()
  if(PELEF_ENABLE_MPI)
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:pelef_mpi_1d>
      $<TARGET_FILE:pelef_mpi_multispecies_1d>
      $<TARGET_FILE:pelef_mpi_h2o2_batch>
      $<TARGET_FILE:pelef_mpi_reactive_transport_1d>
      $<TARGET_FILE:pelef_mpi_reactive_1d>
      $<TARGET_FILE:pelef_mpi_amr_patch_1d>
      $<TARGET_FILE:pelef_mpi_eb_amr_patch_2d>
      $<TARGET_FILE:pelef_mpi_amr_eb_patch_tree_2d>
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      $<TARGET_FILE:pelef_mpi_amr_reactive_1d>
      $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
      $<TARGET_FILE:test_mpi_reactive_integrator_policy>
      $<TARGET_FILE:test_mpi_amr_sparse_integrator_policy>
      $<TARGET_FILE:test_mpi_amr_eb_sparse_integrator_policy>
      $<TARGET_FILE:test_mpi_amr_selected_context_mismatch>
      $<TARGET_FILE:test_mpi_amr_selected_context_3d>
      $<TARGET_FILE:test_mpi_amr_eb_selected_context_mismatch>
      $<TARGET_FILE:test_mpi_amr_eb_selected_checkpoint_presence>
    )
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:test_selected_fixture_mpi_reactive_1d>
      $<TARGET_FILE:test_selected_full_mpi_reactive_1d>
      $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_1d>
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
      $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
      $<TARGET_FILE:test_selected_fixture_mpi_reactive_eb_patch_tree_2d>
      $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
    )
  endif()
  add_test(
    NAME quality_no_executable_stack
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_no_executable_stack.py"
      ${PELEF_STACK_CHECK_TARGETS}
  )
endif()
