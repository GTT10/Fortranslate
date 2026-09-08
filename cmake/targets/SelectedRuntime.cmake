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

