if(PELEF_ENABLE_SPRAY)
  add_test(NAME spray_eb_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_application")
  set_tests_properties(spray_eb_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_executable(test_spray_eb_reactive spray/test_spray_eb_reactive.F90)
  target_link_libraries(test_spray_eb_reactive PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_reactive COMMAND test_spray_eb_reactive)
  set_tests_properties(spray_eb_reactive PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_step spray/test_spray_eb_step.F90)
  target_link_libraries(test_spray_eb_step PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_step COMMAND test_spray_eb_step)
  set_tests_properties(spray_eb_step PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_field spray/test_spray_eb_field.F90)
  target_link_libraries(test_spray_eb_field PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_field COMMAND test_spray_eb_field)
  set_tests_properties(spray_eb_field PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_drift spray/test_spray_eb_drift.F90)
  target_link_libraries(test_spray_eb_drift PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_drift COMMAND test_spray_eb_drift)
  set_tests_properties(spray_eb_drift PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_intersection spray/test_spray_eb_intersection.F90)
  target_link_libraries(test_spray_eb_intersection PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_intersection COMMAND test_spray_eb_intersection)
  set_tests_properties(spray_eb_intersection PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_amr_tree_gas spray/test_spray_amr_tree_gas.F90)
  target_link_libraries(test_spray_amr_tree_gas PRIVATE pelef_spray_tree_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_spray_amr_tree_gas PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME spray_amr_tree_gas COMMAND test_spray_amr_tree_gas)
  set_tests_properties(spray_amr_tree_gas PROPERTIES LABELS "spray;amr" TIMEOUT 30)
  add_executable(test_spray_amr_tree_coupling spray/test_spray_amr_tree_coupling.F90)
  target_link_libraries(test_spray_amr_tree_coupling PRIVATE pelef_spray_tree_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_spray_amr_tree_coupling PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME spray_amr_tree_coupling COMMAND test_spray_amr_tree_coupling)
  set_tests_properties(spray_amr_tree_coupling PROPERTIES LABELS "spray;amr" TIMEOUT 30)
  add_executable(test_spray_amr_tree_partition spray/test_spray_amr_tree_partition.F90)
  target_link_libraries(test_spray_amr_tree_partition PRIVATE pelef_spray_tree_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_spray_amr_tree_partition PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME spray_amr_tree_partition COMMAND test_spray_amr_tree_partition)
  set_tests_properties(spray_amr_tree_partition PROPERTIES LABELS "spray;amr" TIMEOUT 30)
  if(PELEF_ENABLE_MPI)
    add_test(NAME spray_mpi_amr_film_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_film_application"
        --multiport --physical --films)
    set_tests_properties(spray_mpi_amr_film_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_inflow_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_inflow_application"
        --multiport --physical --inflow)
    set_tests_properties(spray_mpi_amr_inflow_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_physical_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_physical_application" --multiport --physical)
    set_tests_properties(spray_mpi_amr_physical_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_rebalance_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_rebalance_application" --multiport --adaptive --rebalance)
    set_tests_properties(spray_mpi_amr_rebalance_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_adaptation
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_adaptation" --multiport --adaptive)
    set_tests_properties(spray_mpi_amr_adaptation PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_ports
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_ports" --multiport)
    set_tests_properties(spray_mpi_amr_ports PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_restart")
    set_tests_properties(spray_mpi_amr_restart PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_application")
    set_tests_properties(spray_mpi_amr_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_executable(test_spray_mpi_amr_balance spray/test_spray_mpi_amr_balance.F90)
    target_link_libraries(test_spray_mpi_amr_balance PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_balance PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film_reactive spray/test_spray_mpi_amr_film_reactive.F90)
    target_link_libraries(test_spray_mpi_amr_film_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film_reactive PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film_checkpoint spray/test_spray_mpi_amr_film_checkpoint.F90)
    target_link_libraries(test_spray_mpi_amr_film_checkpoint PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film_checkpoint PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film_step spray/test_spray_mpi_amr_film_step.F90)
    target_link_libraries(test_spray_mpi_amr_film_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film_step PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film spray/test_spray_mpi_amr_film.F90)
    target_link_libraries(test_spray_mpi_amr_film PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_remap spray/test_spray_mpi_amr_remap.F90)
    target_link_libraries(test_spray_mpi_amr_remap PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_remap PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_reactive spray/test_spray_mpi_amr_reactive.F90)
    target_link_libraries(test_spray_mpi_amr_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_reactive PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_impact spray/test_spray_mpi_amr_impact.F90)
    target_link_libraries(test_spray_mpi_amr_impact PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_impact PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_boundaries spray/test_spray_mpi_amr_boundaries.F90)
    target_link_libraries(test_spray_mpi_amr_boundaries PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_boundaries PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_step spray/test_spray_mpi_amr_step.F90)
    target_link_libraries(test_spray_mpi_amr_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_step PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_gas spray/test_spray_mpi_amr_gas.F90)
    target_link_libraries(test_spray_mpi_amr_gas PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_gas PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_deposition spray/test_spray_mpi_amr_deposition.F90)
    target_link_libraries(test_spray_mpi_amr_deposition PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_deposition PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_sampling spray/test_spray_mpi_amr_sampling.F90)
    target_link_libraries(test_spray_mpi_amr_sampling PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_sampling PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_tree spray/test_spray_mpi_amr_tree.F90)
    target_link_libraries(test_spray_mpi_amr_tree PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_tree PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi spray/test_spray_mpi.F90)
    target_link_libraries(test_spray_mpi PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    foreach(ranks 2 4)
      add_test(NAME spray_mpi_amr_balance_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_balance> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_balance_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_remap_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_remap> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_remap_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_film_reactive_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_checkpoint_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film_checkpoint> ${MPIEXEC_POSTFLAGS}
          "${CMAKE_CURRENT_BINARY_DIR}/film_checkpoint_${ranks}")
      set_tests_properties(spray_mpi_amr_film_checkpoint_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_film_step_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_film_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_reactive_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 120 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_impact_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_impact> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_impact_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_boundaries_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_boundaries> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_boundaries_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_step_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_gas_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_gas> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_gas_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_deposition_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_deposition> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_deposition_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_sampling_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_sampling> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_sampling_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_tree_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_tree> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_tree_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_application_${ranks}
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
          --executable $<TARGET_FILE:pelef_mpi_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
          --ranks ${ranks} --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
          --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_application_${ranks}")
      set_tests_properties(spray_mpi_application_${ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
      add_test(NAME spray_mpi_migration_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_migration_${ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
  endif()
  add_test(NAME spray_acceleration_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --acceleration -9810
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_acceleration")
  set_tests_properties(spray_acceleration_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  foreach(wall_sgs 0 .17)
    add_test(NAME spray_thermal_wall_${wall_sgs} COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --thermal-wall --sgs ${wall_sgs}
      --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_thermal_wall_${wall_sgs}")
    set_tests_properties(spray_thermal_wall_${wall_sgs} PROPERTIES LABELS "spray" TIMEOUT 600)
  endforeach()
  add_test(NAME spray_inflow_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --inflow
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_inflow_restart")
  set_tests_properties(spray_inflow_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  add_test(NAME spray_split_restart COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_spray_etab.py"
    --split --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_split_restart")
  set_tests_properties(spray_split_restart PROPERTIES LABELS "spray" TIMEOUT 900)
  foreach(khrt_variant plain split)
    set(khrt_args)
    if(khrt_variant STREQUAL "split")
      set(khrt_args --split)
    endif()
    add_test(NAME spray_khrt_${khrt_variant}_restart COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_etab.py" --model khrt ${khrt_args}
      --executable "$<TARGET_FILE:pelef_spray_3d>"
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_khrt_${khrt_variant}_restart")
    set_tests_properties(spray_khrt_${khrt_variant}_restart PROPERTIES LABELS "spray" TIMEOUT 900)
  endforeach()
  add_test(NAME spray_etab_restart COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_spray_etab.py"
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_etab_restart")
  set_tests_properties(spray_etab_restart PROPERTIES LABELS "spray" TIMEOUT 900)
  add_test(NAME spray_physical_boundary_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py"
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_physical_boundaries")
  set_tests_properties(spray_physical_boundary_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  add_test(NAME spray_sgs_boundary_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_sgs_boundaries")
  set_tests_properties(spray_sgs_boundary_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  add_executable(test_spray_core spray/test_spray_core.F90)
  target_link_libraries(test_spray_core PRIVATE pelef_spray_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_spray_core PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  foreach(spray_case partition mixture_wall_crossing mixture_wall_rates mixture_wall_exchange mixture_wall_event mixture_terminal mixture_rates mixture_exchange mixture film_evolution wall_event wall_state wall_impact heat_capacity port_states ports surface_cap viscosity density antoine split khrt etab boundaries drag heat evaporation terminal abramzon film_binary injection size_distribution stencil coupling rollback substeps zero invalid
      les_boundary les_algebra les_grid les_zero les_reject)
    add_test(NAME spray_${spray_case} COMMAND test_spray_core ${spray_case})
    set_tests_properties(spray_${spray_case} PROPERTIES LABELS "spray" TIMEOUT 60)
  endforeach()
  add_executable(test_spray_integration spray/test_spray_integration.F90)
  target_link_libraries(test_spray_integration PRIVATE pelef_spray_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_spray_integration PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  foreach(spray_case native coupled rollback checkpoint acceleration materials wall_coupling mixture_wall_coupling mixture_film_grid mixture_grid mixture_reactive)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/spray_${spray_case}")
    add_test(NAME spray_integration_${spray_case} COMMAND test_spray_integration ${spray_case})
    set_tests_properties(spray_integration_${spray_case} PROPERTIES LABELS "spray"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/spray_${spray_case}" TIMEOUT 180)
  endforeach()
  if(PELEF_ENABLE_SUNDIALS)
    add_test(NAME spray_integration_cvode COMMAND test_spray_integration cvode)
    set_tests_properties(spray_integration_cvode PROPERTIES LABELS "spray;cvode" TIMEOUT 180)
  endif()
  add_test(NAME spray_mechanism_portability COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/python/test_mechanism_portability.py")
  set_tests_properties(spray_mechanism_portability PROPERTIES LABELS "spray" TIMEOUT 60)
  add_test(NAME spray_application_cli COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/spray/test_spray_application.py" --exe "$<TARGET_FILE:pelef_spray_3d>")
  set_tests_properties(spray_application_cli PROPERTIES LABELS "spray" TIMEOUT 180)
  if(PELEF_ENABLE_SUNDIALS AND PELEF_ENABLE_CANTERA_REFERENCE)
    add_test(NAME spray_cantera_h2_chemistry COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_cantera.py" --exe "$<TARGET_FILE:pelef_spray_3d>"
      --yaml "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml" --phase ohmech --case h2-chemistry)
    set_tests_properties(spray_cantera_h2_chemistry PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 660)
    if(PELEF_MECHANISM_BUNDLE AND PELEF_MECHANISM_SOURCE)
      file(SHA256 "${PELEF_MECHANISM_SOURCE}" spray_selected_source_sha)
      file(READ "${PROJECT_SOURCE_DIR}/references/spray_methanol.json" spray_reference_json)
      string(JSON spray_methanol_source_sha GET "${spray_reference_json}" files yaml sha256)
      if(spray_selected_source_sha STREQUAL spray_methanol_source_sha)
        foreach(reference_case methanol-chemistry methanol-spray)
          add_test(NAME spray_cantera_${reference_case} COMMAND "${Python3_EXECUTABLE}"
            "${PROJECT_SOURCE_DIR}/tools/check_spray_cantera.py" --exe "$<TARGET_FILE:pelef_spray_3d_selected>"
            --yaml "${PELEF_MECHANISM_SOURCE}" --phase gas --case ${reference_case})
          set_tests_properties(spray_cantera_${reference_case} PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 660)
        endforeach()
        add_test(NAME spray_reactive_cone_restart COMMAND "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/check_spray_cone.py" --exe "$<TARGET_FILE:pelef_spray_3d_selected>"
          --yaml "${PELEF_MECHANISM_SOURCE}")
        set_tests_properties(spray_reactive_cone_restart PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 1200)
        add_test(NAME spray_temporal_refinement COMMAND "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/check_spray_convergence.py" --exe "$<TARGET_FILE:pelef_spray_3d_selected>"
          --yaml "${PELEF_MECHANISM_SOURCE}")
        set_tests_properties(spray_temporal_refinement PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 1800)
      endif()
    endif()
  endif()
endif()
