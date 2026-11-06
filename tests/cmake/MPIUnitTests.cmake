if(PELEF_ENABLE_MPI)
  add_test(NAME application_amr_tree_thermal_slip
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_thermal_slip_3d.py"
      --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" --work "${CMAKE_CURRENT_BINARY_DIR}/thermal_slip_fixed")
  set_tests_properties(application_amr_tree_thermal_slip PROPERTIES TIMEOUT 600 PROCESSORS 3)
  if(TARGET pelef_amr_tree_reactive_3d_selected AND TARGET pelef_mpi_amr_tree_reactive_3d_selected)
    add_test(NAME application_amr_tree_thermal_slip_selected
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_thermal_slip_3d.py"
        --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d_selected>
        --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d_selected>
        --mpiexec "${MPIEXEC_EXECUTABLE}" --work "${CMAKE_CURRENT_BINARY_DIR}/thermal_slip_selected" --selected)
    set_tests_properties(application_amr_tree_thermal_slip_selected PROPERTIES TIMEOUT 600 PROCESSORS 3)
  endif()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME application_amr_tree_mol_default
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_amr_tree_timestep_history.py"
      --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" --work "${CMAKE_CURRENT_BINARY_DIR}" --mol-default)
  set_tests_properties(application_amr_tree_mol_default PROPERTIES TIMEOUT 240 PROCESSORS 2)
  add_test(NAME application_amr_tree_mol_wall_restart
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_amr_tree_timestep_history.py"
      --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" --work "${CMAKE_CURRENT_BINARY_DIR}" --mol --wall-transport)
  set_tests_properties(application_amr_tree_mol_wall_restart PROPERTIES TIMEOUT 240 PROCESSORS 2)
  add_test(NAME application_amr_tree_mol_history
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_amr_tree_timestep_history.py"
      --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" --work "${CMAKE_CURRENT_BINARY_DIR}" --mol)
  set_tests_properties(application_amr_tree_mol_history PROPERTIES TIMEOUT 240 PROCESSORS 2)
  add_test(NAME unit_mpi_amr_tree_mol_stages_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> mol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_tree_mol_stages_3d_2 PROPERTIES TIMEOUT 180 PROCESSORS 2)
  add_test(NAME unit_mpi_amr_tree_rotating_mol_spatial_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_hydro_3d> rotating_mol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_tree_rotating_mol_spatial_3d_2 PROPERTIES TIMEOUT 120 PROCESSORS 2)
  foreach(mol_tree_ranks IN ITEMS 2 4)
    add_test(NAME unit_mpi_amr_tree_mol_spatial_3d_${mol_tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${mol_tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_hydro_3d> mol ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_tree_mol_spatial_3d_${mol_tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${mol_tree_ranks})
  endforeach()
  add_executable(test_mpi_mol_reaction_3d unit/test_mpi_mol_reaction_3d.F90)
  target_link_libraries(test_mpi_mol_reaction_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_mol_reaction_3d_np2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_mol_reaction_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_mol_reaction_3d_np2 PROPERTIES TIMEOUT 120 PROCESSORS 2)
  add_executable(test_mpi_eb_wall_heat_3d unit/test_mpi_eb_wall_heat_3d.F90)
  target_link_libraries(test_mpi_eb_wall_heat_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_wall_heat_3d_np2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_wall_heat_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_wall_heat_3d_np2 PROPERTIES TIMEOUT 120 PROCESSORS 2)
  add_test(NAME mpi_amr_ppm_hierarchy_3d_np8
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_reactive_3d> --ppm-hierarchy-only ${MPIEXEC_POSTFLAGS})
  set_tests_properties(mpi_amr_ppm_hierarchy_3d_np8 PROPERTIES TIMEOUT 240 PROCESSORS 8
    PASS_REGULAR_EXPRESSION "PASS: sparse PPM hierarchy and full split")
  add_test(NAME mpi_amr_ppm_coarse_3d_np8
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_reactive_3d> --ppm-coarse-only ${MPIEXEC_POSTFLAGS})
  set_tests_properties(mpi_amr_ppm_coarse_3d_np8 PROPERTIES TIMEOUT 180 PROCESSORS 8
    PASS_REGULAR_EXPRESSION "PASS: sparse PPM coarse hydro and four-layer halos")
  set(DYNAMIC_AMR_MPI_FLAGS)
  foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND DYNAMIC_AMR_MPI_FLAGS "--mpiexec-preflag=${flag}")
  endforeach()
  foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND DYNAMIC_AMR_MPI_FLAGS "--mpiexec-postflag=${flag}")
  endforeach()
  add_test(NAME mpi_amr_two_level_ppm_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --ppm --riemann-solver pelec --pressure-outlet --reacting --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_two_level_ppm_restart_3d")
  set_tests_properties(mpi_amr_two_level_ppm_restart_3d PROPERTIES TIMEOUT 300 PROCESSORS 4)
  add_test(NAME mpi_amr_two_level_difmag_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --ppm --difmag --riemann-solver pelec --pressure-outlet --reacting --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_two_level_difmag_restart_3d")
  set_tests_properties(mpi_amr_two_level_difmag_restart_3d PROPERTIES TIMEOUT 300 PROCESSORS 4)
  add_test(NAME regression_mpi_dynamic_amr_3d_restart
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_dynamic_amr_3d_restart.py"
      --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 3
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_dynamic_amr_3d_restart")
  add_test(NAME regression_mpi_amr_pressure_outlet_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --pressure-outlet --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_pressure_outlet_restart_3d")
  add_test(NAME regression_mpi_amr_physical_boundary_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_physical_boundary_restart_3d")
  add_test(NAME regression_mpi_amr_fixed_chemistry_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --reacting --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_fixed_chemistry_restart_3d")
  set_tests_properties(regression_mpi_amr_fixed_chemistry_restart_3d PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME regression_mpi_amr_boundary_touching_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --reacting --touching --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_boundary_touching_restart_3d")
  set_tests_properties(regression_mpi_amr_boundary_touching_restart_3d PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME regression_mpi_amr_boundary_touching_dynamic_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --reacting --touching --dynamic-touching --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_boundary_touching_dynamic_restart_3d")
  set_tests_properties(regression_mpi_amr_boundary_touching_dynamic_restart_3d PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME regression_mpi_amr_balanced_boundary_touching_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --reacting --touching --load-balance --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_balanced_boundary_touching_restart_3d")
  set_tests_properties(regression_mpi_amr_balanced_boundary_touching_restart_3d PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME regression_mpi_amr_balanced_boundary_touching_dynamic_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_physical_boundary_restart.py"
      --reacting --touching --load-balance --dynamic-touching --executable "$<TARGET_FILE:pelef_mpi_amr_reactive_3d>"
      --reference-executable "$<TARGET_FILE:pelef_amr_reactive_3d>"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${DYNAMIC_AMR_MPI_FLAGS} --ranks 4 --restart-ranks 2
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_balanced_boundary_touching_dynamic_restart_3d")
  set_tests_properties(regression_mpi_amr_balanced_boundary_touching_dynamic_restart_3d PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  set_tests_properties(regression_mpi_amr_physical_boundary_restart_3d PROPERTIES
    TIMEOUT 600 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  set_tests_properties(regression_mpi_dynamic_amr_3d_restart PROPERTIES
    TIMEOUT 900 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  set(
    MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_3d_transport_restart_0245"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_transport_restart_0245"
  )
  file(MAKE_DIRECTORY
    "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
  )
  add_test(
    NAME regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/checkpoint-np2.log"
      -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
      "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
      fixed_transport_checkpoint_np2 ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    PROPERTIES
      WORKING_DIRECTORY
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
      TIMEOUT 300
  )
  set(MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS)
  foreach(transport_restart_rank IN ITEMS 1 2 4 8)
    set(
      transport_restart_test
      regression_mpi_amr_reactive_3d_transport_restart_0245_np${transport_restart_rank}
    )
    add_test(
      NAME ${transport_restart_test}
      COMMAND "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np${transport_restart_rank}.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        ${transport_restart_rank} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
        "fixed_transport_restart_mpi_np${transport_restart_rank}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${transport_restart_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${transport_restart_rank}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${transport_restart_rank}"
        DEPENDS
          regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
        TIMEOUT 300
    )
    list(APPEND MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS
      ${transport_restart_test})
  endforeach()
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_mpi_amr_reactive_3d_transport_restart_0245_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_reference_0245"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np1"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np2"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np4"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np8"
      --checkpoint
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
      --expected-checkpoint-sha256
        "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256}"
      --expected-coarse-sha256
        "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256}"
      --expected-fine-sha256
        "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256}"
      --expected-schema 4
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --reference-log
        "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/reference.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np1.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np2.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np4.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np8.log"
  )
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_exact
    PROPERTIES DEPENDS
      "regression_amr_reactive_3d_transport_restart_0245_reference;regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2;${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS}"
  )
  add_test(
    NAME regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
      "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
  )
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    PROPERTIES DEPENDS
      "regression_amr_reactive_3d_transport_restart_0245_checkpoint;regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2"
  )

  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/checkpoint-np2.log"
      -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
      selected_transport_checkpoint_np2 ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
      TIMEOUT 600
  )
  set(SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS)
  foreach(selected_transport_restart_rank IN ITEMS 1 2 4 8)
    set(
      selected_transport_restart_test
      regression_selected_mpi_amr_reactive_3d_transport_restart_0245_np${selected_transport_restart_rank}
    )
    add_test(
      NAME ${selected_transport_restart_test}
      COMMAND "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np${selected_transport_restart_rank}.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        ${selected_transport_restart_rank} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
        "selected_transport_restart_mpi_np${selected_transport_restart_rank}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_transport_restart_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_transport_restart_rank}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_transport_restart_rank}"
        DEPENDS
          regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
        TIMEOUT 600
    )
    list(APPEND SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS
      ${selected_transport_restart_test})
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_restart_0245_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_reference_0245"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np8"
      --checkpoint
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
      --expected-checkpoint-sha256
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256}"
      --expected-coarse-sha256
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256}"
      --expected-fine-sha256
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256}"
      --expected-schema 5
      --expected-bundle-sha256
        f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
      --expected-integrator implicit
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --reference-log
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/reference.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np1.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np2.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np4.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np8.log"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_transport_restart_0245_reference;regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2;${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS}"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
      "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint;regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2"
  )
endif()
file(MAKE_DIRECTORY "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_restart_3d_reference
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d_restart/reference.nml"
)
add_test(
  NAME regression_reactive_eb_restart_3d_checkpoint_stop
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_restart_3d_restart
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d_restart/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_reference
  regression_reactive_eb_restart_3d_checkpoint_stop
  regression_reactive_eb_restart_3d_restart
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_restart_3d_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Operator sequence: R-T-H-T-R"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    DEPENDS regression_reactive_eb_restart_3d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_restart_3d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_3d_restart.py"
    --reference
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart_reference.csv"
    --stopped
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart_stopped.csv"
    --restarted
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart.csv"
    --checkpoint
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart.chk"
    --reference-log "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reference.log"
    --stopped-log "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/stopped.log"
    --restarted-log "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/restarted.log"
    --expected-checkpoint-sha256
      "${REACTIVE_EB_RESTART_3D_CHECKPOINT_SHA256}"
    --expected-stopped-sha256
      "${REACTIVE_EB_RESTART_3D_STOPPED_SHA256}"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_check
  PROPERTIES
    DEPENDS
      "regression_reactive_eb_restart_3d_reference;regression_reactive_eb_restart_3d_checkpoint_stop;regression_reactive_eb_restart_3d_restart"
    TIMEOUT 60
)


if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_distribution_3d unit/test_mpi_amr_patch_tree_distribution_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_distribution_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_distribution_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_distribution_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_distribution_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_reflux_budget_3d unit/test_mpi_amr_reflux_budget_3d.F90)
  target_link_libraries(test_mpi_amr_reflux_budget_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_reflux_budget_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_reflux_budget_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_reflux_budget_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_executable(test_mpi_eb_interface_source_3d unit/test_mpi_eb_interface_source_3d.F90)
  target_link_libraries(test_mpi_eb_interface_source_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_interface_source_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_interface_source_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_interface_source_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_executable(test_mpi_amr_patch_tree_cells_3d unit/test_mpi_amr_patch_tree_cells_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_cells_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_cells_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_cells_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_cells_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_ghost_3d unit/test_mpi_amr_patch_tree_ghost_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_ghost_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_ghost_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_ghost_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_ghost_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_restriction_3d unit/test_mpi_amr_patch_tree_restriction_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_restriction_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_restriction_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_restriction_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_restriction_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_flux_peers_3d unit/test_mpi_amr_patch_tree_flux_peers_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_flux_peers_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_flux_peers_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_flux_peers_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_flux_peers_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_flux_restriction_3d unit/test_mpi_amr_patch_tree_flux_restriction_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_flux_restriction_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_flux_restriction_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_flux_restriction_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_flux_restriction_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_thermo_3d unit/test_mpi_amr_patch_tree_thermo_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_thermo_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_thermo_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_thermo_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_thermo_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_hydro_3d unit/test_mpi_amr_patch_tree_hydro_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_hydro_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_amr_patch_tree_sdc_hydro_3d_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_hydro_3d> sdc ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_sdc_hydro_3d_4 PROPERTIES PROCESSORS 4 TIMEOUT 180)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_hydro_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_hydro_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_hydro_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_chemistry_3d unit/test_mpi_amr_patch_tree_chemistry_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_chemistry_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_amr_patch_tree_sdc_reaction_3d_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_chemistry_3d> --sdc-only ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_sdc_reaction_3d_4 PROPERTIES PROCESSORS 4 TIMEOUT 180)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_chemistry_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_chemistry_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_chemistry_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME unit_mpi_amr_patch_tree_pressure_strang_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_chemistry_3d> --pressure-only ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_pressure_strang_3d_2 PROPERTIES
    TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")

endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_transport_3d unit/test_mpi_amr_patch_tree_transport_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_transport_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_amr_patch_tree_constant_sdc_transport_3d_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_transport_3d> ${MPIEXEC_POSTFLAGS} constant_sdc)
  set_tests_properties(unit_mpi_amr_patch_tree_constant_sdc_transport_3d_4 PROPERTIES PROCESSORS 4 TIMEOUT 180)
  add_test(NAME unit_mpi_amr_patch_tree_sdc_transport_3d_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_transport_3d> sdc ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_sdc_transport_3d_4 PROPERTIES PROCESSORS 4 TIMEOUT 180)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_transport_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_transport_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_transport_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_split_3d unit/test_mpi_amr_patch_tree_split_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_split_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_amr_patch_tree_sdc_3d_4
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> sdc ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_sdc_3d_4 PROPERTIES PROCESSORS 4 TIMEOUT 300)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_split_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_split_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_application_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_application_3d.py"
      --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_application_3d")
  set_tests_properties(regression_mpi_amr_tree_application_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
    add_test(NAME regression_mpi_amr_tree_selected_application_3d
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_application_3d.py"
        --selected --serial $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
        --mpi $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
        --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
        --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_selected_application_3d")
    set_tests_properties(regression_mpi_amr_tree_selected_application_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_restart_3d")
  set_tests_properties(regression_mpi_amr_tree_restart_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
  add_test(NAME regression_mpi_amr_tree_mach_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_mach_restart_3d" --scheduled --mach-inlet)
  set_tests_properties(regression_mpi_amr_tree_mach_restart_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
  add_test(NAME regression_mpi_amr_tree_difmag_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_difmag_restart_3d" --subcycling --difmag)
  set_tests_properties(regression_mpi_amr_tree_difmag_restart_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
  add_test(NAME regression_mpi_amr_tree_pressure_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_pressure_restart_3d" --scheduled --pressure-outlet)
  set_tests_properties(regression_mpi_amr_tree_pressure_restart_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_remap_3d unit/test_mpi_amr_patch_tree_remap_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_remap_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_remap_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_remap_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_remap_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_schedule_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --scheduled --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_schedule_3d")
  set_tests_properties(regression_mpi_amr_tree_schedule_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_tagging_3d unit/test_mpi_amr_patch_tree_tagging_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_tagging_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_tagging_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_tagging_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_tagging_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_adaptation_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_adaptation_3d.py"
      --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_adaptation_3d")
  set_tests_properties(regression_mpi_amr_tree_adaptation_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_hydro_level_3d unit/test_mpi_amr_patch_tree_hydro_level_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_hydro_level_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_hydro_level_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_hydro_level_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_hydro_level_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_transport_level_3d unit/test_mpi_amr_patch_tree_transport_level_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_transport_level_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_transport_level_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_transport_level_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_transport_level_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_subcycle_3d unit/test_mpi_amr_patch_tree_subcycle_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_subcycle_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_subcycle_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_subcycle_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_subcycle_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_transport_subcycle_3d unit/test_mpi_amr_patch_tree_transport_subcycle_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_transport_subcycle_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_transport_subcycle_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_transport_subcycle_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_transport_subcycle_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_subcycle_split_3d unit/test_mpi_amr_patch_tree_subcycle_split_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_subcycle_split_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_subcycle_split_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_subcycle_split_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_subcycle_split_3d_${tree_ranks} PROPERTIES
      TIMEOUT 180 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_subcycle_application_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --subcycling --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_subcycle_application_3d")
  set_tests_properties(regression_mpi_amr_tree_subcycle_application_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_subcycle_schedule_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --subcycling --scheduled --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_subcycle_schedule_3d")
  set_tests_properties(regression_mpi_amr_tree_subcycle_schedule_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_subcycle_adaptation_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_adaptation_3d.py"
      --subcycling --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_subcycle_adaptation_3d")
  set_tests_properties(regression_mpi_amr_tree_subcycle_adaptation_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_weighted_adaptation_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_adaptation_3d.py"
      --weighted --subcycling --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_weighted_adaptation_3d")
  set_tests_properties(regression_mpi_amr_tree_weighted_adaptation_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_csv_3d unit/test_mpi_amr_patch_tree_csv_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_csv_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_csv_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_csv_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_csv_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_amr_patch_tree_checkpoint_3d unit/test_mpi_amr_patch_tree_checkpoint_3d.F90)
  target_link_libraries(test_mpi_amr_patch_tree_checkpoint_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_patch_tree_checkpoint_3d_${tree_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_checkpoint_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_patch_tree_checkpoint_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME application_amr_tree_timestep_history
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/application/check_amr_tree_timestep_history.py
      --serial $<TARGET_FILE:pelef_amr_tree_reactive_3d> --mpi $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --mpiexec ${MPIEXEC_EXECUTABLE} --work ${CMAKE_CURRENT_BINARY_DIR}/amr_timestep_history)
  set_tests_properties(application_amr_tree_timestep_history PROPERTIES TIMEOUT 240 PROCESSORS 2 LABELS "amr;restart;mpi")
  add_test(NAME regression_mpi_amr_tree_indexed_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --indexed --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_indexed_restart_3d")
  set_tests_properties(regression_mpi_amr_tree_indexed_restart_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_runtime_balance_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --rebalance --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_runtime_balance_3d")
  set_tests_properties(regression_mpi_amr_tree_runtime_balance_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_recursive_runtime_balance_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_restart_3d.py"
      --rebalance --subcycling --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_recursive_runtime_balance_3d")
  set_tests_properties(regression_mpi_amr_tree_recursive_runtime_balance_3d PROPERTIES TIMEOUT 600 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_dynamic_split_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_dynamic_split_3d.py"
      --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_dynamic_split_3d")
  set_tests_properties(regression_mpi_amr_tree_dynamic_split_3d PROPERTIES TIMEOUT 900 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_amr_tree_adaptive_split_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_amr_tree_dynamic_split_3d.py"
      --adaptive --fixed $<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>
      --selected $<TARGET_FILE:test_selected_elementary_mpi_amr_tree_3d>
      --serial-fixed $<TARGET_FILE:pelef_amr_tree_reactive_3d>
      --serial-selected $<TARGET_FILE:test_selected_elementary_amr_tree_3d>
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --workdir "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_tree_adaptive_split_3d")
  set_tests_properties(regression_mpi_amr_tree_adaptive_split_3d PROPERTIES TIMEOUT 900 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_transfer_3d unit/test_mpi_eb_tree_transfer_3d.F90)
  target_link_libraries(test_mpi_eb_tree_transfer_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_transfer_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transfer_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_transfer_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_halo_3d unit/test_mpi_eb_tree_halo_3d.F90)
  target_link_libraries(test_mpi_eb_tree_halo_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_halo_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_halo_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_halo_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_flux_3d unit/test_mpi_eb_tree_flux_3d.F90)
  target_link_libraries(test_mpi_eb_tree_flux_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_native_adjacent_query_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} native_query)
  set_tests_properties(unit_mpi_eb_native_adjacent_query_3d_3 PROPERTIES
    TIMEOUT 60 PROCESSORS 3 LABELS "mpi;eb;mol;dual_energy;source_funding"
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_flux_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS})
    add_test(NAME unit_mpi_eb_tree_periodic_flux_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> periodic ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_periodic_flux_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
    set_tests_properties(unit_mpi_eb_tree_flux_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME unit_mpi_eb_tree_rotating_flux_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> rotating ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_rotating_flux_3d_3 PROPERTIES
    TIMEOUT 60 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_hydro_3d unit/test_mpi_eb_tree_hydro_3d.F90)
  target_link_libraries(test_mpi_eb_tree_hydro_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_hydro_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_hydro_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  foreach(mol_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_mol_3d_${mol_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${mol_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} mol)
    set_tests_properties(unit_mpi_eb_tree_mol_3d_${mol_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${mol_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
    add_test(NAME unit_mpi_eb_tree_mol_reflux_3d_${mol_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${mol_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} mol_reflux)
    set_tests_properties(unit_mpi_eb_tree_mol_reflux_3d_${mol_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${mol_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_transport_3d unit/test_mpi_eb_tree_transport_3d.F90)
  target_link_libraries(test_mpi_eb_tree_transport_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_transport_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_transport_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME unit_mpi_eb_tree_rotating_timestep_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} rotating_dt)
  set_tests_properties(unit_mpi_eb_tree_rotating_timestep_3d_3 PROPERTIES
    TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_full_3d unit/test_mpi_eb_tree_full_3d.F90)
  add_test(NAME unit_mpi_serial_mol_tree_3d_np2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> mol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_serial_mol_tree_3d_np2 PROPERTIES TIMEOUT 240 PROCESSORS 2)
  add_test(NAME unit_mpi_rotating_mol_tree_3d_np2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> rotating_mol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_rotating_mol_tree_3d_np2 PROPERTIES TIMEOUT 300 PROCESSORS 2)

  target_link_libraries(test_mpi_eb_tree_full_3d PRIVATE pelef_mpi_support)
  foreach(route_mode IN ITEMS route_history route_first routed_replay routed_refusal routed_predictor redist_history flux_replay flux_refusal state_replay state_budget stored_reactor held_auxiliary held_diffusion held_policy held_recursive_policy held_chem_policy held_prepared held_prepared_recursive held_recursive_chemistry held_prepared_serial held_natural_replay held_prepared_parent held_prepared_natural held_scalar_replay held_prepared_scalar_replay held_native_funding held_native_material)
    add_test(NAME unit_mpi_eb_tree_${route_mode}_3d_3
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${route_mode} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_${route_mode}_3d_3 PROPERTIES
      TIMEOUT 180 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  set_tests_properties(unit_mpi_eb_tree_routed_replay_3d_3 unit_mpi_eb_tree_routed_refusal_3d_3
    unit_mpi_eb_tree_routed_predictor_3d_3
    unit_mpi_eb_tree_flux_replay_3d_3 unit_mpi_eb_tree_flux_refusal_3d_3 PROPERTIES TIMEOUT 300)
  set_tests_properties(unit_mpi_eb_tree_state_replay_3d_3 PROPERTIES TIMEOUT 1800)
  set_tests_properties(unit_mpi_eb_tree_held_prepared_3d_3 PROPERTIES TIMEOUT 180)
  set_tests_properties(unit_mpi_eb_tree_held_prepared_recursive_3d_3 PROPERTIES TIMEOUT 240)
  set_tests_properties(unit_mpi_eb_tree_held_recursive_chemistry_3d_3 PROPERTIES TIMEOUT 300)
  set_tests_properties(unit_mpi_eb_tree_held_prepared_serial_3d_3 PROPERTIES TIMEOUT 300)
  set_tests_properties(unit_mpi_eb_tree_held_prepared_parent_3d_3 PROPERTIES TIMEOUT 300)
  # Native prepared parent-window trials rebuild paired spatial operators at
  # all descendant clocks. Observed stack progression reaches the twelve
  # refinements; include serial, owner/idle layouts and atomic refusal checks.
  set_tests_properties(unit_mpi_eb_tree_held_prepared_natural_3d_3 PROPERTIES TIMEOUT 1800)
  # The independent held-energy inventory follows the complete prepared
  # parent-window replay. The same fixture took 2178 seconds under parallel
  # development load; preserve its numerical inputs and budget tolerance.
  set_tests_properties(unit_mpi_eb_tree_held_prepared_scalar_replay_3d_3 PROPERTIES TIMEOUT 3600)
  set_tests_properties(unit_mpi_eb_tree_held_natural_replay_3d_3 unit_mpi_eb_tree_held_scalar_replay_3d_3
    PROPERTIES TIMEOUT 300)
  set_tests_properties(unit_mpi_eb_tree_held_native_funding_3d_3 unit_mpi_eb_tree_held_native_material_3d_3
    PROPERTIES TIMEOUT 1800)
  set_tests_properties(unit_mpi_eb_tree_state_budget_3d_3 PROPERTIES TIMEOUT 300)
  # The coupled held-energy fixture performs native forced chemistry in every
  # cell twice, plus serial and owner/idle layouts; keep its numerical inputs.
  set_tests_properties(unit_mpi_eb_tree_held_chem_policy_3d_3 PROPERTIES TIMEOUT 1200)
  add_test(NAME unit_mpi_eb_tree_recursive_history_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> recursive_history ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_recursive_history_3d_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_interval_history_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> interval_history ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_interval_history_3d_3 PROPERTIES
    TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_rotating_reaction_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> rotating_reaction ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_rotating_reaction_3d_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_sdc_reactor_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} reactor)
  set_tests_properties(unit_mpi_eb_tree_sdc_reactor_3d_2 PROPERTIES
    TIMEOUT 300 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_sdc_redist_level_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} redist_level)
  set_tests_properties(unit_mpi_eb_tree_sdc_redist_level_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_sdc_redist_recursive_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} redist_recursive)
  set_tests_properties(unit_mpi_eb_tree_sdc_redist_recursive_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_sdc_subcycle_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} subcycle)
  set_tests_properties(unit_mpi_eb_tree_sdc_subcycle_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  foreach(tree_ranks IN ITEMS 2)
    add_test(NAME unit_mpi_eb_tree_full_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_full_3d_${tree_ranks} PROPERTIES
      TIMEOUT 180 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME unit_mpi_eb_tree_sdc_state_redist_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} state_redist)
  set_tests_properties(unit_mpi_eb_tree_sdc_state_redist_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_sdc_state_recursive_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} state_recursive)
  set_tests_properties(unit_mpi_eb_tree_sdc_state_recursive_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  foreach(state_redist_order IN ITEMS 2 4)
    foreach(state_redist_mode IN ITEMS state_redist state_recursive)
      set(state_redist_test unit_mpi_eb_tree_sdc_${state_redist_mode}_order${state_redist_order}_3d_2)
      add_test(NAME ${state_redist_test}
        COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} ${state_redist_mode} ${state_redist_order})
      set_tests_properties(${state_redist_test} PROPERTIES
        TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
    endforeach()
  endforeach()
  add_test(NAME unit_mpi_eb_tree_sdc_diffusion_3d_2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} diffusion)
  set_tests_properties(unit_mpi_eb_tree_sdc_diffusion_3d_2 PROPERTIES TIMEOUT 300 PROCESSORS 2)
  add_test(NAME unit_mpi_eb_tree_sdc_seed_aux_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} seed_aux)
  set_tests_properties(unit_mpi_eb_tree_sdc_seed_aux_3d_2 PROPERTIES
    TIMEOUT 300 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_sdc_aux_source_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} aux_source)
  set_tests_properties(unit_mpi_eb_tree_sdc_aux_source_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")

endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_eb_tree_application_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/input.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_application_3d")
  set_tests_properties(regression_mpi_eb_tree_application_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_eb_tree_auto_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/input.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_auto_3d" --automatic)
  set_tests_properties(regression_mpi_eb_tree_auto_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_remap_3d unit/test_mpi_eb_tree_remap_3d.F90)
  target_link_libraries(test_mpi_eb_tree_remap_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_remap_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_remap_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_remap_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME regression_mpi_eb_tree_scheduled_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/scheduled.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_scheduled_3d" --scheduled)
  set_tests_properties(regression_mpi_eb_tree_scheduled_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)

  add_test(NAME regression_mpi_eb_tree_rebalance_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py" --rebalance
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/scheduled.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_rebalance_3d" --scheduled --subcycling)
  set_tests_properties(regression_mpi_eb_tree_rebalance_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_tagging_3d unit/test_mpi_eb_tree_tagging_3d.F90)
  target_link_libraries(test_mpi_eb_tree_tagging_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_tagging_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_tagging_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_tagging_3d_${tree_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME regression_mpi_eb_tree_adaptive_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/adaptive.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_adaptive_3d" --adaptive)
  set_tests_properties(regression_mpi_eb_tree_adaptive_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_subcycle_hydro_3d unit/test_mpi_eb_tree_subcycle_hydro_3d.F90)
  target_link_libraries(test_mpi_eb_tree_subcycle_hydro_3d PRIVATE pelef_mpi_support)
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_subcycle_hydro_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_subcycle_hydro_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_subcycle_hydro_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  foreach(operator IN ITEMS transport full)
    add_executable(test_mpi_eb_tree_subcycle_${operator}_3d unit/test_mpi_eb_tree_subcycle_${operator}_3d.F90)
    target_link_libraries(test_mpi_eb_tree_subcycle_${operator}_3d PRIVATE pelef_mpi_support)
  endforeach()
  foreach(tree_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_eb_tree_subcycle_transport_3d_${tree_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${tree_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_subcycle_transport_3d> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_tree_subcycle_transport_3d_${tree_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${tree_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  foreach(mach_mode IN ITEMS mach mach-common)
    add_test(NAME unit_mpi_eb_tree_${mach_mode}_full_3d_2
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_subcycle_full_3d> ${MPIEXEC_POSTFLAGS} ${mach_mode})
    set_tests_properties(unit_mpi_eb_tree_${mach_mode}_full_3d_2 PROPERTIES
      TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME unit_mpi_eb_tree_subcycle_full_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_subcycle_full_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_subcycle_full_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_mpi_eb_tree_mach_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/input.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_mach_3d" --subcycling --automatic --mach)
  set_tests_properties(regression_mpi_eb_tree_mach_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
  foreach(mode IN ITEMS scheduled automatic)
    if(mode STREQUAL "scheduled")
      set(subcycle_input scheduled.nml)
    else()
      set(subcycle_input input.nml)
    endif()
    add_test(NAME regression_mpi_eb_tree_subcycle_${mode}_3d
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
        --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/${subcycle_input}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/mpi_eb_tree_subcycle_${mode}_3d" --subcycling --${mode})
    set_tests_properties(regression_mpi_eb_tree_subcycle_${mode}_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_csv_3d unit/test_mpi_eb_tree_csv_3d.F90)
  target_link_libraries(test_mpi_eb_tree_csv_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_tree_csv_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_csv_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_csv_3d_3 PROPERTIES TIMEOUT 60 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME unit_mpi_amr_patch_tree_ppm_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} ppm)
  set_tests_properties(unit_mpi_amr_patch_tree_ppm_3d_2 PROPERTIES
    TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME unit_mpi_amr_patch_tree_ppm_subcycle_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_subcycle_3d> ${MPIEXEC_POSTFLAGS} ppm)
  set_tests_properties(unit_mpi_amr_patch_tree_ppm_subcycle_3d_2 PROPERTIES
    TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_amr_patch_tree_split_3d)
  foreach(reflux_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_amr_tree_reflux_recovery_3d_${reflux_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${reflux_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> reflux_budget ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_tree_reflux_recovery_3d_${reflux_ranks} PROPERTIES
      TIMEOUT 120 PROCESSORS ${reflux_ranks})
  endforeach()
  add_test(NAME unit_mpi_amr_patch_tree_mol_conserve_3d_2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> mol_conserve ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_mol_conserve_3d_2 PROPERTIES TIMEOUT 300 PROCESSORS 2)
  add_test(NAME unit_mpi_amr_patch_tree_rotating_mol_3d_2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> rotating_mol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_rotating_mol_3d_2 PROPERTIES TIMEOUT 300 PROCESSORS 2)
  add_test(NAME unit_mpi_amr_patch_tree_rotating_mol_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> rotating_mol ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_patch_tree_rotating_mol_3d_3 PROPERTIES TIMEOUT 300 PROCESSORS 3)
  foreach(variant split subcycle_split)
    add_test(NAME unit_mpi_amr_patch_tree_ppm_${variant}_3d_2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_${variant}_3d> ${MPIEXEC_POSTFLAGS} ppm)
    set_tests_properties(unit_mpi_amr_patch_tree_ppm_${variant}_3d_2 PROPERTIES TIMEOUT 300 PROCESSORS 2)
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  foreach(variant hydro subcycle)
    add_test(NAME unit_mpi_amr_tree_viscosity_${variant}_3d_2
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_${variant}_3d> ${MPIEXEC_POSTFLAGS} viscosity)
    set_tests_properties(unit_mpi_amr_tree_viscosity_${variant}_3d_2 PROPERTIES
      TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME mpi_amr_viscosity_hierarchy_3d_np8
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 8 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_reactive_3d> --viscosity-hierarchy-only ${MPIEXEC_POSTFLAGS})
  set_tests_properties(mpi_amr_viscosity_hierarchy_3d_np8 PROPERTIES TIMEOUT 240 PROCESSORS 8
    PASS_REGULAR_EXPRESSION "PASS: sparse PCM PLM viscosity hierarchy")
endif()

if(TARGET test_mpi_eb_tree_hydro_3d)
  foreach(halo_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_mol_spatial_halo_3d_${halo_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${halo_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} mol_halo)
    set_tests_properties(unit_mpi_eb_tree_mol_spatial_halo_3d_${halo_ranks}
      PROPERTIES PROCESSORS ${halo_ranks} TIMEOUT 300)
  endforeach()
  foreach(allowance_ranks 2 4)
    add_test(NAME unit_mpi_eb_flux_allowance_3d_${allowance_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${allowance_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} flux_budget)
    set_tests_properties(unit_mpi_eb_flux_allowance_3d_${allowance_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${allowance_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME unit_mpi_eb_tree_viscosity_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} viscosity)
  set_tests_properties(unit_mpi_eb_tree_viscosity_3d_3 PROPERTIES
    TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_tree_subcycle_full_3d)
  add_test(NAME unit_mpi_eb_tree_viscosity_subcycle_full_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_subcycle_full_3d> ${MPIEXEC_POSTFLAGS} mach viscosity)
  set_tests_properties(unit_mpi_eb_tree_viscosity_subcycle_full_3d_3 PROPERTIES
    TIMEOUT 300 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  foreach(variant hydro subcycle)
    add_test(NAME unit_mpi_amr_tree_default_${variant}_3d_2
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_${variant}_3d> ${MPIEXEC_POSTFLAGS} default)
    set_tests_properties(unit_mpi_amr_tree_default_${variant}_3d_2 PROPERTIES
      TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME mpi_amr_default_hierarchy_3d_np8
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 8 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_reactive_3d> --default-hierarchy-only ${MPIEXEC_POSTFLAGS})
  set_tests_properties(mpi_amr_default_hierarchy_3d_np8 PROPERTIES TIMEOUT 240 PROCESSORS 8
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME unit_mpi_eb_tree_default_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} default)
  set_tests_properties(unit_mpi_eb_tree_default_3d_3 PROPERTIES TIMEOUT 120 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_tree_hydro_3d)
  add_executable(test_mpi_eb_native_mol_partial_subparent_3d unit/test_mpi_eb_native_mol_partial_subparent_3d.F90)
  target_link_libraries(test_mpi_eb_native_mol_partial_subparent_3d PRIVATE pelef_mpi_support)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_native_mol_partial_subparent_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_mpi_eb_native_mol_partial_subparent_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_partial_subparent_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_native_mol_partial_subparent_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")

  add_executable(test_mpi_eb_native_mol_partial_fine_3d unit/test_mpi_eb_native_mol_partial_fine_3d.F90)
  target_link_libraries(test_mpi_eb_native_mol_partial_fine_3d PRIVATE pelef_mpi_support)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_native_mol_partial_fine_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_mpi_eb_native_mol_partial_fine_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_partial_fine_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_native_mol_partial_fine_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")

  add_executable(test_mpi_eb_native_mol_partial_descendants_3d unit/test_mpi_eb_native_mol_partial_descendants_3d.F90)
  target_link_libraries(test_mpi_eb_native_mol_partial_descendants_3d PRIVATE pelef_mpi_support)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_native_mol_partial_descendants_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_mpi_eb_native_mol_partial_descendants_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_partial_descendants_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_native_mol_partial_descendants_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")

  add_executable(test_mpi_eb_native_mol_partial_recovery_3d unit/test_mpi_eb_native_mol_partial_recovery_3d.F90)
  target_link_libraries(test_mpi_eb_native_mol_partial_recovery_3d PRIVATE pelef_mpi_support)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_native_mol_partial_recovery_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_mpi_eb_native_mol_partial_recovery_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_partial_recovery_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_native_mol_partial_recovery_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 300
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")

  add_executable(test_mpi_eb_native_mol_future_source_3d unit/test_mpi_eb_native_mol_future_source_3d.F90)
  target_link_libraries(test_mpi_eb_native_mol_future_source_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_native_mol_future_source_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_native_mol_future_source_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_funding" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_source_prefix_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} prefix)
  set_tests_properties(unit_mpi_eb_native_mol_source_prefix_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 180
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_prefix_consumer_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} consumer)
  set_tests_properties(unit_mpi_eb_native_mol_prefix_consumer_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 240
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_momentum_consumer_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} consumer momentum)
  set_tests_properties(unit_mpi_eb_native_mol_momentum_consumer_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 240
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_recursive_prefix_retry_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} retry)
  set_tests_properties(unit_mpi_eb_native_mol_recursive_prefix_retry_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_momentum_retry_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} retry momentum)
  set_tests_properties(unit_mpi_eb_native_mol_momentum_retry_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_material_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} replay)
  add_test(NAME unit_mpi_eb_native_mol_receipt_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} receipt)
  set_tests_properties(unit_mpi_eb_native_mol_material_replay_3d_3 unit_mpi_eb_native_mol_receipt_replay_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_funding;replay" PROCESSORS 3 TIMEOUT 1200
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_native_mol_future_source_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()

  add_executable(test_mpi_eb_native_mol_replay_3d unit/test_mpi_eb_native_mol_replay_3d.F90)
  target_link_libraries(test_mpi_eb_native_mol_replay_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_native_mol_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_replay_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_native_mol_replay_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;mol;dual_energy;replay"
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_momentum_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_replay_3d> ${MPIEXEC_POSTFLAGS} momentum)
  set_tests_properties(unit_mpi_eb_native_mol_momentum_replay_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;dual_energy;source_prefix" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_weighted_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_replay_3d> ${MPIEXEC_POSTFLAGS} weighted)
  set_tests_properties(unit_mpi_eb_native_mol_weighted_replay_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;mol;dual_energy;replay"
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  add_test(NAME unit_mpi_eb_native_mol_subparent_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_replay_3d> ${MPIEXEC_POSTFLAGS} subparent)
  set_tests_properties(unit_mpi_eb_native_mol_subparent_replay_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;mol;dual_energy;replay"
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_native_mol_replay_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_executable(test_mpi_eb_held_reset_ledger_3d unit/test_mpi_eb_held_reset_ledger_3d.F90)
  target_link_libraries(test_mpi_eb_held_reset_ledger_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_held_reset_ledger_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_held_reset_ledger_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy")
  add_test(NAME unit_mpi_eb_held_register_closure_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> registers ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_held_register_closure_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;state_redist")
  add_test(NAME unit_mpi_eb_pinned_eos_source_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_pinned_eos_source_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;eos")
  add_test(NAME unit_mpi_eb_mol_temperature_seed_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source corrector_material reactor seed ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_mol_temperature_seed_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;eos")
  foreach(source_kind IN ITEMS material covered_material corrector_material covered_corrector_material
      diffusive_corrector_material)
    add_test(NAME unit_mpi_eb_${source_kind}_source_3d_3
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source ${source_kind} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_${source_kind}_source_3d_3 PROPERTIES
      PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;source_recovery")
  endforeach()
  add_test(NAME unit_mpi_eb_scalar_closure_reset_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source corrector_material disabled total_closure ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_scalar_closure_reset_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;scalar_reset")
  add_test(NAME unit_mpi_eb_reactor_corrector_material_source_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source corrector_material reactor ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_reactor_corrector_material_source_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;source_recovery;reactor_admission;inert")
  add_test(NAME unit_mpi_eb_reactor_covered_corrector_material_source_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source covered_corrector_material reactor ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_reactor_covered_corrector_material_source_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;source_recovery;reactor_admission;inert")
  add_test(NAME unit_mpi_eb_manufactured_caloric_provider_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source caloric_provider reactor ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_manufactured_caloric_provider_3d_3 PROPERTIES
    LABELS "unit;mpi;eb;dual-energy;source;reactor;manufactured" PROCESSORS 3 TIMEOUT 300)
  foreach(source_kind IN ITEMS corrector_material covered_corrector_material diffusive_corrector_material)
    add_test(NAME unit_mpi_eb_kinetic_${source_kind}_source_3d_3
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source ${source_kind} kinetic ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_eb_kinetic_${source_kind}_source_3d_3 PROPERTIES
      PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;source_recovery;kinetics")
  endforeach()
  add_test(NAME unit_mpi_eb_builtin_kinetic_corrector_source_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_held_reset_ledger_3d> eos_source corrector_material builtin ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_builtin_kinetic_corrector_source_3d_3 PROPERTIES
    PROCESSORS 3 TIMEOUT 300 LABELS "mpi;eb;sdc;dual_energy;source_recovery;kinetics;reactor_admission")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_mpi_eb_held_reset_ledger_3d PRIVATE
      -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_executable(test_mpi_eb_godunov_held_3d unit/test_mpi_eb_godunov_held_3d.F90)
  target_link_libraries(test_mpi_eb_godunov_held_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_godunov_held_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_held_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 300)
  add_test(NAME unit_mpi_eb_godunov_held_weight_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} weights 1 1 -1)
  set_tests_properties(unit_mpi_eb_godunov_held_weight_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 300
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_godunov_held_sdc_weight_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} common_weights 1 0 -1)
  add_test(NAME unit_mpi_eb_godunov_held_sdc_weight_reactor_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} common_weights_react 1 1 1)
  set_tests_properties(unit_mpi_eb_godunov_held_sdc_weight_3d_3 unit_mpi_eb_godunov_held_sdc_weight_reactor_3d_3
    PROPERTIES PROCESSORS 3 TIMEOUT 300 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_godunov_held_history_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} history_diff 1 1 1)
  set_tests_properties(unit_mpi_eb_godunov_held_history_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 300
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_godunov_held_clocks_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} clocks_diff 1 0 1)
  set_tests_properties(unit_mpi_eb_godunov_held_clocks_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 300
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_godunov_held_common_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} common 1 0 -1)
  add_test(NAME unit_mpi_eb_godunov_held_selected_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} common 1 1 -1)
  add_test(NAME unit_mpi_eb_godunov_held_recursive_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} recursive 1 0 -1)
  add_test(NAME unit_mpi_eb_godunov_scalar_clock_reset_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} recursive_reset 1 0 -1)
  set_tests_properties(unit_mpi_eb_godunov_scalar_clock_reset_3d_3
    PROPERTIES PROCESSORS 3 TIMEOUT 300 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_godunov_held_recursive_diffusion_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} recursive_diff 1 0 1)
  set_tests_properties(unit_mpi_eb_godunov_held_recursive_3d_3 unit_mpi_eb_godunov_held_recursive_diffusion_3d_3
    PROPERTIES PROCESSORS 3 TIMEOUT 300 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_godunov_held_diffusion_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> ${MPIEXEC_POSTFLAGS} common_diff 1 0 1)
  set_tests_properties(unit_mpi_eb_godunov_held_common_3d_3 unit_mpi_eb_godunov_held_selected_3d_3
    unit_mpi_eb_godunov_held_diffusion_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 300
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  add_test(NAME unit_mpi_eb_tree_rotating_mol_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} rotating_mol)
  set_tests_properties(unit_mpi_eb_tree_rotating_mol_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 180)
  add_test(NAME unit_mpi_eb_tree_ctu_hydro_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} ctu)
  set_tests_properties(unit_mpi_eb_tree_ctu_hydro_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 180)
endif()

if(TARGET test_mpi_eb_tree_subcycle_hydro_3d)
  add_test(NAME unit_mpi_eb_tree_ctu_subcycle_hydro_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_subcycle_hydro_3d> ${MPIEXEC_POSTFLAGS} ctu)
  set_tests_properties(unit_mpi_eb_tree_ctu_subcycle_hydro_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 180)
endif()
if(TARGET test_mpi_eb_tree_subcycle_full_3d)
  foreach(ctu_mode IN ITEMS ctu ctu-common ctu-inflow ctu-inflow-common ctu-periodic ctu-periodic-common ctu-wall ctu-wall-common
      ctu-isothermal ctu-isothermal-common)
    add_test(NAME unit_mpi_eb_tree_${ctu_mode}_full_3d_3
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_subcycle_full_3d> ${MPIEXEC_POSTFLAGS} ${ctu_mode} viscosity)
    set_tests_properties(unit_mpi_eb_tree_${ctu_mode}_full_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 300)
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_eb_tree_held_application
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_held_application")
  set_tests_properties(regression_eb_tree_held_application PROPERTIES PROCESSORS 3 TIMEOUT 600)
  add_test(NAME regression_eb_tree_pressure_application_common
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_pressure_application_common" --paired-pressure)
  set_tests_properties(regression_eb_tree_pressure_application_common
    PROPERTIES PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;restart;application;pressure")
add_test(NAME regression_eb_tree_pressure_application_recursive
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_pressure_application_recursive" --paired-pressure --subcycling)
  set_tests_properties(regression_eb_tree_pressure_application_recursive
    PROPERTIES PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;restart;application;pressure")
  add_test(NAME regression_eb_tree_connected_application_common
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_connected_application_common" --connected-directions)
  set_tests_properties(regression_eb_tree_connected_application_common
    PROPERTIES PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;restart;application;connected")
add_test(NAME regression_eb_tree_connected_application_recursive
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_connected_application_recursive" --connected-directions --subcycling)
  set_tests_properties(regression_eb_tree_connected_application_recursive
    PROPERTIES PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;restart;application;connected")
  foreach(momentum_clock IN ITEMS common recursive)
    set(momentum_clock_args)
    if(momentum_clock STREQUAL "recursive")
      list(APPEND momentum_clock_args --subcycling)
    endif()
    add_test(NAME regression_eb_tree_momentum_application_${momentum_clock}
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
        --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_momentum_application_${momentum_clock}"
        --separate-momentum ${momentum_clock_args})
    set_tests_properties(regression_eb_tree_momentum_application_${momentum_clock}
      PROPERTIES PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;restart;application;momentum")
    add_test(NAME regression_eb_tree_velocity_application_${momentum_clock}
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
        --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_velocity_application_${momentum_clock}"
        --consistent-velocity ${momentum_clock_args})
    set_tests_properties(regression_eb_tree_velocity_application_${momentum_clock}
      PROPERTIES PROCESSORS 3 TIMEOUT 600 LABELS "mpi;eb;restart;application;velocity")
  endforeach()

  add_test(NAME regression_eb_tree_ctu_application
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_eb_tree_ctu_application.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/scheduled.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_ctu_application")
  set_tests_properties(regression_eb_tree_ctu_application PROPERTIES PROCESSORS 3 TIMEOUT 600)
  add_test(NAME regression_eb_tree_isothermal_application
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_eb_tree_ctu_application.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/scheduled.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/eb_tree_isothermal_application" --isothermal)
  set_tests_properties(regression_eb_tree_isothermal_application PROPERTIES PROCESSORS 3 TIMEOUT 600)
endif()

if(TARGET test_mpi_eb_tree_hydro_3d)
  add_test(NAME unit_mpi_eb_tree_ctu_inflow_hydro_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_hydro_3d> ${MPIEXEC_POSTFLAGS} ctu-inflow)
  set_tests_properties(unit_mpi_eb_tree_ctu_inflow_hydro_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 180)
endif()

if(TARGET test_mpi_eb_tree_subcycle_hydro_3d)
  add_test(NAME unit_mpi_eb_tree_ctu_inflow_subcycle_hydro_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_subcycle_hydro_3d> ${MPIEXEC_POSTFLAGS} ctu-inflow)
  set_tests_properties(unit_mpi_eb_tree_ctu_inflow_subcycle_hydro_3d_3 PROPERTIES PROCESSORS 3 TIMEOUT 180)
endif()

if(PELEF_ENABLE_MPI)
  foreach(periodic_step IN ITEMS hydro subcycle_hydro)
    add_test(NAME unit_mpi_eb_tree_periodic_${periodic_step}_3d
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_${periodic_step}_3d> ${MPIEXEC_POSTFLAGS} ctu-periodic)
    set_tests_properties(unit_mpi_eb_tree_periodic_${periodic_step}_3d PROPERTIES
      TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  foreach(wall_step IN ITEMS hydro subcycle_hydro)
    add_test(NAME unit_mpi_eb_tree_wall_${wall_step}_3d
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_${wall_step}_3d> ${MPIEXEC_POSTFLAGS} ctu-wall)
    set_tests_properties(unit_mpi_eb_tree_wall_${wall_step}_3d PROPERTIES
      TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  foreach(thermal_step IN ITEMS transport subcycle_transport)
    add_test(NAME unit_mpi_eb_tree_isothermal_${thermal_step}_3d
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_${thermal_step}_3d> ${MPIEXEC_POSTFLAGS} isothermal)
    set_tests_properties(unit_mpi_eb_tree_isothermal_${thermal_step}_3d PROPERTIES
      TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_sdc_checkpoint_3d unit/test_mpi_eb_tree_sdc_checkpoint_3d.F90)
  target_link_libraries(test_mpi_eb_tree_sdc_checkpoint_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_tree_sdc_checkpoint_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_sdc_checkpoint_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_sdc_checkpoint_3d_2 PROPERTIES TIMEOUT 60 PROCESSORS 2
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_redistribution_3d unit/test_eb_tree_redistribution_3d.F90)
  target_compile_definitions(test_mpi_eb_tree_redistribution_3d PRIVATE PELEF_TEST_MPI_REDIST)
  target_link_libraries(test_mpi_eb_tree_redistribution_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_tree_redistribution_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_redistribution_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_redistribution_3d_2 PROPERTIES
    TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_state_redistribution_3d unit/test_eb_tree_state_redistribution_3d.F90)
  target_compile_definitions(test_mpi_eb_tree_state_redistribution_3d PRIVATE PELEF_TEST_MPI_REDIST)
  target_link_libraries(test_mpi_eb_tree_state_redistribution_3d PRIVATE pelef_mpi_support)
  add_test(NAME unit_mpi_eb_tree_state_redistribution_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_state_redistribution_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_tree_state_redistribution_3d_2 PROPERTIES
    TIMEOUT 120 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  add_test(NAME unit_mpi_eb_tree_sdc_aux_recursive_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} aux_recursive)
  set_tests_properties(unit_mpi_eb_tree_sdc_aux_recursive_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  add_test(NAME unit_mpi_eb_tree_sdc_aux_route_3d_2
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} aux_route)
  set_tests_properties(unit_mpi_eb_tree_sdc_aux_route_3d_2 PROPERTIES
    TIMEOUT 180 PROCESSORS 2 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_tree_transport_3d)
  foreach(transport_halo_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_transport_spatial_halo_3d_${transport_halo_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${transport_halo_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} spatial_halo)
    set_tests_properties(unit_mpi_eb_tree_transport_spatial_halo_3d_${transport_halo_ranks}
      PROPERTIES PROCESSORS ${transport_halo_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_transport_3d)
  foreach(prepared_sdc_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_prepared_mol_sdc_3d_${prepared_sdc_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${prepared_sdc_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} prepared_sdc)
    set_tests_properties(unit_mpi_eb_tree_prepared_mol_sdc_3d_${prepared_sdc_ranks}
      PROPERTIES PROCESSORS ${prepared_sdc_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  foreach(recursive_srd_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_recursive_mol_registers_3d_${recursive_srd_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${recursive_srd_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} recursive_srd)
    set_tests_properties(unit_mpi_eb_tree_recursive_mol_registers_3d_${recursive_srd_ranks}
      PROPERTIES PROCESSORS ${recursive_srd_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_flux_3d)
  foreach(register_close_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_register_close_3d_${register_close_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${register_close_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} register_close)
    set_tests_properties(unit_mpi_eb_tree_register_close_3d_${register_close_ranks}
      PROPERTIES PROCESSORS ${register_close_ranks} TIMEOUT 60)
    add_test(NAME unit_mpi_eb_tree_ancestor_register_close_3d_${register_close_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${register_close_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} gas_ancestor)
    set_tests_properties(unit_mpi_eb_tree_ancestor_register_close_3d_${register_close_ranks}
      PROPERTIES PROCESSORS ${register_close_ranks} TIMEOUT 60 LABELS "mpi;eb;conservation")
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  foreach(nonmol_recursive_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_recursive_nonmol_registers_3d_${nonmol_recursive_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nonmol_recursive_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} nonmol_rec)
    set_tests_properties(unit_mpi_eb_tree_recursive_nonmol_registers_3d_${nonmol_recursive_ranks}
      PROPERTIES PROCESSORS ${nonmol_recursive_ranks} TIMEOUT 300)
  endforeach()
  foreach(nonmol_common_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_common_nonmol_registers_3d_${nonmol_common_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nonmol_common_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} nonmol_srd)
    set_tests_properties(unit_mpi_eb_tree_common_nonmol_registers_3d_${nonmol_common_ranks}
      PROPERTIES PROCESSORS ${nonmol_common_ranks} TIMEOUT 300)
  endforeach()
  foreach(common_srd_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_common_mol_registers_3d_${common_srd_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${common_srd_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} common_srd)
    set_tests_properties(unit_mpi_eb_tree_common_mol_registers_3d_${common_srd_ranks}
      PROPERTIES PROCESSORS ${common_srd_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  foreach(spatial_diffusion_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_spatial_diffusion_callback_3d_${spatial_diffusion_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${spatial_diffusion_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} spatial_diff)
    set_tests_properties(unit_mpi_eb_tree_spatial_diffusion_callback_3d_${spatial_diffusion_ranks}
      PROPERTIES PROCESSORS ${spatial_diffusion_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  foreach(recursive_diffusion_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_recursive_diffusion_callback_3d_${recursive_diffusion_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${recursive_diffusion_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} recursive_diff)
    set_tests_properties(unit_mpi_eb_tree_recursive_diffusion_callback_3d_${recursive_diffusion_ranks}
      PROPERTIES PROCESSORS ${recursive_diffusion_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_transport_3d)
  foreach(multi_sdc_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_multi_mol_registers_3d_${multi_sdc_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${multi_sdc_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} multi_sdc)
    set_tests_properties(unit_mpi_eb_tree_multi_mol_registers_3d_${multi_sdc_ranks}
      PROPERTIES PROCESSORS ${multi_sdc_ranks} TIMEOUT 300)
  endforeach()
  foreach(inert_sdc_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_inert_mol_registers_3d_${inert_sdc_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${inert_sdc_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} inert_iters)
    set_tests_properties(unit_mpi_eb_tree_inert_mol_registers_3d_${inert_sdc_ranks}
      PROPERTIES PROCESSORS ${inert_sdc_ranks} TIMEOUT 300)
  endforeach()
  add_test(NAME unit_mpi_eb_tree_active_ctu_registers_3d_2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} active_ctu)
  set_tests_properties(unit_mpi_eb_tree_active_ctu_registers_3d_2
    PROPERTIES PROCESSORS 2 TIMEOUT 300)
  foreach(godunov_halo_ranks 2 3)
    add_test(NAME unit_mpi_eb_tree_godunov_spatial_halos_3d_${godunov_halo_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${godunov_halo_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} godunov_halo)
    set_tests_properties(unit_mpi_eb_tree_godunov_spatial_halos_3d_${godunov_halo_ranks}
      PROPERTIES PROCESSORS ${godunov_halo_ranks} TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_full_3d)
  foreach(nonmol_callback_mode nonmol_diff nonmol_rd)
    foreach(nonmol_callback_ranks 2 3)
      add_test(NAME unit_mpi_eb_tree_${nonmol_callback_mode}_callback_3d_${nonmol_callback_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nonmol_callback_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} ${nonmol_callback_mode})
      if(nonmol_callback_mode STREQUAL "nonmol_rd")
        set(nonmol_callback_timeout 600)
      else()
        set(nonmol_callback_timeout 300)
      endif()
      set_tests_properties(unit_mpi_eb_tree_${nonmol_callback_mode}_callback_3d_${nonmol_callback_ranks}
        PROPERTIES PROCESSORS ${nonmol_callback_ranks} TIMEOUT ${nonmol_callback_timeout})
    endforeach()
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_transport_3d)
  foreach(physical_srd_ranks 2 3)
    add_test(NAME unit_mpi_state_redist_physical_boundary_3d_${physical_srd_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${physical_srd_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} physical_boundary)
    set_tests_properties(unit_mpi_state_redist_physical_boundary_3d_${physical_srd_ranks}
      PROPERTIES PROCESSORS ${physical_srd_ranks} TIMEOUT 60)
  endforeach()
endif()

if(TARGET test_mpi_eb_tree_transport_3d)
  foreach(constant_eb_mode halo sdc)
    add_test(NAME unit_mpi_eb_tree_constant_${constant_eb_mode}_3d_3
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_transport_3d> ${MPIEXEC_POSTFLAGS} constant_${constant_eb_mode})
    set_tests_properties(unit_mpi_eb_tree_constant_${constant_eb_mode}_3d_3
      PROPERTIES PROCESSORS 3 TIMEOUT 300)
  endforeach()
endif()

if(TARGET test_mpi_amr_patch_tree_transport_level_3d)
  add_test(NAME unit_mpi_amr_patch_tree_constant_level_transport_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_transport_level_3d> ${MPIEXEC_POSTFLAGS} constant)
  set_tests_properties(unit_mpi_amr_patch_tree_constant_level_transport_3d_3 PROPERTIES
    TIMEOUT 120 PROCESSORS 3 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_amr_patch_tree_split_3d)
  add_test(NAME unit_mpi_amr_sdc_auxiliary_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> auxiliary ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_sdc_auxiliary_3d_3 PROPERTIES TIMEOUT 120 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_amr_patch_tree_split_3d)
  add_test(NAME unit_mpi_amr_sdc_recursive_auxiliary_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> recursive_auxiliary ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_amr_sdc_recursive_auxiliary_3d_3 PROPERTIES TIMEOUT 120 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(PELEF_ENABLE_MPI)
  add_test(NAME application_rotating_amr_boundaries
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_rotating_amr_boundaries.py"
      --serial "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_amr_tree_reactive_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/rotating_amr_boundaries")
  set_tests_properties(application_rotating_amr_boundaries PROPERTIES TIMEOUT 600 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  foreach(boundary_ranks IN ITEMS 2 3)
    add_test(NAME unit_mpi_rotating_amr_boundaries_${boundary_ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${boundary_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d> rotating_boundaries ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_rotating_amr_boundaries_${boundary_ranks} PROPERTIES TIMEOUT 300
      PROCESSORS ${boundary_ranks} ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
  add_test(NAME mpi_eb_tree_automatic_timestep_history
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_timestep_history.py"
      --serial "$<TARGET_FILE:pelef_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:pelef_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/automatic_timestep_history")
  set_tests_properties(mpi_eb_tree_automatic_timestep_history PROPERTIES TIMEOUT 300 PROCESSORS 2)
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_source_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> source_replay 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_source_replay_3d_3 PROPERTIES TIMEOUT 300 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_family_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> family_replay 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_family_replay_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_routed_family_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> routed_family_replay 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_routed_family_replay_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_subtree_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> subtree_replay 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_subtree_replay_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_prefix_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> prefix_replay 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_prefix_replay_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_recovery_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> recovery_replay 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_recovery_replay_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_source_stages_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> source_stages 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_source_stages_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_stage_parent_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> stage_parent 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_stage_parent_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_godunov_held_3d)
  add_test(NAME unit_mpi_eb_godunov_source_packets_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_godunov_held_3d> source_packets 1 1 1 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_mpi_eb_godunov_source_packets_3d_3 PROPERTIES TIMEOUT 900 PROCESSORS 3
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET test_mpi_eb_tree_flux_3d)
  add_test(NAME unit_mpi_eb_held_interface_budget_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} held_budget)
  set_tests_properties(unit_mpi_eb_held_interface_budget_3d_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 LABELS "mpi;eb;ueint;conservation")
  add_test(NAME unit_mpi_eb_held_thermal_register_budget_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} held_thermfull)
  set_tests_properties(unit_mpi_eb_held_thermal_register_budget_3d_3 PROPERTIES
    TIMEOUT 60 PROCESSORS 3 LABELS "mpi;eb;ueint;conservation;h2o2_full")
  foreach(gas_thermal_ranks 2 3)
    add_test(NAME unit_mpi_eb_gas_thermal_register_budget_3d_${gas_thermal_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${gas_thermal_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} gas_thermfull)
    set_tests_properties(unit_mpi_eb_gas_thermal_register_budget_3d_${gas_thermal_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${gas_thermal_ranks} LABELS "mpi;eb;conservation;h2o2_full")
  endforeach()
  foreach(gas_fine_ranks 2 3)
    add_test(NAME unit_mpi_eb_gas_fine_eos_budget_3d_${gas_fine_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${gas_fine_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} gas_fine_eos)
    set_tests_properties(unit_mpi_eb_gas_fine_eos_budget_3d_${gas_fine_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${gas_fine_ranks} LABELS "mpi;eb;conservation;h2o2_full")
  endforeach()
  add_test(NAME unit_mpi_eb_held_fine_eos_budget_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} held_fine_eos)
  set_tests_properties(unit_mpi_eb_held_fine_eos_budget_3d_3 PROPERTIES
    TIMEOUT 60 PROCESSORS 3 LABELS "mpi;eb;ueint;conservation;h2o2_full")
  foreach(gas_fine_cov_ranks 2 3)
    add_test(NAME unit_mpi_eb_gas_fine_covered_eos_budget_3d_${gas_fine_cov_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${gas_fine_cov_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} gas_fine_cov)
    set_tests_properties(unit_mpi_eb_gas_fine_covered_eos_budget_3d_${gas_fine_cov_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${gas_fine_cov_ranks} LABELS "mpi;eb;conservation;h2o2_full")
  endforeach()
  add_test(NAME unit_mpi_eb_held_fine_covered_eos_budget_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} held_fine_cov)
  set_tests_properties(unit_mpi_eb_held_fine_covered_eos_budget_3d_3 PROPERTIES
    TIMEOUT 60 PROCESSORS 3 LABELS "mpi;eb;ueint;conservation;h2o2_full")
  foreach(gas_raw_ranks 2 3)
    add_test(NAME unit_mpi_eb_gas_raw_eos_budget_3d_${gas_raw_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${gas_raw_ranks} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} gas_raw_eos)
    set_tests_properties(unit_mpi_eb_gas_raw_eos_budget_3d_${gas_raw_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${gas_raw_ranks} LABELS "mpi;eb;conservation;h2o2_full")
  endforeach()
  add_test(NAME unit_mpi_eb_held_raw_eos_budget_3d_3
    COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_flux_3d> ${MPIEXEC_POSTFLAGS} held_raw_eos)
  set_tests_properties(unit_mpi_eb_held_raw_eos_budget_3d_3 PROPERTIES
    TIMEOUT 60 PROCESSORS 3 LABELS "mpi;eb;ueint;conservation;h2o2_full")
endif()

if(TARGET pelef_spray_tree_runtime AND PELEF_ENABLE_MPI)
  add_executable(test_regular_held_diffusion unit/test_regular_held_diffusion.F90)
  target_link_libraries(test_regular_held_diffusion PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_regular_held_diffusion PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_regular_held_diffusion_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_regular_held_diffusion> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_regular_held_diffusion_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 LABELS "mpi;amr;ueint;mol;transport;les")
endif()

if(TARGET pelef_spray_tree_runtime AND PELEF_ENABLE_MPI)
  add_executable(test_regular_held_boundary_3d unit/test_regular_held_boundary_3d.F90)
  target_link_libraries(test_regular_held_boundary_3d PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_regular_held_boundary_3d PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_regular_held_boundary_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_regular_held_boundary_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_regular_held_boundary_3d_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 LABELS "mpi;amr;ueint;mol;boundary")
endif()

if(TARGET pelef_spray_tree_runtime AND PELEF_ENABLE_MPI)
  add_executable(test_regular_held_godunov_3d unit/test_regular_held_godunov_3d.F90)
  target_link_libraries(test_regular_held_godunov_3d PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_regular_held_godunov_3d PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_regular_held_godunov_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_regular_held_godunov_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_regular_held_godunov_3d_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 LABELS "mpi;amr;ueint;godunov;conservation")
endif()

if(TARGET pelef_spray_tree_runtime AND PELEF_ENABLE_MPI)
  add_executable(test_regular_held_old_reaction_3d unit/test_regular_held_old_reaction_3d.F90)
  target_link_libraries(test_regular_held_old_reaction_3d PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_regular_held_old_reaction_3d PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME unit_regular_held_old_reaction_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_regular_held_old_reaction_3d> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(unit_regular_held_old_reaction_3d_3 PROPERTIES
    TIMEOUT 180 PROCESSORS 3 LABELS "mpi;amr;ueint;mol;chemistry;rollback")
endif()

if(PELEF_ENABLE_MPI AND TARGET test_mpi_eb_native_mol_replay_3d)
  add_test(NAME unit_mpi_eb_native_pressure_replay_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_replay_3d> ${MPIEXEC_POSTFLAGS} pressure)
  add_test(NAME unit_mpi_eb_native_pressure_retry_3d_3
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_native_mol_future_source_3d> ${MPIEXEC_POSTFLAGS} retry pressure)
  set_tests_properties(unit_mpi_eb_native_pressure_replay_3d_3 unit_mpi_eb_native_pressure_retry_3d_3 PROPERTIES
    LABELS "mpi;eb;mol;pressure;rollback" PROCESSORS 3 TIMEOUT 600
    ENVIRONMENT "HWLOC_COMPONENTS=-gl,-x86;OMPI_MCA_hwloc_base_binding_policy=none")
endif()

# Worker 09 proposal: publication failure must agree across all owners.
if(TARGET test_mpi_amr_patch_tree_checkpoint_3d)
  foreach(publication_ranks IN ITEMS 2 3)
    set(publication_work "${CMAKE_CURRENT_BINARY_DIR}/checkpoint-publication-mpi-${publication_ranks}")
    file(MAKE_DIRECTORY "${publication_work}")
    add_test(NAME unit_mpi_amr_tree_checkpoint_publication_${publication_ranks}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" ${publication_ranks}
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_amr_patch_tree_checkpoint_3d>
        publication ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_mpi_amr_tree_checkpoint_publication_${publication_ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${publication_ranks} WORKING_DIRECTORY "${publication_work}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

# Append to tests/cmake/MPIUnitTests.cmake; kept separate from the fixture patch.
if(TARGET test_mpi_amr_patch_tree_split_3d AND Python3_Interpreter_FOUND)
  set(repair_restart_mpi_flags)
  foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND repair_restart_mpi_flags "--preflag=${flag}")
  endforeach()
  foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND repair_restart_mpi_flags "--postflag=${flag}")
  endforeach()
  add_test(NAME unit_mpi_amr_active_repair_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_active_repair_restart.py"
      --exe $<TARGET_FILE:test_mpi_amr_patch_tree_split_3d>
      --work "${CMAKE_CURRENT_BINARY_DIR}/active_repair_restart"
      --mpiexec "${MPIEXEC_EXECUTABLE}" "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      ${repair_restart_mpi_flags})
  set_tests_properties(unit_mpi_amr_active_repair_restart_3d PROPERTIES
    TIMEOUT 180 PROCESSORS 3 LABELS "amr;reflux;restart;conservation")
endif()
