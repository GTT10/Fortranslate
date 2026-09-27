if(PELEF_ENABLE_MPI)
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
endif()

if(PELEF_ENABLE_MPI)
  add_executable(test_mpi_eb_tree_full_3d unit/test_mpi_eb_tree_full_3d.F90)
  target_link_libraries(test_mpi_eb_tree_full_3d PRIVATE pelef_mpi_support)
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
  add_test(NAME unit_mpi_eb_tree_sdc_diffusion_3d_2
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_eb_tree_full_3d> ${MPIEXEC_POSTFLAGS} diffusion)
  set_tests_properties(unit_mpi_eb_tree_sdc_diffusion_3d_2 PROPERTIES TIMEOUT 300 PROCESSORS 2)
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
