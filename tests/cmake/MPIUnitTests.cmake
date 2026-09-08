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

