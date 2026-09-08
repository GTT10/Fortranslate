set(REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_dynamic_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_2d/hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_2d_run
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Completed regrids: [1-9][0-9]*"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_dynamic_2d.py"
    --root
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}/three_level_dynamic_root.csv"
    --middle
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}/three_level_dynamic_middle.csv"
    --finest
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}/three_level_dynamic_finest.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_2d_check
  PROPERTIES
    DEPENDS regression_reactive_eb_amr_three_level_dynamic_2d_run
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}"
)

set(REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_dynamic_restart_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_restart_2d/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_reference
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Completed regrids: [1-9][0-9]*"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart
  PROPERTIES DEPENDS
    regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_dynamic_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_three_level_restart.chk"
    --expected-checkpoint-sha256
      ced92fed867264fd0863b5c85a90702e29647e2494a5d6d69b87a51f194702c5
    --reference
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_reference_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_reference_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_reference_finest.csv"
    --stopped
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_stopped_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_stopped_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_stopped_finest.csv"
    --restarted
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_restarted_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_restarted_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_restarted_finest.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_three_level_dynamic_restart_2d_reference;regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop;regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart"
)

set(REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_restart_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_restart_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_restart_2d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_amr_restart_2d_restart
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_restart_2d/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_restart_2d_reference
  regression_reactive_eb_amr_restart_2d_checkpoint_stop
  regression_reactive_eb_amr_restart_2d_restart
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_amr_restart_2d_restart
  PROPERTIES DEPENDS regression_reactive_eb_amr_restart_2d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_amr_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart.chk"
    --expected-checkpoint-sha256
      40a0781f99a231dadb663365342b92a3e035e6dd195a517a6e45010af87eee17
    --reference
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_reference_coarse_2d.csv"
    --stopped
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_stopped_coarse_2d.csv"
    --expected-stopped-sha256
      2465d2facdd7662091d8e6c1346e8f1a0205ba1b758398309fecb7fe769aa14d
    --restarted
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restarted_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_reference_fine_2d.csv"
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_stopped_fine_2d.csv"
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restarted_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_restart_2d_reference;regression_reactive_eb_amr_restart_2d_checkpoint_stop;regression_reactive_eb_amr_restart_2d_restart"
)

set(REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_chemistry_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_chemistry_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_chemistry_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_chemistry_2d_inert
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_chemistry_2d/inert.nml"
)
add_test(
  NAME regression_reactive_eb_chemistry_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_chemistry_2d/reactive.nml"
)
set_tests_properties(
  regression_reactive_eb_chemistry_2d_reference
  regression_reactive_eb_chemistry_2d_inert
  regression_reactive_eb_chemistry_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_chemistry_2d_parity
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_chemistry_2d.py"
    --reference
      "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_chemistry_reference.csv"
    --reactive
      "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_chemistry.csv"
    --inert
      "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_chemistry_inert.csv"
)
set_tests_properties(
  regression_reactive_eb_chemistry_2d_parity
  PROPERTIES DEPENDS
    "regression_reactive_eb_chemistry_2d_reference;regression_reactive_eb_chemistry_2d_inert;regression_reactive_eb_chemistry_2d_run"
)

if(PELEF_ENABLE_MPI)
  add_executable(
    test_mpi_amr_transport_context_3d
    mpi/test_mpi_amr_transport_context_3d.F90
  )
  target_link_libraries(
    test_mpi_amr_transport_context_3d PRIVATE
    "$<TARGET_NAME_IF_EXISTS:pelef_mpi_support>"
  )
  add_test(
    NAME unit_mpi_amr_transport_context_3d_np2
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_transport_context_3d>
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    unit_mpi_amr_transport_context_3d_np2
    PROPERTIES
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION
        "test_mpi_amr_transport_context_3d: PASS"
      TIMEOUT 60
  )
endif()
