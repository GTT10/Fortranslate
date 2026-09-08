set(
  SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_fixture"
)
set(
  SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_full"
)
set(
  SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_elementary"
)
set(
  SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_invalid"
)
set(
  SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_fixed_rejection"
)
set(
  SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_restart"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_path_alias.nml"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_checkpoint_interval.nml"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_interval.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_restart.nml"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_restart.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_reactive_eb_3d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_eb_3d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_fixture_run_a
  regression_selected_reactive_eb_3d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_3d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_3d.py"
    --input
      "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d.csv"
    --species H2 H
    --nx 4 --ny 4 --nz 4
    --final-time 1.0e-7
    --x-upper 0.004 --y-upper 0.004 --z-upper 0.004
    --plane-axis x --plane-position 0.00163
    --expected-regular 32 --expected-cut 16 --expected-covered 16
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_reactive_eb_3d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_3d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_3d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d.csv"
    "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_3d_fixture_run_a;regression_selected_reactive_eb_3d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_eb_3d_full_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/h2o2_full_transport_chemistry.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_full_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Operator sequence: R-T-H-T-R"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_3d_full_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_3d.py"
    --input
      "${SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY}/selected_full_reactive_eb_3d_transport_chemistry.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 10 --ny 8 --nz 6
    --final-time 2.0e-9
    --x-upper 1.0e-4 --y-upper 1.0e-4 --z-upper 1.0e-4
    --plane-axis x --plane-position 3.95e-5
    --expected-regular 288 --expected-cut 48 --expected-covered 144
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-9
    --minimum-active-density-span 2.0e-2
)
set_tests_properties(
  regression_selected_reactive_eb_3d_full_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_3d_full_run
)

add_test(
  NAME regression_selected_reactive_eb_3d_elementary_run
  COMMAND $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/h2o2_elementary_transport_chemistry.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_elementary_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_3d_elementary_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/reactive_eb_3d_transport_chemistry.csv"
    "${SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY}/selected_elementary_reactive_eb_3d.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_elementary_exact
  PROPERTIES DEPENDS
    "regression_reactive_eb_transport_3d_coupled_run;regression_selected_reactive_eb_3d_elementary_run"
)

add_test(
  NAME regression_selected_reactive_eb_3d_restart_reference
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_checkpoint_stop
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_restart
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_reference
  regression_selected_reactive_eb_3d_restart_checkpoint_stop
  regression_selected_reactive_eb_3d_restart_restart
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_3d_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart continuation: complete"
    DEPENDS regression_selected_reactive_eb_3d_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_mismatch.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}"
    DEPENDS regression_selected_reactive_eb_3d_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_3d_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_3d_restart.py"
    --fixed-reference
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart_reference.csv"
    --reference
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_reference.csv"
    --stopped
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_stopped.csv"
    --restarted
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_restarted.csv"
    --checkpoint
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart.chk"
    --reference-log
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_check
  PROPERTIES
    DEPENDS
      "regression_reactive_eb_restart_3d_reference;regression_selected_reactive_eb_3d_restart_reference;regression_selected_reactive_eb_3d_restart_checkpoint_stop;regression_selected_reactive_eb_3d_restart_restart"
    RESOURCE_LOCK selected_reactive_eb_3d_restart
    TIMEOUT 60
)

foreach(
    selected_reactive_eb_3d_invalid_case IN ITEMS
    unknown_species
    duplicate_species
    nonfinite_composition
    zero_composition
    negative_composition
    temperature
  )
  if(selected_reactive_eb_3d_invalid_case STREQUAL "unknown_species")
    set(selected_reactive_eb_3d_expected "not found exactly once")
    set(selected_reactive_eb_3d_output "unknown")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "duplicate_species")
    set(selected_reactive_eb_3d_expected "duplicate species")
    set(selected_reactive_eb_3d_output "duplicate")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "nonfinite_composition")
    set(selected_reactive_eb_3d_expected "nonfinite mole fraction")
    set(selected_reactive_eb_3d_output "nonfinite")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "zero_composition")
    set(selected_reactive_eb_3d_expected "nonnegative and nonempty")
    set(selected_reactive_eb_3d_output "zero")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "negative_composition")
    set(selected_reactive_eb_3d_expected "nonnegative and nonempty")
    set(selected_reactive_eb_3d_output "negative")
  else()
    set(selected_reactive_eb_3d_expected "outside the mechanism range")
    set(selected_reactive_eb_3d_output "temperature")
  endif()
  add_test(
    NAME
      regression_selected_reactive_eb_3d_${selected_reactive_eb_3d_invalid_case}_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
      --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
      --input
        "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_${selected_reactive_eb_3d_invalid_case}.nml"
      --expected "${selected_reactive_eb_3d_expected}"
      --forbidden-output
        "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_${selected_reactive_eb_3d_output}.csv"
  )
endforeach()

add_test(
  NAME regression_selected_reactive_eb_3d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    --input invalid_path_alias.nml
    --expected "input and output paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_alias.csv"
)
add_test(
  NAME regression_selected_reactive_eb_3d_checkpoint_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    --input invalid_checkpoint_interval.nml
    --expected "input and checkpoint paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_checkpoint_interval.csv"
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    --input invalid_restart.nml
    --expected "input and restart paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_restart.csv"
)
add_test(
  NAME regression_selected_reactive_eb_3d_unsupported_element_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable
      $<TARGET_FILE:test_selected_unsupported_element_reactive_eb_3d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_unsupported_element.nml"
    --expected "chemistry requires supported H/O/N species"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_element.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_unknown_species_rejected
  regression_selected_reactive_eb_3d_duplicate_species_rejected
  regression_selected_reactive_eb_3d_nonfinite_composition_rejected
  regression_selected_reactive_eb_3d_zero_composition_rejected
  regression_selected_reactive_eb_3d_negative_composition_rejected
  regression_selected_reactive_eb_3d_temperature_rejected
  regression_selected_reactive_eb_3d_path_alias_rejected
  regression_selected_reactive_eb_3d_checkpoint_input_alias_rejected
  regression_selected_reactive_eb_3d_restart_input_alias_rejected
  regression_selected_reactive_eb_3d_unsupported_element_rejected
  PROPERTIES
    RESOURCE_LOCK selected_reactive_eb_3d_invalid
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_eb_3d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_eb_3d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture.nml"
    --expected "selected thermo requires explicit opt-in"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d.csv"
)
set_tests_properties(
  regression_fixed_reactive_eb_3d_selected_model_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY}"
)

