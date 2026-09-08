add_executable(
  test_selected_reactive_1d_config
  unit/test_selected_reactive_1d_config.F90
)
target_link_libraries(
  test_selected_reactive_1d_config
  PRIVATE pelef_selected_reactive_1d_runtime
)
add_test(
  NAME unit_selected_reactive_1d_config
  COMMAND test_selected_reactive_1d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/fixture.nml"
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture.nml"
)

set(
  SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_fixture"
)
set(
  SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_full"
)
set(
  SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_invalid"
)
set(
  SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)

add_test(
  NAME regression_selected_reactive_1d_fixture_run
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/fixture.nml"
)
set_tests_properties(
  regression_selected_reactive_1d_fixture_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_1d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_1d.py"
    --input
      "${SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_1d.csv"
    --species H2 H
    --nx 8
    --final-time 1.0e-7
    --require-uniform
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_reactive_1d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_1d_fixture_run
)

add_test(
  NAME regression_selected_reactive_1d_full_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_full_h2o2_1d/uniform.nml"
)
add_test(
  NAME regression_selected_reactive_1d_full_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/h2o2_full.nml"
)
set_tests_properties(
  regression_selected_reactive_1d_full_fixed_run
  regression_selected_reactive_1d_full_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_1d_full_selected_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Species: 10"
)
add_test(
  NAME regression_selected_reactive_1d_full_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_1d.py"
    --input
      "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_1d.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8
    --final-time 2.0e-6
    --require-uniform
)
set_tests_properties(
  regression_selected_reactive_1d_full_structure
  PROPERTIES DEPENDS regression_selected_reactive_1d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_1d_full_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}/reactive_full_h2o2_1d.csv"
    "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_1d.csv"
)
set_tests_properties(
  regression_selected_reactive_1d_full_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_1d_full_fixed_run;regression_selected_reactive_1d_full_selected_run"
)

add_test(
  NAME regression_selected_reactive_1d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_1d.csv"
)
add_test(
  NAME regression_selected_reactive_1d_entropy_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/invalid_entropy_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_entropy_temperature.csv"
)
set_tests_properties(
  regression_selected_reactive_1d_unknown_rejected
  regression_selected_reactive_1d_entropy_temperature_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_1d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_1d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/fixture.nml"
    --expected "Unknown reactive 1D chemistry model"
    --forbidden-output
      "${SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_1d.csv"
)
set_tests_properties(
  regression_fixed_reactive_1d_selected_model_rejected
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)

set(
  SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_fixture"
)
set(
  SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_two_level"
)
set(
  SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_multilevel"
)
set(
  SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_multipatch"
)
set(
  SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_boundary"
)
set(
  SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_invalid"
)
set(
  SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_path_alias.nml"
  "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_amr_reactive_1d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_fixture_run_a
  regression_selected_amr_reactive_1d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Fine level active:  T"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_amr_reactive_1d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d.csv"
    --species H2 H
    --expected-levels 0 1
    --x-upper 0.008
    --final-time 1.0e-7
    --refinement-ratio 2
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-8
)
set_tests_properties(
  regression_selected_amr_reactive_1d_fixture_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_fixture_run_a
)
add_test(
  NAME regression_selected_amr_reactive_1d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d.csv"
    "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d_repeat.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_fixture_run_a;regression_selected_amr_reactive_1d_fixture_run_b"
)

add_test(
  NAME regression_selected_amr_reactive_1d_two_level_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_two_level_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_two_level_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_two_level_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_two_level_fixed_run
  regression_selected_amr_reactive_1d_two_level_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Fine level active:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_two_level_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_two_level.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1
    --x-upper 0.008
    --final-time 1.0e-9
    --refinement-ratio 2
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-12
)
set_tests_properties(
  regression_selected_amr_reactive_1d_two_level_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_two_level_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_two_level_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}/full_h2o2_amr_two_level.csv"
    "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_two_level.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_two_level_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_two_level_fixed_run;regression_selected_amr_reactive_1d_two_level_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multilevel_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multilevel_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multilevel_fixed_run
  regression_selected_amr_reactive_1d_multilevel_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Active AMR levels: 3"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_multilevel.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1 2
    --x-upper 0.012
    --final-time 2.0e-7
    --refinement-ratio 2
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multilevel_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_multilevel_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}/full_h2o2_amr_multilevel.csv"
    "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_multilevel.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multilevel_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_multilevel_fixed_run;regression_selected_amr_reactive_1d_multilevel_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multipatch_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multipatch_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multipatch_fixed_run
  regression_selected_amr_reactive_1d_multipatch_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Active fine patches: 3"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}/selected_full_h2o2_amr_multipatch.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1
    --x-upper 0.012
    --final-time 1.0e-10
    --refinement-ratio 2
    --minimum-finest-regions 3
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multipatch_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_multipatch_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}/full_h2o2_amr_multipatch.csv"
    "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}/selected_full_h2o2_amr_multipatch.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multipatch_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_multipatch_fixed_run;regression_selected_amr_reactive_1d_multipatch_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_boundary_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_boundary_weno_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_boundary_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_boundary_weno_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_boundary_fixed_run
  regression_selected_amr_reactive_1d_boundary_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "AMR hybrid WENO:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_boundary_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}/selected_full_h2o2_amr_boundary_weno.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1 2
    --x-upper 0.012
    --final-time 1.0e-7
    --refinement-ratio 2
    --require-finest-touches-lower
)
set_tests_properties(
  regression_selected_amr_reactive_1d_boundary_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_boundary_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_boundary_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}/full_h2o2_amr_boundary_weno.csv"
    "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}/selected_full_h2o2_amr_boundary_weno.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_boundary_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_boundary_fixed_run;regression_selected_amr_reactive_1d_boundary_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_reactive_1d.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_non_amr_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_non_amr.nml"
    --expected "requires amr_enabled"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_non_amr.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_entropy_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_entropy_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_entropy_temperature.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_checkpoint_interval_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_checkpoint_interval.nml"
    --expected "does not implement checkpoint or restart"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_interval.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_checkpoint_file_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_checkpoint_file.nml"
    --expected "does not implement checkpoint or restart"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_file.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_restart_file_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_restart_file.nml"
    --expected "does not implement checkpoint or restart"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input invalid_path_alias.nml
    --expected "input and output paths must differ"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_alias.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_unknown_rejected
  regression_selected_amr_reactive_1d_non_amr_rejected
  regression_selected_amr_reactive_1d_entropy_temperature_rejected
  regression_selected_amr_reactive_1d_checkpoint_interval_rejected
  regression_selected_amr_reactive_1d_checkpoint_file_rejected
  regression_selected_amr_reactive_1d_restart_file_rejected
  regression_selected_amr_reactive_1d_path_alias_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_amr_reactive_1d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_amr_reactive_1d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture.nml"
    --expected "Unknown reactive 1D chemistry model"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d.csv"
)
set_tests_properties(
  regression_fixed_amr_reactive_1d_selected_model_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)

