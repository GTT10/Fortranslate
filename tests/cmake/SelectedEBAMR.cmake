add_executable(
  test_selected_reactive_eb_amr_2d_config
  unit/test_selected_reactive_eb_amr_2d_config.F90
)
target_link_libraries(
  test_selected_reactive_eb_amr_2d_config
  PRIVATE pelef_selected_reactive_eb_amr_2d_runtime
)
add_test(
  NAME unit_selected_reactive_eb_amr_2d_config
  COMMAND test_selected_reactive_eb_amr_2d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture.nml"
)

set(
  SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_fixture"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_full"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_invalid"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_dynamic_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_three_level_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_dynamic_three_level_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_multipatch_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_path_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_fine_input_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_fine_input_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_checkpoint_input_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_input_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_finest_input_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_finest_input_alias.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_run_a
  regression_selected_reactive_eb_amr_2d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Coarse cut cells: 8"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_coarse_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    --species H2 H
    --nx 8
    --ny 8
    --final-time 1.0e-7
    --expected-regular 32
    --expected-cut 8
    --expected-covered 24
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_fine_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
    --species H2 H
    --nx 12
    --ny 12
    --final-time 1.0e-7
    --expected-regular 84
    --expected-cut 12
    --expected-covered 48
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_coarse_check
  regression_selected_reactive_eb_amr_2d_fixture_fine_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_hierarchy_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_2d.py"
    --coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    --fine
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
    --species H2 H
    --coarse-nx 8
    --coarse-ny 8
    --x-lower 0.0
    --x-upper 0.008
    --y-lower 0.0
    --y-upper 0.008
    --coarse-i-lower 2
    --coarse-i-upper 7
    --coarse-j-lower 2
    --coarse-j-upper 7
    --refinement-ratio 2
    --final-time 1.0e-7
    --coarse-counts 32 8 24
    --fine-counts 84 12 48
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
    --expected-coarse-sha256
      40758bed925456948aa694ee130e766b8b942a2bdfef104c268323a981119c0e
    --expected-fine-sha256
      30728d746ccdbbec2b238e4c738f6e0e54ba9550ec5d9ea609c5d7972551fbfe
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_hierarchy_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_coarse_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse_repeat.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_fine_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_coarse_exact
  regression_selected_reactive_eb_amr_2d_fixture_fine_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_amr_2d_fixture_run_a;regression_selected_reactive_eb_amr_2d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/h2o2_full_fixed.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/h2o2_full_selected.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_fixed_run
  regression_selected_reactive_eb_amr_2d_full_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_coarse_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_coarse.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8
    --ny 8
    --final-time 2.0e-8
    --expected-regular 32
    --expected-cut 8
    --expected-covered 24
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_fine_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_fine.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 12
    --ny 12
    --final-time 2.0e-8
    --expected-regular 84
    --expected-cut 12
    --expected-covered 48
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_coarse_check
  regression_selected_reactive_eb_amr_2d_full_fine_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_hierarchy_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_2d.py"
    --coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_coarse.csv"
    --fine
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_fine.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --coarse-nx 8
    --coarse-ny 8
    --x-lower 0.0
    --x-upper 0.008
    --y-lower 0.0
    --y-upper 0.008
    --coarse-i-lower 2
    --coarse-i-upper 7
    --coarse-j-lower 2
    --coarse-j-upper 7
    --refinement-ratio 2
    --final-time 2.0e-8
    --coarse-counts 32 8 24
    --fine-counts 84 12 48
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
    --expected-coarse-sha256
      $<IF:$<CONFIG:Debug>,bee1bd8baae3606a9ac3483822140c2d0d29324683c9b9deab62c84fe4819d7b,f6e6fe44a629a92c3027f52378a0f892ab8e790a88bdd99dc8fb178aa18605a7>
    --expected-fine-sha256
      $<IF:$<CONFIG:Debug>,be55e526b0cf07fbfff63cd6f874e104599f0d0a019c02ca1145f614acc4079e,a61021c2a77cfd933e4d1af1027ee08b869a4bed4791534080d9aa1cfc069423>
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_hierarchy_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_coarse_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/full_h2o2_reactive_eb_amr_2d_coarse.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_coarse.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_fine_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/full_h2o2_reactive_eb_amr_2d_fine.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_fine.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_coarse_exact
  regression_selected_reactive_eb_amr_2d_full_fine_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_amr_2d_full_fixed_run;regression_selected_reactive_eb_amr_2d_full_selected_run"
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_restart_reference
  regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Completed coarse steps: 2"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_mismatch_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_mismatch_fine.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}"
    DEPENDS regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_2d_restart.py"
    --fixed-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_restart_reference_coarse.csv"
    --fixed-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_restart_reference_fine.csv"
    --reference-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_reference_coarse.csv"
    --reference-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_reference_fine.csv"
    --stopped-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_stopped_coarse.csv"
    --stopped-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_stopped_fine.csv"
    --restarted-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restarted_coarse.csv"
    --restarted-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restarted_fine.csv"
    --expected-reference-coarse-sha256
      5d8e24295aa1331dfdd7a6224fd0809a5df712e2db22c733ac44bed5ffb3f7ac
    --expected-reference-fine-sha256
      cbfef95067e8674239122633711c8e35371ea97e3709f80beeadccbabbc64af2
    --expected-stopped-coarse-sha256
      dda501bd9c27cb090f8bc8f7d06ec55f27e02fdd59e667cff8f89b2edb52c736
    --expected-stopped-fine-sha256
      6f30dd051f9fa46180a539ca62b7dc0f7c237850b9fe1ce0bdab557258bec817
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart.chk"
    --expected-checkpoint-sha256
      6d1a5e35cc001f3dfb7049cba243482bf902acde8da3dd487d4ba93598d784ba
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_restart_reference;regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_restart_restart;regression_selected_reactive_eb_amr_2d_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_restart_reference
  regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_dynamic_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Completed coarse steps: 3"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_mismatch_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_mismatch_fine.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_dynamic_restart_2d.py"
    --fixed-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_restart_reference_coarse.csv"
    --fixed-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_restart_reference_fine.csv"
    --reference-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_reference_coarse.csv"
    --reference-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_reference_fine.csv"
    --stopped-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_stopped_coarse.csv"
    --stopped-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_stopped_fine.csv"
    --restarted-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restarted_coarse.csv"
    --restarted-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restarted_fine.csv"
    --expected-reference-coarse-sha256
      2a67c9b9a618fc36412c0156adec150e3271a5b1ef9ab5ec84e158a8d07c1a87
    --expected-reference-fine-sha256
      66f479a71de2474ae2cad3f126b847238311d5e9a1ae40d0b60d25a1a9b776e9
    --expected-stopped-coarse-sha256
      8c0b0c74beb20ad5915e03578fa16c336b235b1d10d2e3e96b183d261052ba72
    --expected-stopped-fine-sha256
      a0e71a149a0a0991732024bef41f46cdd0b3318946390141ae0a2071f5c5acb2
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart.chk"
    --expected-checkpoint-sha256
      a91050d485886e3109337331f0af2ce0df8f8b8aa99227a2aaa0c69feeb2500d
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_dynamic_restart_reference;regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_dynamic_restart_restart;regression_selected_reactive_eb_amr_2d_dynamic_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_restart
    TIMEOUT 60
)

