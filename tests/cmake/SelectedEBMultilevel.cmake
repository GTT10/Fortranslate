add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_three_level_restart_reference
  regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_three_level_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_three_level_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_three_level_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Completed coarse steps:"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_mismatch_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_mismatch_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_mismatch_finest.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_three_level_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_three_level_restart_2d.py"
    --fixed-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_three_level_reference_root.csv"
    --fixed-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_three_level_reference_middle.csv"
    --fixed-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_three_level_reference_finest.csv"
    --reference-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_reference_root.csv"
    --reference-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_reference_middle.csv"
    --reference-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_reference_finest.csv"
    --stopped-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_stopped_root.csv"
    --stopped-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_stopped_middle.csv"
    --stopped-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_stopped_finest.csv"
    --restarted-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_restarted_root.csv"
    --restarted-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_restarted_middle.csv"
    --restarted-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_restarted_finest.csv"
    --expected-reference-root-sha256 92f9580d35032710b09653e6cfd80d97c133c6df8fd70af0db4f5011714ec980
    --expected-reference-middle-sha256 59e66071d6aaf3fe1fd62de9db72d28163acbb7820c95e9a93f8618ab013c1b9
    --expected-reference-finest-sha256 70062b794483b16de5e7de66a3aeaebfbd914dd9ac1adb7bf3515a9b7a9f4250
    --expected-stopped-root-sha256 53ca2ae393f4586f8632ad5d38fc7d902fc85f6b57c31f8dbbebd778065c2fc4
    --expected-stopped-middle-sha256 38ff58c0bd7c83a87c2cf7752f770c5864d04cb3ed6234ba6fe815424372af11
    --expected-stopped-finest-sha256 46abc4426db3ba608d3bf65f3a5937ef4427f1bc588cfe83f97d117e7a57e3d5
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level.chk"
    --expected-checkpoint-sha256 e8655df9aeb62ece3bcf3e23ee5ed1bcda17a2d834cf214e3ecaaf71f3e952b6
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_three_level_restart_reference;regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_three_level_restart_restart;regression_selected_reactive_eb_amr_2d_three_level_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_three_level_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_three_level_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference
  PROPERTIES
    PASS_REGULAR_EXPRESSION
      "Completed coarse steps:;Completed regrids: 1"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_mismatch_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_mismatch_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_mismatch_finest.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_three_level_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_dynamic_three_level_restart_2d.py"
    --fixed-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_three_level_reference_root.csv"
    --fixed-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_three_level_reference_middle.csv"
    --fixed-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_three_level_reference_finest.csv"
    --reference-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_reference_root.csv"
    --reference-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_reference_middle.csv"
    --reference-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_reference_finest.csv"
    --stopped-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_stopped_root.csv"
    --stopped-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_stopped_middle.csv"
    --stopped-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_stopped_finest.csv"
    --restarted-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_restarted_root.csv"
    --restarted-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_restarted_middle.csv"
    --restarted-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_restarted_finest.csv"
    --expected-reference-root-sha256 91e8ecebce463172ed8c04a1a65589845c9e49c9ae4c1c67e7ae2627d1252ebf
    --expected-reference-middle-sha256 c752ebcd7308a888395b9206bfc944585d88c011bc86d10f7450a8286f5a237e
    --expected-reference-finest-sha256 ceb72783d7a0f531667227862eba3d877af464d46585e86e75cf377c42bef5dd
    --expected-stopped-root-sha256 84dc02136afcef4cb4227d86fb36e8ec7902fee61fd09fbe7a018a55e68bc615
    --expected-stopped-middle-sha256 81b59680ea1347ccded15dc588f763420370c48557f3f83ffa73c01197a7f5b7
    --expected-stopped-finest-sha256 f9d08b869ff4e5cf4b23c0f3193b3a0f2420973e80d0c9ec0b5cb9f56fd78645
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level.chk"
    --expected-checkpoint-sha256 1e8e24e2ead15ac71199cb73146656ffb63be089bafc57eced95fbf419ec1502
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_three_level_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_multipatch_restart_reference
  regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_multipatch_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_multipatch_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_multipatch_restart_reference
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Fine patches: 2"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_mismatch_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_mismatch_fine_patch0001.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_mismatch_fine_patch0002.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_multipatch_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_multipatch_restart_2d.py"
    --fixed-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_multipatch_reference_root.csv"
    --fixed-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_multipatch_reference_fine_patch0001.csv"
    --fixed-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_multipatch_reference_fine_patch0002.csv"
    --reference-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_reference_root.csv"
    --reference-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_reference_fine_patch0001.csv"
    --reference-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_reference_fine_patch0002.csv"
    --stopped-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_stopped_root.csv"
    --stopped-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_stopped_fine_patch0001.csv"
    --stopped-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_stopped_fine_patch0002.csv"
    --restarted-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restarted_root.csv"
    --restarted-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restarted_fine_patch0001.csv"
    --restarted-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restarted_fine_patch0002.csv"
    --expected-reference-root-sha256
      1a2c70bf88bb36466937e0839294ca6363dad34752742e2322a9e49a90507f80
    --expected-reference-patch1-sha256
      e6a2d52e082a9b17ae00325cb5a65c2ee0ee8425646991b0b77e5b17cef42d88
    --expected-reference-patch2-sha256
      ceb4492ca19efbee67ed6e433a787eb48c887193544e1f88acb364d7b7f41b44
    --expected-stopped-root-sha256
      81d4a5ed93947b8790fa19c1ed6df716f8a5405bc84d0da974028d1acb45fdf1
    --expected-stopped-patch1-sha256
      9190c1137c5e41499d421f5fcd97b95d32466eb54a4bbd567a36e21efbe26ad5
    --expected-stopped-patch2-sha256
      2c1833dbd96d1aa3e202ebe87c62085e9d8a8939f77c6e2bc91955b872cda644
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restart.chk"
    --expected-checkpoint-sha256
      72fdd25e457e27ece30739ca52d59536cb3ed247439c14bbd0777e534ad2b368
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_multipatch_restart_reference;regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_multipatch_restart_restart;regression_selected_reactive_eb_amr_2d_multipatch_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_multipatch_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_2d.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_2d_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_embedded_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_temperature.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_temperature_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_path_alias.nml
    --expected "paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_alias_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_checkpoint_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_checkpoint_input_alias.nml
    --expected "checkpoint path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_alias.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_alias_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_output_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_restart_output_alias.nml"
    --expected "restart path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_alias.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_alias_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_checkpoint_restart_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_checkpoint_restart_alias.nml"
    --expected "checkpoint and restart paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_restart.chk"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_restart_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_restart_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_removal_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_three_level_mode.nml"
    --expected "Dynamic three-level EB AMR keeps the finest patch active"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_three_level.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_three_level_fine.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_three_level_finest.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_derived_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_multipatch_derived_alias.nml"
    --expected "multipatch derived paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_multipatch.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_multipatch_fine_patch0001.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_parent_without_dynamic_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_dynamic_parent_mode.nml"
    --expected "Dynamic parent regridding requires dynamic three-level EB AMR"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_dynamic_parent.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_dynamic_parent_fine.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_dynamic_parent_finest.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fine_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_fine_input_alias.nml
    --expected "paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_fine_alias_coarse.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_coarse_fine_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_coarse_fine_alias.nml"
    --expected "fine output must be nonempty and distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_level_alias.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_finest_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_finest_input_alias.nml
    --expected "three-level paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_finest_input_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_finest_input_middle.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_checkpoint_finest_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_checkpoint_finest_alias.nml"
    --expected "checkpoint path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_finest_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_finest_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_finest.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_finest_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_restart_finest_alias.nml"
    --expected "restart path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_finest_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_finest_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_finest.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_unknown_rejected
  regression_selected_reactive_eb_amr_2d_temperature_rejected
  regression_selected_reactive_eb_amr_2d_path_alias_rejected
  regression_selected_reactive_eb_amr_2d_checkpoint_input_alias_rejected
  regression_selected_reactive_eb_amr_2d_restart_output_alias_rejected
  regression_selected_reactive_eb_amr_2d_checkpoint_restart_alias_rejected
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_removal_rejected
  regression_selected_reactive_eb_amr_2d_multipatch_derived_alias_rejected
  regression_selected_reactive_eb_amr_2d_dynamic_parent_without_dynamic_rejected
  regression_selected_reactive_eb_amr_2d_fine_input_alias_rejected
  regression_selected_reactive_eb_amr_2d_coarse_fine_alias_rejected
  regression_selected_reactive_eb_amr_2d_finest_input_alias_rejected
  regression_selected_reactive_eb_amr_2d_checkpoint_finest_alias_rejected
  regression_selected_reactive_eb_amr_2d_restart_finest_alias_rejected
  PROPERTIES
    RESOURCE_LOCK selected_reactive_eb_amr_2d_invalid
    WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_eb_amr_2d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture.nml"
    --expected "Unknown reactive 2D chemistry model"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
)
set_tests_properties(
  regression_fixed_reactive_eb_amr_2d_selected_model_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)

