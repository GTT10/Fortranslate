set(
  SELECTED_CVODE_FIXTURE_A_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_fixture_a"
)
set(
  SELECTED_CVODE_FIXTURE_B_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_fixture_b"
)
set(
  SELECTED_CVODE_SHORT_FINAL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_short_final"
)
set(
  SELECTED_CVODE_H2O2_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_h2o2_full"
)
set(
  SELECTED_CVODE_UNAVAILABLE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_cvode_unavailable"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_CVODE_FIXTURE_A_WORK_DIRECTORY}"
  "${SELECTED_CVODE_FIXTURE_B_WORK_DIRECTORY}"
  "${SELECTED_CVODE_SHORT_FINAL_WORK_DIRECTORY}"
  "${SELECTED_CVODE_H2O2_FULL_WORK_DIRECTORY}"
  "${SELECTED_CVODE_UNAVAILABLE_WORK_DIRECTORY}"
)

if(PELEF_ENABLE_SUNDIALS)
  add_executable(
    test_sundials_constant_volume_reactor
    unit/test_sundials_constant_volume_reactor.F90
  )
  target_link_libraries(
    test_sundials_constant_volume_reactor
    PRIVATE fixture_mechanism pelef_selected_cvode_runtime
  )
  add_test(
    NAME unit_sundials_constant_volume_reactor
    COMMAND test_sundials_constant_volume_reactor
  )
  set_tests_properties(
    unit_sundials_constant_volume_reactor PROPERTIES TIMEOUT 180
  )

  add_executable(
    test_sundials_constant_volume_multi_context
    unit/test_sundials_constant_volume_multi_context.F90
  )
  target_link_libraries(
    test_sundials_constant_volume_multi_context
    PRIVATE fixture_mechanism pelef_selected_cvode_runtime
  )
  add_test(
    NAME unit_sundials_constant_volume_multi_context
    COMMAND test_sundials_constant_volume_multi_context
  )
  set_tests_properties(
    unit_sundials_constant_volume_multi_context PROPERTIES TIMEOUT 180
  )

  add_test(
    NAME regression_selected_cvode_fixture_run_a
    COMMAND $<TARGET_FILE:test_selected_fixture_reactor>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactor/fixture_cvode.nml"
  )
  add_test(
    NAME regression_selected_cvode_fixture_run_b
    COMMAND $<TARGET_FILE:test_selected_fixture_reactor>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactor/fixture_cvode.nml"
  )
  set_tests_properties(
    regression_selected_cvode_fixture_run_a
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_CVODE_FIXTURE_A_WORK_DIRECTORY}"
      PASS_REGULAR_EXPRESSION "CVODE internal steps:"
      TIMEOUT 180
  )
  set_tests_properties(
    regression_selected_cvode_fixture_run_b
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_CVODE_FIXTURE_B_WORK_DIRECTORY}"
      PASS_REGULAR_EXPRESSION "CVODE internal steps:"
      TIMEOUT 180
  )
  add_test(
    NAME regression_selected_cvode_fixture_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_reactor_fixture.py"
      --input
        "${SELECTED_CVODE_FIXTURE_A_WORK_DIRECTORY}/selected_fixture_cvode.csv"
  )
  set_tests_properties(
    regression_selected_cvode_fixture_check
    PROPERTIES DEPENDS regression_selected_cvode_fixture_run_a
  )
  add_test(
    NAME regression_selected_cvode_fixture_deterministic
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_CVODE_FIXTURE_A_WORK_DIRECTORY}/selected_fixture_cvode.csv"
      "${SELECTED_CVODE_FIXTURE_B_WORK_DIRECTORY}/selected_fixture_cvode.csv"
  )
  set_tests_properties(
    regression_selected_cvode_fixture_deterministic
    PROPERTIES DEPENDS
      "regression_selected_cvode_fixture_run_a;regression_selected_cvode_fixture_run_b"
  )

  add_test(
    NAME regression_selected_cvode_short_final_run
    COMMAND $<TARGET_FILE:test_selected_fixture_reactor>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactor/short_final_step_cvode.nml"
  )
  set_tests_properties(
    regression_selected_cvode_short_final_run
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_CVODE_SHORT_FINAL_WORK_DIRECTORY}"
      PASS_REGULAR_EXPRESSION "CVODE internal steps:"
      TIMEOUT 180
  )
  add_test(
    NAME regression_selected_cvode_short_final_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_reactor_fixture.py"
      --input
        "${SELECTED_CVODE_SHORT_FINAL_WORK_DIRECTORY}/selected_short_final_cvode.csv"
      --final-time 1.000005e-7
      --expected-rows 4
  )
  set_tests_properties(
    regression_selected_cvode_short_final_check
    PROPERTIES DEPENDS regression_selected_cvode_short_final_run
  )

  add_test(
    NAME regression_selected_cvode_h2o2_full_run
    COMMAND $<TARGET_FILE:test_selected_full_reactor>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactor/h2o2_full_cvode.nml"
  )
  set_tests_properties(
    regression_selected_cvode_h2o2_full_run
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_CVODE_H2O2_FULL_WORK_DIRECTORY}"
      PASS_REGULAR_EXPRESSION "CVODE Jacobian evaluations:"
      TIMEOUT 300
  )
  add_test(
    NAME regression_selected_cvode_h2o2_full_structure
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_zero_d_h2o2_full.py"
      --input
        "${SELECTED_CVODE_H2O2_FULL_WORK_DIRECTORY}/selected_h2o2_full_cvode.csv"
  )
  set_tests_properties(
    regression_selected_cvode_h2o2_full_structure
    PROPERTIES DEPENDS regression_selected_cvode_h2o2_full_run
  )
else()
  add_test(
    NAME regression_selected_cvode_unavailable_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
      --executable $<TARGET_FILE:test_selected_fixture_reactor>
      --input "${PROJECT_SOURCE_DIR}/cases/selected_reactor/fixture_cvode.nml"
      --expected "CVODE integrator is unavailable"
      --forbidden-output
        "${SELECTED_CVODE_UNAVAILABLE_WORK_DIRECTORY}/selected_fixture_cvode.csv"
  )
  set_tests_properties(
    regression_selected_cvode_unavailable_rejected
    PROPERTIES WORKING_DIRECTORY "${SELECTED_CVODE_UNAVAILABLE_WORK_DIRECTORY}"
  )
endif()

set(
  SELECTED_FIXTURE_A_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactor_fixture_a"
)
set(
  SELECTED_FIXTURE_B_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactor_fixture_b"
)
set(
  SELECTED_FIXTURE_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactor_fixture_invalid"
)
set(
  SELECTED_FIXTURE_SHORT_FINAL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactor_fixture_short_final"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_FIXTURE_A_WORK_DIRECTORY}"
  "${SELECTED_FIXTURE_B_WORK_DIRECTORY}"
  "${SELECTED_FIXTURE_INVALID_WORK_DIRECTORY}"
  "${SELECTED_FIXTURE_SHORT_FINAL_WORK_DIRECTORY}"
)
add_test(
  NAME regression_selected_reactor_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactor>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactor/fixture.nml"
)
add_test(
  NAME regression_selected_reactor_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactor>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactor/fixture.nml"
)
set_tests_properties(
  regression_selected_reactor_fixture_run_a
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_FIXTURE_A_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
set_tests_properties(
  regression_selected_reactor_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_FIXTURE_B_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactor_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactor_fixture.py"
    --input "${SELECTED_FIXTURE_A_WORK_DIRECTORY}/selected_fixture.csv"
)
set_tests_properties(
  regression_selected_reactor_fixture_check
  PROPERTIES DEPENDS regression_selected_reactor_fixture_run_a
)
add_test(
  NAME regression_selected_reactor_fixture_deterministic
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_FIXTURE_A_WORK_DIRECTORY}/selected_fixture.csv"
    "${SELECTED_FIXTURE_B_WORK_DIRECTORY}/selected_fixture.csv"
)
set_tests_properties(
  regression_selected_reactor_fixture_deterministic
  PROPERTIES DEPENDS
    "regression_selected_reactor_fixture_run_a;regression_selected_reactor_fixture_run_b"
)
add_test(
  NAME regression_selected_reactor_short_final_run
  COMMAND $<TARGET_FILE:test_selected_fixture_reactor>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactor/short_final_step.nml"
)
set_tests_properties(
  regression_selected_reactor_short_final_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_FIXTURE_SHORT_FINAL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactor_short_final_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactor_fixture.py"
    --input
      "${SELECTED_FIXTURE_SHORT_FINAL_WORK_DIRECTORY}/selected_short_final.csv"
    --final-time 1.000005e-7
    --expected-rows 4
)
set_tests_properties(
  regression_selected_reactor_short_final_check
  PROPERTIES DEPENDS regression_selected_reactor_short_final_run
)
add_test(
  NAME regression_selected_reactor_unknown_species_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactor>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactor/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_FIXTURE_INVALID_WORK_DIRECTORY}/selected_invalid_unknown.csv"
)
add_test(
  NAME regression_selected_reactor_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactor>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactor/invalid_temperature.nml"
    --expected "outside the mechanism range"
    --forbidden-output
      "${SELECTED_FIXTURE_INVALID_WORK_DIRECTORY}/selected_invalid_temperature.csv"
)
set_tests_properties(
  regression_selected_reactor_unknown_species_rejected
  regression_selected_reactor_temperature_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_FIXTURE_INVALID_WORK_DIRECTORY}"
)

