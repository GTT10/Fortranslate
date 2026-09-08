add_executable(
  test_selected_reactive_2d_config
  unit/test_selected_reactive_2d_config.F90
)
target_link_libraries(
  test_selected_reactive_2d_config
  PRIVATE pelef_selected_reactive_2d_runtime
)
add_test(
  NAME unit_selected_reactive_2d_config
  COMMAND test_selected_reactive_2d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/fixture.nml"
)

set(
  SELECTED_REACTIVE_2D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_fixture"
)
set(
  SELECTED_REACTIVE_2D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_full"
)
set(
  SELECTED_REACTIVE_2D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_full_chemistry_plm"
)
set(
  SELECTED_REACTIVE_2D_WALL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_wall"
)
set(
  SELECTED_REACTIVE_2D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_invalid"
)
set(
  SELECTED_REACTIVE_2D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_2d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_2D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_2D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_2D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_2D_WALL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_2D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)

add_test(
  NAME regression_selected_reactive_2d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_2d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_2d_fixture_run_a
  regression_selected_reactive_2d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_2D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_2d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_2d.py"
    --input
      "${SELECTED_REACTIVE_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_2d.csv"
    --species H2 H
    --nx 4
    --ny 4
    --final-time 1.0e-7
    --require-uniform
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_reactive_2d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_2d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_2d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_2d.csv"
    "${SELECTED_REACTIVE_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_2d_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_2d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_2d_fixture_run_a;regression_selected_reactive_2d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_2d_full_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_full_h2o2_2d/uniform.nml"
)
add_test(
  NAME regression_selected_reactive_2d_full_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/h2o2_full.nml"
)
set_tests_properties(
  regression_selected_reactive_2d_full_fixed_run
  regression_selected_reactive_2d_full_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_2D_FULL_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_2d_full_selected_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Species: 10"
)
add_test(
  NAME regression_selected_reactive_2d_full_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_2d.py"
    --input
      "${SELECTED_REACTIVE_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_2d.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 4
    --ny 4
    --final-time 2.0e-6
    --require-uniform
)
set_tests_properties(
  regression_selected_reactive_2d_full_structure
  PROPERTIES DEPENDS regression_selected_reactive_2d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_2d_full_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_2D_FULL_WORK_DIRECTORY}/reactive_full_h2o2_2d.csv"
    "${SELECTED_REACTIVE_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_2d.csv"
)
set_tests_properties(
  regression_selected_reactive_2d_full_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_2d_full_fixed_run;regression_selected_reactive_2d_full_selected_run"
)

add_test(
  NAME regression_selected_reactive_2d_full_chemistry_plm_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot_2d/hotspot_full_h2o2_plm.nml"
)
add_test(
  NAME regression_selected_reactive_2d_full_chemistry_plm_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/h2o2_full_chemistry_plm.nml"
)
set_tests_properties(
  regression_selected_reactive_2d_full_chemistry_plm_fixed_run
  regression_selected_reactive_2d_full_chemistry_plm_selected_run
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_2D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_2d_full_chemistry_plm_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_2d.py"
    --input
      "${SELECTED_REACTIVE_2D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}/selected_full_h2o2_chemistry_plm_2d.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 4
    --ny 4
    --final-time 1.0e-9
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
set_tests_properties(
  regression_selected_reactive_2d_full_chemistry_plm_structure
  PROPERTIES DEPENDS
    regression_selected_reactive_2d_full_chemistry_plm_selected_run
)
add_test(
  NAME regression_selected_reactive_2d_full_chemistry_plm_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_2D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}/full_h2o2_chemistry_plm_2d.csv"
    "${SELECTED_REACTIVE_2D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}/selected_full_h2o2_chemistry_plm_2d.csv"
)
set_tests_properties(
  regression_selected_reactive_2d_full_chemistry_plm_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_2d_full_chemistry_plm_fixed_run;regression_selected_reactive_2d_full_chemistry_plm_selected_run"
)

add_test(
  NAME regression_selected_reactive_2d_wall_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_boundaries_2d/prescribed_species_wall_full_h2o2.nml"
)
add_test(
  NAME regression_selected_reactive_2d_wall_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/h2o2_full_prescribed_species_wall.nml"
)
set_tests_properties(
  regression_selected_reactive_2d_wall_fixed_run
  regression_selected_reactive_2d_wall_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_2D_WALL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Boundary y: slip_wall slip_wall"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_2d_wall_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_2d.py"
    --input
      "${SELECTED_REACTIVE_2D_WALL_WORK_DIRECTORY}/selected_full_h2o2_prescribed_species_wall_2d.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8
    --ny 8
    --final-time 5.0e-7
)
set_tests_properties(
  regression_selected_reactive_2d_wall_structure
  PROPERTIES DEPENDS regression_selected_reactive_2d_wall_selected_run
)
add_test(
  NAME regression_selected_reactive_2d_wall_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_2D_WALL_WORK_DIRECTORY}/full_h2o2_prescribed_species_wall_2d.csv"
    "${SELECTED_REACTIVE_2D_WALL_WORK_DIRECTORY}/selected_full_h2o2_prescribed_species_wall_2d.csv"
)
set_tests_properties(
  regression_selected_reactive_2d_wall_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_2d_wall_fixed_run;regression_selected_reactive_2d_wall_selected_run"
)

add_test(
  NAME regression_selected_reactive_2d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_2d.csv"
)
set_tests_properties(
  regression_selected_reactive_2d_unknown_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_REACTIVE_2D_INVALID_WORK_DIRECTORY}"
)
foreach(selected_reactive_2d_temperature_case IN ITEMS
    double_hotspot_temperature
    isothermal_ghost_temperature
)
  add_test(
    NAME
      regression_selected_reactive_2d_${selected_reactive_2d_temperature_case}_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
      --executable $<TARGET_FILE:test_selected_fixture_reactive_2d>
      --input
        "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/invalid_${selected_reactive_2d_temperature_case}.nml"
      --expected "temperature is outside the mechanism range"
      --forbidden-output
        "${SELECTED_REACTIVE_2D_INVALID_WORK_DIRECTORY}/invalid_${selected_reactive_2d_temperature_case}_2d.csv"
  )
  set_tests_properties(
    regression_selected_reactive_2d_${selected_reactive_2d_temperature_case}_rejected
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_REACTIVE_2D_INVALID_WORK_DIRECTORY}"
  )
endforeach()
add_test(
  NAME regression_fixed_reactive_2d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_2d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_2d/fixture.nml"
    --expected "Unknown reactive 2D chemistry model"
    --forbidden-output
      "${SELECTED_REACTIVE_2D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_2d.csv"
)
set_tests_properties(
  regression_fixed_reactive_2d_selected_model_rejected
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)

add_executable(
  test_selected_reactive_eb_2d_config
  unit/test_selected_reactive_eb_2d_config.F90
)
target_link_libraries(
  test_selected_reactive_eb_2d_config
  PRIVATE pelef_selected_reactive_eb_2d_runtime
)
add_test(
  NAME unit_selected_reactive_eb_2d_config
  COMMAND test_selected_reactive_eb_2d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/fixture.nml"
)

set(
  SELECTED_REACTIVE_EB_2D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_fixture"
)
set(
  SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_chemistry"
)
set(
  SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_transport"
)
set(
  SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_invalid"
)
set(
  SELECTED_REACTIVE_EB_2D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_2d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_EB_2D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/invalid_path_alias.nml"
  "${SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_reactive_eb_2d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_eb_2d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_fixture_run_a
  regression_selected_reactive_eb_2d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_2D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Cut cells: 4"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_eb_2d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_2d.csv"
    --species H2 H
    --nx 4
    --ny 4
    --final-time 1.0e-7
    --expected-regular 8
    --expected-cut 4
    --expected-covered 4
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
)
set_tests_properties(
  regression_selected_reactive_eb_2d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_2d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_2d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_2d.csv"
    "${SELECTED_REACTIVE_EB_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_2d_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_2d_fixture_run_a;regression_selected_reactive_eb_2d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_chemistry_fixed.nml"
)
add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_chemistry_selected.nml"
)
add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_regular_reference_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_chemistry_regular_selected.nml"
)
add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_inert_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_chemistry_inert_selected.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_chemistry_fixed_run
  regression_selected_reactive_eb_2d_chemistry_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  T"
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_2d_chemistry_regular_reference_run
  regression_selected_reactive_eb_2d_chemistry_inert_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_chemistry.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 4
    --ny 4
    --final-time 2.0e-7
    --expected-regular 8
    --expected-cut 4
    --expected-covered 4
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-8
)
set_tests_properties(
  regression_selected_reactive_eb_2d_chemistry_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_2d_chemistry_selected_run
)
add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}/full_h2o2_reactive_eb_2d_chemistry.csv"
    "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_chemistry.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_chemistry_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_2d_chemistry_fixed_run;regression_selected_reactive_eb_2d_chemistry_selected_run"
)
add_test(
  NAME regression_selected_reactive_eb_2d_chemistry_physics
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_chemistry_2d.py"
    --reference
      "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}/selected_full_h2o2_reactive_2d_chemistry_reference.csv"
    --reactive
      "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_chemistry.csv"
    --inert
      "${SELECTED_REACTIVE_EB_2D_CHEMISTRY_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_chemistry_inert.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_chemistry_physics
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_2d_chemistry_regular_reference_run;regression_selected_reactive_eb_2d_chemistry_inert_run;regression_selected_reactive_eb_2d_chemistry_selected_run"
)

add_test(
  NAME regression_selected_reactive_eb_2d_transport_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_transport_fixed.nml"
)
add_test(
  NAME regression_selected_reactive_eb_2d_transport_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_transport_selected.nml"
)
add_test(
  NAME regression_selected_reactive_eb_2d_transport_reference_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/h2o2_full_transport_reference_selected.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_transport_fixed_run
  regression_selected_reactive_eb_2d_transport_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Geometry: circle"
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_2d_transport_reference_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Geometry: circle"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_2d_transport_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_transport.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8
    --ny 8
    --final-time 2.0e-8
    --expected-regular 48
    --expected-cut 12
    --expected-covered 4
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
set_tests_properties(
  regression_selected_reactive_eb_2d_transport_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_2d_transport_selected_run
)
add_test(
  NAME regression_selected_reactive_eb_2d_transport_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}/full_h2o2_reactive_eb_2d_transport.csv"
    "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_transport.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_transport_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_2d_transport_fixed_run;regression_selected_reactive_eb_2d_transport_selected_run"
)
add_test(
  NAME regression_selected_reactive_eb_2d_transport_physics
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_transport_2d.py"
    --reference
      "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_transport_reference.csv"
    --transport
      "${SELECTED_REACTIVE_EB_2D_TRANSPORT_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_2d_transport.csv"
    --nx 8
    --ny 8
    --final-time 2.0e-8
    --allow-span-increase
)
set_tests_properties(
  regression_selected_reactive_eb_2d_transport_physics
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_2d_transport_reference_run;regression_selected_reactive_eb_2d_transport_selected_run"
)

add_test(
  NAME regression_selected_reactive_eb_2d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_2d.csv"
)
add_test(
  NAME regression_selected_reactive_eb_2d_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/invalid_embedded_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_temperature.csv"
)
add_test(
  NAME regression_selected_reactive_eb_2d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    --input invalid_path_alias.nml
    --expected "input and output paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_alias.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_2d_unknown_rejected
  regression_selected_reactive_eb_2d_temperature_rejected
  regression_selected_reactive_eb_2d_path_alias_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_2D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_eb_2d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_eb_2d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_2d/fixture.nml"
    --expected "Unknown reactive 2D chemistry model"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_2D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_2d.csv"
)
set_tests_properties(
  regression_fixed_reactive_eb_2d_selected_model_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)

