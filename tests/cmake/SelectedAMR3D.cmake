add_executable(
  test_selected_reactive_3d_config
  unit/test_selected_reactive_3d_config.F90
)
target_link_libraries(
  test_selected_reactive_3d_config
  PRIVATE pelef_selected_reactive_3d_runtime
)
add_test(
  NAME unit_selected_reactive_3d_config
  COMMAND test_selected_reactive_3d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/fixture.nml"
)

set(
  SELECTED_REACTIVE_3D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_3d_fixture"
)
set(
  SELECTED_REACTIVE_3D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_3d_invalid"
)
set(
  SELECTED_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_3d_fixed_rejection"
)
set(
  SELECTED_REACTIVE_3D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_3d_full_chemistry_plm"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_3D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}"
)

add_test(
  NAME regression_selected_reactive_3d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_3d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_3d_fixture_run_a
  regression_selected_reactive_3d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_3d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_3d.py"
    --input
      "${SELECTED_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_3d.csv"
    --species H2 H
    --nx 4
    --ny 4
    --nz 4
    --final-time 1.0e-7
    --require-uniform
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_reactive_3d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_3d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_3d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_3d.csv"
    "${SELECTED_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_3d_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_3d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_3d_fixture_run_a;regression_selected_reactive_3d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_3d_full_transport_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/h2o2_full_transport.nml"
)
set_tests_properties(
  regression_selected_reactive_3d_full_transport_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 10"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_3d_full_transport_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_transport_3d.py"
    --control
      "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/reactive_transport_3d_control.csv"
    --transport
      "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/selected_reactive_transport_3d.csv"
    --nx 6
    --ny 6
    --nz 6
    --time 2.0e-6
)
set_tests_properties(
  regression_selected_reactive_3d_full_transport_structure
  PROPERTIES DEPENDS
    "regression_reactive_transport_3d_control_run;regression_selected_reactive_3d_full_transport_run"
)
add_test(
  NAME regression_selected_reactive_3d_full_transport_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/reactive_transport_3d.csv"
    "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/selected_reactive_transport_3d.csv"
)
set_tests_properties(
  regression_selected_reactive_3d_full_transport_exact
  PROPERTIES DEPENDS
    "regression_reactive_transport_3d_run;regression_selected_reactive_3d_full_transport_run"
)

add_test(
  NAME regression_selected_reactive_3d_full_chemistry_plm_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot_3d/hotspot_full_h2o2_plm.nml"
)
add_test(
  NAME regression_selected_reactive_3d_full_chemistry_plm_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/h2o2_full_chemistry_plm.nml"
)
set_tests_properties(
  regression_selected_reactive_3d_full_chemistry_plm_fixed_run
  regression_selected_reactive_3d_full_chemistry_plm_selected_run
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_3D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_3d_full_chemistry_plm_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_3d.py"
    --input
      "${SELECTED_REACTIVE_3D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}/selected_full_h2o2_chemistry_plm_3d.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 4
    --ny 4
    --nz 4
    --final-time 1.0e-9
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
set_tests_properties(
  regression_selected_reactive_3d_full_chemistry_plm_structure
  PROPERTIES DEPENDS
    regression_selected_reactive_3d_full_chemistry_plm_selected_run
)
add_test(
  NAME regression_selected_reactive_3d_full_chemistry_plm_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_3D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}/full_h2o2_chemistry_plm_3d.csv"
    "${SELECTED_REACTIVE_3D_FULL_CHEMISTRY_PLM_WORK_DIRECTORY}/selected_full_h2o2_chemistry_plm_3d.csv"
)
set_tests_properties(
  regression_selected_reactive_3d_full_chemistry_plm_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_3d_full_chemistry_plm_fixed_run;regression_selected_reactive_3d_full_chemistry_plm_selected_run"
)

add_test(
  NAME regression_selected_reactive_3d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_3d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_3d.csv"
)
set_tests_properties(
  regression_selected_reactive_3d_unknown_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_3d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_3d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_3d/fixture.nml"
    --expected "thermo_model"
    --forbidden-output
      "${SELECTED_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_3d.csv"
)
set_tests_properties(
  regression_fixed_reactive_3d_selected_model_rejected
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}"
)

set(
  SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_fixture"
)
set(
  SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_full"
)
set(
  SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_transport"
)
set(
  SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_restart"
)
set(
  SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_invalid"
)
set(
  SELECTED_AMR_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_path_alias.nml"
  "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_checkpoint_input_alias.nml"
  "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_input_alias.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_amr_reactive_3d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/fixture.nml"
)
add_test(
  NAME regression_selected_amr_reactive_3d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_fixture_run_a
  regression_selected_amr_reactive_3d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_amr_reactive_3d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
    --coarse
      "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d_coarse.csv"
    --fine
      "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d_fine.csv"
    --species H2 H
    --nx 4 --ny 4 --nz 4
    --i-lower 2 --i-upper 3
    --j-lower 2 --j-upper 3
    --k-lower 2 --k-upper 3
    --ratio 2
    --final-time 1.0e-7
    --require-uniform
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_amr_reactive_3d_fixture_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_3d_fixture_run_a
)
foreach(selected_amr_3d_level IN ITEMS coarse fine)
  add_test(
    NAME
      regression_selected_amr_reactive_3d_fixture_${selected_amr_3d_level}_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d_${selected_amr_3d_level}.csv"
      "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d_repeat_${selected_amr_3d_level}.csv"
  )
  set_tests_properties(
    regression_selected_amr_reactive_3d_fixture_${selected_amr_3d_level}_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_fixture_run_a;regression_selected_amr_reactive_3d_fixture_run_b"
  )
endforeach()

add_test(
  NAME regression_selected_amr_reactive_3d_full_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_chemistry.nml"
)
add_test(
  NAME regression_selected_amr_reactive_3d_full_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_chemistry.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_full_fixed_run
  regression_selected_amr_reactive_3d_full_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_3d_full_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
    --coarse
      "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_amr_hotspot_chemistry_3d_coarse.csv"
    --fine
      "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_amr_hotspot_chemistry_3d_fine.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8 --ny 8 --nz 8
    --i-lower 3 --i-upper 6
    --j-lower 3 --j-upper 6
    --k-lower 3 --k-upper 6
    --ratio 2
    --final-time 1.0e-9
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-12
)
set_tests_properties(
  regression_selected_amr_reactive_3d_full_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_3d_full_selected_run
)
foreach(selected_amr_3d_level IN ITEMS coarse fine)
  add_test(
    NAME regression_selected_amr_reactive_3d_full_${selected_amr_3d_level}_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/amr_hotspot_chemistry_3d_${selected_amr_3d_level}.csv"
      "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_amr_hotspot_chemistry_3d_${selected_amr_3d_level}.csv"
  )
  set_tests_properties(
    regression_selected_amr_reactive_3d_full_${selected_amr_3d_level}_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_full_fixed_run;regression_selected_amr_reactive_3d_full_selected_run"
  )
endforeach()

add_test(
  NAME regression_selected_amr_reactive_3d_transport_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport.nml"
)
add_test(
  NAME regression_selected_amr_reactive_3d_transport_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_fixed_run
  regression_selected_amr_reactive_3d_transport_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Completed transport fine substeps: 16"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_3d_transport_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
    --coarse
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_amr_hotspot_transport_3d_coarse.csv"
    --fine
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_amr_hotspot_transport_3d_fine.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8 --ny 8 --nz 8
    --i-lower 3 --i-upper 6
    --j-lower 3 --j-upper 6
    --k-lower 3 --k-upper 6
    --ratio 2
    --final-time 1.0e-9
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-12
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_3d_transport_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_3d_transport_exact
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
    --reference-prefix
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/amr_hotspot_transport_3d"
    --candidate-prefix
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_amr_hotspot_transport_3d"
    --final-time 1.0e-9
    --reconstruction characteristic_plm
    --limiter mc
    --expected-coarse-sha256
      $<IF:$<CONFIG:Debug>,5687cb97f834b759b3932bbd3f0b6e30fc5b48e6029e3f47cadfb3ec6a8d74fb,f0713f3a68f48b6537afaafd57209ea8ed4767f9c7d093fe2961bd3a1fa37fae>
    --expected-fine-sha256
      eaf132bbf89974d739c26a35f6b76137f9774438cfe9eebfc12436aa364cff58
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_3d_transport_fixed_run;regression_selected_amr_reactive_3d_transport_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_3d_restart_reference
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart_reference.nml"
)
add_test(
  NAME regression_selected_amr_reactive_3d_checkpoint_stop
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart_checkpoint.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_restart_reference
  regression_selected_amr_reactive_3d_checkpoint_stop
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_amr_reactive_3d_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
add_test(
  NAME regression_selected_amr_reactive_3d_restart_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_restart_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    DEPENDS regression_selected_amr_reactive_3d_checkpoint_stop
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_3d_restart_exact
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
    --reference-prefix
      "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart_reference"
    --candidate-prefix
      "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restarted"
    --checkpoint
      "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
    --final-time 1.0e-9
    --reconstruction characteristic_plm
    --limiter mc
    --expected-schema 3
    --expected-bundle-sha256
      f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
    --expected-integrator implicit
    --expected-checkpoint-sha256
      f44f69a7e1a1657b549b26373feb277840732e5453b0734321cc8c7cdf91c423
    --expected-coarse-sha256
      a03884837909bab719f0783d23ea9150a81f2afc628fee1a01c645c1b47e0ff3
    --expected-fine-sha256
      8d5a992d07b89832624e37af6ddacd510df63633707678693b9053057d4efcbe
)
set_tests_properties(
  regression_selected_amr_reactive_3d_restart_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_3d_restart_reference;regression_selected_amr_reactive_3d_checkpoint_stop;regression_selected_amr_reactive_3d_restart_run"
)

add_test(
  NAME regression_selected_amr_reactive_3d_restart_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_restart.nml"
    --expected "Could not open 3D AMR checkpoint"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart_coarse.csv"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart_fine.csv"
)
add_test(
  NAME regression_selected_amr_reactive_3d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    --input invalid_path_alias.nml
    --expected "input and output paths alias"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_alias_fine.csv"
)
add_test(
  NAME regression_selected_amr_reactive_3d_checkpoint_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    --input invalid_checkpoint_input_alias.nml
    --expected "checkpoint path aliases another file"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_alias_coarse.csv"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_alias_fine.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_restart_rejected
  regression_selected_amr_reactive_3d_path_alias_rejected
  regression_selected_amr_reactive_3d_checkpoint_alias_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_amr_reactive_3d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_amr_reactive_3d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/fixture.nml"
    --expected "thermo_model"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d_coarse.csv"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d_fine.csv"
)
set_tests_properties(
  regression_fixed_amr_reactive_3d_selected_model_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_AMR_REACTIVE_3D_FIXED_REJECTION_WORK_DIRECTORY}"
)

