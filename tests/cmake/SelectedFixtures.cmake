set(PELEF_TEST_SELECTED_FIXTURE_MPI_OPTION)
set(PELEF_TEST_SELECTED_FULL_MPI_OPTION)
if(PELEF_ENABLE_MPI)
  list(
    APPEND PELEF_TEST_SELECTED_FIXTURE_MPI_OPTION
    MPI_REACTIVE_1D_NAME test_selected_fixture_mpi_reactive_1d
    MPI_AMR_REACTIVE_1D_NAME test_selected_fixture_mpi_amr_reactive_1d
    MPI_AMR_REACTIVE_3D_NAME test_selected_fixture_mpi_amr_reactive_3d
    MPI_REACTIVE_EB_PATCH_TREE_2D_NAME
      test_selected_fixture_mpi_reactive_eb_patch_tree_2d
  )
  list(
    APPEND PELEF_TEST_SELECTED_FULL_MPI_OPTION
    MPI_REACTIVE_EB_TREE_3D_NAME test_selected_full_mpi_reactive_eb_tree_3d
    MPI_REACTIVE_1D_NAME test_selected_full_mpi_reactive_1d
    MPI_AMR_REACTIVE_1D_NAME test_selected_full_mpi_amr_reactive_1d
    MPI_AMR_REACTIVE_3D_NAME test_selected_full_mpi_amr_reactive_3d
    MPI_REACTIVE_EB_PATCH_TREE_2D_NAME
      test_selected_full_mpi_reactive_eb_patch_tree_2d
  )
endif()

set(
  PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE
  "${CMAKE_CURRENT_BINARY_DIR}/h2o2_elementary_selected_bundle.json"
)
execute_process(
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compose_h2o2_elementary_bundle.py"
    --kinetics "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_elementary.json"
    --database "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
    --output "${PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE}"
  RESULT_VARIABLE PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE_RESULT
  ERROR_VARIABLE PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE_ERROR
)
if(NOT PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE_RESULT EQUAL 0)
  message(
    FATAL_ERROR
    "failed to compose elementary selected bundle: "
    "${PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE_ERROR}"
  )
endif()
set_property(
  DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/tools/compose_h2o2_elementary_bundle.py"
    "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_elementary.json"
    "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
)

pelef_add_mechanism_bundle(
  NAME fixture
  BUNDLE "${PROJECT_SOURCE_DIR}/tests/fixtures/mechanism_bundle.json"
  PROBE_NAME test_mechanism_bundle
  REACTOR_NAME test_selected_fixture_reactor
  REACTIVE_1D_NAME test_selected_fixture_reactive_1d
  AMR_REACTIVE_1D_NAME test_selected_fixture_amr_reactive_1d
  REACTIVE_2D_NAME test_selected_fixture_reactive_2d
  REACTIVE_EB_2D_NAME test_selected_fixture_reactive_eb_2d
  REACTIVE_EB_AMR_2D_NAME test_selected_fixture_reactive_eb_amr_2d
  REACTIVE_3D_NAME test_selected_fixture_reactive_3d
  AMR_REACTIVE_3D_NAME test_selected_fixture_amr_reactive_3d
  REACTIVE_EB_3D_NAME test_selected_fixture_reactive_eb_3d
  ${PELEF_TEST_SELECTED_FIXTURE_MPI_OPTION}
  ADD_TEST
  REQUIRE_ACTIVITY
)
pelef_add_mechanism_bundle(
  NAME selected_full_fixture
  BUNDLE "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
  SOURCE "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml"
  PROBE_NAME test_selected_full_bundle
  REACTOR_NAME test_selected_full_reactor
  REACTIVE_1D_NAME test_selected_full_reactive_1d
  AMR_REACTIVE_1D_NAME test_selected_full_amr_reactive_1d
  REACTIVE_2D_NAME test_selected_full_reactive_2d
  REACTIVE_EB_2D_NAME test_selected_full_reactive_eb_2d
  REACTIVE_EB_AMR_2D_NAME test_selected_full_reactive_eb_amr_2d
  REACTIVE_3D_NAME test_selected_full_reactive_3d
  AMR_REACTIVE_3D_NAME test_selected_full_amr_reactive_3d
  REACTIVE_EB_3D_NAME test_selected_full_reactive_eb_3d
  REACTIVE_EB_TWO_LEVEL_3D_NAME test_selected_full_reactive_eb_two_level_3d
  REACTIVE_EB_TREE_3D_NAME test_selected_full_reactive_eb_tree_3d
  ${PELEF_TEST_SELECTED_FULL_MPI_OPTION}
  ADD_TEST
)
set(PELEF_TEST_SELECTED_ELEMENTARY_TREE_MPI_OPTION)
if(PELEF_ENABLE_MPI)
  set(PELEF_TEST_SELECTED_ELEMENTARY_TREE_MPI_OPTION
    MPI_AMR_TREE_REACTIVE_3D_NAME test_selected_elementary_mpi_amr_tree_3d)
endif()
pelef_add_mechanism_bundle(
  NAME selected_elementary_fixture
  BUNDLE "${PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE}"
  PROBE_NAME test_selected_elementary_bundle
  AMR_TREE_REACTIVE_3D_NAME test_selected_elementary_amr_tree_3d
  ${PELEF_TEST_SELECTED_ELEMENTARY_TREE_MPI_OPTION}
  REACTIVE_EB_3D_NAME test_selected_elementary_reactive_eb_3d
  ADD_TEST
)
pelef_add_mechanism_bundle(
  NAME selected_unsupported_element_fixture
  BUNDLE "${PROJECT_SOURCE_DIR}/tests/fixtures/unsupported_element_bundle.json"
  PROBE_NAME test_selected_unsupported_element_bundle
  REACTIVE_EB_3D_NAME test_selected_unsupported_element_reactive_eb_3d
)

add_executable(
  test_selected_reactive_eb_3d_config
  unit/test_selected_reactive_eb_3d_config.F90
)
target_link_libraries(
  test_selected_reactive_eb_3d_config
  PRIVATE pelef_selected_reactive_eb_3d_runtime
)
add_test(
  NAME unit_selected_reactive_eb_3d_config
  COMMAND test_selected_reactive_eb_3d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture.nml"
)


add_test(NAME regression_selected_amr_tree_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_tree_3d.py"
    --selected "$<TARGET_FILE:test_selected_elementary_amr_tree_3d>"
    --fixed "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
    --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_tree_3d")
set_tests_properties(regression_selected_amr_tree_3d PROPERTIES TIMEOUT 240)

add_test(NAME regression_amr_tree_restart_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_restart_3d.py"
    --fixed "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
    --selected "$<TARGET_FILE:test_selected_elementary_amr_tree_3d>"
    --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/amr_tree_restart_3d")
set_tests_properties(regression_amr_tree_restart_3d PROPERTIES TIMEOUT 240)

add_test(NAME regression_amr_tree_adaptation_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_adaptation_3d.py"
    --fixed "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
    --selected "$<TARGET_FILE:test_selected_elementary_amr_tree_3d>"
    --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/amr_tree_adaptation_3d")
set_tests_properties(regression_amr_tree_adaptation_3d PROPERTIES TIMEOUT 240)

add_test(NAME regression_amr_tree_subcycle_restart_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_restart_3d.py"
    --fixed "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
    --selected "$<TARGET_FILE:test_selected_elementary_amr_tree_3d>"
    --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/amr_tree_subcycle_restart_3d" --subcycling)
set_tests_properties(regression_amr_tree_subcycle_restart_3d PROPERTIES TIMEOUT 300)

foreach(tree_pressure_mode IN ITEMS common subcycle)
  set(tree_pressure_options --pressure-outlet)
  if(tree_pressure_mode STREQUAL "subcycle")
    list(APPEND tree_pressure_options --subcycling)
  endif()
  add_test(NAME regression_amr_tree_pressure_${tree_pressure_mode}_restart_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_restart_3d.py"
      --fixed "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
      --selected "$<TARGET_FILE:test_selected_elementary_amr_tree_3d>"
      --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/amr_tree_pressure_${tree_pressure_mode}_restart_3d"
      ${tree_pressure_options})
  set_tests_properties(regression_amr_tree_pressure_${tree_pressure_mode}_restart_3d PROPERTIES TIMEOUT 300)
endforeach()


add_test(
  NAME regression_selected_eb_two_level_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_selected_eb_two_level_3d.py"
    --exe $<TARGET_FILE:test_selected_full_reactive_eb_two_level_3d>
    --case "${PROJECT_SOURCE_DIR}/cases/reactive_eb_two_level_3d/selected_h2o2.nml"
    --bundle "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_two_level_3d"
)

set(PELEF_TEST_EB_METHANOL_BUNDLE "" CACHE FILEPATH "Optional pinned FFCM1_Red bundle for EB regression")
if(PELEF_TEST_EB_METHANOL_BUNDLE)
  file(READ "${PELEF_TEST_EB_METHANOL_BUNDLE}" eb_methanol_json)
  file(READ "${PROJECT_SOURCE_DIR}/references/spray_methanol.json" eb_methanol_reference)
  string(JSON eb_methanol_sha GET "${eb_methanol_json}" source sha256)
  string(JSON eb_methanol_expected_sha GET "${eb_methanol_reference}" files yaml sha256)
  if(NOT eb_methanol_sha STREQUAL eb_methanol_expected_sha)
    message(FATAL_ERROR "EB methanol test requires the pinned FFCM1_Red source")
  endif()
pelef_add_mechanism_bundle(
  NAME selected_eb_methanol_fixture
  REACTIVE_EB_3D_NAME test_selected_methanol_eb_3d
  BUNDLE "${PELEF_TEST_EB_METHANOL_BUNDLE}"
  REACTIVE_EB_TWO_LEVEL_3D_NAME test_selected_methanol_eb_two_level_3d
)
add_test(
  NAME regression_selected_methanol_eb_two_level_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_selected_eb_two_level_3d.py"
    --exe $<TARGET_FILE:test_selected_methanol_eb_two_level_3d>
    --case "${PROJECT_SOURCE_DIR}/cases/reactive_eb_two_level_3d/selected_methanol.nml"
    --bundle "${PELEF_TEST_EB_METHANOL_BUNDLE}"
    --species-count 21 --fuel CH3OH
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_methanol_eb_two_level_3d"
)
add_test(
  NAME regression_selected_methanol_eb_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_selected_carbon_eb_3d.py"
    --exe $<TARGET_FILE:test_selected_methanol_eb_3d>
    --case "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/methanol.nml"
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_methanol_eb_3d"
)
endif()

add_test(NAME regression_selected_scheduled_eb_two_level_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_scheduled_eb_two_level_3d.py"
    --exe $<TARGET_FILE:test_selected_full_reactive_eb_two_level_3d>
    --case "${PROJECT_SOURCE_DIR}/cases/reactive_eb_two_level_3d/selected_scheduled.nml"
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_scheduled_eb_two_level_3d")
set_tests_properties(regression_selected_scheduled_eb_two_level_3d PROPERTIES TIMEOUT 180)

add_test(NAME regression_selected_coarse_gap_eb_two_level_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_scheduled_eb_two_level_3d.py"
    --exe $<TARGET_FILE:test_selected_full_reactive_eb_two_level_3d>
    --case "${PROJECT_SOURCE_DIR}/cases/reactive_eb_two_level_3d/selected_scheduled.nml"
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_coarse_gap_eb_two_level_3d" --coarse-gap)
set_tests_properties(regression_selected_coarse_gap_eb_two_level_3d PROPERTIES TIMEOUT 180)

add_test(NAME regression_selected_adaptive_eb_two_level_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_adaptive_eb_two_level_3d.py"
    --exe $<TARGET_FILE:test_selected_full_reactive_eb_two_level_3d>
    --case "${PROJECT_SOURCE_DIR}/cases/reactive_eb_two_level_3d/selected_adaptive.nml"
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_adaptive_eb_two_level_3d")
set_tests_properties(regression_selected_adaptive_eb_two_level_3d PROPERTIES TIMEOUT 180)

add_test(NAME regression_selected_eb_tree_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_eb_tree_application_3d.py"
    --exe "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
    --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/selected.nml"
    --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_3d")
set_tests_properties(regression_selected_eb_tree_3d PROPERTIES TIMEOUT 180)

if(PELEF_ENABLE_MPI)
  add_test(NAME regression_selected_eb_tree_stop_roundoff
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_stop_roundoff.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_stop_roundoff")
  set_tests_properties(regression_selected_eb_tree_stop_roundoff
    PROPERTIES PROCESSORS 3 TIMEOUT 240 LABELS "amr;eb;application;restart")
  add_test(NAME regression_selected_eb_tree_mol_cell_trace
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_mol_cell_trace.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_mol_cell_trace")
  set_tests_properties(regression_selected_eb_tree_mol_cell_trace
    PROPERTIES PROCESSORS 2 TIMEOUT 180 LABELS "amr;eb;application;diagnostics")
  add_test(NAME regression_selected_eb_tree_initial_state
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_initial_state.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_initial_state")
  set_tests_properties(regression_selected_eb_tree_initial_state
    PROPERTIES PROCESSORS 3 TIMEOUT 180 LABELS "amr;eb;application;initialization")
  if(PELEF_ENABLE_SUNDIALS)
    add_test(NAME regression_selected_eb_tree_hot_cvode_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_hot_cvode_restart" --hot-cvode)
    set_tests_properties(regression_selected_eb_tree_hot_cvode_restart
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application")
    add_test(NAME regression_selected_eb_tree_full_cvode_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_full_cvode_restart" --full-cvode)
    set_tests_properties(regression_selected_eb_tree_full_cvode_restart
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application")
    add_test(NAME regression_selected_eb_tree_converged_cvode_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_converged_cvode_restart"
        --hot-cvode --pinned-eos-convergence)
    set_tests_properties(regression_selected_eb_tree_converged_cvode_restart
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;eos")
    add_test(NAME regression_selected_eb_tree_pressure_restart_common
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_pressure_restart_common"
        --hot-cvode --pinned-eos-convergence --paired-pressure)
    set_tests_properties(regression_selected_eb_tree_pressure_restart_common
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;pressure")
add_test(NAME regression_selected_eb_tree_pressure_restart_recursive
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_pressure_restart_recursive"
        --hot-cvode --pinned-eos-convergence --paired-pressure --subcycling)
    set_tests_properties(regression_selected_eb_tree_pressure_restart_recursive
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;pressure")
    add_test(NAME regression_selected_eb_tree_connected_restart_common
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_connected_restart_common"
        --hot-cvode --pinned-eos-convergence --paired-pressure --connected-directions)
    set_tests_properties(regression_selected_eb_tree_connected_restart_common
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;connected")
add_test(NAME regression_selected_eb_tree_connected_restart_recursive
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_connected_restart_recursive"
        --hot-cvode --pinned-eos-convergence --paired-pressure --connected-directions --subcycling)
    set_tests_properties(regression_selected_eb_tree_connected_restart_recursive
      PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;connected")
    foreach(momentum_clock IN ITEMS common recursive)
      set(momentum_clock_args)
      if(momentum_clock STREQUAL "recursive")
        list(APPEND momentum_clock_args --subcycling)
      endif()
      add_test(NAME regression_selected_eb_tree_momentum_restart_${momentum_clock}
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
          --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
          --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
          --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_momentum_restart_${momentum_clock}"
          --hot-cvode --pinned-eos-convergence --separate-momentum ${momentum_clock_args})
      set_tests_properties(regression_selected_eb_tree_momentum_restart_${momentum_clock}
        PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;momentum")
      add_test(NAME regression_selected_eb_tree_velocity_restart_${momentum_clock}
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
          --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
          --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
          --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_velocity_restart_${momentum_clock}"
          --hot-cvode --pinned-eos-convergence --consistent-velocity ${momentum_clock_args})
      set_tests_properties(regression_selected_eb_tree_velocity_restart_${momentum_clock}
        PROPERTIES PROCESSORS 3 TIMEOUT 300 LABELS "chemistry;cvode;amr;eb;restart;application;velocity")
    endforeach()
  endif()
  foreach(reacting_clock IN ITEMS common recursive)
    set(reacting_clock_args)
    if(reacting_clock STREQUAL "recursive")
      list(APPEND reacting_clock_args --subcycling)
    endif()
    add_test(NAME regression_selected_eb_tree_reacting_restart_${reacting_clock}
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_reacting_restart.py"
        --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
        --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_reacting_restart_${reacting_clock}"
        ${reacting_clock_args})
    set_tests_properties(regression_selected_eb_tree_reacting_restart_${reacting_clock}
      PROPERTIES PROCESSORS 3 TIMEOUT 900 LABELS "chemistry;amr;eb;restart;application")
  endforeach()
  add_test(NAME regression_selected_eb_tree_held_application
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_held_application" --selected --chemistry --subcycling)
  set_tests_properties(regression_selected_eb_tree_held_application PROPERTIES PROCESSORS 3 TIMEOUT 600)
  add_test(NAME regression_selected_eb_tree_godunov_held_application
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_godunov_held_application"
      --selected --godunov --chemistry --iterations 2 --state-redist --transport-off)
  set_tests_properties(regression_selected_eb_tree_godunov_held_application PROPERTIES PROCESSORS 3 TIMEOUT 600)
  add_test(NAME regression_selected_eb_tree_godunov_held_recursive_smoke
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_tree_held_app.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_eb_tree_godunov_held_recursive_smoke"
      --selected --godunov --chemistry --iterations 2 --state-redist --transport-off --subcycling --smoke-only)
  set_tests_properties(regression_selected_eb_tree_godunov_held_recursive_smoke PROPERTIES TIMEOUT 180)
  add_test(NAME regression_selected_mpi_eb_tree_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py"
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/selected.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_eb_tree_3d" --selected)
  set_tests_properties(regression_selected_mpi_eb_tree_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)

  add_test(NAME regression_selected_mpi_eb_tree_rebalance_3d
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_mpi_eb_tree_application_3d.py" --rebalance
      --serial "$<TARGET_FILE:test_selected_full_reactive_eb_tree_3d>"
      --mpi "$<TARGET_FILE:test_selected_full_mpi_reactive_eb_tree_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --input "${PROJECT_SOURCE_DIR}/cases/reactive_eb_tree_3d/selected.nml"
      --work "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_eb_tree_rebalance_3d" --selected)
  set_tests_properties(regression_selected_mpi_eb_tree_rebalance_3d PROPERTIES TIMEOUT 180 PROCESSORS 3)
endif()

add_test(NAME regression_amr_tree_mach_restart_3d
  COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_amr_tree_restart_3d.py"
    --fixed "$<TARGET_FILE:pelef_amr_tree_reactive_3d>"
    --selected "$<TARGET_FILE:test_selected_elementary_amr_tree_3d>"
    --case-dir "${PROJECT_SOURCE_DIR}/cases/amr_tree_reactive_3d"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/amr_tree_mach_restart_3d" --mach-inlet --subcycling)
set_tests_properties(regression_amr_tree_mach_restart_3d PROPERTIES TIMEOUT 300)
# Small executable for importing pinned reference states into the production
# tree integrator. External plotfiles are supplied explicitly, not fetched by CI.
add_executable(test_amr_tree_mol_reference
  "${PROJECT_SOURCE_DIR}/tools/reference_probes/amr_tree_mol_probe.F90"
  "${PROJECT_SOURCE_DIR}/src/physics/thermo_database_mod.F90")
set_target_properties(test_amr_tree_mol_reference PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/amr_mol_reference_modules")
target_link_libraries(test_amr_tree_mol_reference PRIVATE pelef_selected_amr_tree_reactive_3d_runtime)

if(PELEF_ENABLE_SUNDIALS)
  set(PELEF_TEST_LIDRYER_EB_MPI_OPTION)
  if(PELEF_ENABLE_MPI)
    set(PELEF_TEST_LIDRYER_EB_MPI_OPTION
      MPI_REACTIVE_EB_TREE_3D_NAME test_lidryer_mpi_reactive_eb_tree_3d)
  endif()
  pelef_add_mechanism_bundle(
    NAME lidryer_reference
    BUNDLE "${PROJECT_SOURCE_DIR}/mechanisms/lidryer.json"
    SOURCE "${PROJECT_SOURCE_DIR}/mechanisms/lidryer.yaml"
    PROBE_NAME test_lidryer_reference_bundle
    REACTIVE_EB_TREE_3D_NAME test_lidryer_reactive_eb_tree_3d
    ${PELEF_TEST_LIDRYER_EB_MPI_OPTION})
  add_executable(test_amr_tree_reacting_reference
    "${PROJECT_SOURCE_DIR}/tools/reference_probes/amr_tree_reacting_probe.F90")
  target_link_libraries(test_amr_tree_reacting_reference PRIVATE
    lidryer_reference_mechanism pelef_selected_amr_tree_reactive_3d_runtime pelef_selected_cvode_runtime)
  if(PELEF_ENABLE_MPI)
    add_executable(test_mpi_amr_tree_reacting_reference
      "${PROJECT_SOURCE_DIR}/tools/reference_probes/amr_tree_reacting_probe.F90")
    target_compile_definitions(test_mpi_amr_tree_reacting_reference PRIVATE PELEF_REFERENCE_MPI)
    target_link_libraries(test_mpi_amr_tree_reacting_reference PRIVATE lidryer_reference_mechanism
      pelef_selected_mpi_amr_tree_reactive_3d_runtime pelef_selected_cvode_runtime MPI::MPI_Fortran)
  endif()
  add_executable(test_eb_tree_reacting_reference
    "${PROJECT_SOURCE_DIR}/tools/reference_probes/eb_tree_reacting_probe.F90")
  target_link_libraries(test_eb_tree_reacting_reference PRIVATE
    lidryer_reference_mechanism pelef_selected_reactive_eb_3d_runtime pelef_selected_cvode_runtime)
  if(PELEF_ENABLE_MPI)
    add_executable(test_eb_tree_initial_reaction unit/test_eb_tree_initial_reaction.F90)
    target_link_libraries(test_eb_tree_initial_reaction PRIVATE lidryer_reference_mechanism
      pelef_selected_mpi_eb_tree_3d_runtime MPI::MPI_Fortran)
    add_test(NAME unit_eb_tree_initial_reaction_3
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_eb_tree_initial_reaction> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_eb_tree_initial_reaction_3 PROPERTIES TIMEOUT 60 PROCESSORS 3)
    add_executable(test_amr_tree_initial_reaction unit/test_amr_tree_initial_reaction.F90)
    target_link_libraries(test_amr_tree_initial_reaction PRIVATE lidryer_reference_mechanism
      pelef_selected_mpi_amr_tree_reactive_3d_runtime MPI::MPI_Fortran)
    add_test(NAME unit_amr_tree_initial_reaction_3
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_amr_tree_initial_reaction> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(unit_amr_tree_initial_reaction_3 PROPERTIES TIMEOUT 60 PROCESSORS 3)
  endif()
endif()
