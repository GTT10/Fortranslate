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
  ${PELEF_TEST_SELECTED_FULL_MPI_OPTION}
  ADD_TEST
)
pelef_add_mechanism_bundle(
  NAME selected_elementary_fixture
  BUNDLE "${PELEF_TEST_ELEMENTARY_SELECTED_BUNDLE}"
  PROBE_NAME test_selected_elementary_bundle
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

