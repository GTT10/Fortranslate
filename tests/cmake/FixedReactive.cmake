

add_executable(
  test_full_h2o2_cfd_support
  unit/test_full_h2o2_cfd_support.F90
)
target_link_libraries(test_full_h2o2_cfd_support PRIVATE pelef_core)
add_test(NAME unit_full_h2o2_cfd_support COMMAND test_full_h2o2_cfd_support)

add_executable(
  test_amr_multilevel_reactive_1d
  regression/test_amr_multilevel_reactive_1d.F90
)
target_link_libraries(test_amr_multilevel_reactive_1d PRIVATE pelef_core)
add_test(
  NAME regression_amr_multilevel_reactive_1d
  COMMAND test_amr_multilevel_reactive_1d
)

add_executable(
  test_amr_multilevel_regrid_1d
  regression/test_amr_multilevel_regrid_1d.F90
)
target_link_libraries(test_amr_multilevel_regrid_1d PRIVATE pelef_core)
add_test(
  NAME regression_amr_multilevel_regrid_1d
  COMMAND test_amr_multilevel_regrid_1d
)

set(REACTIVE_FULL_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/reactive_full_h2o2")
file(MAKE_DIRECTORY "${REACTIVE_FULL_WORK_DIRECTORY}")
add_test(NAME regression_reactive_full_h2o2_1d_run
  COMMAND $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_full_h2o2_1d/uniform.nml")
set_tests_properties(regression_reactive_full_h2o2_1d_run PROPERTIES
  WORKING_DIRECTORY "${REACTIVE_FULL_WORK_DIRECTORY}" TIMEOUT 300)
add_test(NAME regression_reactive_full_h2o2_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_full_h2o2_2d/uniform.nml")
set_tests_properties(regression_reactive_full_h2o2_2d_run PROPERTIES
  WORKING_DIRECTORY "${REACTIVE_FULL_WORK_DIRECTORY}" TIMEOUT 300)
add_test(NAME regression_reactive_full_h2o2_parity
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_full_h2o2.py"
    --one-d "${REACTIVE_FULL_WORK_DIRECTORY}/reactive_full_h2o2_1d.csv"
    --two-d "${REACTIVE_FULL_WORK_DIRECTORY}/reactive_full_h2o2_2d.csv")
set_tests_properties(regression_reactive_full_h2o2_parity PROPERTIES
  DEPENDS "regression_reactive_full_h2o2_1d_run;regression_reactive_full_h2o2_2d_run")

set(REACTIVE_EB_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_circle_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_circle_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_circle_2d/uniform.nml"
)
set_tests_properties(
  regression_reactive_eb_circle_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_circle_2d_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_circle_2d.py"
    --input "${REACTIVE_EB_2D_WORK_DIRECTORY}/reactive_eb_circle_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_circle_2d_structure
  PROPERTIES DEPENDS regression_reactive_eb_circle_2d_run
)

set(REACTIVE_EB_TRANSPORT_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_transport_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_TRANSPORT_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_transport_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_transport_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_transport_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_transport_2d/transport.nml"
)
set_tests_properties(
  regression_reactive_eb_transport_2d_reference
  regression_reactive_eb_transport_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_TRANSPORT_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_transport_2d_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Molecular transport:  T"
)
add_test(
  NAME regression_reactive_eb_transport_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_transport_2d.py"
    --reference
      "${REACTIVE_EB_TRANSPORT_2D_WORK_DIRECTORY}/reactive_eb_transport_reference.csv"
    --transport
      "${REACTIVE_EB_TRANSPORT_2D_WORK_DIRECTORY}/reactive_eb_transport.csv"
)
set_tests_properties(
  regression_reactive_eb_transport_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_transport_2d_reference;regression_reactive_eb_transport_2d_run"
)

set(REACTIVE_EB_AMR_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_2d/uniform.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_2d.py"
    --coarse
      "${REACTIVE_EB_AMR_2D_WORK_DIRECTORY}/reactive_eb_amr_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_2D_WORK_DIRECTORY}/reactive_eb_amr_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_amr_2d_run
)

set(REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_patch_tree_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_patch_tree_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_2d/hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_patch_tree_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
    TIMEOUT 300
    PASS_REGULAR_EXPRESSION "Levels: 4"
)
add_test(
  NAME regression_reactive_eb_patch_tree_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_patch_tree_2d.py"
    --output
      "${REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/reactive_eb_patch_tree_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_patch_tree_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_patch_tree_2d_run
)

set(REACTIVE_EB_PATCH_TREE_BOUNDARY_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_patch_tree_boundary_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_PATCH_TREE_BOUNDARY_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_patch_tree_boundary_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_boundary_2d/hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_patch_tree_boundary_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_PATCH_TREE_BOUNDARY_2D_WORK_DIRECTORY}"
    TIMEOUT 300
    PASS_REGULAR_EXPRESSION "Levels: 4"
)
add_test(
  NAME regression_reactive_eb_patch_tree_boundary_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_patch_tree_2d.py"
    --output
      "${REACTIVE_EB_PATCH_TREE_BOUNDARY_2D_WORK_DIRECTORY}/reactive_eb_patch_tree_boundary_2d.csv"
    --require-x-upper-boundary
    --require-branching
)
set_tests_properties(
  regression_reactive_eb_patch_tree_boundary_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_patch_tree_boundary_2d_run
)

set(REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_patch_tree_restart_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_patch_tree_restart_2d_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --log
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/reference.log"
    --
    $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_patch_tree_restart_2d_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --log
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/stopped.log"
    --
    $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_patch_tree_restart_2d_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --log
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/restarted.log"
    --
    $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_restart_2d/restart.nml"
)
add_test(
  NAME regression_reactive_eb_patch_tree_restart_2d_incompatible
  COMMAND $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_restart_2d/incompatible_restart.nml"
)
set_tests_properties(
  regression_reactive_eb_patch_tree_restart_2d_reference
  regression_reactive_eb_patch_tree_restart_2d_checkpoint_stop
  regression_reactive_eb_patch_tree_restart_2d_restart
  regression_reactive_eb_patch_tree_restart_2d_incompatible
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_patch_tree_restart_2d_restart
  PROPERTIES DEPENDS
    regression_reactive_eb_patch_tree_restart_2d_checkpoint_stop
)
set_tests_properties(
  regression_reactive_eb_patch_tree_restart_2d_incompatible
  PROPERTIES
    DEPENDS regression_reactive_eb_patch_tree_restart_2d_checkpoint_stop
    WILL_FAIL TRUE
)
add_test(
  NAME regression_reactive_eb_patch_tree_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_patch_tree_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/patch_tree_restart.chk"
    --reference
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/patch_tree_restart_reference.csv"
    --stopped
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/patch_tree_restart_stopped.csv"
    --restarted
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/patch_tree_restart_restarted.csv"
    --reference-log
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/reference.log"
    --restarted-logs
      "${REACTIVE_EB_PATCH_TREE_RESTART_2D_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_reactive_eb_patch_tree_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_patch_tree_restart_2d_reference;regression_reactive_eb_patch_tree_restart_2d_checkpoint_stop;regression_reactive_eb_patch_tree_restart_2d_restart"
)

set(REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_transport_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_transport_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_transport_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_transport_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_transport_2d/transport.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_transport_2d_reference
  regression_reactive_eb_amr_transport_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_transport_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_transport_2d.py"
    --reference-coarse
      "${REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY}/reactive_eb_amr_transport_reference_coarse.csv"
    --reference-fine
      "${REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY}/reactive_eb_amr_transport_reference_fine.csv"
    --transport-coarse
      "${REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY}/reactive_eb_amr_transport_coarse.csv"
    --transport-fine
      "${REACTIVE_EB_AMR_TRANSPORT_2D_WORK_DIRECTORY}/reactive_eb_amr_transport_fine.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_transport_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_transport_2d_reference;regression_reactive_eb_amr_transport_2d_run"
)

set(REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_transport_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_transport_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_transport_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_transport_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_transport_2d/transport.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_transport_2d_reference
  regression_reactive_eb_amr_three_level_transport_2d_run
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_three_level_transport_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_transport_2d.py"
    --reference-root
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}/three_level_transport_reference_root.csv"
    --reference-middle
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}/three_level_transport_reference_middle.csv"
    --reference-finest
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}/three_level_transport_reference_finest.csv"
    --transport-root
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}/three_level_transport_root.csv"
    --transport-middle
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}/three_level_transport_middle.csv"
    --transport-finest
      "${REACTIVE_EB_AMR_THREE_LEVEL_TRANSPORT_2D_WORK_DIRECTORY}/three_level_transport_finest.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_transport_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_three_level_transport_2d_reference;regression_reactive_eb_amr_three_level_transport_2d_run"
)

