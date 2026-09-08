set(REACTIVE_EB_AMR_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_physical_boundary_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_physical_boundary_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_physical_boundary_2d/uniform.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_physical_boundary_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_physical_boundary_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_physical_boundary_2d.py"
    --coarse
      "${REACTIVE_EB_AMR_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}/reactive_eb_amr_physical_boundary_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}/reactive_eb_amr_physical_boundary_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_physical_boundary_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_amr_physical_boundary_2d_run
)

set(REACTIVE_EB_AMR_DYNAMIC_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_dynamic_physical_boundary_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_DYNAMIC_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_dynamic_physical_boundary_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_dynamic_physical_boundary_2d/double_hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_dynamic_physical_boundary_2d_run
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_DYNAMIC_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_dynamic_physical_boundary_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_dynamic_physical_boundary_2d.py"
    --coarse
      "${REACTIVE_EB_AMR_DYNAMIC_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}/reactive_eb_amr_dynamic_boundary_coarse_2d.csv"
    --fine-one
      "${REACTIVE_EB_AMR_DYNAMIC_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}/reactive_eb_amr_dynamic_boundary_fine_2d_patch0001.csv"
    --fine-two
      "${REACTIVE_EB_AMR_DYNAMIC_PHYSICAL_BOUNDARY_2D_WORK_DIRECTORY}/reactive_eb_amr_dynamic_boundary_fine_2d_patch0002.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_dynamic_physical_boundary_2d_check
  PROPERTIES
    DEPENDS regression_reactive_eb_amr_dynamic_physical_boundary_2d_run
)

set(REACTIVE_EB_AMR_DYNAMIC_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_dynamic_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_DYNAMIC_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_dynamic_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_2d/dynamic_hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_dynamic_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_DYNAMIC_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_dynamic_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_dynamic_2d.py"
    --coarse
      "${REACTIVE_EB_AMR_DYNAMIC_2D_WORK_DIRECTORY}/reactive_eb_amr_dynamic_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_DYNAMIC_2D_WORK_DIRECTORY}/reactive_eb_amr_dynamic_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_dynamic_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_amr_dynamic_2d_run
)

set(REACTIVE_EB_AMR_MULTIPATCH_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_multipatch_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_MULTIPATCH_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_multipatch_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_multipatch_2d/double_hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_MULTIPATCH_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_multipatch_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_multipatch_2d.py"
    --coarse
      "${REACTIVE_EB_AMR_MULTIPATCH_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_coarse_2d.csv"
    --fine-one
      "${REACTIVE_EB_AMR_MULTIPATCH_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_fine_2d_patch0001.csv"
    --fine-two
      "${REACTIVE_EB_AMR_MULTIPATCH_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_fine_2d_patch0002.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_amr_multipatch_2d_run
)

set(REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_multipatch_transport_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_multipatch_transport_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_multipatch_transport_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_multipatch_transport_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_multipatch_transport_2d/transport.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_transport_2d_reference
  regression_reactive_eb_amr_multipatch_transport_2d_run
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_multipatch_transport_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_multipatch_transport_2d.py"
    --reference-root
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}/multipatch_transport_reference_root.csv"
    --reference-fine
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}/multipatch_transport_reference_fine_patch0001.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}/multipatch_transport_reference_fine_patch0002.csv"
    --transport-root
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}/multipatch_transport_root.csv"
    --transport-fine
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}/multipatch_transport_fine_patch0001.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_TRANSPORT_2D_WORK_DIRECTORY}/multipatch_transport_fine_patch0002.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_transport_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_multipatch_transport_2d_reference;regression_reactive_eb_amr_multipatch_transport_2d_run"
)

set(REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_multipatch_restart_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_multipatch_restart_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_multipatch_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_multipatch_restart_2d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_multipatch_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_amr_multipatch_restart_2d_restart
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_multipatch_restart_2d/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_restart_2d_reference
  regression_reactive_eb_amr_multipatch_restart_2d_checkpoint_stop
  regression_reactive_eb_amr_multipatch_restart_2d_restart
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_restart_2d_restart
  PROPERTIES DEPENDS
    regression_reactive_eb_amr_multipatch_restart_2d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_amr_multipatch_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_multipatch_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart.chk"
    --expected-checkpoint-sha256
      3ce931f5dcd017de4864789f53b57a47e40fe82cc5fa6c87e9292ff7efe1022c
    --reference
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart_reference_coarse_2d.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart_reference_fine_2d_patch0001.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart_reference_fine_2d_patch0002.csv"
    --stopped
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart_stopped_coarse_2d.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart_stopped_fine_2d_patch0001.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restart_stopped_fine_2d_patch0002.csv"
    --restarted
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restarted_coarse_2d.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restarted_fine_2d_patch0001.csv"
      "${REACTIVE_EB_AMR_MULTIPATCH_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_multipatch_restarted_fine_2d_patch0002.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_multipatch_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_multipatch_restart_2d_reference;regression_reactive_eb_amr_multipatch_restart_2d_checkpoint_stop;regression_reactive_eb_amr_multipatch_restart_2d_restart"
)

set(REACTIVE_EB_AMR_DEGRID_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_degrid_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_DEGRID_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_degrid_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_2d/uniform_degrid.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_degrid_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_DEGRID_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_degrid_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_degrid_2d.py"
    --coarse
      "${REACTIVE_EB_AMR_DEGRID_2D_WORK_DIRECTORY}/reactive_eb_amr_degrid_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_DEGRID_2D_WORK_DIRECTORY}/reactive_eb_amr_degrid_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_degrid_2d_check
  PROPERTIES DEPENDS regression_reactive_eb_amr_degrid_2d_run
)

set(REACTIVE_EB_AMR_CHEMISTRY_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_chemistry_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_CHEMISTRY_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_chemistry_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_chemistry_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_chemistry_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_chemistry_2d/amr.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_chemistry_2d_reference
  regression_reactive_eb_amr_chemistry_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_CHEMISTRY_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_chemistry_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_chemistry_2d.py"
    --reference
      "${REACTIVE_EB_AMR_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_amr_chemistry_reference_2d.csv"
    --coarse
      "${REACTIVE_EB_AMR_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_amr_chemistry_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_amr_chemistry_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_chemistry_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_chemistry_2d_reference;regression_reactive_eb_amr_chemistry_2d_run"
)

set(REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_chemistry_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_2d/amr.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_2d_reference
  regression_reactive_eb_amr_three_level_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_three_level_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_2d.py"
    --reference
      "${REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY}/reactive_eb_amr_chemistry_reference_2d.csv"
    --root
      "${REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY}/reactive_eb_amr_three_level_root_2d.csv"
    --middle
      "${REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY}/reactive_eb_amr_three_level_middle_2d.csv"
    --finest
      "${REACTIVE_EB_AMR_THREE_LEVEL_2D_WORK_DIRECTORY}/reactive_eb_amr_three_level_finest_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_three_level_2d_reference;regression_reactive_eb_amr_three_level_2d_run"
)

set(REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_restart_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_restart_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_restart_2d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_restart_2d_restart
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_restart_2d/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_restart_2d_reference
  regression_reactive_eb_amr_three_level_restart_2d_checkpoint_stop
  regression_reactive_eb_amr_three_level_restart_2d_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_restart_2d_restart
  PROPERTIES DEPENDS
    regression_reactive_eb_amr_three_level_restart_2d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_amr_three_level_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart.chk"
    --reference
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_reference_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_reference_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_reference_finest.csv"
    --stopped
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_stopped_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_stopped_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_stopped_finest.csv"
    --restarted
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_restarted_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_restarted_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_RESTART_2D_WORK_DIRECTORY}/three_level_restart_restarted_finest.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_three_level_restart_2d_reference;regression_reactive_eb_amr_three_level_restart_2d_checkpoint_stop;regression_reactive_eb_amr_three_level_restart_2d_restart"
)

