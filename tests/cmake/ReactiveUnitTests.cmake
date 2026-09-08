add_test(
  NAME regression_sod_pelec_run
  COMMAND $<TARGET_FILE:pelef> "${PROJECT_SOURCE_DIR}/cases/sod/sod_pelec.nml"
)
set_tests_properties(
  regression_sod_pelec_run
  PROPERTIES
    WORKING_DIRECTORY "${SOD_PELEC_WORK_DIRECTORY}"
)

add_test(
  NAME regression_sod_pelec_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_sod.py"
    --input "${SOD_PELEC_WORK_DIRECTORY}/sod_pelec.csv"
    --nx 400
    --time 0.2
    --gamma 1.4
    --density-l1-max 2.0e-3
    --pressure-l1-max 1.5e-3
    --mass-error-max 2.0e-12
    --energy-error-max 2.0e-12
    --momentum-error-max 2.0e-12
)
set_tests_properties(
  regression_sod_pelec_compare
  PROPERTIES
    DEPENDS regression_sod_pelec_run
)

set(SOD_PELEC_PLM_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sod_pelec_plm")
file(MAKE_DIRECTORY "${SOD_PELEC_PLM_WORK_DIRECTORY}")

add_test(
  NAME regression_sod_pelec_plm_run
  COMMAND
    $<TARGET_FILE:pelef>
    "${PROJECT_SOURCE_DIR}/cases/sod/sod_pelec_plm.nml"
)
set_tests_properties(
  regression_sod_pelec_plm_run
  PROPERTIES
    WORKING_DIRECTORY "${SOD_PELEC_PLM_WORK_DIRECTORY}"
)

add_test(
  NAME regression_sod_pelec_plm_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_sod.py"
    --input "${SOD_PELEC_PLM_WORK_DIRECTORY}/sod_pelec_plm.csv"
    --nx 400
    --time 0.2
    --gamma 1.4
    --density-l1-max 1.6e-3
    --pressure-l1-max 1.0e-3
    --mass-error-max 2.0e-12
    --energy-error-max 2.0e-12
    --momentum-error-max 2.0e-12
)
set_tests_properties(
  regression_sod_pelec_plm_compare
  PROPERTIES
    DEPENDS regression_sod_pelec_plm_run
)

set(SHU_OSHER_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/shu_osher")
file(MAKE_DIRECTORY "${SHU_OSHER_WORK_DIRECTORY}")

add_test(
  NAME regression_shu_osher_run
  COMMAND
    $<TARGET_FILE:pelef>
    "${PROJECT_SOURCE_DIR}/cases/shu_osher/shu_osher.nml"
)
set_tests_properties(
  regression_shu_osher_run
  PROPERTIES
    WORKING_DIRECTORY "${SHU_OSHER_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_shu_osher_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_shu_osher.py"
    --input "${SHU_OSHER_WORK_DIRECTORY}/shu_osher.csv"
)
set_tests_properties(
  regression_shu_osher_compare
  PROPERTIES
    DEPENDS regression_shu_osher_run
)

set(
  SHU_OSHER_PELEC_PLM_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/shu_osher_pelec_plm"
)
file(MAKE_DIRECTORY "${SHU_OSHER_PELEC_PLM_WORK_DIRECTORY}")

add_test(
  NAME regression_shu_osher_pelec_plm_run
  COMMAND
    $<TARGET_FILE:pelef>
    "${PROJECT_SOURCE_DIR}/cases/shu_osher/shu_osher_pelec_plm.nml"
)
set_tests_properties(
  regression_shu_osher_pelec_plm_run
  PROPERTIES
    WORKING_DIRECTORY "${SHU_OSHER_PELEC_PLM_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_shu_osher_pelec_plm_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_shu_osher.py"
    --input
    "${SHU_OSHER_PELEC_PLM_WORK_DIRECTORY}/shu_osher_pelec_plm.csv"
    --density-squared-reference 112.78740201441508
    --density-moment-reference -7.071310814334509
    --signature-relative-tolerance 5.0e-6
)
set_tests_properties(
  regression_shu_osher_pelec_plm_compare
  PROPERTIES
    DEPENDS regression_shu_osher_pelec_plm_run
)

set(SEDOV_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sedov")
file(MAKE_DIRECTORY "${SEDOV_WORK_DIRECTORY}")

add_test(
  NAME regression_sedov_run
  COMMAND
    $<TARGET_FILE:pelef>
    "${PROJECT_SOURCE_DIR}/cases/sedov/sedov.nml"
)
set_tests_properties(
  regression_sedov_run
  PROPERTIES
    WORKING_DIRECTORY "${SEDOV_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_sedov_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_sedov.py"
    --input "${SEDOV_WORK_DIRECTORY}/sedov.csv"
    --density-squared-reference 1.4513811919127926
    --pressure-integral-reference 2.0606157390942492
    --shock-radius-reference 0.13187500000000002
    --signature-relative-tolerance 5.0e-6
)
set_tests_properties(
  regression_sedov_compare
  PROPERTIES
    DEPENDS regression_sedov_run
)

add_executable(
  test_directional_flux_2d
  unit/test_directional_flux_2d.F90
)
target_link_libraries(test_directional_flux_2d PRIVATE pelef_core)
add_test(NAME unit_directional_flux_2d COMMAND test_directional_flux_2d)

add_executable(test_mesh_3d unit/test_mesh_3d.F90)
target_link_libraries(test_mesh_3d PRIVATE pelef_core)
add_test(NAME unit_mesh_3d COMMAND test_mesh_3d)

add_executable(test_eb_geometry_3d unit/test_eb_geometry_3d.F90)
target_link_libraries(test_eb_geometry_3d PRIVATE pelef_core)
add_test(NAME unit_eb_geometry_3d COMMAND test_eb_geometry_3d)

add_executable(
  test_reactive_eb_cfl_3d
  unit/test_reactive_eb_cfl_3d.F90
)
target_link_libraries(test_reactive_eb_cfl_3d PRIVATE pelef_core)
add_test(NAME unit_reactive_eb_cfl_3d COMMAND test_reactive_eb_cfl_3d)

add_executable(
  test_reactive_eb_3d_config
  unit/test_reactive_eb_3d_config.F90
)
target_link_libraries(test_reactive_eb_3d_config PRIVATE pelef_core)
add_test(NAME unit_reactive_eb_3d_config COMMAND test_reactive_eb_3d_config)

add_executable(
  test_reactive_eb_3d_checkpoint
  unit/test_reactive_eb_3d_checkpoint.F90
)
target_link_libraries(test_reactive_eb_3d_checkpoint PRIVATE pelef_core)
add_test(
  NAME unit_reactive_eb_3d_checkpoint
  COMMAND test_reactive_eb_3d_checkpoint
)

add_executable(
  test_eb_reactive_wall_flux_3d
  unit/test_eb_reactive_wall_flux_3d.F90
)
target_link_libraries(test_eb_reactive_wall_flux_3d PRIVATE pelef_core)
add_test(
  NAME unit_eb_reactive_wall_flux_3d
  COMMAND test_eb_reactive_wall_flux_3d
)

add_executable(
  test_eb_reactive_redistribution_3d
  unit/test_eb_reactive_redistribution_3d.F90
)
target_link_libraries(test_eb_reactive_redistribution_3d PRIVATE pelef_core)
add_test(
  NAME unit_eb_reactive_redistribution_3d
  COMMAND test_eb_reactive_redistribution_3d
)

add_executable(
  test_reactive_eb_hydro_3d
  regression/test_reactive_eb_hydro_3d.F90
)
target_link_libraries(test_reactive_eb_hydro_3d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_eb_hydro_3d
  COMMAND test_reactive_eb_hydro_3d
)

add_executable(
  test_reactive_eb_chemistry_3d
  unit/test_reactive_eb_chemistry_3d.F90
)
target_link_libraries(test_reactive_eb_chemistry_3d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_eb_chemistry_3d
  COMMAND test_reactive_eb_chemistry_3d
)

add_executable(
  test_reactive_eb_transport_3d
  unit/test_reactive_eb_transport_3d.F90
)
target_link_libraries(test_reactive_eb_transport_3d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_eb_transport_3d
  COMMAND test_reactive_eb_transport_3d
)

set(REACTIVE_EB_3D_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_3d")
file(MAKE_DIRECTORY "${REACTIVE_EB_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_3d_state_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/density_sheet_state.nml"
)
set_tests_properties(
  regression_reactive_eb_3d_state_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Redistribution: state_redist"
)
add_test(
  NAME regression_reactive_eb_3d_flux_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/density_sheet_flux.nml"
)
set_tests_properties(
  regression_reactive_eb_3d_flux_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Redistribution: flux_redist"
)
add_test(
  NAME regression_reactive_eb_3d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_3d.py"
    --state "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_state.csv"
    --flux "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_flux.csv"
)
set_tests_properties(
  regression_reactive_eb_3d_check
  PROPERTIES
    DEPENDS
      "regression_reactive_eb_3d_state_run;regression_reactive_eb_3d_flux_run"
)

add_test(
  NAME regression_reactive_eb_chemistry_3d_inert_run
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_chemistry_inert.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/chemistry_inert.nml"
)
set_tests_properties(
  regression_reactive_eb_chemistry_3d_inert_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  F"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_chemistry_3d_reactive_run
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_chemistry_reactive.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/chemistry_reactive.nml"
)
set_tests_properties(
  regression_reactive_eb_chemistry_3d_reactive_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_chemistry_3d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_chemistry_3d.py"
    --baseline "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_state.csv"
    --inert "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_chemistry_inert.csv"
    --reactive "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_chemistry_reactive.csv"
    --inert-log
      "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_chemistry_inert.log"
    --reactive-log
      "${REACTIVE_EB_3D_WORK_DIRECTORY}/reactive_eb_3d_chemistry_reactive.log"
)
set(REACTIVE_EB_CHEMISTRY_3D_CHECK_DEPENDS
  regression_reactive_eb_3d_state_run
  regression_reactive_eb_chemistry_3d_inert_run
  regression_reactive_eb_chemistry_3d_reactive_run
)
set_tests_properties(
  regression_reactive_eb_chemistry_3d_check
  PROPERTIES
    DEPENDS
      "${REACTIVE_EB_CHEMISTRY_3D_CHECK_DEPENDS}"
    TIMEOUT 60
)

set(REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_transport_3d")
file(MAKE_DIRECTORY "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_transport_3d_control_run
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/control.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/transport_control.nml"
)
add_test(
  NAME regression_reactive_eb_transport_3d_run
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/transport.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/transport_only.nml"
)
add_test(
  NAME regression_reactive_eb_transport_3d_coupled_run
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/coupled.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d/transport_chemistry.nml"
)
set_tests_properties(
  regression_reactive_eb_transport_3d_control_run
  regression_reactive_eb_transport_3d_run
  regression_reactive_eb_transport_3d_coupled_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_transport_3d_control_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Operator sequence: H"
)
set_tests_properties(
  regression_reactive_eb_transport_3d_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Operator sequence: T-H-T"
)
set_tests_properties(
  regression_reactive_eb_transport_3d_coupled_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Operator sequence: R-T-H-T-R"
)
add_test(
  NAME regression_reactive_eb_transport_3d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_transport_3d.py"
    --control
      "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/reactive_eb_3d_transport_control.csv"
    --transport
      "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/reactive_eb_3d_transport_only.csv"
    --coupled
      "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/reactive_eb_3d_transport_chemistry.csv"
    --control-log "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/control.log"
    --transport-log "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/transport.log"
    --coupled-log "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/coupled.log"
)
set_tests_properties(
  regression_reactive_eb_transport_3d_check
  PROPERTIES
    DEPENDS
      "regression_reactive_eb_transport_3d_control_run;regression_reactive_eb_transport_3d_run;regression_reactive_eb_transport_3d_coupled_run"
    TIMEOUT 60
)

set(REACTIVE_EB_RESTART_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_restart_3d")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(REACTIVE_EB_RESTART_3D_CHECKPOINT_SHA256
    "2acd59da30152435b832b02782ac7edfa21c6bd26e11e793279209f662e2bcbb")
  set(REACTIVE_EB_RESTART_3D_STOPPED_SHA256
    "11606d8180fd593f65a417ccdc7e5eb5419dc0b751c8ecc2ce00341879e2111c")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(REACTIVE_EB_RESTART_3D_CHECKPOINT_SHA256
    "91f8bacb835413af5c291cf75a839f93e923d5e124ac4a4f9a8c90114c57b835")
  set(REACTIVE_EB_RESTART_3D_STOPPED_SHA256
    "ee4d8f656c37b454a81e20e8f5fedc86119d138950fb372fc90003378a973e35")
else()
  message(FATAL_ERROR
    "3D EB restart parity requires CMAKE_BUILD_TYPE=Debug or Release")
endif()

set(
  AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_reactive_3d_transport_restart_0245"
)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256
    "6a1b1ba7bece94e4f59e0da95fb6dc8ddc5bc5d356f3d88f56d37e184e0b573d")
  set(AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256
    "0b0955e45d7d5e6296bdd69270cd78ed02eba4f814a3789d7218d861fd0eac6f")
  set(AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256
    "45d51c690e8f3b1d20bb224deae24cbbb5ec4dc0c648a9b1367d6963acaeb37b")
  set(SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256
    "6fa40282c6a12598bf7147efe2253540e49a66abb5e907aac8b5226443fdf342")
  set(SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256
    "34d781bd9faef30a0fe5f859300fac459380403d67b00795f1c3da95e294cb39")
  set(SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256
    "e4c3f6fe27ed9c8e62d33912c42adc6bdef4ee1c3e8cd0b6aebfdf34c42956ed")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256
    "15c2e5a3a087a38fe34ece1fa575bbc7857ec7c651ae9ee84f3fd8d52615afc5")
  set(AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256
    "c81eb517ab4111f66eafa09b70d1f9543251e206d62143cfe02e50ca74e4f441")
  set(AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256
    "45d51c690e8f3b1d20bb224deae24cbbb5ec4dc0c648a9b1367d6963acaeb37b")
  set(SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256
    "52526217eb88e869167b8ef8b958f3cf6176b5346746a5d9ff1d1757ca526df3")
  set(SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256
    "21470b7ac704274637266990f8bd77bbeacdba5c2808c2943fdc17fb106ee35a")
  set(SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256
    "e4c3f6fe27ed9c8e62d33912c42adc6bdef4ee1c3e8cd0b6aebfdf34c42956ed")
else()
  message(FATAL_ERROR
    "0.245 transport restart parity requires CMAKE_BUILD_TYPE=Debug or Release")
endif()
file(MAKE_DIRECTORY "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_reactive_3d_transport_restart_0245_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_reference.nml"
)
add_test(
  NAME regression_amr_reactive_3d_transport_restart_0245_checkpoint
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/checkpoint.log"
    -- $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
)
add_test(
  NAME regression_amr_reactive_3d_transport_restart_0245_run
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart.log"
    -- $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
)
set_tests_properties(
  regression_amr_reactive_3d_transport_restart_0245_reference
  regression_amr_reactive_3d_transport_restart_0245_checkpoint
  regression_amr_reactive_3d_transport_restart_0245_run
  PROPERTIES
    WORKING_DIRECTORY
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_amr_reactive_3d_transport_restart_0245_checkpoint
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
set_tests_properties(
  regression_amr_reactive_3d_transport_restart_0245_run
  PROPERTIES
    DEPENDS regression_amr_reactive_3d_transport_restart_0245_checkpoint
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
)
add_test(
  NAME regression_amr_reactive_3d_transport_restart_0245_exact
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
    --reference-prefix
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_reference_0245"
    --candidate-prefix
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_restarted_0245"
    --checkpoint
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
    --expected-checkpoint-sha256
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256}"
    --expected-coarse-sha256
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256}"
    --expected-fine-sha256
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256}"
    --expected-schema 4
    --final-time 1.0e-9
    --reconstruction characteristic_plm
    --limiter mc
    --reference-log
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/reference.log"
    --candidate-log
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart.log"
)
set_tests_properties(
  regression_amr_reactive_3d_transport_restart_0245_exact
  PROPERTIES DEPENDS
    "regression_amr_reactive_3d_transport_restart_0245_reference;regression_amr_reactive_3d_transport_restart_0245_checkpoint;regression_amr_reactive_3d_transport_restart_0245_run"
)

set(
  SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_3d_transport_restart_0245"
)
file(MAKE_DIRECTORY
  "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}")
add_test(
  NAME regression_selected_amr_reactive_3d_transport_restart_0245_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_reference.nml"
)
add_test(
  NAME regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/checkpoint.log"
    -- $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
)
add_test(
  NAME regression_selected_amr_reactive_3d_transport_restart_0245_run
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart.log"
    -- $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_restart_0245_reference
  regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint
  regression_selected_amr_reactive_3d_transport_restart_0245_run
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
    TIMEOUT 600
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_restart_0245_run
  PROPERTIES
    DEPENDS regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
)
add_test(
  NAME regression_selected_amr_reactive_3d_transport_restart_0245_exact
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
    --reference-prefix
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_reference_0245"
    --candidate-prefix
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_restarted_0245"
    --checkpoint
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
    --expected-checkpoint-sha256
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_CHECKPOINT_SHA256}"
    --expected-coarse-sha256
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_COARSE_SHA256}"
    --expected-fine-sha256
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_FINE_SHA256}"
    --expected-schema 5
    --expected-bundle-sha256
      f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
    --expected-integrator implicit
    --final-time 1.0e-9
    --reconstruction characteristic_plm
    --limiter mc
    --reference-log
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/reference.log"
    --candidate-log
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart.log"
)
set_tests_properties(
  regression_selected_amr_reactive_3d_transport_restart_0245_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_3d_transport_restart_0245_reference;regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint;regression_selected_amr_reactive_3d_transport_restart_0245_run"
)

