add_executable(
  test_directional_flux_3d
  unit/test_directional_flux_3d.F90
)
target_link_libraries(test_directional_flux_3d PRIVATE pelef_core)
add_test(NAME unit_directional_flux_3d COMMAND test_directional_flux_3d)

add_executable(
  test_euler_3d_dimensional_reduction
  unit/test_euler_3d_dimensional_reduction.F90
)
target_link_libraries(test_euler_3d_dimensional_reduction PRIVATE pelef_core)
add_test(
  NAME unit_euler_3d_dimensional_reduction
  COMMAND test_euler_3d_dimensional_reduction
)

add_executable(
  test_reactive_euler_3d_dimensional_reduction
  unit/test_reactive_euler_3d_dimensional_reduction.F90
)
target_link_libraries(
  test_reactive_euler_3d_dimensional_reduction PRIVATE pelef_core
)
add_test(
  NAME unit_reactive_euler_3d_dimensional_reduction
  COMMAND test_reactive_euler_3d_dimensional_reduction
)

add_executable(
  test_reactive_chemistry_3d
  unit/test_reactive_chemistry_3d.F90
)
target_link_libraries(test_reactive_chemistry_3d PRIVATE pelef_core)
add_test(NAME unit_reactive_chemistry_3d COMMAND test_reactive_chemistry_3d)

add_executable(
  test_reactive_transport_3d
  unit/test_reactive_transport_3d.F90
)
target_link_libraries(test_reactive_transport_3d PRIVATE pelef_core)
add_test(NAME unit_reactive_transport_3d COMMAND test_reactive_transport_3d)
set_tests_properties(unit_reactive_transport_3d PROPERTIES TIMEOUT 300)

add_executable(
  test_amr_reactive_transport_3d
  unit/test_amr_reactive_transport_3d.F90
)
target_link_libraries(test_amr_reactive_transport_3d PRIVATE pelef_core)
add_test(
  NAME unit_amr_reactive_transport_3d
  COMMAND test_amr_reactive_transport_3d
)
set_tests_properties(unit_amr_reactive_transport_3d PROPERTIES TIMEOUT 300)

add_executable(
  test_amr_hierarchy_3d
  unit/test_amr_hierarchy_3d.F90
)
target_link_libraries(test_amr_hierarchy_3d PRIVATE pelef_core)
add_test(NAME unit_amr_hierarchy_3d COMMAND test_amr_hierarchy_3d)

add_executable(
  test_amr_reactive_3d_config
  unit/test_amr_reactive_3d_config.F90
)
target_link_libraries(test_amr_reactive_3d_config PRIVATE pelef_core)
add_test(NAME unit_amr_reactive_3d_config COMMAND test_amr_reactive_3d_config)

add_executable(
  test_amr_reactive_chemistry_3d
  unit/test_amr_reactive_chemistry_3d.F90
)
target_link_libraries(test_amr_reactive_chemistry_3d PRIVATE pelef_core)
add_test(
  NAME unit_amr_reactive_chemistry_3d
  COMMAND test_amr_reactive_chemistry_3d
)
set_tests_properties(unit_amr_reactive_chemistry_3d PROPERTIES TIMEOUT 300)

add_executable(
  test_amr_reactive_3d_checkpoint
  unit/test_amr_reactive_3d_checkpoint.F90
)
target_link_libraries(test_amr_reactive_3d_checkpoint PRIVATE pelef_core)
add_test(
  NAME unit_amr_reactive_3d_checkpoint
  COMMAND test_amr_reactive_3d_checkpoint
)

add_executable(
  test_amr_reactive_3d
  regression/test_amr_reactive_3d.F90
)
target_link_libraries(test_amr_reactive_3d PRIVATE pelef_core)
add_test(NAME regression_amr_reactive_3d COMMAND test_amr_reactive_3d)
set_tests_properties(regression_amr_reactive_3d PROPERTIES TIMEOUT 300)

set(AMR_REACTIVE_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_reactive_3d")
file(MAKE_DIRECTORY "${AMR_REACTIVE_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_reactive_3d_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave.nml"
)
set_tests_properties(
  regression_amr_reactive_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "static two-level 3D reactive AMR"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
    --coarse
      "${AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d_coarse.csv"
    --fine
      "${AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d_fine.csv"
    --nx 8
    --ny 8
    --nz 8
    --i-lower 3
    --i-upper 6
    --j-lower 3
    --j-upper 6
    --k-lower 3
    --k-upper 6
    --ratio 2
    --time 1.0e-6
)
set_tests_properties(
  regression_amr_reactive_3d_check
  PROPERTIES DEPENDS regression_amr_reactive_3d_run
)

set(AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_reactive_chemistry_3d")
file(MAKE_DIRECTORY "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_reactive_chemistry_3d_control_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_control.nml"
)
set_tests_properties(
  regression_amr_reactive_chemistry_3d_control_run
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  F"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_chemistry_3d_chemistry_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_chemistry.nml"
)
set_tests_properties(
  regression_amr_reactive_chemistry_3d_chemistry_run
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  T"
    TIMEOUT 300
)
foreach(amr_reactive_chemistry_3d_mode IN ITEMS control chemistry)
  add_test(
    NAME regression_amr_reactive_chemistry_3d_${amr_reactive_chemistry_3d_mode}_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
      --coarse
        "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/amr_hotspot_${amr_reactive_chemistry_3d_mode}_3d_coarse.csv"
      --fine
        "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/amr_hotspot_${amr_reactive_chemistry_3d_mode}_3d_fine.csv"
      --nx 8 --ny 8 --nz 8
      --i-lower 3 --i-upper 6
      --j-lower 3 --j-upper 6
      --k-lower 3 --k-upper 6
      --ratio 2 --time 1.0e-9
  )
  set_tests_properties(
    regression_amr_reactive_chemistry_3d_${amr_reactive_chemistry_3d_mode}_check
    PROPERTIES DEPENDS
      regression_amr_reactive_chemistry_3d_${amr_reactive_chemistry_3d_mode}_run
  )
endforeach()
add_test(
  NAME regression_amr_reactive_chemistry_3d_activity
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_chemistry_3d.py"
    --control-prefix
      "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/amr_hotspot_control_3d"
    --reactive-prefix
      "${AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/amr_hotspot_chemistry_3d"
    --nx 8 --ny 8 --nz 8
    --i-lower 3 --i-upper 6
    --j-lower 3 --j-upper 6
    --k-lower 3 --k-upper 6
    --ratio 2
)
set_tests_properties(
  regression_amr_reactive_chemistry_3d_activity
  PROPERTIES DEPENDS
    "regression_amr_reactive_chemistry_3d_control_run;regression_amr_reactive_chemistry_3d_chemistry_run"
)

add_test(
  NAME regression_amr_reactive_3d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml"
)
set_tests_properties(
  regression_amr_reactive_3d_checkpoint_stop
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_3d_restart
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_restart.nml"
)
set_tests_properties(
  regression_amr_reactive_3d_restart
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    DEPENDS regression_amr_reactive_3d_checkpoint_stop
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_3d_restart_compare
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_amr_reactive_3d_restart.py"
    --reference-coarse
      "${AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d_coarse.csv"
    --reference-fine
      "${AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d_fine.csv"
    --restart-coarse
      "${AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d_restart_coarse.csv"
    --restart-fine
      "${AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d_restart_fine.csv"
)
set_tests_properties(
  regression_amr_reactive_3d_restart_compare
  PROPERTIES DEPENDS
    "regression_amr_reactive_3d_run;regression_amr_reactive_3d_restart"
)

set(AMR_REACTIVE_PLM_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_reactive_plm_3d")
file(MAKE_DIRECTORY "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_reactive_plm_3d_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm.nml"
)
set_tests_properties(
  regression_amr_reactive_plm_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_plm_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
    --coarse
      "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d_coarse.csv"
    --fine
      "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d_fine.csv"
    --nx 8
    --ny 8
    --nz 8
    --i-lower 3
    --i-upper 6
    --j-lower 3
    --j-upper 6
    --k-lower 3
    --k-upper 6
    --ratio 2
    --time 1.0e-6
)
set_tests_properties(
  regression_amr_reactive_plm_3d_check
  PROPERTIES DEPENDS regression_amr_reactive_plm_3d_run
)
add_test(
  NAME regression_amr_reactive_plm_3d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm_checkpoint.nml"
)
set_tests_properties(
  regression_amr_reactive_plm_3d_checkpoint_stop
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_plm_3d_restart
  COMMAND $<TARGET_FILE:pelef_amr_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm_restart.nml"
)
set_tests_properties(
  regression_amr_reactive_plm_3d_restart
  PROPERTIES
    WORKING_DIRECTORY "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    DEPENDS regression_amr_reactive_plm_3d_checkpoint_stop
    TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_plm_3d_restart_compare
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_amr_reactive_3d_restart.py"
    --reference-coarse
      "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d_coarse.csv"
    --reference-fine
      "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d_fine.csv"
    --restart-coarse
      "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d_restart_coarse.csv"
    --restart-fine
      "${AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d_restart_fine.csv"
)
set_tests_properties(
  regression_amr_reactive_plm_3d_restart_compare
  PROPERTIES DEPENDS
    "regression_amr_reactive_plm_3d_run;regression_amr_reactive_plm_3d_restart"
)

add_executable(
  test_reactive_transport_3d_regression
  regression/test_reactive_transport_3d.F90
)
target_link_libraries(test_reactive_transport_3d_regression PRIVATE pelef_core)
add_test(
  NAME regression_reactive_transport_3d
  COMMAND test_reactive_transport_3d_regression
)
set_tests_properties(regression_reactive_transport_3d PROPERTIES TIMEOUT 300)

add_executable(
  test_entropy_wave_3d
  regression/test_entropy_wave_3d.F90
)
target_link_libraries(test_entropy_wave_3d PRIVATE pelef_core)
add_test(NAME regression_entropy_wave_3d COMMAND test_entropy_wave_3d)
set_tests_properties(regression_entropy_wave_3d PROPERTIES TIMEOUT 300)

add_executable(
  test_reactive_entropy_wave_3d
  regression/test_reactive_entropy_wave_3d.F90
)
target_link_libraries(test_reactive_entropy_wave_3d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_entropy_wave_3d
  COMMAND test_reactive_entropy_wave_3d
)
set_tests_properties(
  regression_reactive_entropy_wave_3d PROPERTIES TIMEOUT 300
)

add_executable(
  test_reactive_entropy_wave_plm_3d
  regression/test_reactive_entropy_wave_plm_3d.F90
)
target_link_libraries(test_reactive_entropy_wave_plm_3d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_entropy_wave_plm_3d
  COMMAND test_reactive_entropy_wave_plm_3d
)
set_tests_properties(
  regression_reactive_entropy_wave_plm_3d PROPERTIES TIMEOUT 300
)

set(REACTIVE_ENTROPY_WAVE_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_entropy_wave_3d")
file(MAKE_DIRECTORY "${REACTIVE_ENTROPY_WAVE_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_entropy_wave_3d_run
  COMMAND $<TARGET_FILE:pelef_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_entropy_wave_3d/entropy_wave.nml"
)
set_tests_properties(
  regression_reactive_entropy_wave_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_ENTROPY_WAVE_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION
      "PeleF [0-9]+\\.[0-9]+\\.[0-9]+ 3D reactive flow"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_entropy_wave_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_entropy_wave_3d.py"
    --input
      "${REACTIVE_ENTROPY_WAVE_3D_WORK_DIRECTORY}/reactive_entropy_wave_3d.csv"
    --nx 12
    --ny 12
    --nz 12
    --time 1.0e-6
    --thermo-model full_h2o2
    --maximum-l1 0.007
)
set_tests_properties(
  regression_reactive_entropy_wave_3d_check
  PROPERTIES DEPENDS regression_reactive_entropy_wave_3d_run
)

set(REACTIVE_ENTROPY_WAVE_PLM_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_entropy_wave_plm_3d")
file(MAKE_DIRECTORY "${REACTIVE_ENTROPY_WAVE_PLM_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_entropy_wave_plm_3d_run
  COMMAND $<TARGET_FILE:pelef_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_entropy_wave_3d/entropy_wave_plm.nml"
)
set_tests_properties(
  regression_reactive_entropy_wave_plm_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_ENTROPY_WAVE_PLM_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_entropy_wave_plm_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_entropy_wave_3d.py"
    --input
      "${REACTIVE_ENTROPY_WAVE_PLM_3D_WORK_DIRECTORY}/reactive_entropy_wave_plm_3d.csv"
    --nx 12
    --ny 12
    --nz 12
    --time 1.0e-6
    --thermo-model full_h2o2
    --maximum-l1 0.002
)
set_tests_properties(
  regression_reactive_entropy_wave_plm_3d_check
  PROPERTIES DEPENDS regression_reactive_entropy_wave_plm_3d_run
)

set(REACTIVE_HOTSPOT_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_hotspot_3d")
file(MAKE_DIRECTORY "${REACTIVE_HOTSPOT_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_hotspot_3d_run
  COMMAND $<TARGET_FILE:pelef_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot_3d/hotspot.nml"
)
set_tests_properties(
  regression_reactive_hotspot_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_HOTSPOT_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_hotspot_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot_3d.py"
    --input "${REACTIVE_HOTSPOT_3D_WORK_DIRECTORY}/reactive_hotspot_3d.csv"
    --nx 6
    --ny 6
    --nz 6
    --time 1.0e-6
)
set_tests_properties(
  regression_reactive_hotspot_3d_check
  PROPERTIES DEPENDS regression_reactive_hotspot_3d_run
)

set(REACTIVE_TRANSPORT_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_transport_3d")
file(MAKE_DIRECTORY "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_transport_3d_control_run
  COMMAND $<TARGET_FILE:pelef_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_transport_3d/control.nml"
)
set_tests_properties(
  regression_reactive_transport_3d_control_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Molecular transport:  F"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_transport_3d_run
  COMMAND $<TARGET_FILE:pelef_reactive_3d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_transport_3d/transport.nml"
)
set_tests_properties(
  regression_reactive_transport_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Molecular transport:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_transport_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_transport_3d.py"
    --control
      "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/reactive_transport_3d_control.csv"
    --transport
      "${REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/reactive_transport_3d.csv"
    --nx 6
    --ny 6
    --nz 6
    --time 2.0e-6
)
set_tests_properties(
  regression_reactive_transport_3d_check
  PROPERTIES DEPENDS
    "regression_reactive_transport_3d_control_run;regression_reactive_transport_3d_run"
)

set(ENTROPY_WAVE_3D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/entropy_wave_3d")
file(MAKE_DIRECTORY "${ENTROPY_WAVE_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_entropy_wave_3d_run
  COMMAND $<TARGET_FILE:pelef3d>
    "${PROJECT_SOURCE_DIR}/cases/entropy_wave_3d/entropy_wave.nml"
)
set_tests_properties(
  regression_entropy_wave_3d_run
  PROPERTIES
    WORKING_DIRECTORY "${ENTROPY_WAVE_3D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "PeleF [0-9]+\\.[0-9]+\\.[0-9]+ 3D"
    TIMEOUT 300
)
add_test(
  NAME regression_entropy_wave_3d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_entropy_wave_3d.py"
    --input "${ENTROPY_WAVE_3D_WORK_DIRECTORY}/entropy_wave_3d.csv"
    --nx 16
    --ny 16
    --nz 16
    --time 0.05
    --maximum-l1 0.04
)
set_tests_properties(
  regression_entropy_wave_3d_check
  PROPERTIES DEPENDS regression_entropy_wave_3d_run
)

add_executable(
  test_ctu_transverse_correction
  unit/test_ctu_transverse_correction.F90
)
target_link_libraries(test_ctu_transverse_correction PRIVATE pelef_core)
add_test(
  NAME unit_ctu_transverse_correction
  COMMAND test_ctu_transverse_correction
)

add_executable(
  test_isentropic_vortex_2d
  regression/test_isentropic_vortex_2d.F90
)
target_link_libraries(test_isentropic_vortex_2d PRIVATE pelef_core)
add_test(
  NAME regression_isentropic_vortex_2d_convergence
  COMMAND test_isentropic_vortex_2d
)
set_tests_properties(
  regression_isentropic_vortex_2d_convergence
  PROPERTIES
    TIMEOUT 180
)

set(VORTEX_2D_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/isentropic_vortex_2d")
file(MAKE_DIRECTORY "${VORTEX_2D_WORK_DIRECTORY}")

add_test(
  NAME regression_isentropic_vortex_2d_run
  COMMAND
    $<TARGET_FILE:pelef2d>
    "${PROJECT_SOURCE_DIR}/cases/isentropic_vortex/vortex.nml"
)
set_tests_properties(
  regression_isentropic_vortex_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${VORTEX_2D_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_isentropic_vortex_2d_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_isentropic_vortex.py"
    --input "${VORTEX_2D_WORK_DIRECTORY}/isentropic_vortex.csv"
)
set_tests_properties(
  regression_isentropic_vortex_2d_compare
  PROPERTIES
    DEPENDS regression_isentropic_vortex_2d_run
)

add_executable(
  test_ctu_dimensional_reduction
  unit/test_ctu_dimensional_reduction.F90
)
target_link_libraries(test_ctu_dimensional_reduction PRIVATE pelef_core)
add_test(
  NAME unit_ctu_dimensional_reduction
  COMMAND test_ctu_dimensional_reduction
)

add_executable(test_multispecies_state unit/test_multispecies_state.F90)
target_link_libraries(test_multispecies_state PRIVATE pelef_core)
add_test(NAME unit_multispecies_state COMMAND test_multispecies_state)

add_executable(test_multispecies_flux unit/test_multispecies_flux.F90)
target_link_libraries(test_multispecies_flux PRIVATE pelef_core)
add_test(NAME unit_multispecies_flux COMMAND test_multispecies_flux)

add_executable(
  test_multispecies_reconstruction
  unit/test_multispecies_reconstruction.F90
)
target_link_libraries(test_multispecies_reconstruction PRIVATE pelef_core)
add_test(
  NAME unit_multispecies_reconstruction
  COMMAND test_multispecies_reconstruction
)

add_executable(
  test_multispecies_ctu_dimensional_reduction
  unit/test_multispecies_ctu_dimensional_reduction.F90
)
target_link_libraries(test_multispecies_ctu_dimensional_reduction PRIVATE pelef_core)
add_test(
  NAME unit_multispecies_ctu_dimensional_reduction
  COMMAND test_multispecies_ctu_dimensional_reduction
)

add_executable(
  test_multispecies_entropy_wave
  regression/test_multispecies_entropy_wave.F90
)
target_link_libraries(test_multispecies_entropy_wave PRIVATE pelef_core)
add_test(
  NAME regression_multispecies_entropy_wave
  COMMAND test_multispecies_entropy_wave
)

add_executable(
  test_multispecies_advection_2d
  regression/test_multispecies_advection_2d.F90
)
target_link_libraries(test_multispecies_advection_2d PRIVATE pelef_core)
add_test(
  NAME regression_multispecies_advection_2d
  COMMAND test_multispecies_advection_2d
)
set_tests_properties(
  regression_multispecies_advection_2d
  PROPERTIES TIMEOUT 180
)

set(MULTISPEC_SOD_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/multispec_sod")
file(MAKE_DIRECTORY "${MULTISPEC_SOD_WORK_DIRECTORY}")

add_test(
  NAME regression_multispec_sod_run
  COMMAND
    $<TARGET_FILE:pelef_ms>
    "${PROJECT_SOURCE_DIR}/cases/multispec_sod/multispec_sod.nml"
)
set_tests_properties(
  regression_multispec_sod_run
  PROPERTIES
    WORKING_DIRECTORY "${MULTISPEC_SOD_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_multispec_sod_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_multispec_sod.py"
    --input "${MULTISPEC_SOD_WORK_DIRECTORY}/multispec_sod.csv"
)
set_tests_properties(
  regression_multispec_sod_compare
  PROPERTIES DEPENDS regression_multispec_sod_run
)

add_executable(test_nasa7_thermo unit/test_nasa7_thermo.F90)
target_link_libraries(test_nasa7_thermo PRIVATE pelef_core)
add_test(NAME unit_nasa7_thermo COMMAND test_nasa7_thermo)

add_executable(test_mixture_thermo unit/test_mixture_thermo.F90)
target_link_libraries(test_mixture_thermo PRIVATE pelef_core)
add_test(NAME unit_mixture_thermo COMMAND test_mixture_thermo)

add_executable(
  test_isomerization_reactor
  unit/test_isomerization_reactor.F90
)
target_link_libraries(test_isomerization_reactor PRIVATE pelef_core)
add_test(NAME unit_isomerization_reactor COMMAND test_isomerization_reactor)

set(ZERO_D_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/zero_d_isomerization")
file(MAKE_DIRECTORY "${ZERO_D_WORK_DIRECTORY}")

add_test(
  NAME regression_zero_d_isomerization_run
  COMMAND
    $<TARGET_FILE:pelef0d>
    "${PROJECT_SOURCE_DIR}/cases/zero_d_isomerization/reactor.nml"
)
set_tests_properties(
  regression_zero_d_isomerization_run
  PROPERTIES
    WORKING_DIRECTORY "${ZERO_D_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_zero_d_isomerization_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_zero_d_isomerization.py"
    --input "${ZERO_D_WORK_DIRECTORY}/zero_d_isomerization.csv"
)
set_tests_properties(
  regression_zero_d_isomerization_compare
  PROPERTIES DEPENDS regression_zero_d_isomerization_run
)

add_executable(
  test_elementary_kinetics
  unit/test_elementary_kinetics.F90
)
target_link_libraries(test_elementary_kinetics PRIVATE pelef_core)
add_test(NAME unit_elementary_kinetics COMMAND test_elementary_kinetics)

add_executable(
  test_constant_volume_h2o2
  unit/test_constant_volume_h2o2.F90
)
target_link_libraries(test_constant_volume_h2o2 PRIVATE pelef_core)
add_test(NAME unit_constant_volume_h2o2 COMMAND test_constant_volume_h2o2)
set_tests_properties(unit_constant_volume_h2o2 PROPERTIES TIMEOUT 180)

add_test(
  NAME generated_h2o2_mechanism_clean
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_generated_mechanism.py"
    --generator "${PROJECT_SOURCE_DIR}/tools/generate_elementary_mechanism.py"
    --input "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_elementary.json"
    --committed
    "${PROJECT_SOURCE_DIR}/src/generated/h2o2_elementary_mechanism_mod.F90"
)

set(H2O2_ZERO_D_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/zero_d_h2o2")
file(MAKE_DIRECTORY "${H2O2_ZERO_D_WORK_DIRECTORY}")

add_test(
  NAME regression_zero_d_h2o2_run
  COMMAND
    $<TARGET_FILE:pelef0d_h2o2>
    "${PROJECT_SOURCE_DIR}/cases/zero_d_h2o2/reactor.nml"
)
set_tests_properties(
  regression_zero_d_h2o2_run
  PROPERTIES
    WORKING_DIRECTORY "${H2O2_ZERO_D_WORK_DIRECTORY}"
    TIMEOUT 180
)

add_test(
  NAME regression_zero_d_h2o2_structure
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_zero_d_h2o2.py"
    --input "${H2O2_ZERO_D_WORK_DIRECTORY}/zero_d_h2o2.csv"
)
set_tests_properties(
  regression_zero_d_h2o2_structure
  PROPERTIES DEPENDS regression_zero_d_h2o2_run
)

if(PELEF_ENABLE_CANTERA_REFERENCE)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import cantera"
    RESULT_VARIABLE PELEF_CANTERA_IMPORT_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(NOT PELEF_CANTERA_IMPORT_RESULT EQUAL 0)
    message(FATAL_ERROR "PELEF_ENABLE_CANTERA_REFERENCE requires Python Cantera")
  endif()
  add_test(
    NAME regression_zero_d_h2o2_cantera
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_h2o2_cantera.py"
      --input "${H2O2_ZERO_D_WORK_DIRECTORY}/zero_d_h2o2.csv"
      --mechanism
      "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_elementary_cantera.yaml"
  )
  set_tests_properties(
    regression_zero_d_h2o2_cantera
    PROPERTIES DEPENDS regression_zero_d_h2o2_run TIMEOUT 180
  )
endif()

