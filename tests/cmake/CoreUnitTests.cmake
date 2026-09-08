add_executable(test_eos unit/test_eos.F90)
target_link_libraries(test_eos PRIVATE pelef_core)
add_test(NAME unit_eos COMMAND test_eos)

add_executable(test_rusanov unit/test_rusanov.F90)
target_link_libraries(test_rusanov PRIVATE pelef_core)
add_test(NAME unit_rusanov COMMAND test_rusanov)

add_executable(test_pelec_riemann unit/test_pelec_riemann.F90)
target_link_libraries(test_pelec_riemann PRIVATE pelef_core)
add_test(NAME unit_pelec_riemann COMMAND test_pelec_riemann)

add_executable(test_slope_limiter unit/test_slope_limiter.F90)
target_link_libraries(test_slope_limiter PRIVATE pelef_core)
add_test(NAME unit_slope_limiter COMMAND test_slope_limiter)

add_executable(test_eb_geometry_2d unit/test_eb_geometry_2d.F90)
target_link_libraries(test_eb_geometry_2d PRIVATE pelef_core)
add_test(NAME unit_eb_geometry_2d COMMAND test_eb_geometry_2d)

add_executable(
  test_amr_eb_hierarchy_2d
  unit/test_amr_eb_hierarchy_2d.F90
)
target_link_libraries(test_amr_eb_hierarchy_2d PRIVATE pelef_core)
add_test(NAME unit_amr_eb_hierarchy_2d COMMAND test_amr_eb_hierarchy_2d)

add_executable(
  test_amr_eb_multilevel_2d
  unit/test_amr_eb_multilevel_2d.F90
)
target_link_libraries(test_amr_eb_multilevel_2d PRIVATE pelef_core)
add_test(NAME unit_amr_eb_multilevel_2d COMMAND test_amr_eb_multilevel_2d)

add_executable(
  test_amr_eb_flux_register_2d
  unit/test_amr_eb_flux_register_2d.F90
)
target_link_libraries(test_amr_eb_flux_register_2d PRIVATE pelef_core)
add_test(NAME unit_amr_eb_flux_register_2d COMMAND test_amr_eb_flux_register_2d)

add_executable(
  test_amr_eb_reactive_2d
  unit/test_amr_eb_reactive_2d.F90
)
target_link_libraries(test_amr_eb_reactive_2d PRIVATE pelef_core)
add_test(NAME unit_amr_eb_reactive_2d COMMAND test_amr_eb_reactive_2d)

add_executable(
  test_amr_eb_regrid_2d
  unit/test_amr_eb_regrid_2d.F90
)
target_link_libraries(test_amr_eb_regrid_2d PRIVATE pelef_core)
add_test(NAME unit_amr_eb_regrid_2d COMMAND test_amr_eb_regrid_2d)

add_executable(
  test_amr_eb_multipatch_2d
  unit/test_amr_eb_multipatch_2d.F90
)
target_link_libraries(test_amr_eb_multipatch_2d PRIVATE pelef_core)
add_test(NAME unit_amr_eb_multipatch_2d COMMAND test_amr_eb_multipatch_2d)

add_executable(
  test_reactive_eb_amr_2d_driver
  unit/test_reactive_eb_amr_2d_driver.F90
)
target_link_libraries(test_reactive_eb_amr_2d_driver PRIVATE pelef_core)
add_test(
  NAME unit_reactive_eb_amr_2d_driver
  COMMAND test_reactive_eb_amr_2d_driver
)

add_executable(
  test_eb_reactive_wall_flux_2d
  unit/test_eb_reactive_wall_flux_2d.F90
)
target_link_libraries(test_eb_reactive_wall_flux_2d PRIVATE pelef_core)
add_test(
  NAME unit_eb_reactive_wall_flux_2d
  COMMAND test_eb_reactive_wall_flux_2d
)

add_executable(
  test_eb_reactive_flux_divergence_2d
  unit/test_eb_reactive_flux_divergence_2d.F90
)
target_link_libraries(test_eb_reactive_flux_divergence_2d PRIVATE pelef_core)
add_test(
  NAME unit_eb_reactive_flux_divergence_2d
  COMMAND test_eb_reactive_flux_divergence_2d
)

add_executable(
  test_eb_reactive_redistribution_2d
  unit/test_eb_reactive_redistribution_2d.F90
)
target_link_libraries(test_eb_reactive_redistribution_2d PRIVATE pelef_core)
add_test(
  NAME unit_eb_reactive_redistribution_2d
  COMMAND test_eb_reactive_redistribution_2d
)

add_executable(
  test_eb_reactive_hydro_2d
  unit/test_eb_reactive_hydro_2d.F90
)
target_link_libraries(test_eb_reactive_hydro_2d PRIVATE pelef_core)
add_test(
  NAME unit_eb_reactive_hydro_2d
  COMMAND test_eb_reactive_hydro_2d
)

add_executable(
  test_reactive_eb_2d_driver
  unit/test_reactive_eb_2d_driver.F90
)
target_link_libraries(test_reactive_eb_2d_driver PRIVATE pelef_core)
add_test(
  NAME unit_reactive_eb_2d_driver
  COMMAND test_reactive_eb_2d_driver
)

add_executable(
  test_reactive_eb_transport_2d
  unit/test_reactive_eb_transport_2d.F90
)
target_link_libraries(test_reactive_eb_transport_2d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_eb_transport_2d
  COMMAND test_reactive_eb_transport_2d
)

add_executable(
  test_amr_eb_transport_2d
  unit/test_amr_eb_transport_2d.F90
)
target_link_libraries(test_amr_eb_transport_2d PRIVATE pelef_core)
add_test(
  NAME unit_amr_eb_transport_2d
  COMMAND test_amr_eb_transport_2d
)

add_executable(
  test_amr_eb_multilevel_transport_2d
  unit/test_amr_eb_multilevel_transport_2d.F90
)
target_link_libraries(test_amr_eb_multilevel_transport_2d PRIVATE pelef_core)
add_test(
  NAME unit_amr_eb_multilevel_transport_2d
  COMMAND test_amr_eb_multilevel_transport_2d
)

add_executable(test_amr_hierarchy_1d unit/test_amr_hierarchy_1d.F90)
target_link_libraries(test_amr_hierarchy_1d PRIVATE pelef_core)
add_test(NAME unit_amr_hierarchy_1d COMMAND test_amr_hierarchy_1d)

add_executable(
  test_amr_multilevel_hierarchy_1d
  unit/test_amr_multilevel_hierarchy_1d.F90
)
target_link_libraries(test_amr_multilevel_hierarchy_1d PRIVATE pelef_core)
add_test(
  NAME unit_amr_multilevel_hierarchy_1d
  COMMAND test_amr_multilevel_hierarchy_1d
)

add_executable(
  test_amr_multipatch_1d
  unit/test_amr_multipatch_1d.F90
)
target_link_libraries(test_amr_multipatch_1d PRIVATE pelef_core)
add_test(NAME unit_amr_multipatch_1d COMMAND test_amr_multipatch_1d)

add_executable(
  test_amr_patch_tree_1d
  unit/test_amr_patch_tree_1d.F90
)
target_link_libraries(test_amr_patch_tree_1d PRIVATE pelef_core)
add_test(NAME unit_amr_patch_tree_1d COMMAND test_amr_patch_tree_1d)

add_executable(
  test_amr_patch_tree_reactive_1d
  regression/test_amr_patch_tree_reactive_1d.F90
)
target_link_libraries(test_amr_patch_tree_reactive_1d PRIVATE pelef_core)
add_test(
  NAME regression_amr_patch_tree_reactive_1d
  COMMAND test_amr_patch_tree_reactive_1d
)
set_tests_properties(
  regression_amr_patch_tree_reactive_1d
  PROPERTIES TIMEOUT 300
)
add_test(
  NAME regression_amr_patch_tree_reactive_1d_fixed_checkpoint_sha256
  COMMAND "${CMAKE_COMMAND}" -E sha256sum
    "${CMAKE_CURRENT_BINARY_DIR}/amr_patch_tree_reactive.chk"
)
set_tests_properties(
  regression_amr_patch_tree_reactive_1d_fixed_checkpoint_sha256
  PROPERTIES
    DEPENDS regression_amr_patch_tree_reactive_1d
    PASS_REGULAR_EXPRESSION
      "1337b54d2761f2e3a2d25e264e9d6d073d1934756b206da16634d3140d933ed6"
)

add_executable(
  test_selected_amr_patch_tree_checkpoint_1d
  unit/test_selected_amr_patch_tree_checkpoint_1d.F90
)
target_link_libraries(
  test_selected_amr_patch_tree_checkpoint_1d PRIVATE pelef_core
)
add_test(
  NAME unit_selected_amr_patch_tree_checkpoint_1d
  COMMAND test_selected_amr_patch_tree_checkpoint_1d
)

add_executable(
  test_amr_multipatch_reactive_1d
  regression/test_amr_multipatch_reactive_1d.F90
)
target_link_libraries(test_amr_multipatch_reactive_1d PRIVATE pelef_core)
add_test(
  NAME regression_amr_multipatch_reactive_1d
  COMMAND test_amr_multipatch_reactive_1d
)
set_tests_properties(
  regression_amr_multipatch_reactive_1d
  PROPERTIES TIMEOUT 300
)

add_executable(
  test_amr_multipatch_reactive_physics_1d
  regression/test_amr_multipatch_reactive_physics_1d.F90
)
target_link_libraries(
  test_amr_multipatch_reactive_physics_1d PRIVATE pelef_core
)
add_test(
  NAME regression_amr_multipatch_reactive_physics_1d
  COMMAND test_amr_multipatch_reactive_physics_1d
)
set_tests_properties(
  regression_amr_multipatch_reactive_physics_1d
  PROPERTIES TIMEOUT 300
)

set(AMR_MULTIPATCH_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_multipatch_entropy_wave")
file(MAKE_DIRECTORY "${AMR_MULTIPATCH_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_multipatch_entropy_wave_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/amr_multipatch_entropy_wave/entropy_wave.nml"
)
set_tests_properties(
  regression_amr_multipatch_entropy_wave_run
  PROPERTIES WORKING_DIRECTORY "${AMR_MULTIPATCH_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_multipatch_entropy_wave_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_multipatch_output.py"
    --input
    "${AMR_MULTIPATCH_WORK_DIRECTORY}/amr_multipatch_entropy_wave.csv"
    --minimum-patches 2
)
set_tests_properties(
  regression_amr_multipatch_entropy_wave_check
  PROPERTIES DEPENDS regression_amr_multipatch_entropy_wave_run
)

add_executable(
  test_amr_multipatch_dynamic_1d
  regression/test_amr_multipatch_dynamic_1d.F90
)
target_link_libraries(test_amr_multipatch_dynamic_1d PRIVATE pelef_core)
add_test(
  NAME regression_amr_multipatch_dynamic_1d
  COMMAND test_amr_multipatch_dynamic_1d
)
set_tests_properties(
  regression_amr_multipatch_dynamic_1d
  PROPERTIES TIMEOUT 300
)

add_executable(test_amr_regrid_1d unit/test_amr_regrid_1d.F90)
target_link_libraries(test_amr_regrid_1d PRIVATE pelef_core)
add_test(NAME unit_amr_regrid_1d COMMAND test_amr_regrid_1d)

add_executable(
  test_amr_reactive_1d
  regression/test_amr_reactive_1d.F90
)
target_link_libraries(test_amr_reactive_1d PRIVATE pelef_core)
add_test(NAME regression_amr_reactive_1d COMMAND test_amr_reactive_1d)
set_tests_properties(regression_amr_reactive_1d PROPERTIES TIMEOUT 240)

add_executable(
  test_amr_reactive_transport_1d
  regression/test_amr_reactive_transport_1d.F90
)
target_link_libraries(test_amr_reactive_transport_1d PRIVATE pelef_core)
add_test(
  NAME regression_amr_reactive_transport_1d
  COMMAND test_amr_reactive_transport_1d
)
set_tests_properties(
  regression_amr_reactive_transport_1d
  PROPERTIES TIMEOUT 300
)

add_executable(
  test_amr_reactive_plm
  regression/test_amr_reactive_plm.F90
)
target_link_libraries(test_amr_reactive_plm PRIVATE pelef_core)
add_test(NAME regression_amr_reactive_plm COMMAND test_amr_reactive_plm)
set_tests_properties(regression_amr_reactive_plm PROPERTIES TIMEOUT 240)

set(AMR_REACTIVE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_reactive_hotspot")
file(MAKE_DIRECTORY "${AMR_REACTIVE_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_reactive_hotspot_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/amr_reactive_hotspot/hotspot.nml"
)
set_tests_properties(
  regression_amr_reactive_hotspot_run
  PROPERTIES WORKING_DIRECTORY "${AMR_REACTIVE_WORK_DIRECTORY}" TIMEOUT 300
)
add_test(
  NAME regression_amr_reactive_hotspot_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_hotspot.py"
    --input "${AMR_REACTIVE_WORK_DIRECTORY}/amr_reactive_hotspot.csv"
)
set_tests_properties(
  regression_amr_reactive_hotspot_check
  PROPERTIES DEPENDS regression_amr_reactive_hotspot_run
)

set(AMR_MULTILEVEL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_multilevel_reactive_hotspot")
file(MAKE_DIRECTORY "${AMR_MULTILEVEL_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_multilevel_hotspot_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/amr_multilevel_reactive_hotspot/hotspot.nml"
)
set_tests_properties(
  regression_amr_multilevel_hotspot_run
  PROPERTIES WORKING_DIRECTORY "${AMR_MULTILEVEL_WORK_DIRECTORY}" TIMEOUT 300
)
add_test(
  NAME regression_amr_multilevel_hotspot_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_multilevel_hotspot.py"
    --input
    "${AMR_MULTILEVEL_WORK_DIRECTORY}/amr_multilevel_reactive_hotspot.csv"
)
set_tests_properties(
  regression_amr_multilevel_hotspot_check
  PROPERTIES DEPENDS regression_amr_multilevel_hotspot_run
)

set(AMR_MULTILEVEL_PPM_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_multilevel_characteristic_ppm_hotspot")
file(MAKE_DIRECTORY "${AMR_MULTILEVEL_PPM_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_multilevel_characteristic_ppm_hotspot_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/amr_multilevel_reactive_hotspot/hotspot_characteristic_ppm.nml"
)
set_tests_properties(
  regression_amr_multilevel_characteristic_ppm_hotspot_run
  PROPERTIES WORKING_DIRECTORY "${AMR_MULTILEVEL_PPM_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_multilevel_characteristic_ppm_hotspot_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_multilevel_hotspot.py"
    --input
    "${AMR_MULTILEVEL_PPM_WORK_DIRECTORY}/amr_multilevel_characteristic_ppm_hotspot.csv"
)
set_tests_properties(
  regression_amr_multilevel_characteristic_ppm_hotspot_check
  PROPERTIES DEPENDS
    regression_amr_multilevel_characteristic_ppm_hotspot_run
)

set(AMR_MULTILEVEL_WENO_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_multilevel_weno5z_hotspot")
file(MAKE_DIRECTORY "${AMR_MULTILEVEL_WENO_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_multilevel_weno5z_hotspot_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/amr_multilevel_reactive_hotspot/hotspot_weno5z.nml"
)
set_tests_properties(
  regression_amr_multilevel_weno5z_hotspot_run
  PROPERTIES WORKING_DIRECTORY "${AMR_MULTILEVEL_WENO_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_multilevel_weno5z_hotspot_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_multilevel_hotspot.py"
    --input
    "${AMR_MULTILEVEL_WENO_WORK_DIRECTORY}/amr_multilevel_weno5z_hotspot.csv"
)
set_tests_properties(
  regression_amr_multilevel_weno5z_hotspot_check
  PROPERTIES DEPENDS regression_amr_multilevel_weno5z_hotspot_run
)

foreach(WENO_VARIANT IN ITEMS weno7z weno3z)
  set(WENO_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/amr_multilevel_${WENO_VARIANT}_hotspot")
  file(MAKE_DIRECTORY "${WENO_WORK_DIRECTORY}")
  add_test(
    NAME regression_amr_multilevel_${WENO_VARIANT}_hotspot_run
    COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
      "${PROJECT_SOURCE_DIR}/cases/amr_multilevel_reactive_hotspot/hotspot_${WENO_VARIANT}.nml"
  )
  set_tests_properties(
    regression_amr_multilevel_${WENO_VARIANT}_hotspot_run
    PROPERTIES WORKING_DIRECTORY "${WENO_WORK_DIRECTORY}"
      TIMEOUT 300
  )
  add_test(
    NAME regression_amr_multilevel_${WENO_VARIANT}_hotspot_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_amr_multilevel_hotspot.py"
      --input
      "${WENO_WORK_DIRECTORY}/amr_multilevel_${WENO_VARIANT}_hotspot.csv"
  )
  set_tests_properties(
    regression_amr_multilevel_${WENO_VARIANT}_hotspot_check
    PROPERTIES DEPENDS
      regression_amr_multilevel_${WENO_VARIANT}_hotspot_run
  )
endforeach()

set(AMR_BOUNDARY_WENO_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/amr_multilevel_boundary_weno7z_hotspot")
file(MAKE_DIRECTORY "${AMR_BOUNDARY_WENO_WORK_DIRECTORY}")
add_test(
  NAME regression_amr_multilevel_boundary_weno7z_hotspot_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/amr_multilevel_reactive_hotspot/hotspot_boundary_weno7z.nml"
)
set_tests_properties(
  regression_amr_multilevel_boundary_weno7z_hotspot_run
  PROPERTIES WORKING_DIRECTORY "${AMR_BOUNDARY_WENO_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_amr_multilevel_boundary_weno7z_hotspot_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_amr_multilevel_hotspot.py"
    --input
    "${AMR_BOUNDARY_WENO_WORK_DIRECTORY}/amr_multilevel_boundary_weno7z_hotspot.csv"
)
set_tests_properties(
  regression_amr_multilevel_boundary_weno7z_hotspot_check
  PROPERTIES DEPENDS regression_amr_multilevel_boundary_weno7z_hotspot_run
)

add_executable(test_plm_reconstruction unit/test_plm_reconstruction.F90)
target_link_libraries(test_plm_reconstruction PRIVATE pelef_core)
add_test(NAME unit_plm_reconstruction COMMAND test_plm_reconstruction)

add_executable(
  test_pelec_plm_tracing
  unit/test_pelec_plm_tracing.F90
)
target_link_libraries(test_pelec_plm_tracing PRIVATE pelef_core)
add_test(NAME unit_pelec_plm_tracing COMMAND test_pelec_plm_tracing)

add_executable(
  test_pelec_fourth_order_flattening
  unit/test_pelec_fourth_order_flattening.F90
)
target_link_libraries(test_pelec_fourth_order_flattening PRIVATE pelef_core)
add_test(
  NAME unit_pelec_fourth_order_flattening
  COMMAND test_pelec_fourth_order_flattening
)

add_executable(test_entropy_wave regression/test_entropy_wave.F90)
target_link_libraries(test_entropy_wave PRIVATE pelef_core)
add_test(NAME regression_entropy_wave_convergence COMMAND test_entropy_wave)

add_executable(
  test_entropy_wave_pelec
  regression/test_entropy_wave_pelec.F90
)
target_link_libraries(test_entropy_wave_pelec PRIVATE pelef_core)
add_test(
  NAME regression_entropy_wave_pelec_convergence
  COMMAND test_entropy_wave_pelec
)

add_executable(
  test_entropy_wave_pelec_plm
  regression/test_entropy_wave_pelec_plm.F90
)
target_link_libraries(test_entropy_wave_pelec_plm PRIVATE pelef_core)
add_test(
  NAME regression_entropy_wave_pelec_plm_convergence
  COMMAND test_entropy_wave_pelec_plm
)

add_executable(
  test_entropy_wave_pelec_plm4
  regression/test_entropy_wave_pelec_plm4.F90
)
target_link_libraries(test_entropy_wave_pelec_plm4 PRIVATE pelef_core)
add_test(
  NAME regression_entropy_wave_pelec_plm4_convergence
  COMMAND test_entropy_wave_pelec_plm4
)

set(SOD_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sod")
file(MAKE_DIRECTORY "${SOD_WORK_DIRECTORY}")

add_test(
  NAME regression_sod_run
  COMMAND $<TARGET_FILE:pelef> "${PROJECT_SOURCE_DIR}/cases/sod/sod.nml"
)
set_tests_properties(
  regression_sod_run
  PROPERTIES
    WORKING_DIRECTORY "${SOD_WORK_DIRECTORY}"
)

add_test(
  NAME regression_sod_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_sod.py"
    --input "${SOD_WORK_DIRECTORY}/sod.csv"
    --nx 400
    --time 0.2
    --gamma 1.4
    --density-l1-max 3.0e-2
    --pressure-l1-max 3.0e-2
    --mass-error-max 2.0e-12
    --energy-error-max 2.0e-12
    --momentum-error-max 2.0e-12
)
set_tests_properties(
  regression_sod_compare
  PROPERTIES
    DEPENDS regression_sod_run
)

set(SOD_PLM_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sod_plm")
file(MAKE_DIRECTORY "${SOD_PLM_WORK_DIRECTORY}")

add_test(
  NAME regression_sod_plm_run
  COMMAND $<TARGET_FILE:pelef> "${PROJECT_SOURCE_DIR}/cases/sod/sod_plm.nml"
)
set_tests_properties(
  regression_sod_plm_run
  PROPERTIES
    WORKING_DIRECTORY "${SOD_PLM_WORK_DIRECTORY}"
)

add_test(
  NAME regression_sod_plm_compare
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/compare_sod.py"
    --input "${SOD_PLM_WORK_DIRECTORY}/sod_plm.csv"
    --nx 400
    --time 0.2
    --gamma 1.4
    --density-l1-max 4.0e-3
    --pressure-l1-max 3.0e-3
    --mass-error-max 2.0e-12
    --energy-error-max 2.0e-12
    --momentum-error-max 2.0e-12
)
set_tests_properties(
  regression_sod_plm_compare
  PROPERTIES
    DEPENDS regression_sod_plm_run
)

set(SOD_PELEC_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sod_pelec")
file(MAKE_DIRECTORY "${SOD_PELEC_WORK_DIRECTORY}")

