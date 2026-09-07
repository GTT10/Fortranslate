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

if(PELEF_ENABLE_MPI)
  set(
    MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_3d_transport_restart_0245"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_transport_restart_0245"
  )
  file(MAKE_DIRECTORY
    "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
  )
  add_test(
    NAME regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/checkpoint-np2.log"
      -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
      "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
      fixed_transport_checkpoint_np2 ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    PROPERTIES
      WORKING_DIRECTORY
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
      TIMEOUT 300
  )
  set(MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS)
  foreach(transport_restart_rank IN ITEMS 1 2 4 8)
    set(
      transport_restart_test
      regression_mpi_amr_reactive_3d_transport_restart_0245_np${transport_restart_rank}
    )
    add_test(
      NAME ${transport_restart_test}
      COMMAND "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np${transport_restart_rank}.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        ${transport_restart_rank} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
        "fixed_transport_restart_mpi_np${transport_restart_rank}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${transport_restart_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${transport_restart_rank}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${transport_restart_rank}"
        DEPENDS
          regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
        TIMEOUT 300
    )
    list(APPEND MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS
      ${transport_restart_test})
  endforeach()
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_mpi_amr_reactive_3d_transport_restart_0245_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_reference_0245"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np1"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np2"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np4"
      --candidate-prefix
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/fixed_transport_restart_mpi_np8"
      --checkpoint
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
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
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np1.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np2.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np4.log"
      --candidate-log
        "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np8.log"
  )
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_exact
    PROPERTIES DEPENDS
      "regression_amr_reactive_3d_transport_restart_0245_reference;regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2;${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS}"
  )
  add_test(
    NAME regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
      "${MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/amr_hotspot_transport_restart_0245_fixed.chk"
  )
  set_tests_properties(
    regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    PROPERTIES DEPENDS
      "regression_amr_reactive_3d_transport_restart_0245_checkpoint;regression_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2"
  )

  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/checkpoint-np2.log"
      -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
      selected_transport_checkpoint_np2 ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
      TIMEOUT 600
  )
  set(SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS)
  foreach(selected_transport_restart_rank IN ITEMS 1 2 4 8)
    set(
      selected_transport_restart_test
      regression_selected_mpi_amr_reactive_3d_transport_restart_0245_np${selected_transport_restart_rank}
    )
    add_test(
      NAME ${selected_transport_restart_test}
      COMMAND "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np${selected_transport_restart_rank}.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        ${selected_transport_restart_rank} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
        "selected_transport_restart_mpi_np${selected_transport_restart_rank}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_transport_restart_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_transport_restart_rank}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_transport_restart_rank}"
        DEPENDS
          regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2
        TIMEOUT 600
    )
    list(APPEND SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS
      ${selected_transport_restart_test})
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_restart_0245_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_reference_0245"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_transport_restart_mpi_np8"
      --checkpoint
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
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
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np1.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np2.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np4.log"
      --candidate-log
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/restart-np8.log"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_transport_restart_0245_reference;regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2;${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_TESTS}"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
      "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RESTART_0245_WORK_DIRECTORY}/selected_h2o2_full_transport_restart_0245.chk"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_transport_restart_0245_checkpoint;regression_selected_mpi_amr_reactive_3d_transport_restart_0245_checkpoint_np2"
  )
endif()
file(MAKE_DIRECTORY "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_restart_3d_reference
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d_restart/reference.nml"
)
add_test(
  NAME regression_reactive_eb_restart_3d_checkpoint_stop
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_restart_3d_restart
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:pelef_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_3d_restart/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_reference
  regression_reactive_eb_restart_3d_checkpoint_stop
  regression_reactive_eb_restart_3d_restart
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_restart_3d_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Operator sequence: R-T-H-T-R"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    DEPENDS regression_reactive_eb_restart_3d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_restart_3d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_3d_restart.py"
    --reference
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart_reference.csv"
    --stopped
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart_stopped.csv"
    --restarted
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart.csv"
    --checkpoint
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart.chk"
    --reference-log "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reference.log"
    --stopped-log "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/stopped.log"
    --restarted-log "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/restarted.log"
    --expected-checkpoint-sha256
      "${REACTIVE_EB_RESTART_3D_CHECKPOINT_SHA256}"
    --expected-stopped-sha256
      "${REACTIVE_EB_RESTART_3D_STOPPED_SHA256}"
)
set_tests_properties(
  regression_reactive_eb_restart_3d_check
  PROPERTIES
    DEPENDS
      "regression_reactive_eb_restart_3d_reference;regression_reactive_eb_restart_3d_checkpoint_stop;regression_reactive_eb_restart_3d_restart"
    TIMEOUT 60
)

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

add_executable(test_reactive_general_eos unit/test_reactive_general_eos.F90)
target_link_libraries(test_reactive_general_eos PRIVATE pelef_core)
add_test(NAME unit_reactive_general_eos COMMAND test_reactive_general_eos)

add_executable(test_reactive_hllc unit/test_reactive_hllc.F90)
target_link_libraries(test_reactive_hllc PRIVATE pelef_core)
add_test(NAME unit_reactive_hllc COMMAND test_reactive_hllc)

add_executable(test_reactive_pelec_riemann unit/test_reactive_pelec_riemann.F90)
target_link_libraries(test_reactive_pelec_riemann PRIVATE pelef_core)
add_test(NAME unit_reactive_pelec_riemann COMMAND test_reactive_pelec_riemann)

add_executable(test_reactive_ppm unit/test_reactive_ppm.F90)
target_link_libraries(test_reactive_ppm PRIVATE pelef_core)
add_test(NAME unit_reactive_ppm COMMAND test_reactive_ppm)

add_executable(test_weno_reconstruction unit/test_weno_reconstruction.F90)
target_link_libraries(test_weno_reconstruction PRIVATE pelef_core)
add_test(NAME unit_weno_reconstruction COMMAND test_weno_reconstruction)

add_executable(
  test_reactive_ppm_convergence
  regression/test_reactive_ppm_convergence.F90
)
target_link_libraries(test_reactive_ppm_convergence PRIVATE pelef_core)
add_test(
  NAME regression_reactive_ppm_convergence
  COMMAND test_reactive_ppm_convergence
)
set_tests_properties(
  regression_reactive_ppm_convergence
  PROPERTIES TIMEOUT 240
)

add_executable(
  test_reactive_characteristic_ppm_convergence
  regression/test_reactive_characteristic_ppm_convergence.F90
)
target_link_libraries(
  test_reactive_characteristic_ppm_convergence PRIVATE pelef_core
)
add_test(
  NAME regression_reactive_characteristic_ppm_convergence
  COMMAND test_reactive_characteristic_ppm_convergence
)
set_tests_properties(
  regression_reactive_characteristic_ppm_convergence
  PROPERTIES TIMEOUT 240
)

add_executable(
  test_reactive_ppm_shock
  regression/test_reactive_ppm_shock.F90
)
target_link_libraries(test_reactive_ppm_shock PRIVATE pelef_core)
add_test(
  NAME regression_reactive_ppm_shock
  COMMAND test_reactive_ppm_shock
)
set_tests_properties(
  regression_reactive_ppm_shock
  PROPERTIES TIMEOUT 240
)

add_executable(
  test_reactive_uniform_reduction
  unit/test_reactive_uniform_reduction.F90
)
target_link_libraries(test_reactive_uniform_reduction PRIVATE pelef_core)
add_test(
  NAME unit_reactive_uniform_reduction
  COMMAND test_reactive_uniform_reduction
)
set_tests_properties(unit_reactive_uniform_reduction PROPERTIES TIMEOUT 180)

add_executable(
  test_reactive_entropy_wave
  regression/test_reactive_entropy_wave.F90
)
target_link_libraries(test_reactive_entropy_wave PRIVATE pelef_core)
add_test(
  NAME regression_reactive_entropy_wave
  COMMAND test_reactive_entropy_wave
)
set_tests_properties(regression_reactive_entropy_wave PROPERTIES TIMEOUT 180)

add_executable(
  test_reactive_composition_wave
  regression/test_reactive_composition_wave.F90
)
target_link_libraries(test_reactive_composition_wave PRIVATE pelef_core)
add_test(
  NAME regression_reactive_composition_wave
  COMMAND test_reactive_composition_wave
)
set_tests_properties(
  regression_reactive_composition_wave
  PROPERTIES TIMEOUT 180
)

set(
  REACTIVE_COMPOSITION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_composition_wave"
)
file(MAKE_DIRECTORY "${REACTIVE_COMPOSITION_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_composition_wave_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_composition_wave/composition_wave.nml"
)
set_tests_properties(
  regression_reactive_composition_wave_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_COMPOSITION_WORK_DIRECTORY}"
    TIMEOUT 180
)
add_test(
  NAME regression_reactive_composition_wave_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_composition_wave.py"
    --input
    "${REACTIVE_COMPOSITION_WORK_DIRECTORY}/reactive_composition_wave.csv"
)
set_tests_properties(
  regression_reactive_composition_wave_check
  PROPERTIES DEPENDS regression_reactive_composition_wave_run
)

add_executable(
  test_reactive_material_contact
  regression/test_reactive_material_contact.F90
)
target_link_libraries(test_reactive_material_contact PRIVATE pelef_core)
add_test(
  NAME regression_reactive_material_contact
  COMMAND test_reactive_material_contact
)
set_tests_properties(
  regression_reactive_material_contact
  PROPERTIES TIMEOUT 180
)

set(REACTIVE_HOTSPOT_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/reactive_hotspot")
file(MAKE_DIRECTORY "${REACTIVE_HOTSPOT_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_hotspot_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot/hotspot.nml"
)
set_tests_properties(
  regression_reactive_hotspot_run
  PROPERTIES WORKING_DIRECTORY "${REACTIVE_HOTSPOT_WORK_DIRECTORY}" TIMEOUT 180
)
add_test(
  NAME regression_reactive_hotspot_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot.py"
    --input "${REACTIVE_HOTSPOT_WORK_DIRECTORY}/reactive_hotspot.csv"
)
set_tests_properties(
  regression_reactive_hotspot_check
  PROPERTIES DEPENDS regression_reactive_hotspot_run
)

set(
  REACTIVE_HOTSPOT_HLLC_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_hotspot_hllc"
)
file(MAKE_DIRECTORY "${REACTIVE_HOTSPOT_HLLC_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_hotspot_hllc_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot/hotspot_hllc.nml"
)
set_tests_properties(
  regression_reactive_hotspot_hllc_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_HOTSPOT_HLLC_WORK_DIRECTORY}"
    TIMEOUT 180
)
add_test(
  NAME regression_reactive_hotspot_hllc_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot.py"
    --input
    "${REACTIVE_HOTSPOT_HLLC_WORK_DIRECTORY}/reactive_hotspot_hllc.csv"
)
set_tests_properties(
  regression_reactive_hotspot_hllc_check
  PROPERTIES DEPENDS regression_reactive_hotspot_hllc_run
)

set(
  REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_hotspot_characteristic_ppm"
)
file(MAKE_DIRECTORY "${REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_hotspot_characteristic_ppm_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot/hotspot_characteristic_ppm.nml"
)
set_tests_properties(
  regression_reactive_hotspot_characteristic_ppm_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_WORK_DIRECTORY}"
    TIMEOUT 180
)
add_test(
  NAME regression_reactive_hotspot_characteristic_ppm_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot.py"
    --input
    "${REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_WORK_DIRECTORY}/reactive_hotspot_characteristic_ppm.csv"
)
set_tests_properties(
  regression_reactive_hotspot_characteristic_ppm_check
  PROPERTIES DEPENDS regression_reactive_hotspot_characteristic_ppm_run
)

add_executable(
  test_reactive_hotspot_reference
  regression/test_reactive_hotspot_reference.F90
)
target_link_libraries(test_reactive_hotspot_reference PRIVATE pelef_core)
add_test(
  NAME regression_reactive_hotspot_reference
  COMMAND test_reactive_hotspot_reference
)
set_tests_properties(regression_reactive_hotspot_reference PROPERTIES TIMEOUT 240)

add_executable(
  test_reactive_directional_flux_2d
  unit/test_reactive_directional_flux_2d.F90
)
target_link_libraries(test_reactive_directional_flux_2d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_directional_flux_2d
  COMMAND test_reactive_directional_flux_2d
)

add_executable(
  test_reactive_directional_flux_3d
  unit/test_reactive_directional_flux_3d.F90
)
target_link_libraries(test_reactive_directional_flux_3d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_directional_flux_3d
  COMMAND test_reactive_directional_flux_3d
)

add_executable(
  test_reactive_ctu_dimensional_reduction
  unit/test_reactive_ctu_dimensional_reduction.F90
)
target_link_libraries(test_reactive_ctu_dimensional_reduction PRIVATE pelef_core)
add_test(
  NAME unit_reactive_ctu_dimensional_reduction
  COMMAND test_reactive_ctu_dimensional_reduction
)
set_tests_properties(
  unit_reactive_ctu_dimensional_reduction
  PROPERTIES TIMEOUT 180
)

add_executable(
  test_reactive_characteristic_ppm_cell_2d
  unit/test_reactive_characteristic_ppm_cell_2d.F90
)
target_link_libraries(
  test_reactive_characteristic_ppm_cell_2d PRIVATE pelef_core
)
add_test(
  NAME unit_reactive_characteristic_ppm_cell_2d
  COMMAND test_reactive_characteristic_ppm_cell_2d
)

add_executable(
  test_reactive_diagonal_wave_2d
  regression/test_reactive_diagonal_wave_2d.F90
)
target_link_libraries(test_reactive_diagonal_wave_2d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_diagonal_wave_2d
  COMMAND test_reactive_diagonal_wave_2d
)
set_tests_properties(
  regression_reactive_diagonal_wave_2d
  PROPERTIES TIMEOUT 480
)

add_executable(
  test_reactive_characteristic_ppm_2d
  regression/test_reactive_characteristic_ppm_2d.F90
)
target_link_libraries(test_reactive_characteristic_ppm_2d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_characteristic_ppm_2d
  COMMAND test_reactive_characteristic_ppm_2d
)
set_tests_properties(
  regression_reactive_characteristic_ppm_2d
  PROPERTIES TIMEOUT 480
)

add_executable(
  test_reactive_material_contact_2d
  regression/test_reactive_material_contact_2d.F90
)
target_link_libraries(test_reactive_material_contact_2d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_material_contact_2d
  COMMAND test_reactive_material_contact_2d
)
set_tests_properties(
  regression_reactive_material_contact_2d
  PROPERTIES TIMEOUT 480
)

add_executable(
  test_reactive_vortex_2d
  regression/test_reactive_vortex_2d.F90
)
target_link_libraries(test_reactive_vortex_2d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_vortex_2d
  COMMAND test_reactive_vortex_2d
)
set_tests_properties(
  regression_reactive_vortex_2d
  PROPERTIES TIMEOUT 300
)

set(
  REACTIVE_HOTSPOT_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_hotspot_2d"
)
file(MAKE_DIRECTORY "${REACTIVE_HOTSPOT_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_hotspot_2d_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot_2d/hotspot.nml"
)
set_tests_properties(
  regression_reactive_hotspot_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_HOTSPOT_2D_WORK_DIRECTORY}"
    TIMEOUT 480
)
add_test(
  NAME regression_reactive_hotspot_2d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot_2d.py"
    --input
    "${REACTIVE_HOTSPOT_2D_WORK_DIRECTORY}/reactive_hotspot_2d.csv"
    --nx 24
    --ny 24
)
set_tests_properties(
  regression_reactive_hotspot_2d_check
  PROPERTIES DEPENDS regression_reactive_hotspot_2d_run
)

set(
  REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_hotspot_characteristic_ppm_2d"
)
file(MAKE_DIRECTORY "${REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_hotspot_characteristic_ppm_2d_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_hotspot_2d/hotspot_characteristic_ppm.nml"
)
set_tests_properties(
  regression_reactive_hotspot_characteristic_ppm_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_2D_WORK_DIRECTORY}"
    TIMEOUT 480
)
add_test(
  NAME regression_reactive_hotspot_characteristic_ppm_2d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot_2d.py"
    --input
    "${REACTIVE_HOTSPOT_CHARACTERISTIC_PPM_2D_WORK_DIRECTORY}/reactive_hotspot_characteristic_ppm_2d.csv"
    --nx 24
    --ny 24
)
set_tests_properties(
  regression_reactive_hotspot_characteristic_ppm_2d_check
  PROPERTIES DEPENDS regression_reactive_hotspot_characteristic_ppm_2d_run
)

add_executable(
  test_reactive_ppm_shock_2d
  regression/test_reactive_ppm_shock_2d.F90
)
target_link_libraries(test_reactive_ppm_shock_2d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_ppm_shock_2d
  COMMAND test_reactive_ppm_shock_2d
)
set_tests_properties(
  regression_reactive_ppm_shock_2d
  PROPERTIES TIMEOUT 300
)


add_executable(
  test_mixture_transport
  unit/test_mixture_transport.F90
)
target_link_libraries(test_mixture_transport PRIVATE pelef_core)
add_test(
  NAME unit_mixture_transport
  COMMAND test_mixture_transport
)

add_executable(
  test_reactive_diffusive_flux
  unit/test_reactive_diffusive_flux.F90
)
target_link_libraries(test_reactive_diffusive_flux PRIVATE pelef_core)
add_test(
  NAME unit_reactive_diffusive_flux
  COMMAND test_reactive_diffusive_flux
)

add_executable(
  test_reactive_transport_1d
  regression/test_reactive_transport_1d.F90
)
target_link_libraries(test_reactive_transport_1d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_transport_1d
  COMMAND test_reactive_transport_1d
)
set_tests_properties(
  regression_reactive_transport_1d
  PROPERTIES TIMEOUT 300
)

set(
  TRANSPORT_PROBE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/transport_probe"
)
file(MAKE_DIRECTORY "${TRANSPORT_PROBE_WORK_DIRECTORY}")
add_test(
  NAME regression_transport_probe_run
  COMMAND
    $<TARGET_FILE:pelef_transport_probe>
    "${TRANSPORT_PROBE_WORK_DIRECTORY}/transport_probe.csv"
)
set_tests_properties(
  regression_transport_probe_run
  PROPERTIES WORKING_DIRECTORY "${TRANSPORT_PROBE_WORK_DIRECTORY}"
)
if(PELEF_ENABLE_CANTERA_REFERENCE)
  add_test(
    NAME regression_transport_cantera
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_transport_cantera.py"
      --input
      "${TRANSPORT_PROBE_WORK_DIRECTORY}/transport_probe.csv"
  )
  set_tests_properties(
    regression_transport_cantera
    PROPERTIES DEPENDS regression_transport_probe_run
  )
endif()

set(
  REACTIVE_TRANSPORT_1D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_transport_1d"
)
file(MAKE_DIRECTORY "${REACTIVE_TRANSPORT_1D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_transport_1d_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_transport_1d/transport_pulse.nml"
)
set_tests_properties(
  regression_reactive_transport_1d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_TRANSPORT_1D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_transport_1d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_transport_1d.py"
    --input
    "${REACTIVE_TRANSPORT_1D_WORK_DIRECTORY}/reactive_transport_pulse.csv"
    --nx 96
)
set_tests_properties(
  regression_reactive_transport_1d_check
  PROPERTIES DEPENDS regression_reactive_transport_1d_run
)


add_executable(
  test_reactive_transport_2d
  unit/test_reactive_transport_2d.F90
)
target_link_libraries(test_reactive_transport_2d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_transport_2d
  COMMAND test_reactive_transport_2d
)

add_executable(
  test_reactive_transport_2d_regression
  regression/test_reactive_transport_2d.F90
)
target_link_libraries(test_reactive_transport_2d_regression PRIVATE pelef_core)
add_test(
  NAME regression_reactive_transport_2d
  COMMAND test_reactive_transport_2d_regression
)
set_tests_properties(
  regression_reactive_transport_2d
  PROPERTIES TIMEOUT 300
)

set(
  REACTIVE_TRANSPORT_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_transport_2d"
)
file(MAKE_DIRECTORY "${REACTIVE_TRANSPORT_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_transport_2d_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_transport_2d/transport_hotspot.nml"
)
set_tests_properties(
  regression_reactive_transport_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_TRANSPORT_2D_WORK_DIRECTORY}"
    TIMEOUT 480
)
add_test(
  NAME regression_reactive_transport_2d_check
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_hotspot_2d.py"
    --input
    "${REACTIVE_TRANSPORT_2D_WORK_DIRECTORY}/reactive_transport_hotspot_2d.csv"
    --nx 20
    --ny 20
)
set_tests_properties(
  regression_reactive_transport_2d_check
  PROPERTIES DEPENDS regression_reactive_transport_2d_run
)

add_executable(
  test_reactive_boundary_2d
  unit/test_reactive_boundary_2d.F90
)
target_link_libraries(test_reactive_boundary_2d PRIVATE pelef_core)
add_test(
  NAME unit_reactive_boundary_2d
  COMMAND test_reactive_boundary_2d
)

add_executable(
  test_reactive_physical_boundaries_2d
  regression/test_reactive_physical_boundaries_2d.F90
)
target_link_libraries(test_reactive_physical_boundaries_2d PRIVATE pelef_core)
add_test(
  NAME regression_reactive_physical_boundaries_2d
  COMMAND test_reactive_physical_boundaries_2d
)
set_tests_properties(
  regression_reactive_physical_boundaries_2d
  PROPERTIES TIMEOUT 300
)

set(
  REACTIVE_PRESCRIBED_SPECIES_WALL_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_prescribed_species_wall_2d"
)
file(MAKE_DIRECTORY "${REACTIVE_PRESCRIBED_SPECIES_WALL_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_prescribed_species_wall_2d_run
  COMMAND
    $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_boundaries_2d/prescribed_species_wall.nml"
)
set_tests_properties(
  regression_reactive_prescribed_species_wall_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_PRESCRIBED_SPECIES_WALL_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)


add_executable(test_pressure_dependent_kinetics unit/test_pressure_dependent_kinetics.F90)
target_link_libraries(test_pressure_dependent_kinetics PRIVATE pelef_core)
add_test(NAME unit_pressure_dependent_kinetics COMMAND test_pressure_dependent_kinetics)

add_executable(test_full_h2o2_jacobian unit/test_full_h2o2_jacobian.F90)
target_link_libraries(test_full_h2o2_jacobian PRIVATE pelef_core)
add_test(NAME unit_full_h2o2_jacobian COMMAND test_full_h2o2_jacobian)

add_executable(test_implicit_h2o2_reactor unit/test_implicit_h2o2_reactor.F90)
target_link_libraries(test_implicit_h2o2_reactor PRIVATE pelef_core)
add_test(NAME unit_implicit_h2o2_reactor COMMAND test_implicit_h2o2_reactor)
set_tests_properties(unit_implicit_h2o2_reactor PROPERTIES TIMEOUT 180)

add_executable(
  test_selected_reactor_config
  unit/test_selected_reactor_config.F90
)
target_link_libraries(
  test_selected_reactor_config PRIVATE pelef_selected_runtime_core
)
add_test(
  NAME unit_selected_reactor_config
  COMMAND test_selected_reactor_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactor/fixture.nml"
)

add_test(
  NAME unit_generate_mechanism_bundle
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/python/test_generate_mechanism_bundle.py"
)

add_test(
  NAME unit_ingest_cantera_mechanism
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/unit/test_ingest_cantera_mechanism.py"
)

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

if(PELEF_ENABLE_MPI)
  set(
    SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_1d_fixture"
  )
  set(
    SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_1d_full"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
    "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}"
  )

  set(SELECTED_MPI_REACTIVE_1D_FIXTURE_RUN_TESTS)
  set(SELECTED_MPI_REACTIVE_1D_FULL_RUN_TESTS)
  foreach(selected_mpi_ranks IN ITEMS 1 2 4)
    set(
      fixture_test_name
      regression_selected_mpi_reactive_1d_fixture_np${selected_mpi_ranks}
    )
    add_test(
      NAME ${fixture_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_fixture_mpi_reactive_1d>
        "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np${selected_mpi_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${fixture_test_name}
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_ranks}
        PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
        TIMEOUT 300
    )
    list(APPEND SELECTED_MPI_REACTIVE_1D_FIXTURE_RUN_TESTS ${fixture_test_name})

    set(
      full_test_name
      regression_selected_mpi_reactive_1d_full_np${selected_mpi_ranks}
    )
    add_test(
      NAME ${full_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_1d>
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np${selected_mpi_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${full_test_name}
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_ranks}
        PASS_REGULAR_EXPRESSION "Chemistry integrator: implicit"
        TIMEOUT 300
    )
    list(APPEND SELECTED_MPI_REACTIVE_1D_FULL_RUN_TESTS ${full_test_name})
  endforeach()

  add_test(
    NAME regression_selected_mpi_reactive_1d_fixture_rank_parity
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_reactive_1d.py"
      "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np1.csv"
      "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np2.csv"
      "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np4.csv"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_1d_fixture_rank_parity
    PROPERTIES DEPENDS "${SELECTED_MPI_REACTIVE_1D_FIXTURE_RUN_TESTS}"
  )

  add_test(
    NAME regression_selected_mpi_reactive_1d_full_rank_parity
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_reactive_1d.py"
      "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np1.csv"
      "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np2.csv"
      "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np4.csv"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_1d_full_rank_parity
    PROPERTIES DEPENDS "${SELECTED_MPI_REACTIVE_1D_FULL_RUN_TESTS}"
  )

  set(
    SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_1d"
  )
  file(MAKE_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}")
  set(SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS)

  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_fixed_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_amr_reactive_1d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/fixed_reference.nml"
      "fixed_reference.csv"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_fixed_np1
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "MPI ranks: 1"
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
    regression_selected_mpi_amr_reactive_1d_fixed_np1
  )

  set(SELECTED_MPI_AMR_REACTIVE_1D_REFERENCE_TESTS)
  foreach(selected_mpi_amr_ranks IN ITEMS 1 2 4)
    set(
      selected_mpi_amr_reference_test
      regression_selected_mpi_amr_reactive_1d_reference_np${selected_mpi_amr_ranks}
    )
    add_test(
      NAME ${selected_mpi_amr_reference_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/selected_reference.nml"
        "selected_reference_np${selected_mpi_amr_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_amr_reference_test}
      PROPERTIES
        WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_ranks}"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_REFERENCE_TESTS
      ${selected_mpi_amr_reference_test}
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
      ${selected_mpi_amr_reference_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/checkpoint_stop.nml"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
    regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
  )

  set(SELECTED_MPI_AMR_REACTIVE_1D_RESTART_TESTS)
  foreach(selected_mpi_amr_restart_ranks IN ITEMS 2 4)
    set(
      selected_mpi_amr_restart_test
      regression_selected_mpi_amr_reactive_1d_restart_np${selected_mpi_amr_restart_ranks}
    )
    add_test(
      NAME ${selected_mpi_amr_restart_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/restart.nml"
        "selected_restart_np${selected_mpi_amr_restart_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_amr_restart_test}
      PROPERTIES
        WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_restart_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_restart_ranks}"
        DEPENDS regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_RESTART_TESTS
      ${selected_mpi_amr_restart_test}
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
      ${selected_mpi_amr_restart_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_restart_check
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_sparse_amr_restart_1d.py"
      --checkpoint
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_sparse_mpi_amr.chk"
      --fixed-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/fixed_reference.csv"
      --selected-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_reference_np1.csv"
      --selected-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_reference_np2.csv"
      --selected-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_reference_np4.csv"
      --restarted
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_restart_np2.csv"
      --restarted
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_restart_np4.csv"
      --stopped
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_checkpoint_stop.csv"
      --expected-bundle-sha256
        f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
      --expected-checkpoint-sha256
        0e5cc002ac5cceb3fa7732f8fe46382f140e038d19cac17664230ce8c500c2fd
      --final-time 1.0e-9
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_restart_check
    PROPERTIES DEPENDS "${SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS}"
  )

  foreach(
      selected_mpi_amr_mismatch
      IN ITEMS
        sha256 integrator composition species baseline density state_species
        near_floor_species state_closure temperature geometry
    )
    set(
      mismatch_work_directory
      "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/mismatch_${selected_mpi_amr_mismatch}"
    )
    file(MAKE_DIRECTORY "${mismatch_work_directory}")
    configure_file(
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/restart.nml"
      "${mismatch_work_directory}/restart.nml"
      COPYONLY
    )
    set(
      mutate_test_name
      regression_selected_mpi_amr_reactive_1d_mutate_${selected_mpi_amr_mismatch}
    )
    add_test(
      NAME ${mutate_test_name}
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/mutate_selected_mpi_amr_checkpoint_1d.py"
        --input
          "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_sparse_mpi_amr.chk"
        --output "${mismatch_work_directory}/selected_sparse_mpi_amr.chk"
        --field "${selected_mpi_amr_mismatch}"
    )
    set_tests_properties(
      ${mutate_test_name}
      PROPERTIES DEPENDS
        regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
    )

    set(
      reject_test_name
      regression_selected_mpi_amr_reactive_1d_reject_${selected_mpi_amr_mismatch}
    )
    add_test(
      NAME ${reject_test_name}
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
        --executable $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
        --input "${mismatch_work_directory}/restart.nml"
        --expected "Sparse MPI AMR restart read failed"
        --forbidden-output
          "${mismatch_work_directory}/selected_restarted.csv"
    )
    set_tests_properties(
      ${reject_test_name}
      PROPERTIES
        WORKING_DIRECTORY "${mismatch_work_directory}"
        DEPENDS ${mutate_test_name}
        TIMEOUT 60
      )
  endforeach()

  set(
    SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY
    "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/invalid"
  )
  file(MAKE_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}")
  configure_file(
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_path_alias.nml"
    "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
    COPYONLY
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_path_alias_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
      --executable $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
      --input invalid_path_alias.nml
      --expected "Selected sparse MPI AMR input and output paths alias"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_mpi_amr_alias.csv"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_path_alias_rejected
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      TIMEOUT 60
  )

  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_eb_patch_tree_2d"
  )
  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/invalid"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}"
  )

  set(SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS)
  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_fixed_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_fixed.nml"
      fixed_np1.csv
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_fixed_np1
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Chemistry level advances:"
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS
    regression_selected_mpi_reactive_eb_patch_tree_2d_fixed_np1
  )

  set(
    FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/fixed_restart"
  )
  file(
    MAKE_DIRECTORY
    "${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
  )
  add_test(
    NAME regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_restart_2d/checkpoint_stop.nml"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    PROPERTIES
      WORKING_DIRECTORY
        "${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
      TIMEOUT 300
  )
  set(FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      list(
        APPEND FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS
        --expected-sha256
          6b817c871fe03eec3c96601e77c38d2e14e0505b96f9a590f4044b539ce85b3a
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
      list(
        APPEND FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS
        --expected-sha256
          d3ed831ad70a488f7303bbc296b07f63f503c9f02c04a0ddff81d0bb8e7f1de1
      )
    endif()
  endif()
  add_test(
    NAME regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_check
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint.py"
      --checkpoint
        "${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/patch_tree_restart.chk"
      ${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS}
  )
  set_tests_properties(
    regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_check
    PROPERTIES DEPENDS
      regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
  )

  foreach(selected_mpi_eb_tree_ranks IN ITEMS 1 2 4)
    set(
      selected_mpi_eb_tree_test
      regression_selected_mpi_reactive_eb_patch_tree_2d_np${selected_mpi_eb_tree_ranks}
    )
    add_test(
      NAME ${selected_mpi_eb_tree_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_eb_tree_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
        "selected_np${selected_mpi_eb_tree_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_eb_tree_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_eb_tree_ranks}
        PASS_REGULAR_EXPRESSION "Chemistry integrator: implicit"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS
      ${selected_mpi_eb_tree_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_exact
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d.py"
      --reference
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/fixed_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_np2.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_np4.csv"
      --species H2 H O O2 OH H2O HO2 H2O2 AR N2
      --expected-levels 0 1 2 3
      --final-time 1.0e-12
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_exact
    PROPERTIES DEPENDS "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS}"
  )

  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_restart"
  )
  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/invalid"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY}"
  )
  set(SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS)

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_restart_fixed_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_fixed_reference.nml"
      fixed_reference_np1.csv
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_fixed_np1
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Completed root steps: 2"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_fixed_np1
  )

  foreach(selected_mpi_eb_restart_ranks IN ITEMS 1 2 4)
    set(
      selected_mpi_eb_restart_reference_test
      regression_selected_mpi_reactive_eb_patch_tree_2d_restart_selected_np${selected_mpi_eb_restart_ranks}
    )
    add_test(
      NAME ${selected_mpi_eb_restart_reference_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_eb_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_selected_reference.nml"
        "selected_reference_np${selected_mpi_eb_restart_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_eb_restart_reference_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_eb_restart_ranks}
        PASS_REGULAR_EXPRESSION "Completed root steps: 2"
        RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
      ${selected_mpi_eb_restart_reference_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_checkpoint_stop.nml"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
    regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
  )

  foreach(selected_mpi_eb_restart_ranks IN ITEMS 2 4)
    set(
      selected_mpi_eb_restart_test
      regression_selected_mpi_reactive_eb_patch_tree_2d_restart_np${selected_mpi_eb_restart_ranks}
    )
    add_test(
      NAME ${selected_mpi_eb_restart_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_eb_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml"
        "restarted_np${selected_mpi_eb_restart_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_eb_restart_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_eb_restart_ranks}
        PASS_REGULAR_EXPRESSION "Restarted:  T"
        DEPENDS
          regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
        RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
      ${selected_mpi_eb_restart_test}
    )
  endforeach()

  set(SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_HASH_ARGUMENTS)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
     CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_HASH_ARGUMENTS
      --expected-checkpoint-sha256
        2edd758ca1ff888735f5719de8ed9969b36f60d45bb0b3c523e47c2e35510f4d
      --expected-final-sha256
        5394fd007283d746ec4507832ccce8a7013fd4fb9b5c1b10d50d7c91e216975a
      --expected-stopped-sha256
        bba321c7ace389a6e0d95af0d1a07edd55c9cde870daee514e9b6c1c49e5fc3c
    )
  endif()
  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_restart_check
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d_restart.py"
      --checkpoint
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_mpi_eb_restart.chk"
      --reference
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/fixed_reference_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_reference_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_reference_np2.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_reference_np4.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/restarted_np2.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/restarted_np4.csv"
      --stopped
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_mpi_eb_restart_stopped.csv"
      --species H2 H O O2 OH H2O HO2 H2O2 AR N2
      --expected-levels 0 1 2 3
      --final-time 1.0e-12
      --bundle-sha256
        f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
      --integrator implicit
      --composition 0.29570 1.0e-5 1.0e-5 0.14784 1.0e-5 0.0 0.0 0.0 0.0
        0.55643
      ${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_HASH_ARGUMENTS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_check
    PROPERTIES
      DEPENDS "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS}"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
  )

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_restart_failures
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d_restart_failures.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 1
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml"
      --checkpoint
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_mpi_eb_restart.chk"
      --work-directory
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY}/corruptions"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_failures
    PROPERTIES
      DEPENDS
        regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 300
  )

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_aliases
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d_aliases.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 1
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --reference-input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_selected_reference.nml"
      --checkpoint-input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_checkpoint_stop.nml"
      --restart-input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml"
      --work-directory
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY}/aliases"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_aliases
    PROPERTIES
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 180
  )

  add_test(
    NAME regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_fixture_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/fixture_selected.nml"
      "fixture_selected_np1.csv"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_np1
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
      TIMEOUT 300
  )
  add_test(
    NAME regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_valid
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d.py"
      --reference
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/fixture_selected_np1.csv"
      --species H2 H
      --expected-levels 0 1 2 3
      --final-time 1.0e-12
  )
  set_tests_properties(
    regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_valid
    PROPERTIES DEPENDS
      regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_np1
  )

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_model_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_fixed.nml"
      --expected "Selected sparse MPI EB requires chemistry_model='selected'"
      --forbidden-output
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}/selected_mpi_eb_tree_fixed.csv"
  )
  add_test(
    NAME regression_fixed_mpi_reactive_eb_patch_tree_2d_selected_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
      --expected "Unknown reactive 2D chemistry model"
      --forbidden-output
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}/selected_mpi_eb_tree_selected.csv"
  )
  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_alias_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
      --expected "Selected sparse MPI EB input and output paths alias"
      --forbidden-output
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}/selected_mpi_eb_tree_selected.csv"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_model_rejected
    regression_fixed_mpi_reactive_eb_patch_tree_2d_selected_rejected
    regression_selected_mpi_reactive_eb_patch_tree_2d_alias_rejected
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      TIMEOUT 60
  )
endif()

add_executable(
  test_selected_reactive_1d_config
  unit/test_selected_reactive_1d_config.F90
)
target_link_libraries(
  test_selected_reactive_1d_config
  PRIVATE pelef_selected_reactive_1d_runtime
)
add_test(
  NAME unit_selected_reactive_1d_config
  COMMAND test_selected_reactive_1d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/fixture.nml"
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture.nml"
)

set(
  SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_fixture"
)
set(
  SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_full"
)
set(
  SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_invalid"
)
set(
  SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_1d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)

add_test(
  NAME regression_selected_reactive_1d_fixture_run
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/fixture.nml"
)
set_tests_properties(
  regression_selected_reactive_1d_fixture_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_1d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_1d.py"
    --input
      "${SELECTED_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_1d.csv"
    --species H2 H
    --nx 8
    --final-time 1.0e-7
    --require-uniform
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_reactive_1d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_1d_fixture_run
)

add_test(
  NAME regression_selected_reactive_1d_full_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_full_h2o2_1d/uniform.nml"
)
add_test(
  NAME regression_selected_reactive_1d_full_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/h2o2_full.nml"
)
set_tests_properties(
  regression_selected_reactive_1d_full_fixed_run
  regression_selected_reactive_1d_full_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_1d_full_selected_run
  PROPERTIES PASS_REGULAR_EXPRESSION "Species: 10"
)
add_test(
  NAME regression_selected_reactive_1d_full_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_1d.py"
    --input
      "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_1d.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8
    --final-time 2.0e-6
    --require-uniform
)
set_tests_properties(
  regression_selected_reactive_1d_full_structure
  PROPERTIES DEPENDS regression_selected_reactive_1d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_1d_full_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}/reactive_full_h2o2_1d.csv"
    "${SELECTED_REACTIVE_1D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_1d.csv"
)
set_tests_properties(
  regression_selected_reactive_1d_full_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_1d_full_fixed_run;regression_selected_reactive_1d_full_selected_run"
)

add_test(
  NAME regression_selected_reactive_1d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_1d.csv"
)
add_test(
  NAME regression_selected_reactive_1d_entropy_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/invalid_entropy_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_entropy_temperature.csv"
)
set_tests_properties(
  regression_selected_reactive_1d_unknown_rejected
  regression_selected_reactive_1d_entropy_temperature_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_1d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_1d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_1d/fixture.nml"
    --expected "Unknown reactive 1D chemistry model"
    --forbidden-output
      "${SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_1d.csv"
)
set_tests_properties(
  regression_fixed_reactive_1d_selected_model_rejected
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)

set(
  SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_fixture"
)
set(
  SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_two_level"
)
set(
  SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_multilevel"
)
set(
  SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_multipatch"
)
set(
  SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_boundary"
)
set(
  SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_invalid"
)
set(
  SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_amr_reactive_1d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_path_alias.nml"
  "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_amr_reactive_1d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_fixture_run_a
  regression_selected_amr_reactive_1d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Fine level active:  T"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_amr_reactive_1d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d.csv"
    --species H2 H
    --expected-levels 0 1
    --x-upper 0.008
    --final-time 1.0e-7
    --refinement-ratio 2
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-8
)
set_tests_properties(
  regression_selected_amr_reactive_1d_fixture_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_fixture_run_a
)
add_test(
  NAME regression_selected_amr_reactive_1d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d.csv"
    "${SELECTED_AMR_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d_repeat.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_fixture_run_a;regression_selected_amr_reactive_1d_fixture_run_b"
)

add_test(
  NAME regression_selected_amr_reactive_1d_two_level_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_two_level_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_two_level_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_two_level_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_two_level_fixed_run
  regression_selected_amr_reactive_1d_two_level_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Fine level active:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_two_level_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_two_level.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1
    --x-upper 0.008
    --final-time 1.0e-9
    --refinement-ratio 2
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-12
)
set_tests_properties(
  regression_selected_amr_reactive_1d_two_level_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_two_level_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_two_level_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}/full_h2o2_amr_two_level.csv"
    "${SELECTED_AMR_REACTIVE_1D_TWO_LEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_two_level.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_two_level_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_two_level_fixed_run;regression_selected_amr_reactive_1d_two_level_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multilevel_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multilevel_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multilevel_fixed_run
  regression_selected_amr_reactive_1d_multilevel_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Active AMR levels: 3"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_multilevel.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1 2
    --x-upper 0.012
    --final-time 2.0e-7
    --refinement-ratio 2
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multilevel_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_multilevel_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_multilevel_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}/full_h2o2_amr_multilevel.csv"
    "${SELECTED_AMR_REACTIVE_1D_MULTILEVEL_WORK_DIRECTORY}/selected_full_h2o2_amr_multilevel.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multilevel_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_multilevel_fixed_run;regression_selected_amr_reactive_1d_multilevel_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multipatch_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_multipatch_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multipatch_fixed_run
  regression_selected_amr_reactive_1d_multipatch_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Active fine patches: 3"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}/selected_full_h2o2_amr_multipatch.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1
    --x-upper 0.012
    --final-time 1.0e-10
    --refinement-ratio 2
    --minimum-finest-regions 3
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multipatch_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_multipatch_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_multipatch_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}/full_h2o2_amr_multipatch.csv"
    "${SELECTED_AMR_REACTIVE_1D_MULTIPATCH_WORK_DIRECTORY}/selected_full_h2o2_amr_multipatch.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_multipatch_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_multipatch_fixed_run;regression_selected_amr_reactive_1d_multipatch_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_boundary_fixed_run
  COMMAND $<TARGET_FILE:pelef_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_boundary_weno_fixed.nml"
)
add_test(
  NAME regression_selected_amr_reactive_1d_boundary_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/h2o2_full_boundary_weno_selected.nml"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_boundary_fixed_run
  regression_selected_amr_reactive_1d_boundary_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "AMR hybrid WENO:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_amr_reactive_1d_boundary_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_1d.py"
    --input
      "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}/selected_full_h2o2_amr_boundary_weno.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --expected-levels 0 1 2
    --x-upper 0.012
    --final-time 1.0e-7
    --refinement-ratio 2
    --require-finest-touches-lower
)
set_tests_properties(
  regression_selected_amr_reactive_1d_boundary_check
  PROPERTIES DEPENDS regression_selected_amr_reactive_1d_boundary_selected_run
)
add_test(
  NAME regression_selected_amr_reactive_1d_boundary_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}/full_h2o2_amr_boundary_weno.csv"
    "${SELECTED_AMR_REACTIVE_1D_BOUNDARY_WORK_DIRECTORY}/selected_full_h2o2_amr_boundary_weno.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_boundary_exact
  PROPERTIES DEPENDS
    "regression_selected_amr_reactive_1d_boundary_fixed_run;regression_selected_amr_reactive_1d_boundary_selected_run"
)

add_test(
  NAME regression_selected_amr_reactive_1d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_reactive_1d.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_non_amr_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_non_amr.nml"
    --expected "requires amr_enabled"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_non_amr.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_entropy_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_entropy_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_entropy_temperature.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_checkpoint_interval_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_checkpoint_interval.nml"
    --expected "does not implement checkpoint or restart"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_interval.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_checkpoint_file_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_checkpoint_file.nml"
    --expected "does not implement checkpoint or restart"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_file.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_restart_file_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_restart_file.nml"
    --expected "does not implement checkpoint or restart"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart.csv"
)
add_test(
  NAME regression_selected_amr_reactive_1d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    --input invalid_path_alias.nml
    --expected "input and output paths must differ"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_alias.csv"
)
set_tests_properties(
  regression_selected_amr_reactive_1d_unknown_rejected
  regression_selected_amr_reactive_1d_non_amr_rejected
  regression_selected_amr_reactive_1d_entropy_temperature_rejected
  regression_selected_amr_reactive_1d_checkpoint_interval_rejected
  regression_selected_amr_reactive_1d_checkpoint_file_rejected
  regression_selected_amr_reactive_1d_restart_file_rejected
  regression_selected_amr_reactive_1d_path_alias_rejected
  PROPERTIES WORKING_DIRECTORY "${SELECTED_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_amr_reactive_1d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_amr_reactive_1d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/fixture.nml"
    --expected "Unknown reactive 1D chemistry model"
    --forbidden-output
      "${SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_amr_reactive_1d.csv"
)
set_tests_properties(
  regression_fixed_amr_reactive_1d_selected_model_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_AMR_REACTIVE_1D_FIXED_REJECTION_WORK_DIRECTORY}"
)

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

add_executable(
  test_selected_reactive_eb_amr_2d_config
  unit/test_selected_reactive_eb_amr_2d_config.F90
)
target_link_libraries(
  test_selected_reactive_eb_amr_2d_config
  PRIVATE pelef_selected_reactive_eb_amr_2d_runtime
)
add_test(
  NAME unit_selected_reactive_eb_amr_2d_config
  COMMAND test_selected_reactive_eb_amr_2d_config
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture.nml"
)

set(
  SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_fixture"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_full"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_invalid"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_dynamic_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_three_level_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_dynamic_three_level_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_multipatch_restart"
)
set(
  SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_amr_2d_fixed_rejection"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_path_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_fine_input_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_fine_input_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_checkpoint_input_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_input_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_finest_input_alias.nml"
  "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_finest_input_alias.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_run_a
  regression_selected_reactive_eb_amr_2d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Coarse cut cells: 8"
    TIMEOUT 180
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_coarse_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    --species H2 H
    --nx 8
    --ny 8
    --final-time 1.0e-7
    --expected-regular 32
    --expected-cut 8
    --expected-covered 24
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_fine_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
    --species H2 H
    --nx 12
    --ny 12
    --final-time 1.0e-7
    --expected-regular 84
    --expected-cut 12
    --expected-covered 48
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_coarse_check
  regression_selected_reactive_eb_amr_2d_fixture_fine_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_hierarchy_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_2d.py"
    --coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    --fine
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
    --species H2 H
    --coarse-nx 8
    --coarse-ny 8
    --x-lower 0.0
    --x-upper 0.008
    --y-lower 0.0
    --y-upper 0.008
    --coarse-i-lower 2
    --coarse-i-upper 7
    --coarse-j-lower 2
    --coarse-j-upper 7
    --refinement-ratio 2
    --final-time 1.0e-7
    --coarse-counts 32 8 24
    --fine-counts 84 12 48
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-5
    --expected-coarse-sha256
      40758bed925456948aa694ee130e766b8b942a2bdfef104c268323a981119c0e
    --expected-fine-sha256
      30728d746ccdbbec2b238e4c738f6e0e54ba9550ec5d9ea609c5d7972551fbfe
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_hierarchy_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_coarse_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse_repeat.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fixture_fine_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_fixture_coarse_exact
  regression_selected_reactive_eb_amr_2d_fixture_fine_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_amr_2d_fixture_run_a;regression_selected_reactive_eb_amr_2d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_fixed_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/h2o2_full_fixed.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_selected_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/h2o2_full_selected.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_fixed_run
  regression_selected_reactive_eb_amr_2d_full_selected_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry:  T"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_coarse_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_coarse.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 8
    --ny 8
    --final-time 2.0e-8
    --expected-regular 32
    --expected-cut 8
    --expected-covered 24
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_fine_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_2d.py"
    --input
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_fine.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 12
    --ny 12
    --final-time 2.0e-8
    --expected-regular 84
    --expected-cut 12
    --expected-covered 48
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_coarse_check
  regression_selected_reactive_eb_amr_2d_full_fine_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_hierarchy_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_2d.py"
    --coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_coarse.csv"
    --fine
      "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_fine.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --coarse-nx 8
    --coarse-ny 8
    --x-lower 0.0
    --x-upper 0.008
    --y-lower 0.0
    --y-upper 0.008
    --coarse-i-lower 2
    --coarse-i-upper 7
    --coarse-j-lower 2
    --coarse-j-upper 7
    --refinement-ratio 2
    --final-time 2.0e-8
    --coarse-counts 32 8 24
    --fine-counts 84 12 48
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-10
    --expected-coarse-sha256
      $<IF:$<CONFIG:Debug>,bee1bd8baae3606a9ac3483822140c2d0d29324683c9b9deab62c84fe4819d7b,f6e6fe44a629a92c3027f52378a0f892ab8e790a88bdd99dc8fb178aa18605a7>
    --expected-fine-sha256
      $<IF:$<CONFIG:Debug>,be55e526b0cf07fbfff63cd6f874e104599f0d0a019c02ca1145f614acc4079e,a61021c2a77cfd933e4d1af1027ee08b869a4bed4791534080d9aa1cfc069423>
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_hierarchy_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_amr_2d_full_selected_run
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_coarse_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/full_h2o2_reactive_eb_amr_2d_coarse.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_coarse.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_full_fine_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/full_h2o2_reactive_eb_amr_2d_fine.csv"
    "${SELECTED_REACTIVE_EB_AMR_2D_FULL_WORK_DIRECTORY}/selected_full_h2o2_reactive_eb_amr_2d_fine.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_full_coarse_exact
  regression_selected_reactive_eb_amr_2d_full_fine_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_amr_2d_full_fixed_run;regression_selected_reactive_eb_amr_2d_full_selected_run"
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_restart_reference
  regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Completed coarse steps: 2"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_mismatch_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_mismatch_fine.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}"
    DEPENDS regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_2d_restart.py"
    --fixed-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_restart_reference_coarse.csv"
    --fixed-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_restart_reference_fine.csv"
    --reference-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_reference_coarse.csv"
    --reference-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_reference_fine.csv"
    --stopped-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_stopped_coarse.csv"
    --stopped-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart_stopped_fine.csv"
    --restarted-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restarted_coarse.csv"
    --restarted-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restarted_fine.csv"
    --expected-reference-coarse-sha256
      5d8e24295aa1331dfdd7a6224fd0809a5df712e2db22c733ac44bed5ffb3f7ac
    --expected-reference-fine-sha256
      cbfef95067e8674239122633711c8e35371ea97e3709f80beeadccbabbc64af2
    --expected-stopped-coarse-sha256
      dda501bd9c27cb090f8bc8f7d06ec55f27e02fdd59e667cff8f89b2edb52c736
    --expected-stopped-fine-sha256
      6f30dd051f9fa46180a539ca62b7dc0f7c237850b9fe1ce0bdab557258bec817
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_restart.chk"
    --expected-checkpoint-sha256
      6d1a5e35cc001f3dfb7049cba243482bf902acde8da3dd487d4ba93598d784ba
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_restart_reference;regression_selected_reactive_eb_amr_2d_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_restart_restart;regression_selected_reactive_eb_amr_2d_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_restart_reference
  regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_dynamic_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Completed coarse steps: 3"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_mismatch_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_mismatch_fine.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_dynamic_restart_2d.py"
    --fixed-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_restart_reference_coarse.csv"
    --fixed-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_restart_reference_fine.csv"
    --reference-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_reference_coarse.csv"
    --reference-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_reference_fine.csv"
    --stopped-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_stopped_coarse.csv"
    --stopped-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart_stopped_fine.csv"
    --restarted-coarse
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restarted_coarse.csv"
    --restarted-fine
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restarted_fine.csv"
    --expected-reference-coarse-sha256
      2a67c9b9a618fc36412c0156adec150e3271a5b1ef9ab5ec84e158a8d07c1a87
    --expected-reference-fine-sha256
      66f479a71de2474ae2cad3f126b847238311d5e9a1ae40d0b60d25a1a9b776e9
    --expected-stopped-coarse-sha256
      8c0b0c74beb20ad5915e03578fa16c336b235b1d10d2e3e96b183d261052ba72
    --expected-stopped-fine-sha256
      a0e71a149a0a0991732024bef41f46cdd0b3318946390141ae0a2071f5c5acb2
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_restart.chk"
    --expected-checkpoint-sha256
      a91050d485886e3109337331f0af2ce0df8f8b8aa99227a2aaa0c69feeb2500d
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_dynamic_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_dynamic_restart_reference;regression_selected_reactive_eb_amr_2d_dynamic_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_dynamic_restart_restart;regression_selected_reactive_eb_amr_2d_dynamic_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_three_level_restart_reference
  regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_three_level_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_three_level_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_three_level_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Completed coarse steps:"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_three_level_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_mismatch_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_mismatch_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_mismatch_finest.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_three_level_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_three_level_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_three_level_restart_2d.py"
    --fixed-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_three_level_reference_root.csv"
    --fixed-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_three_level_reference_middle.csv"
    --fixed-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_three_level_reference_finest.csv"
    --reference-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_reference_root.csv"
    --reference-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_reference_middle.csv"
    --reference-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_reference_finest.csv"
    --stopped-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_stopped_root.csv"
    --stopped-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_stopped_middle.csv"
    --stopped-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_stopped_finest.csv"
    --restarted-root
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_restarted_root.csv"
    --restarted-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_restarted_middle.csv"
    --restarted-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level_restarted_finest.csv"
    --expected-reference-root-sha256 92f9580d35032710b09653e6cfd80d97c133c6df8fd70af0db4f5011714ec980
    --expected-reference-middle-sha256 59e66071d6aaf3fe1fd62de9db72d28163acbb7820c95e9a93f8618ab013c1b9
    --expected-reference-finest-sha256 70062b794483b16de5e7de66a3aeaebfbd914dd9ac1adb7bf3515a9b7a9f4250
    --expected-stopped-root-sha256 53ca2ae393f4586f8632ad5d38fc7d902fc85f6b57c31f8dbbebd778065c2fc4
    --expected-stopped-middle-sha256 38ff58c0bd7c83a87c2cf7752f770c5864d04cb3ed6234ba6fe815424372af11
    --expected-stopped-finest-sha256 46abc4426db3ba608d3bf65f3a5937ef4427f1bc588cfe83f97d117e7a57e3d5
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_three_level.chk"
    --expected-checkpoint-sha256 e8655df9aeb62ece3bcf3e23ee5ed1bcda17a2d834cf214e3ecaaf71f3e952b6
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_three_level_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_three_level_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_three_level_restart_reference;regression_selected_reactive_eb_amr_2d_three_level_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_three_level_restart_restart;regression_selected_reactive_eb_amr_2d_three_level_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_three_level_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_three_level_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference
  PROPERTIES
    PASS_REGULAR_EXPRESSION
      "Completed coarse steps:;Completed regrids: 1"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_dynamic_three_level_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_mismatch_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_mismatch_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_mismatch_finest.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_three_level_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_dynamic_three_level_restart_2d.py"
    --fixed-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_three_level_reference_root.csv"
    --fixed-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_three_level_reference_middle.csv"
    --fixed-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_dynamic_three_level_reference_finest.csv"
    --reference-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_reference_root.csv"
    --reference-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_reference_middle.csv"
    --reference-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_reference_finest.csv"
    --stopped-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_stopped_root.csv"
    --stopped-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_stopped_middle.csv"
    --stopped-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_stopped_finest.csv"
    --restarted-root
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_restarted_root.csv"
    --restarted-middle
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_restarted_middle.csv"
    --restarted-finest
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level_restarted_finest.csv"
    --expected-reference-root-sha256 91e8ecebce463172ed8c04a1a65589845c9e49c9ae4c1c67e7ae2627d1252ebf
    --expected-reference-middle-sha256 c752ebcd7308a888395b9206bfc944585d88c011bc86d10f7450a8286f5a237e
    --expected-reference-finest-sha256 ceb72783d7a0f531667227862eba3d877af464d46585e86e75cf377c42bef5dd
    --expected-stopped-root-sha256 84dc02136afcef4cb4227d86fb36e8ec7902fee61fd09fbe7a018a55e68bc615
    --expected-stopped-middle-sha256 81b59680ea1347ccded15dc588f763420370c48557f3f83ffa73c01197a7f5b7
    --expected-stopped-finest-sha256 f9d08b869ff4e5cf4b23c0f3193b3a0f2420973e80d0c9ec0b5cb9f56fd78645
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_dynamic_three_level.chk"
    --expected-checkpoint-sha256 1e8e24e2ead15ac71199cb73146656ffb63be089bafc57eced95fbf419ec1502
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_DYNAMIC_THREE_LEVEL_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_reference;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_restart;regression_selected_reactive_eb_amr_2d_dynamic_three_level_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_dynamic_three_level_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_reference.log"
    -- $<TARGET_FILE:pelef_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/fixed_reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_reference
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_restart
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_multipatch_restart_reference
  regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
  regression_selected_reactive_eb_amr_2d_multipatch_restart_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_multipatch_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference
  regression_selected_reactive_eb_amr_2d_multipatch_restart_reference
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Fine patches: 2"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart source:"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_full_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d_multipatch_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_mismatch_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_mismatch_fine_patch0001.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_mismatch_fine_patch0002.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}"
    DEPENDS
      regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_amr_2d_multipatch_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_amr_multipatch_restart_2d.py"
    --fixed-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_multipatch_reference_root.csv"
    --fixed-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_multipatch_reference_fine_patch0001.csv"
    --fixed-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/fixed_selected_reactive_eb_amr_multipatch_reference_fine_patch0002.csv"
    --reference-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_reference_root.csv"
    --reference-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_reference_fine_patch0001.csv"
    --reference-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_reference_fine_patch0002.csv"
    --stopped-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_stopped_root.csv"
    --stopped-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_stopped_fine_patch0001.csv"
    --stopped-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_stopped_fine_patch0002.csv"
    --restarted-root
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restarted_root.csv"
    --restarted-patch1
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restarted_fine_patch0001.csv"
    --restarted-patch2
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restarted_fine_patch0002.csv"
    --expected-reference-root-sha256
      1a2c70bf88bb36466937e0839294ca6363dad34752742e2322a9e49a90507f80
    --expected-reference-patch1-sha256
      e6a2d52e082a9b17ae00325cb5a65c2ee0ee8425646991b0b77e5b17cef42d88
    --expected-reference-patch2-sha256
      ceb4492ca19efbee67ed6e433a787eb48c887193544e1f88acb364d7b7f41b44
    --expected-stopped-root-sha256
      81d4a5ed93947b8790fa19c1ed6df716f8a5405bc84d0da974028d1acb45fdf1
    --expected-stopped-patch1-sha256
      9190c1137c5e41499d421f5fcd97b95d32466eb54a4bbd567a36e21efbe26ad5
    --expected-stopped-patch2-sha256
      2c1833dbd96d1aa3e202ebe87c62085e9d8a8939f77c6e2bc91955b872cda644
    --checkpoint
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/selected_reactive_eb_amr_multipatch_restart.chk"
    --expected-checkpoint-sha256
      72fdd25e457e27ece30739ca52d59536cb3ed247439c14bbd0777e534ad2b368
    --reference-log
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_AMR_2D_MULTIPATCH_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_multipatch_restart_check
  PROPERTIES
    DEPENDS
      "regression_selected_reactive_eb_amr_2d_multipatch_restart_fixed_reference;regression_selected_reactive_eb_amr_2d_multipatch_restart_reference;regression_selected_reactive_eb_amr_2d_multipatch_restart_checkpoint_stop;regression_selected_reactive_eb_amr_2d_multipatch_restart_restart;regression_selected_reactive_eb_amr_2d_multipatch_restart_composition_mismatch_rejected"
    RESOURCE_LOCK selected_reactive_eb_amr_2d_multipatch_restart
    TIMEOUT 60
)

add_test(
  NAME regression_selected_reactive_eb_amr_2d_unknown_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_unknown_species.nml"
    --expected "not found exactly once"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_2d.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_2d_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_temperature_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_embedded_temperature.nml"
    --expected "temperature is outside the mechanism range"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_temperature.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_temperature_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_path_alias.nml
    --expected "paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_alias_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_checkpoint_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_checkpoint_input_alias.nml
    --expected "checkpoint path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_alias.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_alias_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_output_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_restart_output_alias.nml"
    --expected "restart path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_alias.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_alias_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_checkpoint_restart_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_checkpoint_restart_alias.nml"
    --expected "checkpoint and restart paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_restart.chk"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_restart_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_restart_fine.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_three_level_removal_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_three_level_mode.nml"
    --expected "Dynamic three-level EB AMR keeps the finest patch active"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_three_level.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_three_level_fine.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_three_level_finest.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_multipatch_derived_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_multipatch_derived_alias.nml"
    --expected "multipatch derived paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_multipatch.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_multipatch_fine_patch0001.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_dynamic_parent_without_dynamic_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_dynamic_parent_mode.nml"
    --expected "Dynamic parent regridding requires dynamic three-level EB AMR"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_dynamic_parent.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_dynamic_parent_fine.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_dynamic_parent_finest.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_fine_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_fine_input_alias.nml
    --expected "paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_fine_alias_coarse.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_coarse_fine_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_coarse_fine_alias.nml"
    --expected "fine output must be nonempty and distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_level_alias.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_finest_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input invalid_finest_input_alias.nml
    --expected "three-level paths must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_finest_input_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_finest_input_middle.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_checkpoint_finest_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_checkpoint_finest_alias.nml"
    --expected "checkpoint path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_finest_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_finest_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_checkpoint_finest.csv"
)
add_test(
  NAME regression_selected_reactive_eb_amr_2d_restart_finest_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/invalid_restart_finest_alias.nml"
    --expected "restart path must be distinct"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_finest_root.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_finest_middle.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_amr_restart_finest.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_amr_2d_unknown_rejected
  regression_selected_reactive_eb_amr_2d_temperature_rejected
  regression_selected_reactive_eb_amr_2d_path_alias_rejected
  regression_selected_reactive_eb_amr_2d_checkpoint_input_alias_rejected
  regression_selected_reactive_eb_amr_2d_restart_output_alias_rejected
  regression_selected_reactive_eb_amr_2d_checkpoint_restart_alias_rejected
  regression_selected_reactive_eb_amr_2d_dynamic_three_level_removal_rejected
  regression_selected_reactive_eb_amr_2d_multipatch_derived_alias_rejected
  regression_selected_reactive_eb_amr_2d_dynamic_parent_without_dynamic_rejected
  regression_selected_reactive_eb_amr_2d_fine_input_alias_rejected
  regression_selected_reactive_eb_amr_2d_coarse_fine_alias_rejected
  regression_selected_reactive_eb_amr_2d_finest_input_alias_rejected
  regression_selected_reactive_eb_amr_2d_checkpoint_finest_alias_rejected
  regression_selected_reactive_eb_amr_2d_restart_finest_alias_rejected
  PROPERTIES
    RESOURCE_LOCK selected_reactive_eb_amr_2d_invalid
    WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_AMR_2D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_eb_amr_2d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_amr_2d/fixture.nml"
    --expected "Unknown reactive 2D chemistry model"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_coarse.csv"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_amr_2d_fine.csv"
)
set_tests_properties(
  regression_fixed_reactive_eb_amr_2d_selected_model_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_AMR_2D_FIXED_REJECTION_WORK_DIRECTORY}"
)

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

if(PELEF_ENABLE_MPI)
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_fixture"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_full"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_transport"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_restart"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_invalid"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
  )
  configure_file(
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_path_alias.nml"
    "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
    COPYONLY
  )
  configure_file(
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_checkpoint_input_alias.nml"
    "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_input_alias.nml"
    COPYONLY
  )

  set(SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_RUN_TESTS)
  foreach(selected_mpi_amr_3d_ranks IN ITEMS 1 2 4)
    set(
      fixture_test_name
      regression_selected_mpi_amr_reactive_3d_fixture_np${selected_mpi_amr_3d_ranks}
    )
    add_test(
      NAME ${fixture_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/fixture.nml"
        "selected_fixture_mpi_np${selected_mpi_amr_3d_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${fixture_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_3d_ranks}"
        TIMEOUT 180
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_RUN_TESTS
      ${fixture_test_name}
    )
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_fixture_np4
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_fixture_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
      --coarse
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np4_coarse.csv"
      --fine
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np4_fine.csv"
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
    regression_selected_mpi_amr_reactive_3d_fixture_check
    PROPERTIES DEPENDS regression_selected_mpi_amr_reactive_3d_fixture_np4
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_fixture_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np4"
      --final-time 1.0e-7
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_fixture_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_fixture_run_a;${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_RUN_TESTS}"
  )

  set(SELECTED_MPI_AMR_REACTIVE_3D_FULL_RUN_TESTS)
  foreach(selected_mpi_amr_3d_ranks IN ITEMS 1 2 4 8)
    set(
      full_test_name
      regression_selected_mpi_amr_reactive_3d_full_np${selected_mpi_amr_3d_ranks}
    )
    add_test(
      NAME ${full_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_chemistry.nml"
        "selected_full_mpi_np${selected_mpi_amr_3d_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${full_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_3d_ranks}"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_3D_FULL_RUN_TESTS ${full_test_name}
    )
  endforeach()
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_full_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
      --coarse
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np8_coarse.csv"
      --fine
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np8_fine.csv"
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
    regression_selected_mpi_amr_reactive_3d_full_check
    PROPERTIES DEPENDS regression_selected_mpi_amr_reactive_3d_full_np8
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_full_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_amr_hotspot_chemistry_3d"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np8"
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --expected-coarse-sha256
        3313a9a6eff5d6175625e4bbe3613a782cdd794314fd9fbdd53ddaa309fc9516
      --expected-fine-sha256
        4ff0785c3462a49648fd694d090eb23bffd5af05a43a33360cd187417bd8a31d
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_full_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_full_selected_run;${SELECTED_MPI_AMR_REACTIVE_3D_FULL_RUN_TESTS}"
  )

  set(SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RUN_TESTS)
  foreach(selected_mpi_amr_3d_transport_ranks IN ITEMS 1 2 4 8)
    set(
      transport_test_name
      regression_selected_mpi_amr_reactive_3d_transport_np${selected_mpi_amr_3d_transport_ranks}
    )
    add_test(
      NAME ${transport_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_transport_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport.nml"
        "selected_transport_mpi_np${selected_mpi_amr_3d_transport_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${transport_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_transport_ranks}
        PASS_REGULAR_EXPRESSION
          "Completed transport fine substeps: 16"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RUN_TESTS
      ${transport_test_name}
    )
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
      --coarse
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np8_coarse.csv"
      --fine
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np8_fine.csv"
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
    regression_selected_mpi_amr_reactive_3d_transport_check
    PROPERTIES DEPENDS regression_selected_mpi_amr_reactive_3d_transport_np8
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_amr_hotspot_transport_3d"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np8"
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --expected-coarse-sha256
        $<IF:$<CONFIG:Debug>,5687cb97f834b759b3932bbd3f0b6e30fc5b48e6029e3f47cadfb3ec6a8d74fb,f0713f3a68f48b6537afaafd57209ea8ed4767f9c7d093fe2961bd3a1fa37fae>
      --expected-fine-sha256
        eaf132bbf89974d739c26a35f6b76137f9774438cfe9eebfc12436aa364cff58
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_transport_selected_run;${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RUN_TESTS}"
  )

  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_checkpoint_np2
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart_checkpoint.nml"
      selected_full_mpi_checkpoint_np2
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_checkpoint_np2
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
      TIMEOUT 300
  )
  set(SELECTED_MPI_AMR_REACTIVE_3D_RESTART_TESTS)
  foreach(selected_mpi_amr_3d_restart_ranks IN ITEMS 1 2 4 8)
    set(
      restart_test_name
      regression_selected_mpi_amr_reactive_3d_restart_np${selected_mpi_amr_3d_restart_ranks}
    )
    add_test(
      NAME ${restart_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart.nml"
        "selected_full_mpi_restart_np${selected_mpi_amr_3d_restart_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${restart_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_restart_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_3d_restart_ranks}"
        DEPENDS regression_selected_mpi_amr_reactive_3d_checkpoint_np2
        TIMEOUT 300
    )
    list(APPEND SELECTED_MPI_AMR_REACTIVE_3D_RESTART_TESTS ${restart_test_name})
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_restart_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_restart_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart_reference"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np8"
      --checkpoint
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
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
    regression_selected_mpi_amr_reactive_3d_restart_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_restart_reference;regression_selected_mpi_amr_reactive_3d_checkpoint_np2;${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_TESTS}"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_checkpoint_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
      "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_checkpoint_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_checkpoint_stop;regression_selected_mpi_amr_reactive_3d_checkpoint_np2"
  )

  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_restart_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_restart.nml"
      --expected "Could not open 3D AMR checkpoint"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart_coarse.csv"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart_fine.csv"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_path_alias_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      --argument invalid_path_alias.nml
      --expected "input and output paths alias"
      --preserved-input invalid_path_alias.nml
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_alias_fine.csv"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_checkpoint_alias_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      --argument invalid_checkpoint_input_alias.nml
      --expected "checkpoint path aliases another file"
      --preserved-input invalid_checkpoint_input_alias.nml
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_alias_coarse.csv"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_alias_fine.csv"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_restart_rejected
    regression_selected_mpi_amr_reactive_3d_path_alias_rejected
    regression_selected_mpi_amr_reactive_3d_checkpoint_alias_rejected
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
      PROCESSORS 2
      TIMEOUT 60
  )
endif()

set(
  SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_fixture"
)
set(
  SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_full"
)
set(
  SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_elementary"
)
set(
  SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_invalid"
)
set(
  SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_fixed_rejection"
)
set(
  SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_reactive_eb_3d_restart"
)
file(
  MAKE_DIRECTORY
  "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY}"
  "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_path_alias.nml"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_checkpoint_interval.nml"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_interval.nml"
  COPYONLY
)
configure_file(
  "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_restart.nml"
  "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_restart.nml"
  COPYONLY
)

add_test(
  NAME regression_selected_reactive_eb_3d_fixture_run_a
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture.nml"
)
add_test(
  NAME regression_selected_reactive_eb_3d_fixture_run_b
  COMMAND $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture_repeat.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_fixture_run_a
  regression_selected_reactive_eb_3d_fixture_run_b
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 2"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_3d_fixture_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_3d.py"
    --input
      "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d.csv"
    --species H2 H
    --nx 4 --ny 4 --nz 4
    --final-time 1.0e-7
    --x-upper 0.004 --y-upper 0.004 --z-upper 0.004
    --plane-axis x --plane-position 0.00163
    --expected-regular 32 --expected-cut 16 --expected-covered 16
    --activity-species H2
    --initial-mass-fraction 0.8888888888888889
    --minimum-change 1.0e-6
)
set_tests_properties(
  regression_selected_reactive_eb_3d_fixture_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_3d_fixture_run_a
)
add_test(
  NAME regression_selected_reactive_eb_3d_fixture_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d.csv"
    "${SELECTED_REACTIVE_EB_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d_repeat.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_fixture_exact
  PROPERTIES DEPENDS
    "regression_selected_reactive_eb_3d_fixture_run_a;regression_selected_reactive_eb_3d_fixture_run_b"
)

add_test(
  NAME regression_selected_reactive_eb_3d_full_run
  COMMAND $<TARGET_FILE:test_selected_full_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/h2o2_full_transport_chemistry.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_full_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Operator sequence: R-T-H-T-R"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_3d_full_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_3d.py"
    --input
      "${SELECTED_REACTIVE_EB_3D_FULL_WORK_DIRECTORY}/selected_full_reactive_eb_3d_transport_chemistry.csv"
    --species H2 H O O2 OH H2O HO2 H2O2 AR N2
    --nx 10 --ny 8 --nz 6
    --final-time 2.0e-9
    --x-upper 1.0e-4 --y-upper 1.0e-4 --z-upper 1.0e-4
    --plane-axis x --plane-position 3.95e-5
    --expected-regular 288 --expected-cut 48 --expected-covered 144
    --activity-species H2
    --initial-mass-fraction 0.028502723329253223
    --minimum-change 1.0e-9
    --minimum-active-density-span 2.0e-2
)
set_tests_properties(
  regression_selected_reactive_eb_3d_full_check
  PROPERTIES DEPENDS regression_selected_reactive_eb_3d_full_run
)

add_test(
  NAME regression_selected_reactive_eb_3d_elementary_run
  COMMAND $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/h2o2_elementary_transport_chemistry.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_elementary_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_reactive_eb_3d_elementary_exact
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${REACTIVE_EB_TRANSPORT_3D_WORK_DIRECTORY}/reactive_eb_3d_transport_chemistry.csv"
    "${SELECTED_REACTIVE_EB_3D_ELEMENTARY_WORK_DIRECTORY}/selected_elementary_reactive_eb_3d.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_elementary_exact
  PROPERTIES DEPENDS
    "regression_reactive_eb_transport_3d_coupled_run;regression_selected_reactive_eb_3d_elementary_run"
)

add_test(
  NAME regression_selected_reactive_eb_3d_restart_reference
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/reference.log"
    -- $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/reference.nml"
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_checkpoint_stop
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/stopped.log"
    -- $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/checkpoint_stop.nml"
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_restart
  COMMAND
    "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
    --output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/restarted.log"
    -- $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/restart.nml"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_reference
  regression_selected_reactive_eb_3d_restart_checkpoint_stop
  regression_selected_reactive_eb_3d_restart_restart
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}"
    FAIL_REGULAR_EXPRESSION
      "ERROR STOP;Program received signal;Floating-point exception"
    RESOURCE_LOCK selected_reactive_eb_3d_restart
    TIMEOUT 300
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_reference
  PROPERTIES PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_checkpoint_stop
  PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_restart
  PROPERTIES
    PASS_REGULAR_EXPRESSION "Restart continuation: complete"
    DEPENDS regression_selected_reactive_eb_3d_restart_checkpoint_stop
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_composition_mismatch_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d_restart/mismatched_composition.nml"
    --expected "selected checkpoint composition mismatch"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_mismatch.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_composition_mismatch_rejected
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}"
    DEPENDS regression_selected_reactive_eb_3d_restart_checkpoint_stop
    RESOURCE_LOCK selected_reactive_eb_3d_restart
    TIMEOUT 60
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_selected_reactive_eb_3d_restart.py"
    --fixed-reference
      "${REACTIVE_EB_RESTART_3D_WORK_DIRECTORY}/reactive_eb_3d_restart_reference.csv"
    --reference
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_reference.csv"
    --stopped
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_stopped.csv"
    --restarted
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart_restarted.csv"
    --checkpoint
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/selected_reactive_eb_3d_restart.chk"
    --reference-log
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/reference.log"
    --stopped-log
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/stopped.log"
    --restarted-log
      "${SELECTED_REACTIVE_EB_3D_RESTART_WORK_DIRECTORY}/restarted.log"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_restart_check
  PROPERTIES
    DEPENDS
      "regression_reactive_eb_restart_3d_reference;regression_selected_reactive_eb_3d_restart_reference;regression_selected_reactive_eb_3d_restart_checkpoint_stop;regression_selected_reactive_eb_3d_restart_restart"
    RESOURCE_LOCK selected_reactive_eb_3d_restart
    TIMEOUT 60
)

foreach(
    selected_reactive_eb_3d_invalid_case IN ITEMS
    unknown_species
    duplicate_species
    nonfinite_composition
    zero_composition
    negative_composition
    temperature
  )
  if(selected_reactive_eb_3d_invalid_case STREQUAL "unknown_species")
    set(selected_reactive_eb_3d_expected "not found exactly once")
    set(selected_reactive_eb_3d_output "unknown")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "duplicate_species")
    set(selected_reactive_eb_3d_expected "duplicate species")
    set(selected_reactive_eb_3d_output "duplicate")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "nonfinite_composition")
    set(selected_reactive_eb_3d_expected "nonfinite mole fraction")
    set(selected_reactive_eb_3d_output "nonfinite")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "zero_composition")
    set(selected_reactive_eb_3d_expected "nonnegative and nonempty")
    set(selected_reactive_eb_3d_output "zero")
  elseif(selected_reactive_eb_3d_invalid_case STREQUAL "negative_composition")
    set(selected_reactive_eb_3d_expected "nonnegative and nonempty")
    set(selected_reactive_eb_3d_output "negative")
  else()
    set(selected_reactive_eb_3d_expected "outside the mechanism range")
    set(selected_reactive_eb_3d_output "temperature")
  endif()
  add_test(
    NAME
      regression_selected_reactive_eb_3d_${selected_reactive_eb_3d_invalid_case}_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
      --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
      --input
        "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_${selected_reactive_eb_3d_invalid_case}.nml"
      --expected "${selected_reactive_eb_3d_expected}"
      --forbidden-output
        "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_${selected_reactive_eb_3d_output}.csv"
  )
endforeach()

add_test(
  NAME regression_selected_reactive_eb_3d_path_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    --input invalid_path_alias.nml
    --expected "input and output paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_alias.csv"
)
add_test(
  NAME regression_selected_reactive_eb_3d_checkpoint_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    --input invalid_checkpoint_interval.nml
    --expected "input and checkpoint paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_checkpoint_interval.csv"
)
add_test(
  NAME regression_selected_reactive_eb_3d_restart_input_alias_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    --input invalid_restart.nml
    --expected "input and restart paths must differ"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_restart.csv"
)
add_test(
  NAME regression_selected_reactive_eb_3d_unsupported_element_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable
      $<TARGET_FILE:test_selected_unsupported_element_reactive_eb_3d>
    --input
      "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/invalid_unsupported_element.nml"
    --expected "chemistry requires supported H/O/N species"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}/invalid_selected_reactive_eb_3d_element.csv"
)
set_tests_properties(
  regression_selected_reactive_eb_3d_unknown_species_rejected
  regression_selected_reactive_eb_3d_duplicate_species_rejected
  regression_selected_reactive_eb_3d_nonfinite_composition_rejected
  regression_selected_reactive_eb_3d_zero_composition_rejected
  regression_selected_reactive_eb_3d_negative_composition_rejected
  regression_selected_reactive_eb_3d_temperature_rejected
  regression_selected_reactive_eb_3d_path_alias_rejected
  regression_selected_reactive_eb_3d_checkpoint_input_alias_rejected
  regression_selected_reactive_eb_3d_restart_input_alias_rejected
  regression_selected_reactive_eb_3d_unsupported_element_rejected
  PROPERTIES
    RESOURCE_LOCK selected_reactive_eb_3d_invalid
    WORKING_DIRECTORY "${SELECTED_REACTIVE_EB_3D_INVALID_WORK_DIRECTORY}"
)
add_test(
  NAME regression_fixed_reactive_eb_3d_selected_model_rejected
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
    --executable $<TARGET_FILE:pelef_reactive_eb_3d>
    --input "${PROJECT_SOURCE_DIR}/cases/selected_reactive_eb_3d/fixture.nml"
    --expected "selected thermo requires explicit opt-in"
    --forbidden-output
      "${SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY}/selected_fixture_reactive_eb_3d.csv"
)
set_tests_properties(
  regression_fixed_reactive_eb_3d_selected_model_rejected
  PROPERTIES WORKING_DIRECTORY
    "${SELECTED_REACTIVE_EB_3D_FIXED_REJECTION_WORK_DIRECTORY}"
)

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

add_test(
  NAME unit_mechanism_source_hash_rejection
  COMMAND
    "${CMAKE_COMMAND}" -P
    "${PROJECT_SOURCE_DIR}/tests/cmake/check_mechanism_source_rejection.cmake"
)
add_test(
  NAME unit_mechanism_bundle_contract_rejection
  COMMAND
    "${CMAKE_COMMAND}" -P
    "${PROJECT_SOURCE_DIR}/tests/cmake/check_mechanism_bundle_contract_rejection.cmake"
)

add_test(
  NAME generated_h2o2_full_mechanism_clean
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_generated_mechanism.py"
    --generator "${PROJECT_SOURCE_DIR}/tools/generate_elementary_mechanism.py"
    --input "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
    --committed "${PROJECT_SOURCE_DIR}/src/generated/h2o2_full_mechanism_mod.F90"
    --source "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml")

set(H2O2_FULL_WORK_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/zero_d_h2o2_full")
file(MAKE_DIRECTORY "${H2O2_FULL_WORK_DIRECTORY}")
add_test(NAME regression_zero_d_h2o2_full_run
  COMMAND $<TARGET_FILE:pelef0d_h2o2_full>
    "${PROJECT_SOURCE_DIR}/cases/zero_d_h2o2_full/reactor.nml")
set_tests_properties(regression_zero_d_h2o2_full_run PROPERTIES
  WORKING_DIRECTORY "${H2O2_FULL_WORK_DIRECTORY}" TIMEOUT 300)
add_test(NAME regression_zero_d_h2o2_full_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_zero_d_h2o2_full.py"
    --input "${H2O2_FULL_WORK_DIRECTORY}/zero_d_h2o2_full.csv")
set_tests_properties(regression_zero_d_h2o2_full_structure PROPERTIES
  DEPENDS regression_zero_d_h2o2_full_run)

set(
  SELECTED_H2O2_FULL_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/selected_zero_d_h2o2_full"
)
file(MAKE_DIRECTORY "${SELECTED_H2O2_FULL_WORK_DIRECTORY}")
add_test(
  NAME regression_selected_zero_d_h2o2_full_run
  COMMAND $<TARGET_FILE:test_selected_full_reactor>
    "${PROJECT_SOURCE_DIR}/cases/selected_reactor/h2o2_full.nml"
)
set_tests_properties(
  regression_selected_zero_d_h2o2_full_run
  PROPERTIES
    WORKING_DIRECTORY "${SELECTED_H2O2_FULL_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Species: 10"
    TIMEOUT 300
)
add_test(
  NAME regression_selected_zero_d_h2o2_full_structure
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_zero_d_h2o2_full.py"
    --input
      "${SELECTED_H2O2_FULL_WORK_DIRECTORY}/selected_h2o2_full.csv"
)
set_tests_properties(
  regression_selected_zero_d_h2o2_full_structure
  PROPERTIES DEPENDS regression_selected_zero_d_h2o2_full_run
)
add_test(
  NAME regression_selected_zero_d_h2o2_full_exact_parity
  COMMAND "${CMAKE_COMMAND}" -E compare_files
    "${H2O2_FULL_WORK_DIRECTORY}/zero_d_h2o2_full.csv"
    "${SELECTED_H2O2_FULL_WORK_DIRECTORY}/selected_h2o2_full.csv"
)
set_tests_properties(
  regression_selected_zero_d_h2o2_full_exact_parity
  PROPERTIES DEPENDS
    "regression_zero_d_h2o2_full_run;regression_selected_zero_d_h2o2_full_run"
)
if(PELEF_ENABLE_CANTERA_REFERENCE)
  set(
    H2O2_FULL_INGEST_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/generated_h2o2_full_ingest"
  )
  file(MAKE_DIRECTORY "${H2O2_FULL_INGEST_WORK_DIRECTORY}")
  add_test(
    NAME generated_h2o2_full_ingest
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/ingest_cantera_mechanism.py"
      --input "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml"
      --phase ohmech
      --output "${H2O2_FULL_INGEST_WORK_DIRECTORY}/h2o2_full.json"
      --module-name h2o2_full_mechanism_mod
      --loader-name load_h2o2_full_mechanism
      --kernel-name h2o2_full_production_rates
      --jacobian-name h2o2_full_mass_fraction_jacobian
      --thermo-loader-name load_h2o2_full_thermo_data
      --transport-loader-name load_h2o2_full_transport_data
      --symbol-prefix h2o2_full
      --description
        "Full 10-species, 29-reaction H2/O2 bundle ingested from pinned Cantera YAML"
  )
  add_test(
    NAME generated_h2o2_full_ingest_clean
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
      "${H2O2_FULL_INGEST_WORK_DIRECTORY}/h2o2_full.json"
  )
  set_tests_properties(
    generated_h2o2_full_ingest_clean
    PROPERTIES DEPENDS generated_h2o2_full_ingest
  )
  add_test(NAME regression_zero_d_h2o2_full_cantera
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_h2o2_full_cantera.py"
      --input "${H2O2_FULL_WORK_DIRECTORY}/zero_d_h2o2_full.csv"
      --mechanism "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml"
      --phase ohmech)
  set_tests_properties(regression_zero_d_h2o2_full_cantera PROPERTIES
    DEPENDS regression_zero_d_h2o2_full_run TIMEOUT 300)
  if(PELEF_ENABLE_SUNDIALS)
    add_test(
      NAME regression_selected_cvode_h2o2_full_cantera
      COMMAND "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_h2o2_full_cantera.py"
        --input
          "${SELECTED_CVODE_H2O2_FULL_WORK_DIRECTORY}/selected_h2o2_full_cvode.csv"
        --mechanism "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml"
        --phase ohmech
    )
    set_tests_properties(
      regression_selected_cvode_h2o2_full_cantera
      PROPERTIES
        DEPENDS regression_selected_cvode_h2o2_full_run
        TIMEOUT 300
    )
  endif()
endif()



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

set(REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_dynamic_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_2d/hotspot.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_2d_run
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Completed regrids: [1-9][0-9]*"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_dynamic_2d.py"
    --root
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}/three_level_dynamic_root.csv"
    --middle
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}/three_level_dynamic_middle.csv"
    --finest
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}/three_level_dynamic_finest.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_2d_check
  PROPERTIES
    DEPENDS regression_reactive_eb_amr_three_level_dynamic_2d_run
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_2D_WORK_DIRECTORY}"
)

set(REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_three_level_dynamic_restart_2d")
file(MAKE_DIRECTORY
  "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_three_level_dynamic_restart_2d/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_reference
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart
  PROPERTIES
    WORKING_DIRECTORY
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}"
    PASS_REGULAR_EXPRESSION "Completed regrids: [1-9][0-9]*"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart
  PROPERTIES DEPENDS
    regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_amr_three_level_dynamic_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_three_level_dynamic_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_three_level_restart.chk"
    --expected-checkpoint-sha256
      ced92fed867264fd0863b5c85a90702e29647e2494a5d6d69b87a51f194702c5
    --reference
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_reference_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_reference_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_reference_finest.csv"
    --stopped
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_stopped_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_stopped_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_stopped_finest.csv"
    --restarted
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_restarted_root.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_restarted_middle.csv"
      "${REACTIVE_EB_AMR_THREE_LEVEL_DYNAMIC_RESTART_2D_WORK_DIRECTORY}/dynamic_restart_restarted_finest.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_three_level_dynamic_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_three_level_dynamic_restart_2d_reference;regression_reactive_eb_amr_three_level_dynamic_restart_2d_checkpoint_stop;regression_reactive_eb_amr_three_level_dynamic_restart_2d_restart"
)

set(REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_amr_restart_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_amr_restart_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_restart_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_amr_restart_2d_checkpoint_stop
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_restart_2d/checkpoint_stop.nml"
)
add_test(
  NAME regression_reactive_eb_amr_restart_2d_restart
  COMMAND $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_amr_restart_2d/restart.nml"
)
set_tests_properties(
  regression_reactive_eb_amr_restart_2d_reference
  regression_reactive_eb_amr_restart_2d_checkpoint_stop
  regression_reactive_eb_amr_restart_2d_restart
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
set_tests_properties(
  regression_reactive_eb_amr_restart_2d_restart
  PROPERTIES DEPENDS regression_reactive_eb_amr_restart_2d_checkpoint_stop
)
add_test(
  NAME regression_reactive_eb_amr_restart_2d_check
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_amr_restart_2d.py"
    --checkpoint
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart.chk"
    --expected-checkpoint-sha256
      40a0781f99a231dadb663365342b92a3e035e6dd195a517a6e45010af87eee17
    --reference
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_reference_coarse_2d.csv"
    --stopped
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_stopped_coarse_2d.csv"
    --expected-stopped-sha256
      2465d2facdd7662091d8e6c1346e8f1a0205ba1b758398309fecb7fe769aa14d
    --restarted
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restarted_coarse_2d.csv"
    --fine
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_reference_fine_2d.csv"
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restart_stopped_fine_2d.csv"
      "${REACTIVE_EB_AMR_RESTART_2D_WORK_DIRECTORY}/reactive_eb_amr_restarted_fine_2d.csv"
)
set_tests_properties(
  regression_reactive_eb_amr_restart_2d_check
  PROPERTIES DEPENDS
    "regression_reactive_eb_amr_restart_2d_reference;regression_reactive_eb_amr_restart_2d_checkpoint_stop;regression_reactive_eb_amr_restart_2d_restart"
)

set(REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/reactive_eb_chemistry_2d")
file(MAKE_DIRECTORY "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}")
add_test(
  NAME regression_reactive_eb_chemistry_2d_reference
  COMMAND $<TARGET_FILE:pelef_reactive_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_chemistry_2d/reference.nml"
)
add_test(
  NAME regression_reactive_eb_chemistry_2d_inert
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_chemistry_2d/inert.nml"
)
add_test(
  NAME regression_reactive_eb_chemistry_2d_run
  COMMAND $<TARGET_FILE:pelef_reactive_eb_2d>
    "${PROJECT_SOURCE_DIR}/cases/reactive_eb_chemistry_2d/reactive.nml"
)
set_tests_properties(
  regression_reactive_eb_chemistry_2d_reference
  regression_reactive_eb_chemistry_2d_inert
  regression_reactive_eb_chemistry_2d_run
  PROPERTIES
    WORKING_DIRECTORY "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}"
    TIMEOUT 300
)
add_test(
  NAME regression_reactive_eb_chemistry_2d_parity
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_reactive_eb_chemistry_2d.py"
    --reference
      "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_chemistry_reference.csv"
    --reactive
      "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_chemistry.csv"
    --inert
      "${REACTIVE_EB_CHEMISTRY_2D_WORK_DIRECTORY}/reactive_eb_chemistry_inert.csv"
)
set_tests_properties(
  regression_reactive_eb_chemistry_2d_parity
  PROPERTIES DEPENDS
    "regression_reactive_eb_chemistry_2d_reference;regression_reactive_eb_chemistry_2d_inert;regression_reactive_eb_chemistry_2d_run"
)

if(PELEF_ENABLE_MPI)
  add_executable(
    test_mpi_amr_transport_context_3d
    mpi/test_mpi_amr_transport_context_3d.F90
  )
  target_link_libraries(
    test_mpi_amr_transport_context_3d PRIVATE
    "$<TARGET_NAME_IF_EXISTS:pelef_mpi_support>"
  )
  add_test(
    NAME unit_mpi_amr_transport_context_3d_np2
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_mpi_amr_transport_context_3d>
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    unit_mpi_amr_transport_context_3d_np2
    PROPERTIES
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION
        "test_mpi_amr_transport_context_3d: PASS"
      TIMEOUT 60
  )
endif()
