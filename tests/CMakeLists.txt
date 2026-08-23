find_package(Python3 REQUIRED COMPONENTS Interpreter)

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

add_test(
  NAME generated_h2o2_full_mechanism_clean
  COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_generated_mechanism.py"
    --generator "${PROJECT_SOURCE_DIR}/tools/generate_elementary_mechanism.py"
    --input "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json"
    --committed "${PROJECT_SOURCE_DIR}/src/generated/h2o2_full_mechanism_mod.F90")

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
if(PELEF_ENABLE_CANTERA_REFERENCE)
  add_test(NAME regression_zero_d_h2o2_full_cantera
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_h2o2_full_cantera.py"
      --input "${H2O2_FULL_WORK_DIRECTORY}/zero_d_h2o2_full.csv"
      --mechanism h2o2.yaml --phase ohmech)
  set_tests_properties(regression_zero_d_h2o2_full_cantera PROPERTIES
    DEPENDS regression_zero_d_h2o2_full_run TIMEOUT 300)
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
