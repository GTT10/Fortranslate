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

