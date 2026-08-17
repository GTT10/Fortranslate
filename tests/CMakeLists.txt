find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_eos unit/test_eos.F90)
target_link_libraries(test_eos PRIVATE pelef_core)
add_test(NAME unit_eos COMMAND test_eos)

add_executable(test_rusanov unit/test_rusanov.F90)
target_link_libraries(test_rusanov PRIVATE pelef_core)
add_test(NAME unit_rusanov COMMAND test_rusanov)

add_executable(test_slope_limiter unit/test_slope_limiter.F90)
target_link_libraries(test_slope_limiter PRIVATE pelef_core)
add_test(NAME unit_slope_limiter COMMAND test_slope_limiter)

add_executable(test_plm_reconstruction unit/test_plm_reconstruction.F90)
target_link_libraries(test_plm_reconstruction PRIVATE pelef_core)
add_test(NAME unit_plm_reconstruction COMMAND test_plm_reconstruction)

add_executable(test_entropy_wave regression/test_entropy_wave.F90)
target_link_libraries(test_entropy_wave PRIVATE pelef_core)
add_test(NAME regression_entropy_wave_convergence COMMAND test_entropy_wave)

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
