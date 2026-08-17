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
