find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_executable(test_eos unit/test_eos.F90)
target_link_libraries(test_eos PRIVATE pelef_core)
add_test(NAME unit_eos COMMAND test_eos)

add_executable(test_rusanov unit/test_rusanov.F90)
target_link_libraries(test_rusanov PRIVATE pelef_core)
add_test(NAME unit_rusanov COMMAND test_rusanov)

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
