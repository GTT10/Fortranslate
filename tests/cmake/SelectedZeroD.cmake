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

