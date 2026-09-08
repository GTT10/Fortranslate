set(full_bundle "${CMAKE_CURRENT_LIST_DIR}/../fixtures/mechanism_bundle.json")
set(
  legacy_bundle
  "${CMAKE_CURRENT_LIST_DIR}/../fixtures/mechanism_bundle_legacy.json"
)
set(contract_script "${CMAKE_CURRENT_LIST_DIR}/check_mechanism_bundle_contract.cmake")

execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    "-DPELEF_TEST_MECHANISM_BUNDLE=${full_bundle}"
    -P "${contract_script}"
  RESULT_VARIABLE accepted_result
  OUTPUT_VARIABLE accepted_output
  ERROR_VARIABLE accepted_error
)
if(NOT accepted_result EQUAL 0)
  message(
    FATAL_ERROR
    "a complete thermo/transport bundle was rejected: "
    "${accepted_output}${accepted_error}"
  )
endif()

execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    "-DPELEF_TEST_MECHANISM_BUNDLE=${legacy_bundle}"
    -P "${contract_script}"
  RESULT_VARIABLE legacy_result
  OUTPUT_VARIABLE legacy_output
  ERROR_VARIABLE legacy_error
)
if(legacy_result EQUAL 0)
  message(FATAL_ERROR "a legacy kinetics-only bundle reached the build helper")
endif()
set(legacy_combined "${legacy_output}\n${legacy_error}")
string(
  FIND "${legacy_combined}"
  "selected mechanism build contract requires explicit jacobian_name"
  legacy_message_position
)
if(legacy_message_position EQUAL -1)
  message(FATAL_ERROR "legacy bundle rejection lacked the required diagnostic")
endif()
string(FIND "${legacy_combined}" "generator-only" generator_only_position)
if(generator_only_position EQUAL -1)
  message(FATAL_ERROR "legacy bundle rejection did not state its supported route")
endif()

message(STATUS "mechanism bundle build-contract rejection: PASS")
