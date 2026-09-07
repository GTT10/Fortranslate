set(source "${CMAKE_CURRENT_LIST_DIR}/../fixtures/mechanism_bundle.json")
file(SHA256 "${source}" expected_sha256)
execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    "-DPELEF_MECHANISM_SOURCE=${source}"
    "-DPELEF_MECHANISM_SHA256=${expected_sha256}"
    "-DPELEF_MECHANISM_SOURCE_FILE=mechanism_bundle.json"
    -P "${CMAKE_CURRENT_LIST_DIR}/../../cmake/CheckMechanismSource.cmake"
  RESULT_VARIABLE accepted_result
  OUTPUT_VARIABLE accepted_output
  ERROR_VARIABLE accepted_error
)
if(NOT accepted_result EQUAL 0)
  message(
    FATAL_ERROR
    "the correct mechanism source was rejected: "
    "${accepted_output}${accepted_error}"
  )
endif()

execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    "-DPELEF_MECHANISM_SOURCE=${source}"
    "-DPELEF_MECHANISM_SHA256=0000000000000000000000000000000000000000000000000000000000000000"
    -P "${CMAKE_CURRENT_LIST_DIR}/../../cmake/CheckMechanismSource.cmake"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error
)

if(result EQUAL 0)
  message(FATAL_ERROR "an incorrect mechanism source SHA-256 was accepted")
endif()
set(combined "${output}\n${error}")
string(FIND "${combined}" "mechanism source SHA-256 mismatch" message_position)
if(message_position EQUAL -1)
  message(FATAL_ERROR "source rejection did not contain the required diagnostic")
endif()

execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    "-DPELEF_MECHANISM_SOURCE=${source}"
    "-DPELEF_MECHANISM_SHA256=${expected_sha256}"
    "-DPELEF_MECHANISM_SOURCE_FILE=wrong-name.json"
    -P "${CMAKE_CURRENT_LIST_DIR}/../../cmake/CheckMechanismSource.cmake"
  RESULT_VARIABLE basename_result
  OUTPUT_VARIABLE basename_output
  ERROR_VARIABLE basename_error
)
if(basename_result EQUAL 0)
  message(FATAL_ERROR "an incorrect mechanism source basename was accepted")
endif()
set(basename_combined "${basename_output}\n${basename_error}")
string(FIND "${basename_combined}" "mechanism source basename mismatch" basename_position)
if(basename_position EQUAL -1)
  message(FATAL_ERROR "basename rejection lacked the required diagnostic")
endif()

message(STATUS "mechanism source rejection contract: PASS")
