if(NOT DEFINED TEST_BINARY_DIR)
  message(FATAL_ERROR "TEST_BINARY_DIR is required")
endif()

foreach(rejected_precision single extended)
  execute_process(
    COMMAND "${CMAKE_COMMAND}"
      -DTEST_PRECISION=${rejected_precision}
      -DTEST_BINARY_DIR=${TEST_BINARY_DIR}
      -P
      "${CMAKE_CURRENT_LIST_DIR}/check_sundials_precision_contract.cmake"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
  )
  if(result EQUAL 0)
    message(
      FATAL_ERROR
      "SUNDIALS ${rejected_precision}-precision contract was accepted"
    )
  endif()
  set(combined_output "${output}\n${error}")
  if(NOT combined_output MATCHES "requires a double-precision SUNDIALS ABI")
    message(
      FATAL_ERROR
      "SUNDIALS precision rejection did not report the expected diagnostic:\n"
      "${combined_output}"
    )
  endif()
endforeach()
