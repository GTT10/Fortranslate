if(NOT DEFINED TEST_BUILD_DIR OR TEST_BUILD_DIR STREQUAL "")
  message(FATAL_ERROR "TEST_BUILD_DIR is required")
endif()
if(NOT DEFINED TEST_INSTALL_PREFIX OR TEST_INSTALL_PREFIX STREQUAL "")
  message(FATAL_ERROR "TEST_INSTALL_PREFIX is required")
endif()
if(NOT DEFINED TEST_CONFIG OR TEST_CONFIG STREQUAL "")
  message(FATAL_ERROR "TEST_CONFIG is required")
endif()
if(NOT IS_ABSOLUTE "${TEST_BUILD_DIR}" OR
    NOT IS_ABSOLUTE "${TEST_INSTALL_PREFIX}")
  message(FATAL_ERROR "test build and install paths must be absolute")
endif()

set(test_build_dir "${TEST_BUILD_DIR}")
set(test_install_prefix "${TEST_INSTALL_PREFIX}")
cmake_path(NORMAL_PATH test_build_dir)
cmake_path(NORMAL_PATH test_install_prefix)
cmake_path(
  IS_PREFIX test_build_dir "${test_install_prefix}" NORMALIZE
  install_prefix_is_private
)
if(NOT install_prefix_is_private OR
    test_install_prefix STREQUAL test_build_dir)
  message(FATAL_ERROR "refusing to clear a non-private install prefix")
endif()

file(REMOVE_RECURSE "${test_install_prefix}")
file(MAKE_DIRECTORY "${test_install_prefix}")
execute_process(
  COMMAND
    "${CMAKE_COMMAND}" --install "${test_build_dir}"
    --prefix "${test_install_prefix}"
    --config "${TEST_CONFIG}"
  RESULT_VARIABLE install_status
  OUTPUT_VARIABLE install_output
  ERROR_VARIABLE install_error
)
if(NOT install_output STREQUAL "")
  message("${install_output}")
endif()
if(NOT install_status EQUAL 0)
  message(FATAL_ERROR
    "clean-prefix install failed (${install_status}): ${install_error}")
endif()
