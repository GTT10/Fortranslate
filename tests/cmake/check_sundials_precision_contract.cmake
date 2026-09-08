if(NOT DEFINED TEST_PRECISION)
  message(FATAL_ERROR "TEST_PRECISION must be double, single, or extended")
endif()
if(NOT DEFINED TEST_BINARY_DIR)
  message(FATAL_ERROR "TEST_BINARY_DIR is required")
endif()

set(test_directory "${TEST_BINARY_DIR}/sundials_precision_${TEST_PRECISION}")
set(config_header "${test_directory}/sundials_config.h")
file(MAKE_DIRECTORY "${test_directory}")
if(TEST_PRECISION STREQUAL "double")
  file(WRITE "${config_header}" "#define SUNDIALS_DOUBLE_PRECISION 1\n")
elseif(TEST_PRECISION STREQUAL "single")
  file(WRITE "${config_header}" "#define SUNDIALS_SINGLE_PRECISION 1\n")
elseif(TEST_PRECISION STREQUAL "extended")
  file(WRITE "${config_header}" "#define SUNDIALS_EXTENDED_PRECISION 1\n")
else()
  message(FATAL_ERROR "TEST_PRECISION must be double, single, or extended")
endif()

include(
  "${CMAKE_CURRENT_LIST_DIR}/../../cmake/RequireSundialsDoublePrecision.cmake"
)
pelef_require_sundials_double_precision("${config_header}")
