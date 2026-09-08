if(NOT DEFINED PELEF_TEST_MECHANISM_BUNDLE)
  message(FATAL_ERROR "PELEF_TEST_MECHANISM_BUNDLE is required")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/MechanismBundle.cmake")
_pelef_read_mechanism_bundle_build_contract(
  "${PELEF_TEST_MECHANISM_BUNDLE}"
  test_bundle
)

if(NOT test_bundle_nspecies GREATER 0 OR
    NOT test_bundle_nreactions GREATER 0)
  message(FATAL_ERROR "mechanism bundle build contract returned invalid sizes")
endif()

message(STATUS "mechanism bundle build contract: PASS")
