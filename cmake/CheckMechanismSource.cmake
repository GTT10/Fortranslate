if(NOT DEFINED PELEF_MECHANISM_SOURCE)
  message(FATAL_ERROR "PELEF_MECHANISM_SOURCE is required")
endif()
if(NOT DEFINED PELEF_MECHANISM_SHA256)
  message(FATAL_ERROR "PELEF_MECHANISM_SHA256 is required")
endif()
if(NOT EXISTS "${PELEF_MECHANISM_SOURCE}")
  message(FATAL_ERROR "mechanism source does not exist: ${PELEF_MECHANISM_SOURCE}")
endif()

file(SHA256 "${PELEF_MECHANISM_SOURCE}" actual_sha256)
if(NOT actual_sha256 STREQUAL PELEF_MECHANISM_SHA256)
  message(
    FATAL_ERROR
    "mechanism source SHA-256 mismatch: expected ${PELEF_MECHANISM_SHA256}, "
    "found ${actual_sha256}"
  )
endif()

if(DEFINED PELEF_MECHANISM_SOURCE_FILE AND
   NOT PELEF_MECHANISM_SOURCE_FILE STREQUAL "")
  get_filename_component(actual_source_file "${PELEF_MECHANISM_SOURCE}" NAME)
  if(NOT actual_source_file STREQUAL PELEF_MECHANISM_SOURCE_FILE)
    message(
      FATAL_ERROR
      "mechanism source basename mismatch: expected "
      "${PELEF_MECHANISM_SOURCE_FILE}, found ${actual_source_file}"
    )
  endif()
endif()
