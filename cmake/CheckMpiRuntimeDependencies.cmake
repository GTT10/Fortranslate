if(NOT DEFINED PELEF_MPI_EXECUTABLE OR NOT EXISTS "${PELEF_MPI_EXECUTABLE}")
  message(FATAL_ERROR "MPI runtime audit requires an existing executable")
endif()

file(
  GET_RUNTIME_DEPENDENCIES
  EXECUTABLES "${PELEF_MPI_EXECUTABLE}"
  RESOLVED_DEPENDENCIES_VAR resolved_dependencies
  UNRESOLVED_DEPENDENCIES_VAR unresolved_dependencies
  CONFLICTING_DEPENDENCIES_PREFIX conflicting_dependencies
)
if(unresolved_dependencies)
  message(
    FATAL_ERROR
    "unresolved runtime dependencies for ${PELEF_MPI_EXECUTABLE}: "
    "${unresolved_dependencies}"
  )
endif()
if(conflicting_dependencies_FILENAMES)
  message(
    FATAL_ERROR
    "conflicting runtime dependencies for ${PELEF_MPI_EXECUTABLE}: "
    "${conflicting_dependencies_FILENAMES}"
  )
endif()

set(mpi_abi_versions)
foreach(runtime_dependency IN LISTS resolved_dependencies)
  get_filename_component(runtime_name "${runtime_dependency}" NAME)
  if(runtime_name MATCHES "^libmpi.*\\.so\\.([0-9]+)")
    list(APPEND mpi_abi_versions "${CMAKE_MATCH_1}")
  elseif(runtime_name MATCHES "^libmpi.*\\.([0-9]+)\\.dylib")
    list(APPEND mpi_abi_versions "${CMAKE_MATCH_1}")
  endif()
endforeach()
list(REMOVE_DUPLICATES mpi_abi_versions)
list(LENGTH mpi_abi_versions mpi_abi_version_count)
if(mpi_abi_version_count GREATER 1)
  message(
    FATAL_ERROR
    "mixed MPI runtime ABI versions for ${PELEF_MPI_EXECUTABLE}: "
    "${mpi_abi_versions}"
  )
endif()
if(mpi_abi_version_count EQUAL 1)
  message(
    STATUS
    "MPI runtime dependency audit: PASS (ABI ${mpi_abi_versions})"
  )
else()
  message(
    STATUS
    "MPI runtime dependency audit: PASS (static or unversioned MPI)"
  )
endif()
