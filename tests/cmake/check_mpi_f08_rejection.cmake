execute_process(
  COMMAND
    "${CMAKE_COMMAND}" -DTEST_MPI_F08=OFF -P
    "${CMAKE_CURRENT_LIST_DIR}/check_mpi_f08_contract.cmake"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error
)

if(result EQUAL 0)
  message(FATAL_ERROR "an unavailable mpi_f08 module was accepted")
endif()

set(combined "${output}\n${error}")
string(
  FIND "${combined}"
  "requires a usable MPI Fortran 2008 mpi_f08 module"
  message_position
)
if(message_position EQUAL -1)
  message(FATAL_ERROR "mpi_f08 rejection did not contain the required diagnostic")
endif()

message(STATUS "mpi_f08 rejection contract: PASS")
