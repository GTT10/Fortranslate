execute_process(
  COMMAND
    "${CMAKE_COMMAND}"
    -DTEST_MPI_F08=ON
    -DTEST_MPI_FORTRAN_COMPILER=${CMAKE_COMMAND}
    -DTEST_MPIEXEC_EXECUTABLE=/opt/foreign-mpi/bin/mpiexec
    -P "${CMAKE_CURRENT_LIST_DIR}/check_mpi_f08_contract.cmake"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error
)

if(result EQUAL 0)
  message(FATAL_ERROR "MPI wrapper/launcher mismatch was accepted")
endif()
set(combined "${output}\n${error}")
string(
  FIND "${combined}"
  "must resolve to the same"
  message_position
)
if(message_position EQUAL -1)
  message(FATAL_ERROR "MPI launcher rejection diagnostic was missing")
endif()

message(STATUS "MPI wrapper/launcher rejection contract: PASS")
