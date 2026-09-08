include_guard(GLOBAL)

function(pelef_require_mpi_f08)
  if(NOT MPI_Fortran_HAVE_F08_MODULE)
    message(
      FATAL_ERROR
      "PeleF requires a usable MPI Fortran 2008 mpi_f08 module. "
      "The selected MPI wrapper '${MPI_Fortran_COMPILER}' is not compatible "
      "with CMAKE_Fortran_COMPILER='${CMAKE_Fortran_COMPILER}'. Set "
      "MPI_Fortran_COMPILER to a wrapper from the same compiler toolchain."
    )
  endif()

  if(NOT MPI_Fortran_COMPILER OR NOT MPIEXEC_EXECUTABLE)
    message(
      FATAL_ERROR
      "PeleF requires explicit MPI_Fortran_COMPILER and MPIEXEC_EXECUTABLE "
      "paths from one MPI installation."
    )
  endif()
  get_filename_component(
    mpi_fortran_compiler_real "${MPI_Fortran_COMPILER}" REALPATH
  )
  get_filename_component(mpiexec_real "${MPIEXEC_EXECUTABLE}" REALPATH)
  get_filename_component(
    mpi_fortran_compiler_directory "${mpi_fortran_compiler_real}" DIRECTORY
  )
  get_filename_component(mpiexec_directory "${mpiexec_real}" DIRECTORY)
  if(NOT mpi_fortran_compiler_directory STREQUAL mpiexec_directory)
    message(
      FATAL_ERROR
      "MPI_Fortran_COMPILER='${MPI_Fortran_COMPILER}' and "
      "MPIEXEC_EXECUTABLE='${MPIEXEC_EXECUTABLE}' must resolve to the same "
      "MPI installation directory."
    )
  endif()
endfunction()

function(pelef_require_coherent_mpi_runtime target_name)
  if(NOT TARGET "${target_name}")
    message(FATAL_ERROR "unknown MPI runtime-audit target: ${target_name}")
  endif()
  add_custom_command(
    TARGET "${target_name}"
    POST_BUILD
    COMMAND
      "${CMAKE_COMMAND}"
      "-DPELEF_MPI_EXECUTABLE=$<TARGET_FILE:${target_name}>"
      -P "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/CheckMpiRuntimeDependencies.cmake"
    VERBATIM
  )
endfunction()
