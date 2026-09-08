cmake_minimum_required(VERSION 3.23)

project(
  PeleF
  VERSION 0.246.0
  DESCRIPTION "Modern Fortran reimplementation of selected PeleC algorithms"
  LANGUAGES Fortran
)

option(PELEF_ENABLE_SPRAY "Build experimental conservative spray/LES applications" OFF)
option(PELEF_ENABLE_TESTS "Build PeleF tests" ON)
option(PELEF_ENABLE_CANTERA_REFERENCE "Enable live Cantera parity tests" OFF)
option(PELEF_ENABLE_MPI "Build MPI verification targets" OFF)
option(
  PELEF_ENABLE_SUNDIALS
  "Build the optional official-Fortran CVODE selected-reactor backend"
  OFF
)
set(
  PELEF_MECHANISM_BUNDLE "" CACHE FILEPATH
  "Optional normalized JSON mechanism bundle to compile and probe"
)
set(
  PELEF_MECHANISM_SOURCE "" CACHE FILEPATH
  "Optional source YAML whose SHA-256 must match the selected bundle"
)
if(PELEF_MECHANISM_SOURCE AND NOT PELEF_MECHANISM_BUNDLE)
  message(FATAL_ERROR "PELEF_MECHANISM_SOURCE requires PELEF_MECHANISM_BUNDLE")
endif()

if(PELEF_ENABLE_TESTS OR PELEF_MECHANISM_BUNDLE)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

include(cmake/RequireMpiF08.cmake)
include(cmake/MechanismBundle.cmake)
include(cmake/RequireSundialsDoublePrecision.cmake)
if(PELEF_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  pelef_require_mpi_f08()
endif()
if(PELEF_ENABLE_SUNDIALS)
  find_package(
    SUNDIALS 7.2.0 EXACT CONFIG REQUIRED
    COMPONENTS
      fcvode_mod
      fnvecserial_mod
      fsunmatrixdense_mod
      fsunlinsoldense_mod
  )
  foreach(
      sundials_target
      SUNDIALS::fcvode_mod_static
      SUNDIALS::fnvecserial_mod_static
      SUNDIALS::fsunmatrixdense_mod_static
      SUNDIALS::fsunlinsoldense_mod_static
    )
    if(NOT TARGET "${sundials_target}")
      message(
        FATAL_ERROR
        "SUNDIALS official Fortran target is unavailable: ${sundials_target}"
      )
    endif()
  endforeach()
  get_target_property(
    sundials_include_directories
    SUNDIALS::fcvode_mod_static
    INTERFACE_INCLUDE_DIRECTORIES
  )
  set(sundials_config_header "")
  foreach(sundials_include_directory IN LISTS sundials_include_directories)
    if(
        EXISTS
          "${sundials_include_directory}/sundials/sundials_config.h"
      )
      set(
        sundials_config_header
        "${sundials_include_directory}/sundials/sundials_config.h"
      )
      break()
    endif()
  endforeach()
  if(NOT sundials_config_header)
    message(
      FATAL_ERROR
      "SUNDIALS official precision header is unavailable from "
      "SUNDIALS::fcvode_mod_static"
    )
  endif()
  pelef_require_sundials_double_precision("${sundials_config_header}")
endif()
if(PELEF_ENABLE_TESTS)
  include(CTest)
  enable_testing()
endif()

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include("${PROJECT_SOURCE_DIR}/cmake/targets/Core.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/targets/SelectedRuntime.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/targets/Spray.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/targets/Applications.cmake")
if(PELEF_ENABLE_TESTS)
  add_subdirectory(tests)
endif()


include("${PROJECT_SOURCE_DIR}/cmake/targets/MPI.cmake")
include("${PROJECT_SOURCE_DIR}/cmake/targets/InstalledTests.cmake")
