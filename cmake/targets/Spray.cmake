# Optional, isolated Lagrangian spray and deviatoric SGS runtime.
# Default OFF preserves the existing fixed/selected installation contracts.
if(PELEF_ENABLE_SPRAY)
  function(pelef_add_spray_runtime target base)
    add_library(${target} STATIC
      src/spray/spray_parcel_mod.F90
      src/spray/spray_coupling_3d_mod.F90
      src/spray/spray_reactive_3d_mod.F90
      src/spray/spray_checkpoint_mod.F90
      src/spray/spray_application_mod.F90
      src/les/les_smagorinsky_3d_mod.F90
    )
    set_target_properties(${target} PROPERTIES
      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${target}_modules")
    target_include_directories(${target} PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/${target}_modules>")
    target_link_libraries(${target} PUBLIC ${base})
    if(PELEF_ENABLE_SUNDIALS)
      target_sources(${target} PRIVATE src/chemistry/sundials_constant_volume_reactor_mod.F90)
      target_compile_definitions(${target} PRIVATE PELEF_SPRAY_CVODE=1)
      target_link_libraries(${target} PRIVATE SUNDIALS::fcvode_mod_static
        SUNDIALS::fnvecserial_mod_static SUNDIALS::fsunmatrixdense_mod_static SUNDIALS::fsunlinsoldense_mod_static)
    endif()
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${target} PRIVATE -Wall -Wextra -Wimplicit-interface)
      if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
      endif()
    endif()
  endfunction()
  pelef_add_spray_runtime(pelef_spray_runtime pelef_core)
  if(PELEF_MECHANISM_BUNDLE)
    pelef_add_spray_runtime(pelef_selected_spray_runtime pelef_selected_reactive_3d_runtime)
  endif()
  file(SHA256 "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_full.json" PELEF_SPRAY_FIXED_SHA256)
  configure_file(app/pelef_spray_3d.F90.in "${CMAKE_CURRENT_BINARY_DIR}/pelef_spray_3d.F90" @ONLY)
  add_executable(pelef_spray_3d "${CMAKE_CURRENT_BINARY_DIR}/pelef_spray_3d.F90")
  target_link_libraries(pelef_spray_3d PRIVATE pelef_spray_runtime)
  install(TARGETS pelef_spray_3d RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
