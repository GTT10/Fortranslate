if(PELEF_ENABLE_SPRAY)
  add_executable(test_spray_core spray/test_spray_core.F90)
  target_link_libraries(test_spray_core PRIVATE pelef_spray_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_spray_core PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  foreach(spray_case drag heat evaporation terminal abramzon film_binary injection size_distribution stencil coupling rollback substeps zero invalid
      les_algebra les_grid les_zero les_reject)
    add_test(NAME spray_${spray_case} COMMAND test_spray_core ${spray_case})
    set_tests_properties(spray_${spray_case} PROPERTIES LABELS "spray" TIMEOUT 60)
  endforeach()
  add_executable(test_spray_integration spray/test_spray_integration.F90)
  target_link_libraries(test_spray_integration PRIVATE pelef_spray_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_spray_integration PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  foreach(spray_case native coupled rollback checkpoint)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/spray_${spray_case}")
    add_test(NAME spray_integration_${spray_case} COMMAND test_spray_integration ${spray_case})
    set_tests_properties(spray_integration_${spray_case} PROPERTIES LABELS "spray"
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/spray_${spray_case}" TIMEOUT 180)
  endforeach()
  if(PELEF_ENABLE_SUNDIALS)
    add_test(NAME spray_integration_cvode COMMAND test_spray_integration cvode)
    set_tests_properties(spray_integration_cvode PROPERTIES LABELS "spray;cvode" TIMEOUT 180)
  endif()
  add_test(NAME spray_mechanism_portability COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/python/test_mechanism_portability.py")
  set_tests_properties(spray_mechanism_portability PROPERTIES LABELS "spray" TIMEOUT 60)
  add_test(NAME spray_application_cli COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/spray/test_spray_application.py" --exe "$<TARGET_FILE:pelef_spray_3d>")
  set_tests_properties(spray_application_cli PROPERTIES LABELS "spray" TIMEOUT 180)
  if(PELEF_ENABLE_SUNDIALS AND PELEF_ENABLE_CANTERA_REFERENCE)
    add_test(NAME spray_cantera_h2_chemistry COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_cantera.py" --exe "$<TARGET_FILE:pelef_spray_3d>"
      --yaml "${PROJECT_SOURCE_DIR}/mechanisms/h2o2_cantera.yaml" --phase ohmech --case h2-chemistry)
    set_tests_properties(spray_cantera_h2_chemistry PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 660)
    if(PELEF_MECHANISM_BUNDLE AND PELEF_MECHANISM_SOURCE)
      file(SHA256 "${PELEF_MECHANISM_SOURCE}" spray_selected_source_sha)
      file(READ "${PROJECT_SOURCE_DIR}/references/spray_methanol.json" spray_reference_json)
      string(JSON spray_methanol_source_sha GET "${spray_reference_json}" files yaml sha256)
      if(spray_selected_source_sha STREQUAL spray_methanol_source_sha)
        foreach(reference_case methanol-chemistry methanol-spray)
          add_test(NAME spray_cantera_${reference_case} COMMAND "${Python3_EXECUTABLE}"
            "${PROJECT_SOURCE_DIR}/tools/check_spray_cantera.py" --exe "$<TARGET_FILE:pelef_spray_3d_selected>"
            --yaml "${PELEF_MECHANISM_SOURCE}" --phase gas --case ${reference_case})
          set_tests_properties(spray_cantera_${reference_case} PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 660)
        endforeach()
        add_test(NAME spray_reactive_cone_restart COMMAND "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/check_spray_cone.py" --exe "$<TARGET_FILE:pelef_spray_3d_selected>"
          --yaml "${PELEF_MECHANISM_SOURCE}")
        set_tests_properties(spray_reactive_cone_restart PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 1200)
        add_test(NAME spray_temporal_refinement COMMAND "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/check_spray_convergence.py" --exe "$<TARGET_FILE:pelef_spray_3d_selected>"
          --yaml "${PELEF_MECHANISM_SOURCE}")
        set_tests_properties(spray_temporal_refinement PROPERTIES LABELS "spray;cvode;cantera" TIMEOUT 1800)
      endif()
    endif()
  endif()
endif()
