if(PELEF_ENABLE_SPRAY)
  add_test(NAME spray_eb_difmag_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_difmag_application" --force --restart --difmag .037)
  set_tests_properties(spray_eb_difmag_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_rate_files_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/spray/test_spray_application.py" --exe "$<TARGET_FILE:pelef_spray_3d>"
    SprayApplication.test_external_rate_files SprayApplication.test_piecewise_injection_profile)
  set_tests_properties(spray_rate_files_application PROPERTIES LABELS "spray;application" TIMEOUT 120)
  add_test(NAME spray_eb_rate_files_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_rate_files_application" --properties --ports --rate --rate-files)
  set_tests_properties(spray_eb_rate_files_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_force_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_force_application" --force --restart)
  set_tests_properties(spray_eb_force_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_splash_outflow_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_splash_outflow_application" --impacts splash --outflow --restart)
  set_tests_properties(spray_eb_splash_outflow_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_periodic_film_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_periodic_film_application" --periodic --impacts deposit --restart)
  set_tests_properties(spray_eb_periodic_film_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_outflow_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_outflow_application" --outflow --restart)
  set_tests_properties(spray_eb_outflow_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_port_models_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_port_models_application" --properties --ports --rate --models)
  set_tests_properties(spray_eb_port_models_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_port_profiles_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_port_profiles_application" --properties --ports --rate)
  set_tests_properties(spray_eb_port_profiles_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_ports_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_ports_application" --properties --ports)
  set_tests_properties(spray_eb_ports_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_rate_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_rate_application" --rate)
  set_tests_properties(spray_eb_rate_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_properties_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_properties_application" --properties)
  set_tests_properties(spray_eb_properties_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_nozzle_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray_nozzle.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_nozzle_application")
  set_tests_properties(spray_eb_nozzle_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_khrt_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_khrt_application" --restart --khrt)
  set_tests_properties(spray_eb_khrt_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_split_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_split_application" --scheduled --restart --split)
  set_tests_properties(spray_eb_split_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  foreach(impact_mode deposit splash)
    add_test(NAME spray_eb_app_${impact_mode} COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
      --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_app_${impact_mode}" --scheduled --restart --impacts ${impact_mode})
    set_tests_properties(spray_eb_app_${impact_mode} PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  endforeach()
  add_executable(test_spray_eb_impact_reactive spray/test_spray_eb_impact_reactive.F90)
  target_link_libraries(test_spray_eb_impact_reactive PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_impact_reactive COMMAND test_spray_eb_impact_reactive)
  set_tests_properties(spray_eb_impact_reactive PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_impact_coupling spray/test_spray_eb_impact_coupling.F90)
  target_link_libraries(test_spray_eb_impact_coupling PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_impact_coupling COMMAND test_spray_eb_impact_coupling)
  set_tests_properties(spray_eb_impact_coupling PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_impact spray/test_spray_eb_impact.F90)
  target_link_libraries(test_spray_eb_impact PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_impact COMMAND test_spray_eb_impact)
  set_tests_properties(spray_eb_impact PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_tree_stencil spray/test_spray_eb_tree_stencil.F90)
  target_link_libraries(test_spray_eb_tree_stencil PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_tree_stencil COMMAND test_spray_eb_tree_stencil)
  set_tests_properties(spray_eb_tree_stencil PROPERTIES LABELS "spray;eb;amr" TIMEOUT 30)
  add_test(NAME spray_eb_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_restart" --scheduled --restart)
  set_tests_properties(spray_eb_restart PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_schedule COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_schedule" --scheduled)
  set_tests_properties(spray_eb_schedule PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_test(NAME spray_eb_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_application")
  set_tests_properties(spray_eb_application PROPERTIES LABELS "spray;eb;application" TIMEOUT 60)
  add_executable(test_spray_eb_reactive spray/test_spray_eb_reactive.F90)
  target_link_libraries(test_spray_eb_reactive PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_reactive COMMAND test_spray_eb_reactive)
  set_tests_properties(spray_eb_reactive PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_step spray/test_spray_eb_step.F90)
  target_link_libraries(test_spray_eb_step PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_step COMMAND test_spray_eb_step)
  set_tests_properties(spray_eb_step PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_field spray/test_spray_eb_field.F90)
  target_link_libraries(test_spray_eb_field PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_field COMMAND test_spray_eb_field)
  set_tests_properties(spray_eb_field PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_drift spray/test_spray_eb_drift.F90)
  target_link_libraries(test_spray_eb_drift PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_drift COMMAND test_spray_eb_drift)
  set_tests_properties(spray_eb_drift PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_eb_intersection spray/test_spray_eb_intersection.F90)
  target_link_libraries(test_spray_eb_intersection PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_eb_intersection COMMAND test_spray_eb_intersection)
  set_tests_properties(spray_eb_intersection PROPERTIES LABELS "spray;eb" TIMEOUT 30)
  add_executable(test_spray_amr_tree_gas spray/test_spray_amr_tree_gas.F90)
  target_link_libraries(test_spray_amr_tree_gas PRIVATE pelef_spray_tree_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_spray_amr_tree_gas PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME spray_amr_tree_gas COMMAND test_spray_amr_tree_gas)
  set_tests_properties(spray_amr_tree_gas PROPERTIES LABELS "spray;amr" TIMEOUT 30)
  add_executable(test_spray_amr_tree_coupling spray/test_spray_amr_tree_coupling.F90)
  target_link_libraries(test_spray_amr_tree_coupling PRIVATE pelef_spray_tree_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_spray_amr_tree_coupling PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME spray_amr_tree_coupling COMMAND test_spray_amr_tree_coupling)
  set_tests_properties(spray_amr_tree_coupling PROPERTIES LABELS "spray;amr" TIMEOUT 30)
  add_executable(test_spray_amr_tree_partition spray/test_spray_amr_tree_partition.F90)
  target_link_libraries(test_spray_amr_tree_partition PRIVATE pelef_spray_tree_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_spray_amr_tree_partition PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  add_test(NAME spray_amr_tree_partition COMMAND test_spray_amr_tree_partition)
  set_tests_properties(spray_amr_tree_partition PROPERTIES LABELS "spray;amr" TIMEOUT 30)
  if(PELEF_ENABLE_MPI)
    add_test(NAME spray_mpi_amr_film_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_film_application"
        --multiport --physical --films)
    set_tests_properties(spray_mpi_amr_film_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_inflow_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_inflow_application"
        --multiport --physical --inflow)
    set_tests_properties(spray_mpi_amr_inflow_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_physical_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_physical_application" --multiport --physical)
    set_tests_properties(spray_mpi_amr_physical_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_rebalance_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_rebalance_application" --multiport --adaptive --rebalance)
    set_tests_properties(spray_mpi_amr_rebalance_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_adaptation
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_adaptation" --multiport --adaptive)
    set_tests_properties(spray_mpi_amr_adaptation PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_ports
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_ports" --multiport)
    set_tests_properties(spray_mpi_amr_ports PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_split_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_split_application" --split)
    set_tests_properties(spray_mpi_amr_split_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_discharge_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_discharge_application"
        --multiport --rates --discharge)
    set_tests_properties(spray_mpi_amr_discharge_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_rate_files_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_rate_files_application"
        --multiport --rates --rate-files)
    set_tests_properties(spray_mpi_amr_rate_files_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_rates_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_rates_application" --multiport --rates)
    set_tests_properties(spray_mpi_amr_rates_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_sizes_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_sizes_application" --multiport --sizes)
    set_tests_properties(spray_mpi_amr_sizes_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_properties_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_properties_application" --properties)
    set_tests_properties(spray_mpi_amr_properties_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_force_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_force_application" --force)
    set_tests_properties(spray_mpi_amr_force_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_pressure_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_pressure_application"
        --multiport --physical --les --pressure)
    set_tests_properties(spray_mpi_amr_pressure_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_application" --les --dynamic)
    set_tests_properties(spray_mpi_amr_dynamic_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_subcycle
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_subcycle" --les --dynamic --subcycling)
    set_tests_properties(spray_mpi_amr_dynamic_subcycle PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_adaptation
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_adaptation"
        --les --dynamic --multiport --adaptive --rebalance)
    set_tests_properties(spray_mpi_amr_dynamic_adaptation PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_ppm_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_ppm_restart"
        --ppm --subcycling)
    set_tests_properties(spray_mpi_amr_ppm_restart PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_difmag_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_difmag_restart"
        --ppm --subcycling --difmag)
    set_tests_properties(spray_mpi_amr_difmag_restart PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_open
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_open"
        --les --dynamic --dynamic-open)
    set_tests_properties(spray_mpi_amr_dynamic_open PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_mach
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_mach"
        --les --dynamic --dynamic-open --mach-inlet)
    set_tests_properties(spray_mpi_amr_dynamic_mach PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_pressure
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_pressure"
        --les --dynamic --dynamic-open --dynamic-pressure)
    set_tests_properties(spray_mpi_amr_dynamic_pressure PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_dynamic_walls
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_dynamic_walls"
        --les --dynamic --dynamic-walls)
    set_tests_properties(spray_mpi_amr_dynamic_walls PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 300 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_les_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe "$<TARGET_FILE:pelef_mpi_amr_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_les_application" --les)
    set_tests_properties(spray_mpi_amr_les_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray_restart.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_restart")
    set_tests_properties(spray_mpi_amr_restart PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_test(NAME spray_mpi_amr_application
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_amr_spray.py"
        --exe $<TARGET_FILE:pelef_mpi_amr_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_amr_application")
    set_tests_properties(spray_mpi_amr_application PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 180 PROCESSORS 4)
    add_executable(test_spray_mpi_amr_balance spray/test_spray_mpi_amr_balance.F90)
    target_link_libraries(test_spray_mpi_amr_balance PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_balance PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film_reactive spray/test_spray_mpi_amr_film_reactive.F90)
    target_link_libraries(test_spray_mpi_amr_film_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film_reactive PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film_checkpoint spray/test_spray_mpi_amr_film_checkpoint.F90)
    target_link_libraries(test_spray_mpi_amr_film_checkpoint PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film_checkpoint PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_film_step spray/test_spray_mpi_eb_film_step.F90)
    target_link_libraries(test_spray_mpi_eb_film_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_eb_film_step PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film_step spray/test_spray_mpi_amr_film_step.F90)
    target_link_libraries(test_spray_mpi_amr_film_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film_step PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_film spray/test_spray_mpi_eb_film.F90)
    target_link_libraries(test_spray_mpi_eb_film PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_eb_film PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_film spray/test_spray_mpi_amr_film.F90)
    target_link_libraries(test_spray_mpi_amr_film PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_film PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_remap spray/test_spray_mpi_amr_remap.F90)
    target_link_libraries(test_spray_mpi_amr_remap PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_remap PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_reactive spray/test_spray_mpi_amr_reactive.F90)
    add_test(NAME spray_mpi_amr_split_2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_spray_mpi_amr_reactive> split ${MPIEXEC_POSTFLAGS})
    set_tests_properties(spray_mpi_amr_split_2 PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 120 PROCESSORS 2)
    add_test(NAME spray_mpi_amr_force_2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_spray_mpi_amr_reactive> force ${MPIEXEC_POSTFLAGS})
    set_tests_properties(spray_mpi_amr_force_2 PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 120 PROCESSORS 2)
    foreach(les_ranks 2 4)
      add_test(NAME spray_mpi_amr_les_coupled_${les_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${les_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_reactive> les ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_les_coupled_${les_ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 120 PROCESSORS ${les_ranks})
    endforeach()
    target_link_libraries(test_spray_mpi_amr_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_reactive PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_impact spray/test_spray_mpi_amr_impact.F90)
    target_link_libraries(test_spray_mpi_amr_impact PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_impact PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_boundaries spray/test_spray_mpi_amr_boundaries.F90)
    target_link_libraries(test_spray_mpi_amr_boundaries PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_boundaries PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_step spray/test_spray_mpi_amr_step.F90)
    target_link_libraries(test_spray_mpi_amr_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_step PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_gas spray/test_spray_mpi_amr_gas.F90)
    target_link_libraries(test_spray_mpi_amr_gas PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_gas PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_deposition spray/test_spray_mpi_amr_deposition.F90)
    target_link_libraries(test_spray_mpi_amr_deposition PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_deposition PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_amr_sampling spray/test_spray_mpi_amr_sampling.F90)
    target_link_libraries(test_spray_mpi_amr_sampling PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_sampling PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_impact spray/test_spray_mpi_eb_impact.F90)
    target_link_libraries(test_spray_mpi_eb_impact PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_eb_impact PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_impact_reactive spray/test_spray_mpi_eb_impact_reactive.F90)
    target_link_libraries(test_spray_mpi_eb_impact_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_eb_impact_reactive PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_film_reactive spray/test_spray_mpi_eb_film_reactive.F90)
    target_link_libraries(test_spray_mpi_eb_film_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_eb_film_reactive PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_reactive spray/test_spray_mpi_eb_reactive.F90)
    target_link_libraries(test_spray_mpi_eb_reactive PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_checkpoint spray/test_spray_mpi_eb_checkpoint.F90)
    target_link_libraries(test_spray_mpi_eb_checkpoint PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_balance spray/test_spray_mpi_eb_balance.F90)
    target_link_libraries(test_spray_mpi_eb_balance PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_remap spray/test_spray_mpi_eb_remap.F90)
    target_link_libraries(test_spray_mpi_eb_remap PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_impact_step spray/test_spray_mpi_eb_impact_step.F90)
    target_link_libraries(test_spray_mpi_eb_impact_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_eb_impact_step PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_spray_mpi_eb_step spray/test_spray_mpi_eb_step.F90)
    target_link_libraries(test_spray_mpi_eb_step PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_deposit spray/test_spray_mpi_eb_deposit.F90)
    target_link_libraries(test_spray_mpi_eb_deposit PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_sample spray/test_spray_mpi_eb_sample.F90)
    target_link_libraries(test_spray_mpi_eb_sample PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_eb_tree spray/test_spray_mpi_eb_tree.F90)
    target_link_libraries(test_spray_mpi_eb_tree PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_amr_tree spray/test_spray_mpi_amr_tree.F90)
    target_link_libraries(test_spray_mpi_amr_tree PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(test_spray_mpi_amr_tree PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
    add_executable(test_les_mpi_dynamic spray/test_les_mpi_dynamic.F90)
    target_link_libraries(test_les_mpi_dynamic PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    foreach(dynamic_ranks 2 4)
      add_test(NAME spray_les_mpi_dynamic_${dynamic_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${dynamic_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_les_mpi_dynamic> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_les_mpi_dynamic_${dynamic_ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 60 PROCESSORS ${dynamic_ranks})
    endforeach()
    add_test(NAME spray_mpi_difmag_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2 --ppm --difmag .037
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_difmag_restart")
    set_tests_properties(spray_mpi_difmag_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
    add_test(NAME spray_mpi_ppm_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2 --ppm
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_ppm_restart")
    set_tests_properties(spray_mpi_ppm_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
    add_test(NAME spray_mpi_dynamic_mach_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable $<TARGET_FILE:pelef_mpi_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
        --ranks 2 --les-model dynamic --mach-inlet --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_dynamic_mach_restart")
    set_tests_properties(spray_mpi_dynamic_mach_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
    add_executable(test_les_mpi_boundary spray/test_les_mpi_boundary.F90)
    target_link_libraries(test_les_mpi_boundary PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    foreach(boundary_ranks 2 4)
      add_test(NAME spray_les_mpi_boundary_${boundary_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${boundary_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_les_mpi_boundary> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_les_mpi_boundary_${boundary_ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 60 PROCESSORS ${boundary_ranks})
    endforeach()
    add_executable(test_les_mpi_filter spray/test_les_mpi_filter.F90)
    target_link_libraries(test_les_mpi_filter PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    foreach(filter_ranks 2 4)
      add_test(NAME spray_les_mpi_filter_${filter_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${filter_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_les_mpi_filter> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_les_mpi_filter_${filter_ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 60 PROCESSORS ${filter_ranks})
    endforeach()
    add_executable(test_les_mpi_models spray/test_les_mpi_models.F90)
    target_link_libraries(test_les_mpi_models PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    foreach(les_ranks 2 4)
      add_test(NAME spray_les_mpi_models_${les_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${les_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_les_mpi_models> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_les_mpi_models_${les_ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 60 PROCESSORS ${les_ranks})
    endforeach()
    add_executable(test_spray_rate_file spray/test_spray_rate_file.F90)
    target_link_libraries(test_spray_rate_file PRIVATE pelef_spray_runtime)
    add_test(NAME spray_rate_file COMMAND test_spray_rate_file "${CMAKE_CURRENT_BINARY_DIR}/rate_file_fixture.dat")
    set_tests_properties(spray_rate_file PROPERTIES LABELS "spray" TIMEOUT 20)
    add_executable(test_spray_mpi spray/test_spray_mpi.F90)
    target_link_libraries(test_spray_mpi PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    add_executable(test_spray_mpi_ppm spray/test_spray_mpi_ppm.F90)
    target_link_libraries(test_spray_mpi_ppm PRIVATE pelef_spray_runtime MPI::MPI_Fortran)
    add_test(NAME spray_mpi_ppm_4 COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_spray_mpi_ppm> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(spray_mpi_ppm_4 PROPERTIES LABELS "spray;mpi" TIMEOUT 120 PROCESSORS 4
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")

    foreach(ranks 2 4)
      add_test(NAME spray_mpi_amr_balance_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_balance> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_balance_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_remap_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_remap> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_remap_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_film_reactive_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_checkpoint_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film_checkpoint> ${MPIEXEC_POSTFLAGS}
          "${CMAKE_CURRENT_BINARY_DIR}/film_checkpoint_${ranks}")
      set_tests_properties(spray_mpi_amr_film_checkpoint_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_film_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_film_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_film_step_${ranks} PROPERTIES LABELS "spray;mpi;amr;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_film_step_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_film_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_film> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_film_${ranks} PROPERTIES LABELS "spray;mpi;amr;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_film_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_film> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_film_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_reactive_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 120 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_impact_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_impact> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_impact_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_boundaries_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_boundaries> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_boundaries_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_step_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_gas_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_gas> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_gas_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_deposition_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_deposition> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_deposition_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_sampling_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_sampling> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_sampling_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_impact_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_impact> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_impact_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_impact_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_impact_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_impact_reactive_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_film_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_film_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_film_reactive_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_reactive_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_reactive> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_reactive_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_balance_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_balance> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_balance_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_remap_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_remap> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_remap_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_impact_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_impact_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_impact_step_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_wall_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_step> ${MPIEXEC_POSTFLAGS} wall)
      set_tests_properties(spray_mpi_eb_wall_step_${ranks} PROPERTIES TIMEOUT 120 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_periodic_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_step> ${MPIEXEC_POSTFLAGS} periodic)
      set_tests_properties(spray_mpi_eb_periodic_step_${ranks} PROPERTIES
        LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_step_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_step> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_step_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_deposit_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_deposit> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_deposit_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      foreach(exchange IN ITEMS sample deposit)
        add_test(NAME spray_mpi_eb_periodic_${exchange}_${ranks}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_spray_mpi_eb_${exchange}> ${MPIEXEC_POSTFLAGS} periodic)
        set_tests_properties(spray_mpi_eb_periodic_${exchange}_${ranks} PROPERTIES
          LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      endforeach()
      add_test(NAME spray_mpi_eb_sample_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_sample> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_sample_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_eb_tree_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_eb_tree> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_eb_tree_${ranks} PROPERTIES LABELS "spray;mpi;eb" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_amr_tree_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi_amr_tree> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_amr_tree_${ranks} PROPERTIES LABELS "spray;mpi;amr" TIMEOUT 60 PROCESSORS ${ranks})
      add_test(NAME spray_mpi_application_${ranks}
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
          --executable $<TARGET_FILE:pelef_mpi_spray_3d> --mpiexec "${MPIEXEC_EXECUTABLE}"
          --ranks ${ranks} --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
          --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_application_${ranks}")
      set_tests_properties(spray_mpi_application_${ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
      add_test(NAME spray_mpi_migration_${ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_spray_mpi> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_mpi_migration_${ranks} PROPERTIES LABELS "spray;mpi" TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
  endif()
  add_test(NAME spray_acceleration_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --acceleration -9810
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_acceleration")
  set_tests_properties(spray_acceleration_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  foreach(wall_sgs 0 .17)
    add_test(NAME spray_thermal_wall_${wall_sgs} COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --thermal-wall --sgs ${wall_sgs}
      --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_thermal_wall_${wall_sgs}")
    set_tests_properties(spray_thermal_wall_${wall_sgs} PROPERTIES LABELS "spray" TIMEOUT 600)
  endforeach()
  add_test(NAME spray_inflow_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --inflow
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_inflow_restart")
  set_tests_properties(spray_inflow_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  add_test(NAME spray_split_restart COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_spray_etab.py"
    --split --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_split_restart")
  set_tests_properties(spray_split_restart PROPERTIES LABELS "spray" TIMEOUT 900)
  foreach(khrt_variant plain split)
    set(khrt_args)
    if(khrt_variant STREQUAL "split")
      set(khrt_args --split)
    endif()
    add_test(NAME spray_khrt_${khrt_variant}_restart COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_etab.py" --model khrt ${khrt_args}
      --executable "$<TARGET_FILE:pelef_spray_3d>"
      --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_khrt_${khrt_variant}_restart")
    set_tests_properties(spray_khrt_${khrt_variant}_restart PROPERTIES LABELS "spray" TIMEOUT 900)
  endforeach()
  add_test(NAME spray_etab_restart COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/check_spray_etab.py"
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_etab_restart")
  set_tests_properties(spray_etab_restart PROPERTIES LABELS "spray" TIMEOUT 900)
  add_test(NAME spray_physical_boundary_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py"
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_physical_boundaries")
  set_tests_properties(spray_physical_boundary_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  add_test(NAME spray_sgs_boundary_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_sgs_boundaries")
  set_tests_properties(spray_sgs_boundary_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  if(PELEF_ENABLE_MPI)
    foreach(les_model wale vreman dynamic)
      add_test(NAME spray_mpi_${les_model}_restart
        COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
          --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2
          --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
          --les-model ${les_model} --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_${les_model}_restart")
      set_tests_properties(spray_mpi_${les_model}_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
    endforeach()
  endif()
  foreach(les_model wale vreman dynamic)
    add_test(NAME spray_${les_model}_boundary_restart COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17 --les-model ${les_model}
      --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_${les_model}_boundaries")
    set_tests_properties(spray_${les_model}_boundary_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  endforeach()
  add_test(NAME spray_difmag_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --ppm --inflow --mach-inlet --difmag .05
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_difmag_restart")
  set_tests_properties(spray_difmag_restart PROPERTIES LABELS "spray" TIMEOUT 300)
  add_test(NAME spray_ppm_mach_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --ppm --inflow --mach-inlet
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_ppm_mach_restart")
  set_tests_properties(spray_ppm_mach_restart PROPERTIES LABELS "spray" TIMEOUT 300)
  add_test(NAME spray_dynamic_mach_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --inflow --mach-inlet --sgs .17
    --les-model dynamic --dynamic-physical
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_dynamic_mach_restart")
  set_tests_properties(spray_dynamic_mach_restart PROPERTIES LABELS "spray;les" TIMEOUT 600)
  add_test(NAME spray_dynamic_wall_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17 --les-model dynamic --dynamic-physical --thermal-wall
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_dynamic_wall_restart")
  set_tests_properties(spray_dynamic_wall_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  add_test(NAME spray_dynamic_harmonic_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17 --les-model dynamic --les-harmonic
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_dynamic_harmonic_restart")
  set_tests_properties(spray_dynamic_harmonic_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  if(PELEF_ENABLE_MPI)
    add_test(NAME spray_mpi_dynamic_harmonic_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --les-model dynamic --les-harmonic --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_dynamic_harmonic_restart")
    set_tests_properties(spray_mpi_dynamic_harmonic_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
    add_test(NAME spray_mpi_dynamic_wall_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --les-model dynamic --dynamic-walls --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_dynamic_wall_restart")
    set_tests_properties(spray_mpi_dynamic_wall_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
  endif()
  add_test(NAME spray_isotropic_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17 --les-model wale --les-isotropic
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_isotropic_restart")
  set_tests_properties(spray_isotropic_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  if(PELEF_ENABLE_MPI)
    add_test(NAME spray_mpi_isotropic_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --les-model wale --les-isotropic --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_isotropic_restart")
    set_tests_properties(spray_mpi_isotropic_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
  endif()
  add_test(NAME spray_pressure_outlet_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --sgs .17 --les-model wale --les-isotropic --pressure-outlet
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_pressure_outlet_restart")
  set_tests_properties(spray_pressure_outlet_restart PROPERTIES LABELS "spray" TIMEOUT 600)
  if(PELEF_ENABLE_MPI)
    add_test(NAME spray_mpi_pressure_outlet_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --les-model wale --les-isotropic --pressure-outlet
        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_pressure_outlet_restart")
    set_tests_properties(spray_mpi_pressure_outlet_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
  endif()
  if(PELEF_ENABLE_MPI)
    add_executable(test_les_mpi_amr_level spray/test_les_mpi_amr_level.F90)
    target_link_libraries(test_les_mpi_amr_level PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_test(NAME spray_les_mpi_amr_level_2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_les_mpi_amr_level> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(spray_les_mpi_amr_level_2 PROPERTIES LABELS "spray;amr;mpi" TIMEOUT 120 PROCESSORS 2)
    add_executable(test_les_mpi_amr_subcycle spray/test_les_mpi_amr_subcycle.F90)
    target_link_libraries(test_les_mpi_amr_subcycle PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    add_test(NAME spray_les_mpi_amr_subcycle_2
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_les_mpi_amr_subcycle> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(spray_les_mpi_amr_subcycle_2 PROPERTIES LABELS "spray;amr;mpi" TIMEOUT 180 PROCESSORS 2)
    add_executable(test_les_mpi_amr_tree spray/test_les_mpi_amr_tree.F90)
    target_link_libraries(test_les_mpi_amr_tree PRIVATE pelef_spray_tree_runtime MPI::MPI_Fortran)
    foreach(les_ranks 2 4)
      add_test(NAME spray_les_mpi_amr_tree_${les_ranks}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${les_ranks} ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_les_mpi_amr_tree> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(spray_les_mpi_amr_tree_${les_ranks} PROPERTIES LABELS "spray;amr;mpi" TIMEOUT 180 PROCESSORS ${les_ranks})
    endforeach()
  endif()
  add_executable(test_les_amr_tree spray/test_les_amr_tree.F90)
  add_executable(test_les_amr_subcycle spray/test_les_amr_subcycle.F90)
  target_link_libraries(test_les_amr_subcycle PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_les_amr_subcycle COMMAND test_les_amr_subcycle)
  set_tests_properties(spray_les_amr_subcycle PROPERTIES LABELS "spray;amr" TIMEOUT 180)
  add_executable(test_les_dynamic_boundaries spray/test_les_dynamic_boundaries.F90)
  target_link_libraries(test_les_dynamic_boundaries PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_les_dynamic_boundaries COMMAND test_les_dynamic_boundaries)
  set_tests_properties(spray_les_dynamic_boundaries PROPERTIES LABELS "spray;les" TIMEOUT 120)
  add_executable(test_les_amr_level spray/test_les_amr_level.F90)
  target_link_libraries(test_les_amr_level PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_les_amr_level COMMAND test_les_amr_level)
  set_tests_properties(spray_les_amr_level PROPERTIES LABELS "spray;amr" TIMEOUT 60)
  target_link_libraries(test_les_amr_tree PRIVATE pelef_spray_tree_runtime)
  add_test(NAME spray_les_amr_tree COMMAND test_les_amr_tree)
  set_tests_properties(spray_les_amr_tree PROPERTIES LABELS "spray;amr" TIMEOUT 120)
  add_executable(test_spray_core spray/test_spray_core.F90)
  target_link_libraries(test_spray_core PRIVATE pelef_spray_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_spray_core PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  foreach(spray_case partition mixture_wall_crossing mixture_wall_rates mixture_wall_exchange mixture_wall_event mixture_terminal mixture_rates mixture_exchange mixture film_evolution wall_event wall_state wall_impact heat_capacity port_states ports surface_cap viscosity density antoine split khrt etab boundaries drag heat evaporation terminal abramzon film_binary injection size_distribution stencil coupling rollback substeps zero invalid
      vapor_cp vapor_skin liquid_mw many_to_one pelemp_cap watson_properties les_boundary les_face_fluxes les_pressure_outlet les_isotropic les_face_mean les_algebra les_filter les_filter_approximations les_dynamic_step les_dynamic_flux les_dynamic_ghosted les_dynamic_field les_dynamic_algebra les_alternative_algebra les_wale les_vreman les_grid les_zero les_reject)
    add_test(NAME spray_${spray_case} COMMAND test_spray_core ${spray_case})
    set_tests_properties(spray_${spray_case} PROPERTIES LABELS "spray" TIMEOUT 60)
  endforeach()
  add_executable(test_spray_integration spray/test_spray_integration.F90)
  target_link_libraries(test_spray_integration PRIVATE pelef_spray_runtime)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_spray_integration PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
  endif()
  foreach(spray_case native coupled coupled_wale coupled_vreman coupled_dynamic rollback checkpoint acceleration materials wall_coupling mixture_wall_coupling mixture_film_grid mixture_grid mixture_reactive)
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

if(PELEF_ENABLE_SPRAY AND PELEF_ENABLE_MPI)
  add_test(NAME spray_mpi_rate_files_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
    --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2
    --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
    --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_rate_files_application" --rate-files)
  set_tests_properties(spray_mpi_rate_files_application PROPERTIES LABELS "spray;mpi;application" TIMEOUT 180 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_rate_files_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_nozzle.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_rate_files_application" --rate --rate-files)
  set_tests_properties(spray_mpi_eb_rate_files_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_properties_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_nozzle.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_properties_application" --rate --properties)
  set_tests_properties(spray_mpi_eb_properties_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 180 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_rate_velocity_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_nozzle.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_rate_velocity_application" --rate)
  set_tests_properties(spray_mpi_eb_rate_velocity_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 180 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_nozzle_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_nozzle.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_nozzle_application")
  set_tests_properties(spray_mpi_eb_nozzle_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 180 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_split_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_split_application" --split)
  set_tests_properties(spray_mpi_eb_split_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_khrt_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_khrt_application" --khrt)
  set_tests_properties(spray_mpi_eb_khrt_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_reflected_impact_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_reflected_impact_application"
    --reflect-domain)
  set_tests_properties(spray_mpi_eb_reflected_impact_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_periodic_impact_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_periodic_impact_application"
    --periodic)
  set_tests_properties(spray_mpi_eb_periodic_impact_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_wall_khrt_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_wall_khrt_application"
    --khrt --domain-wall)
  set_tests_properties(spray_mpi_eb_wall_khrt_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_periodic_khrt_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_periodic_khrt_application"
    --khrt --periodic)
  set_tests_properties(spray_mpi_eb_periodic_khrt_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_pcm_wall_khrt_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_pcm_wall_khrt_application"
    --khrt --domain-wall --pcm)
  set_tests_properties(spray_mpi_eb_pcm_wall_khrt_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_force_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_force_application" --force)
  set_tests_properties(spray_mpi_eb_force_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_specular_outflow_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_specular_outflow_application"
    --specular-outflow)
  set_tests_properties(spray_mpi_eb_specular_outflow_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_outflow_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_outflow_application" --outflow)
  set_tests_properties(spray_mpi_eb_outflow_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_impact_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_impact.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_impact_application")
  set_tests_properties(spray_mpi_eb_impact_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_application")
  set_tests_properties(spray_mpi_eb_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 180 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_schedule COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_schedule" --scheduled)
  set_tests_properties(spray_mpi_eb_schedule PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 180 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_restart" --scheduled --restart)
  set_tests_properties(spray_mpi_eb_restart PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_mach COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_mach" --restart --mach)
  set_tests_properties(spray_mpi_eb_mach PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_difmag COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_difmag" --restart --mach --difmag .037)
  set_tests_properties(spray_mpi_eb_difmag PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_pelemp_cap COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_pelemp_cap" --restart --pelemp-cap)
  set_tests_properties(spray_mpi_eb_pelemp_cap PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_ports COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_ports" --scheduled --restart --regrid --ports)
  set_tests_properties(spray_mpi_eb_ports PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_profiles COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_profiles" --scheduled --restart --regrid --ports --profiles)
  set_tests_properties(spray_mpi_eb_profiles PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_regrid_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_regrid_application" --scheduled --restart --regrid)
  set_tests_properties(spray_mpi_eb_regrid_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
  add_test(NAME spray_mpi_eb_rebalance_application COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
    --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_rebalance_application" --scheduled --restart --rebalance)
  set_tests_properties(spray_mpi_eb_rebalance_application PROPERTIES LABELS "spray;mpi;eb;application" TIMEOUT 240 PROCESSORS 4)
endif()

if(PELEF_ENABLE_SPRAY AND PELEF_ENABLE_MPI)
  add_test(NAME spray_mpi_eb_particle_checkpoint COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_particle_checkpoint.py"
    --exe "$<TARGET_FILE:test_spray_mpi_eb_checkpoint>" --mpiexec "${MPIEXEC_EXECUTABLE}"
    --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_particle_checkpoint")
  set_tests_properties(spray_mpi_eb_particle_checkpoint PROPERTIES LABELS "spray;mpi;eb;restart" TIMEOUT 150 PROCESSORS 4)
endif()

if(TARGET test_spray_mpi_ppm)
  add_test(NAME spray_mpi_plm_viscosity_4 COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4
    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_spray_mpi_ppm> ${MPIEXEC_POSTFLAGS} plm)
  set_tests_properties(spray_mpi_plm_viscosity_4 PROPERTIES LABELS "spray;mpi" TIMEOUT 120 PROCESSORS 4
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
endif()

if(TARGET pelef_spray_3d)
  add_test(NAME spray_plm_difmag_restart COMMAND "${Python3_EXECUTABLE}"
    "${PROJECT_SOURCE_DIR}/tools/check_spray_boundaries.py" --inflow --mach-inlet --difmag .05
    --executable "$<TARGET_FILE:pelef_spray_3d>" --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_plm_difmag_restart")
  set_tests_properties(spray_plm_difmag_restart PROPERTIES LABELS "spray" TIMEOUT 300)
endif()
if(TARGET pelef_mpi_spray_3d)
    add_test(NAME spray_mpi_plm_difmag_restart
      COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/spray/test_mpi_spray_application.py"
        --executable "$<TARGET_FILE:pelef_mpi_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}" --ranks 2 --difmag .037
        --input "${PROJECT_SOURCE_DIR}/cases/mpi_spray_3d/periodic_mixture.nml"
        --work-dir "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_plm_difmag_restart")
    set_tests_properties(spray_mpi_plm_difmag_restart PROPERTIES LABELS "spray;mpi" TIMEOUT 200 PROCESSORS 4)
endif()

if(TARGET test_spray_mpi_eb_reactive)
  foreach(ctu_mode IN ITEMS ctu ctu-subcycle ctu-inflow ctu-inflow-subcycle ctu-periodic ctu-periodic-subcycle ctu-wall ctu-wall-subcycle
      ctu-isothermal ctu-isothermal-subcycle)
    add_test(NAME spray_mpi_eb_reactive_${ctu_mode}
      COMMAND "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 4 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_spray_mpi_eb_reactive> ${MPIEXEC_POSTFLAGS} ${ctu_mode})
    set_tests_properties(spray_mpi_eb_reactive_${ctu_mode} PROPERTIES
      TIMEOUT 180 PROCESSORS 4 ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  endforeach()
endif()

if(PELEF_ENABLE_SPRAY AND PELEF_ENABLE_MPI)
  add_test(NAME spray_mpi_eb_ctu_restart
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
      --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_ctu_restart"
      --ctu --restart)
  set_tests_properties(spray_mpi_eb_ctu_restart PROPERTIES TIMEOUT 600 PROCESSORS 4 LABELS "spray;mpi;eb")
  add_test(NAME spray_mpi_eb_isothermal_restart
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_mpi_eb_spray.py"
      --exe "$<TARGET_FILE:pelef_mpi_eb_spray_3d>" --mpiexec "${MPIEXEC_EXECUTABLE}"
      --source "${PROJECT_SOURCE_DIR}" --work "${CMAKE_CURRENT_BINARY_DIR}/spray_mpi_eb_isothermal_restart"
      --isothermal --restart --difmag .037)
  set_tests_properties(spray_mpi_eb_isothermal_restart PROPERTIES TIMEOUT 600 PROCESSORS 4 LABELS "spray;mpi;eb")
endif()

if(TARGET test_spray_eb_reactive)
  add_test(NAME spray_eb_reactive_ctu COMMAND test_spray_eb_reactive ctu)
endif()

if(PELEF_ENABLE_SPRAY)
  add_test(NAME spray_eb_ctu_restart
    COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tests/application/check_eb_spray.py"
      --exe "$<TARGET_FILE:pelef_eb_spray_3d>" --source "${PROJECT_SOURCE_DIR}"
      --work "${CMAKE_CURRENT_BINARY_DIR}/spray_eb_ctu_restart" --ctu --restart)
  set_tests_properties(spray_eb_ctu_restart PROPERTIES TIMEOUT 180 LABELS "spray;eb")
endif()
