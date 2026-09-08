if(PELEF_ENABLE_MPI)
  add_library(
    pelef_mpi_support
    src/parallel/mpi_domain_1d_mod.F90
    src/parallel/mpi_amr_reactive_3d_mod.F90
    src/parallel/mpi_amr_sparse_reactive_3d_mod.F90
    src/parallel/mpi_amr_patch_1d_mod.F90
    src/parallel/mpi_amr_sparse_patch_1d_mod.F90
    src/parallel/mpi_amr_eb_patch_tree_2d_mod.F90
    src/parallel/mpi_amr_eb_patch_tree_io_2d_mod.F90
    src/driver/mpi_reactive_eb_patch_tree_2d_application_mod.F90
    src/parallel/mpi_amr_eb_patch_2d_mod.F90
    src/parallel/mpi_amr_eb_io_2d_mod.F90
    src/parallel/mpi_reactive_transport_1d_mod.F90
    src/parallel/mpi_reactive_1d_mod.F90
    src/driver/mpi_reactive_1d_application_mod.F90
    src/driver/mpi_amr_reactive_1d_application_mod.F90
    src/driver/mpi_amr_reactive_3d_application_mod.F90
  )
  target_link_libraries(pelef_mpi_support PUBLIC pelef_core MPI::MPI_Fortran)
  set_target_properties(pelef_mpi_support PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mpi_modules")
  target_include_directories(pelef_mpi_support PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mpi_modules>")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(pelef_mpi_support PRIVATE -Wall -Wextra -Wpedantic -Wimplicit-interface -Wconversion-extra)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      target_compile_options(pelef_mpi_support PRIVATE -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow)
    endif()
  endif()
  add_executable(pelef_mpi_1d app/pelef_mpi_1d.F90)
  target_link_libraries(pelef_mpi_1d PRIVATE pelef_mpi_support)
  add_executable(pelef_mpi_multispecies_1d app/pelef_mpi_multispecies_1d.F90)
  target_link_libraries(pelef_mpi_multispecies_1d PRIVATE pelef_mpi_support)
  add_executable(pelef_mpi_h2o2_batch app/pelef_mpi_h2o2_batch.F90)
  target_link_libraries(pelef_mpi_h2o2_batch PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_reactive_transport_1d
    app/pelef_mpi_reactive_transport_1d.F90
  )
  target_link_libraries(
    pelef_mpi_reactive_transport_1d PRIVATE pelef_mpi_support
  )
  add_executable(pelef_mpi_reactive_1d app/pelef_mpi_reactive_1d.F90)
  target_link_libraries(pelef_mpi_reactive_1d PRIVATE pelef_mpi_support)
  pelef_require_coherent_mpi_runtime(pelef_mpi_reactive_1d)
  add_executable(
    pelef_mpi_amr_patch_1d
    app/pelef_mpi_amr_patch_1d.F90
  )
  target_link_libraries(pelef_mpi_amr_patch_1d PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_eb_amr_patch_2d
    app/pelef_mpi_eb_amr_patch_2d.F90
  )
  target_link_libraries(
    pelef_mpi_eb_amr_patch_2d PRIVATE pelef_mpi_support
  )
  add_executable(
    pelef_mpi_amr_eb_patch_tree_2d
    app/pelef_mpi_amr_eb_patch_tree_2d.F90
  )
  target_link_libraries(
    pelef_mpi_amr_eb_patch_tree_2d PRIVATE pelef_mpi_support
  )
  add_executable(
    pelef_mpi_reactive_eb_patch_tree_2d
    app/pelef_mpi_reactive_eb_patch_tree_2d.F90
  )
  target_link_libraries(
    pelef_mpi_reactive_eb_patch_tree_2d PRIVATE pelef_mpi_support
  )
  add_executable(
    pelef_mpi_amr_reactive_1d
    app/pelef_mpi_amr_reactive_1d.F90
  )
  target_link_libraries(pelef_mpi_amr_reactive_1d PRIVATE pelef_mpi_support)
  add_executable(
    pelef_mpi_amr_reactive_3d
    app/pelef_mpi_amr_reactive_3d.F90
  )
  target_link_libraries(pelef_mpi_amr_reactive_3d PRIVATE pelef_mpi_support)
  if(PELEF_ENABLE_TESTS)
    add_executable(
      test_mpi_reactive_integrator_policy
      tests/mpi/test_mpi_reactive_integrator_policy.F90
    )
    target_link_libraries(
      test_mpi_reactive_integrator_policy PRIVATE pelef_mpi_support
    )
    foreach(mpi_integrator_policy_ranks IN ITEMS 1 2 4)
      add_test(
        NAME
          unit_mpi_reactive_integrator_policy_np${mpi_integrator_policy_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_integrator_policy_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_reactive_integrator_policy>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_reactive_integrator_policy_np${mpi_integrator_policy_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_integrator_policy_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_reactive_integrator_policy: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_sparse_integrator_policy
      tests/mpi/test_mpi_amr_sparse_integrator_policy.F90
    )
    target_link_libraries(
      test_mpi_amr_sparse_integrator_policy PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_integrator_policy_ranks IN ITEMS 1 2 4)
      add_test(
        NAME
          unit_mpi_amr_sparse_integrator_policy_np${mpi_amr_integrator_policy_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_integrator_policy_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_sparse_integrator_policy>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_amr_sparse_integrator_policy_np${mpi_amr_integrator_policy_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_integrator_policy_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_sparse_integrator_policy: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_eb_sparse_integrator_policy
      tests/mpi/test_mpi_amr_eb_sparse_integrator_policy.F90
    )
    target_link_libraries(
      test_mpi_amr_eb_sparse_integrator_policy PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_eb_integrator_policy_ranks IN ITEMS 1 2 4)
      add_test(
        NAME
          unit_mpi_amr_eb_sparse_integrator_policy_np${mpi_amr_eb_integrator_policy_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_eb_integrator_policy_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_eb_sparse_integrator_policy>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_amr_eb_sparse_integrator_policy_np${mpi_amr_eb_integrator_policy_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_eb_integrator_policy_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_eb_sparse_integrator_policy: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_selected_context_mismatch
      tests/mpi/test_mpi_amr_selected_context_mismatch.F90
    )
    target_link_libraries(
      test_mpi_amr_selected_context_mismatch PRIVATE pelef_mpi_support
    )
    set(
      MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_amr_selected_context_mismatch"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_selected_context_mismatch
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
        --launcher "${MPIEXEC_EXECUTABLE}"
        "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --ranks 2
        --executable $<TARGET_FILE:test_mpi_amr_selected_context_mismatch>
        --expected "Selected sparse MPI AMR rank context mismatch"
        --forbidden-output
          "${MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}/mpi_amr_context_mismatch.csv"
    )
    set_tests_properties(
      regression_mpi_amr_selected_context_mismatch
      PROPERTIES
        WORKING_DIRECTORY
          "${MPI_AMR_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        TIMEOUT 60
    )

    add_executable(
      test_mpi_amr_selected_context_3d
      tests/mpi/test_mpi_amr_selected_context_3d.F90
    )
    target_link_libraries(
      test_mpi_amr_selected_context_3d PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_selected_context_3d_ranks IN ITEMS 2 4)
      add_test(
        NAME
          unit_mpi_amr_selected_context_3d_np${mpi_amr_selected_context_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_selected_context_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_selected_context_3d>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        unit_mpi_amr_selected_context_3d_np${mpi_amr_selected_context_3d_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_selected_context_3d_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_selected_context_3d: PASS"
          TIMEOUT 60
      )
    endforeach()

    add_executable(
      test_mpi_amr_eb_selected_context_mismatch
      tests/mpi/test_mpi_amr_eb_selected_context_mismatch.F90
    )
    target_link_libraries(
      test_mpi_amr_eb_selected_context_mismatch PRIVATE pelef_mpi_support
    )
    set(
      MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/tests/mpi_amr_eb_selected_context_mismatch"
    )
    file(
      MAKE_DIRECTORY
      "${MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}"
    )
    add_test(
      NAME regression_mpi_amr_eb_selected_context_mismatch
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
        --launcher "${MPIEXEC_EXECUTABLE}"
        "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
        --ranks 2
        --executable $<TARGET_FILE:test_mpi_amr_eb_selected_context_mismatch>
        --expected "Selected sparse MPI EB rank context mismatch"
        --forbidden-output
          "${MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}/mpi_amr_eb_context_mismatch.csv"
    )
    set_tests_properties(
      regression_mpi_amr_eb_selected_context_mismatch
      PROPERTIES
        WORKING_DIRECTORY
          "${MPI_AMR_EB_SELECTED_CONTEXT_MISMATCH_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        TIMEOUT 60
    )

    add_executable(
      test_mpi_amr_eb_selected_checkpoint_presence
      tests/mpi/test_mpi_amr_eb_selected_checkpoint_presence.F90
    )
    target_link_libraries(
      test_mpi_amr_eb_selected_checkpoint_presence PRIVATE pelef_mpi_support
    )
    add_test(
      NAME unit_mpi_amr_eb_selected_checkpoint_presence_np2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_mpi_amr_eb_selected_checkpoint_presence>
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      unit_mpi_amr_eb_selected_checkpoint_presence_np2
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION
          "test_mpi_amr_eb_selected_checkpoint_presence: PASS"
        TIMEOUT 60
    )

    set(
      SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/tests/selected_mpi_reactive_1d_full"
    )
    file(MAKE_DIRECTORY "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_reactive_1d_fixed_np1_reference
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_reactive_1d>
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/fixed_np1.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_mpi_reactive_1d_fixed_np1_reference
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 1
        PASS_REGULAR_EXPRESSION "MPI ranks: 1"
        TIMEOUT 300
    )
    add_test(
      NAME regression_selected_mpi_reactive_1d_full_fixed_parity
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_reactive_1d.py"
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/fixed_np1.csv"
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np1.csv"
    )
    set_tests_properties(
      regression_selected_mpi_reactive_1d_full_fixed_parity
      PROPERTIES
        DEPENDS
          "regression_mpi_reactive_1d_fixed_np1_reference;regression_selected_mpi_reactive_1d_full_np1"
    )

    add_executable(
      test_mpi_amr_reactive_3d
      tests/mpi/test_mpi_amr_reactive_3d.F90
    )
    target_link_libraries(test_mpi_amr_reactive_3d PRIVATE pelef_mpi_support)
    foreach(mpi_amr_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_reactive_3d> ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_3d_ranks}
          PASS_REGULAR_EXPRESSION "test_mpi_amr_reactive_3d: PASS"
          TIMEOUT 300
      )
    endforeach()

    add_executable(
      test_mpi_amr_reactive_chemistry_3d
      tests/mpi/test_mpi_amr_reactive_chemistry_3d.F90
    )
    target_link_libraries(
      test_mpi_amr_reactive_chemistry_3d PRIVATE pelef_mpi_support
    )
    foreach(mpi_amr_chemistry_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_chemistry_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:test_mpi_amr_reactive_chemistry_3d>
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_3d_ranks}
        PROPERTIES
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_chemistry_3d_ranks}
          PASS_REGULAR_EXPRESSION
            "test_mpi_amr_reactive_chemistry_3d: PASS"
          TIMEOUT 300
      )
    endforeach()

    set(
      MPI_AMR_REACTIVE_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "static two-level 3D reactive AMR"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_3D_RUN_TESTS "")
    foreach(mpi_amr_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave.nml"
          "mpi_amr_np${mpi_amr_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_3d_ranks}
          PASS_REGULAR_EXPRESSION "distributed-slab static two-level 3D"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_3D_RUN_TESTS
        regression_mpi_amr_reactive_3d_np${mpi_amr_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np1_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np1_fine.csv"
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --time 1.0e-6
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_3d_np1
    )
    add_test(
      NAME regression_mpi_amr_reactive_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_np8"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_3d_serial_reference;${MPI_AMR_REACTIVE_3D_RUN_TESTS}"
    )

    set(
      MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_transport_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_transport_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "Completed transport fine substeps: 16"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_TRANSPORT_3D_RUN_TESTS "")
    foreach(mpi_amr_transport_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_transport_3d_np${mpi_amr_transport_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_transport_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport.nml"
          "mpi_amr_transport_np${mpi_amr_transport_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_transport_3d_np${mpi_amr_transport_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_transport_3d_ranks}
          PASS_REGULAR_EXPRESSION "Completed transport fine substeps: 16"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_TRANSPORT_3D_RUN_TESTS
        regression_mpi_amr_reactive_transport_3d_np${mpi_amr_transport_3d_ranks}
      )
    endforeach()
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_np8
      PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
    )
    add_test(
      NAME regression_mpi_amr_reactive_transport_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np8_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np8_fine.csv"
        --species H2 H O O2 OH H2O HO2 H2O2 AR N2
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --final-time 1.0e-9
        --activity-species H2
        --initial-mass-fraction 0.028502723329253223
        --minimum-change 1.0e-12
    )
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_transport_3d_np8
    )
    add_test(
      NAME regression_mpi_amr_reactive_transport_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/amr_hotspot_transport_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_TRANSPORT_3D_WORK_DIRECTORY}/mpi_amr_transport_np8"
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --expected-coarse-sha256
          $<IF:$<CONFIG:Debug>,5687cb97f834b759b3932bbd3f0b6e30fc5b48e6029e3f47cadfb3ec6a8d74fb,f0713f3a68f48b6537afaafd57209ea8ed4767f9c7d093fe2961bd3a1fa37fae>
        --expected-fine-sha256
          eaf132bbf89974d739c26a35f6b76137f9774438cfe9eebfc12436aa364cff58
    )
    set_tests_properties(
      regression_mpi_amr_reactive_transport_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_transport_3d_serial_reference;${MPI_AMR_REACTIVE_TRANSPORT_3D_RUN_TESTS}"
    )

    set(
      MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_chemistry_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_chemistry_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_chemistry.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_chemistry_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "Chemistry:  T"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_CHEMISTRY_3D_RUN_TESTS "")
    foreach(mpi_amr_chemistry_app_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_app_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_chemistry_app_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_chemistry.nml"
          "mpi_amr_chemistry_np${mpi_amr_chemistry_app_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_app_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_chemistry_app_ranks}
          PASS_REGULAR_EXPRESSION "Chemistry:  T"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_CHEMISTRY_3D_RUN_TESTS
        regression_mpi_amr_reactive_chemistry_3d_np${mpi_amr_chemistry_app_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_chemistry_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np8_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np8_fine.csv"
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --time 1.0e-9
    )
    set_tests_properties(
      regression_mpi_amr_reactive_chemistry_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_chemistry_3d_np8
    )
    add_test(
      NAME regression_mpi_amr_reactive_chemistry_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/amr_hotspot_chemistry_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_CHEMISTRY_3D_WORK_DIRECTORY}/mpi_amr_chemistry_np8"
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --expected-coarse-sha256
          3313a9a6eff5d6175625e4bbe3613a782cdd794314fd9fbdd53ddaa309fc9516
        --expected-fine-sha256
          4ff0785c3462a49648fd694d090eb23bffd5af05a43a33360cd187417bd8a31d
    )
    set_tests_properties(
      regression_mpi_amr_reactive_chemistry_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_chemistry_3d_serial_reference;${MPI_AMR_REACTIVE_CHEMISTRY_3D_RUN_TESTS}"
    )

    add_test(
      NAME regression_mpi_amr_reactive_3d_checkpoint_np2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_checkpoint.nml"
        "mpi_amr_checkpoint_np2"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_3D_RESTART_TESTS "")
    foreach(mpi_amr_3d_ranks IN ITEMS 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_3d_restart_np${mpi_amr_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_restart.nml"
          "mpi_amr_restart_np${mpi_amr_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_3d_restart_np${mpi_amr_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_3d_ranks}
          PASS_REGULAR_EXPRESSION "Restored coarse steps: 4"
          DEPENDS regression_mpi_amr_reactive_3d_checkpoint_np2
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_3D_RESTART_TESTS
        regression_mpi_amr_reactive_3d_restart_np${mpi_amr_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_3d_restart_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_restart_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/mpi_amr_restart_np8"
        --checkpoint
          "${MPI_AMR_REACTIVE_3D_WORK_DIRECTORY}/amr_reactive_3d.chk"
        --final-time 1.0e-6
    )
    set_tests_properties(
      regression_mpi_amr_reactive_3d_restart_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_3d_serial_reference;${MPI_AMR_REACTIVE_3D_RESTART_TESTS}"
    )

    set(
      MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY
      "${CMAKE_CURRENT_BINARY_DIR}/mpi_amr_reactive_plm_3d"
    )
    file(MAKE_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}")
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_serial_reference
      COMMAND
        $<TARGET_FILE:pelef_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm.nml"
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_serial_reference
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
        PASS_REGULAR_EXPRESSION "Reconstruction: characteristic_plm"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_PLM_3D_RUN_TESTS "")
    foreach(mpi_amr_plm_3d_ranks IN ITEMS 1 2 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_plm_3d_np${mpi_amr_plm_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_plm_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm.nml"
          "mpi_amr_plm_np${mpi_amr_plm_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_plm_3d_np${mpi_amr_plm_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_plm_3d_ranks}
          PASS_REGULAR_EXPRESSION "rank-local sparse distributed-slab"
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_PLM_3D_RUN_TESTS
        regression_mpi_amr_reactive_plm_3d_np${mpi_amr_plm_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_check
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_amr_reactive_3d.py"
        --coarse
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np1_coarse.csv"
        --fine
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np1_fine.csv"
        --nx 8 --ny 8 --nz 8
        --i-lower 3 --i-upper 6
        --j-lower 3 --j-upper 6
        --k-lower 3 --k-upper 6
        --ratio 2 --time 1.0e-6
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_check
      PROPERTIES DEPENDS regression_mpi_amr_reactive_plm_3d_np1
    )
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np1"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np2"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_np8"
        --expected-coarse-sha256
          864e9cabe8f5c2bee19f97e9279bdada7647305171955c6c767cc83e4db15c21
        --expected-fine-sha256
          8f764e38491e9cf62f8c0f17eac2470c4705f7d80bfe7e686db4ba54de2f58d7
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_plm_3d_serial_reference;${MPI_AMR_REACTIVE_PLM_3D_RUN_TESTS}"
    )
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_checkpoint_np2
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm_checkpoint.nml"
        "mpi_amr_plm_checkpoint_np2"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        TIMEOUT 300
    )
    set(MPI_AMR_REACTIVE_PLM_3D_RESTART_TESTS "")
    foreach(mpi_amr_plm_3d_ranks IN ITEMS 4 8)
      add_test(
        NAME regression_mpi_amr_reactive_plm_3d_restart_np${mpi_amr_plm_3d_ranks}
        COMMAND
          "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          "${mpi_amr_plm_3d_ranks}" ${MPIEXEC_PREFLAGS}
          $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_entropy_wave_plm_restart.nml"
          "mpi_amr_plm_restart_np${mpi_amr_plm_3d_ranks}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        regression_mpi_amr_reactive_plm_3d_restart_np${mpi_amr_plm_3d_ranks}
        PROPERTIES
          WORKING_DIRECTORY "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${mpi_amr_plm_3d_ranks}
          PASS_REGULAR_EXPRESSION "Restored coarse steps: 4"
          DEPENDS regression_mpi_amr_reactive_plm_3d_checkpoint_np2
          TIMEOUT 300
      )
      list(
        APPEND MPI_AMR_REACTIVE_PLM_3D_RESTART_TESTS
        regression_mpi_amr_reactive_plm_3d_restart_np${mpi_amr_plm_3d_ranks}
      )
    endforeach()
    add_test(
      NAME regression_mpi_amr_reactive_plm_3d_restart_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_restart_np4"
        --candidate-prefix
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/mpi_amr_plm_restart_np8"
        --checkpoint
          "${MPI_AMR_REACTIVE_PLM_3D_WORK_DIRECTORY}/amr_reactive_plm_3d.chk"
        --final-time 1.0e-6
        --reconstruction characteristic_plm
        --limiter mc
        --expected-checkpoint-sha256
          8f51be7f192dd61ab0e1f4eccf665870c2103205d13ea562fb52d32e981dbf9d
        --expected-coarse-sha256
          864e9cabe8f5c2bee19f97e9279bdada7647305171955c6c767cc83e4db15c21
        --expected-fine-sha256
          8f764e38491e9cf62f8c0f17eac2470c4705f7d80bfe7e686db4ba54de2f58d7
    )
    set_tests_properties(
      regression_mpi_amr_reactive_plm_3d_restart_exact
      PROPERTIES DEPENDS
        "regression_mpi_amr_reactive_plm_3d_serial_reference;${MPI_AMR_REACTIVE_PLM_3D_RESTART_TESTS}"
    )
  endif()
  install(
    TARGETS
      pelef_mpi_1d
      pelef_mpi_multispecies_1d
      pelef_mpi_h2o2_batch
      pelef_mpi_reactive_transport_1d
      pelef_mpi_reactive_1d
      pelef_mpi_amr_patch_1d
      pelef_mpi_eb_amr_patch_2d
      pelef_mpi_amr_eb_patch_tree_2d
      pelef_mpi_reactive_eb_patch_tree_2d
      pelef_mpi_amr_reactive_1d
      pelef_mpi_amr_reactive_3d
    RUNTIME DESTINATION bin
  )
endif()

