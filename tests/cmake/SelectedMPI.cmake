if(PELEF_ENABLE_MPI)
  set(
    SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_1d_fixture"
  )
  set(
    SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_1d_full"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}"
    "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}"
  )

  set(SELECTED_MPI_REACTIVE_1D_FIXTURE_RUN_TESTS)
  set(SELECTED_MPI_REACTIVE_1D_FULL_RUN_TESTS)
  foreach(selected_mpi_ranks IN ITEMS 1 2 4)
    set(
      fixture_test_name
      regression_selected_mpi_reactive_1d_fixture_np${selected_mpi_ranks}
    )
    add_test(
      NAME ${fixture_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_fixture_mpi_reactive_1d>
        "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np${selected_mpi_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${fixture_test_name}
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_ranks}
        PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
        TIMEOUT 300
    )
    list(APPEND SELECTED_MPI_REACTIVE_1D_FIXTURE_RUN_TESTS ${fixture_test_name})

    set(
      full_test_name
      regression_selected_mpi_reactive_1d_full_np${selected_mpi_ranks}
    )
    add_test(
      NAME ${full_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_1d>
        "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np${selected_mpi_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${full_test_name}
      PROPERTIES
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_ranks}
        PASS_REGULAR_EXPRESSION "Chemistry integrator: implicit"
        TIMEOUT 300
    )
    list(APPEND SELECTED_MPI_REACTIVE_1D_FULL_RUN_TESTS ${full_test_name})
  endforeach()

  add_test(
    NAME regression_selected_mpi_reactive_1d_fixture_rank_parity
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_reactive_1d.py"
      "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np1.csv"
      "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np2.csv"
      "${SELECTED_MPI_REACTIVE_1D_FIXTURE_WORK_DIRECTORY}/fixture_np4.csv"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_1d_fixture_rank_parity
    PROPERTIES DEPENDS "${SELECTED_MPI_REACTIVE_1D_FIXTURE_RUN_TESTS}"
  )

  add_test(
    NAME regression_selected_mpi_reactive_1d_full_rank_parity
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_reactive_1d.py"
      "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np1.csv"
      "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np2.csv"
      "${SELECTED_MPI_REACTIVE_1D_FULL_WORK_DIRECTORY}/full_np4.csv"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_1d_full_rank_parity
    PROPERTIES DEPENDS "${SELECTED_MPI_REACTIVE_1D_FULL_RUN_TESTS}"
  )

  set(
    SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_1d"
  )
  file(MAKE_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}")
  set(SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS)

  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_fixed_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_amr_reactive_1d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/fixed_reference.nml"
      "fixed_reference.csv"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_fixed_np1
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "MPI ranks: 1"
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
    regression_selected_mpi_amr_reactive_1d_fixed_np1
  )

  set(SELECTED_MPI_AMR_REACTIVE_1D_REFERENCE_TESTS)
  foreach(selected_mpi_amr_ranks IN ITEMS 1 2 4)
    set(
      selected_mpi_amr_reference_test
      regression_selected_mpi_amr_reactive_1d_reference_np${selected_mpi_amr_ranks}
    )
    add_test(
      NAME ${selected_mpi_amr_reference_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/selected_reference.nml"
        "selected_reference_np${selected_mpi_amr_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_amr_reference_test}
      PROPERTIES
        WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_ranks}"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_REFERENCE_TESTS
      ${selected_mpi_amr_reference_test}
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
      ${selected_mpi_amr_reference_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/checkpoint_stop.nml"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
    regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
  )

  set(SELECTED_MPI_AMR_REACTIVE_1D_RESTART_TESTS)
  foreach(selected_mpi_amr_restart_ranks IN ITEMS 2 4)
    set(
      selected_mpi_amr_restart_test
      regression_selected_mpi_amr_reactive_1d_restart_np${selected_mpi_amr_restart_ranks}
    )
    add_test(
      NAME ${selected_mpi_amr_restart_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/restart.nml"
        "selected_restart_np${selected_mpi_amr_restart_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_amr_restart_test}
      PROPERTIES
        WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_restart_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_restart_ranks}"
        DEPENDS regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_RESTART_TESTS
      ${selected_mpi_amr_restart_test}
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS
      ${selected_mpi_amr_restart_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_restart_check
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_sparse_amr_restart_1d.py"
      --checkpoint
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_sparse_mpi_amr.chk"
      --fixed-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/fixed_reference.csv"
      --selected-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_reference_np1.csv"
      --selected-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_reference_np2.csv"
      --selected-reference
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_reference_np4.csv"
      --restarted
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_restart_np2.csv"
      --restarted
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_restart_np4.csv"
      --stopped
        "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_checkpoint_stop.csv"
      --expected-bundle-sha256
        f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
      --expected-checkpoint-sha256
        0e5cc002ac5cceb3fa7732f8fe46382f140e038d19cac17664230ce8c500c2fd
      --final-time 1.0e-9
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_restart_check
    PROPERTIES DEPENDS "${SELECTED_MPI_AMR_REACTIVE_1D_RUN_TESTS}"
  )

  foreach(
      selected_mpi_amr_mismatch
      IN ITEMS
        sha256 integrator composition species baseline density state_species
        near_floor_species state_closure temperature geometry
    )
    set(
      mismatch_work_directory
      "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/mismatch_${selected_mpi_amr_mismatch}"
    )
    file(MAKE_DIRECTORY "${mismatch_work_directory}")
    configure_file(
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_sparse_amr_restart_1d/restart.nml"
      "${mismatch_work_directory}/restart.nml"
      COPYONLY
    )
    set(
      mutate_test_name
      regression_selected_mpi_amr_reactive_1d_mutate_${selected_mpi_amr_mismatch}
    )
    add_test(
      NAME ${mutate_test_name}
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/mutate_selected_mpi_amr_checkpoint_1d.py"
        --input
          "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/selected_sparse_mpi_amr.chk"
        --output "${mismatch_work_directory}/selected_sparse_mpi_amr.chk"
        --field "${selected_mpi_amr_mismatch}"
    )
    set_tests_properties(
      ${mutate_test_name}
      PROPERTIES DEPENDS
        regression_selected_mpi_amr_reactive_1d_checkpoint_stop_np1
    )

    set(
      reject_test_name
      regression_selected_mpi_amr_reactive_1d_reject_${selected_mpi_amr_mismatch}
    )
    add_test(
      NAME ${reject_test_name}
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
        --executable $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
        --input "${mismatch_work_directory}/restart.nml"
        --expected "Sparse MPI AMR restart read failed"
        --forbidden-output
          "${mismatch_work_directory}/selected_restarted.csv"
    )
    set_tests_properties(
      ${reject_test_name}
      PROPERTIES
        WORKING_DIRECTORY "${mismatch_work_directory}"
        DEPENDS ${mutate_test_name}
        TIMEOUT 60
      )
  endforeach()

  set(
    SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY
    "${SELECTED_MPI_AMR_REACTIVE_1D_WORK_DIRECTORY}/invalid"
  )
  file(MAKE_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}")
  configure_file(
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_1d/invalid_path_alias.nml"
    "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
    COPYONLY
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_1d_path_alias_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_expected_failure.py"
      --executable $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
      --input invalid_path_alias.nml
      --expected "Selected sparse MPI AMR input and output paths alias"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}/invalid_selected_mpi_amr_alias.csv"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_1d_path_alias_rejected
    PROPERTIES
      WORKING_DIRECTORY "${SELECTED_MPI_AMR_REACTIVE_1D_INVALID_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      TIMEOUT 60
  )

  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_reactive_eb_patch_tree_2d"
  )
  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/invalid"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}"
  )

  set(SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS)
  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_fixed_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_fixed.nml"
      fixed_np1.csv
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_fixed_np1
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Chemistry level advances:"
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS
    regression_selected_mpi_reactive_eb_patch_tree_2d_fixed_np1
  )

  set(
    FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/fixed_restart"
  )
  file(
    MAKE_DIRECTORY
    "${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
  )
  add_test(
    NAME regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/reactive_eb_patch_tree_restart_2d/checkpoint_stop.nml"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    PROPERTIES
      WORKING_DIRECTORY
        "${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
      TIMEOUT 300
  )
  set(FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      list(
        APPEND FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS
        --expected-sha256
          6b817c871fe03eec3c96601e77c38d2e14e0505b96f9a590f4044b539ce85b3a
      )
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
      list(
        APPEND FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS
        --expected-sha256
          d3ed831ad70a488f7303bbc296b07f63f503c9f02c04a0ddff81d0bb8e7f1de1
      )
    endif()
  endif()
  add_test(
    NAME regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_check
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint.py"
      --checkpoint
        "${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/patch_tree_restart.chk"
      ${FIXED_MPI_REACTIVE_EB_PATCH_TREE_2D_HASH_ARGUMENTS}
  )
  set_tests_properties(
    regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_check
    PROPERTIES DEPENDS
      regression_fixed_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
  )

  foreach(selected_mpi_eb_tree_ranks IN ITEMS 1 2 4)
    set(
      selected_mpi_eb_tree_test
      regression_selected_mpi_reactive_eb_patch_tree_2d_np${selected_mpi_eb_tree_ranks}
    )
    add_test(
      NAME ${selected_mpi_eb_tree_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_eb_tree_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
        "selected_np${selected_mpi_eb_tree_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_eb_tree_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_eb_tree_ranks}
        PASS_REGULAR_EXPRESSION "Chemistry integrator: implicit"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS
      ${selected_mpi_eb_tree_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_exact
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d.py"
      --reference
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/fixed_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_np2.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_np4.csv"
      --species H2 H O O2 OH H2O HO2 H2O2 AR N2
      --expected-levels 0 1 2 3
      --final-time 1.0e-12
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_exact
    PROPERTIES DEPENDS "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RUN_TESTS}"
  )

  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/selected_restart"
  )
  set(
    SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/invalid"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
    "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY}"
  )
  set(SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS)

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_restart_fixed_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_fixed_reference.nml"
      fixed_reference_np1.csv
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_fixed_np1
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Completed root steps: 2"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_fixed_np1
  )

  foreach(selected_mpi_eb_restart_ranks IN ITEMS 1 2 4)
    set(
      selected_mpi_eb_restart_reference_test
      regression_selected_mpi_reactive_eb_patch_tree_2d_restart_selected_np${selected_mpi_eb_restart_ranks}
    )
    add_test(
      NAME ${selected_mpi_eb_restart_reference_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_eb_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_selected_reference.nml"
        "selected_reference_np${selected_mpi_eb_restart_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_eb_restart_reference_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_eb_restart_ranks}
        PASS_REGULAR_EXPRESSION "Completed root steps: 2"
        RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
      ${selected_mpi_eb_restart_reference_test}
    )
  endforeach()

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_checkpoint_stop.nml"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint:  T"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 300
  )
  list(
    APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
    regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
  )

  foreach(selected_mpi_eb_restart_ranks IN ITEMS 2 4)
    set(
      selected_mpi_eb_restart_test
      regression_selected_mpi_reactive_eb_patch_tree_2d_restart_np${selected_mpi_eb_restart_ranks}
    )
    add_test(
      NAME ${selected_mpi_eb_restart_test}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_eb_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml"
        "restarted_np${selected_mpi_eb_restart_ranks}.csv"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${selected_mpi_eb_restart_test}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_eb_restart_ranks}
        PASS_REGULAR_EXPRESSION "Restarted:  T"
        DEPENDS
          regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
        RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS
      ${selected_mpi_eb_restart_test}
    )
  endforeach()

  set(SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_HASH_ARGUMENTS)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
     CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
    list(
      APPEND SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_HASH_ARGUMENTS
      --expected-checkpoint-sha256
        2edd758ca1ff888735f5719de8ed9969b36f60d45bb0b3c523e47c2e35510f4d
      --expected-final-sha256
        5394fd007283d746ec4507832ccce8a7013fd4fb9b5c1b10d50d7c91e216975a
      --expected-stopped-sha256
        bba321c7ace389a6e0d95af0d1a07edd55c9cde870daee514e9b6c1c49e5fc3c
    )
  endif()
  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_restart_check
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d_restart.py"
      --checkpoint
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_mpi_eb_restart.chk"
      --reference
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/fixed_reference_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_reference_np1.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_reference_np2.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_reference_np4.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/restarted_np2.csv"
      --candidate
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/restarted_np4.csv"
      --stopped
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_mpi_eb_restart_stopped.csv"
      --species H2 H O O2 OH H2O HO2 H2O2 AR N2
      --expected-levels 0 1 2 3
      --final-time 1.0e-12
      --bundle-sha256
        f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
      --integrator implicit
      --composition 0.29570 1.0e-5 1.0e-5 0.14784 1.0e-5 0.0 0.0 0.0 0.0
        0.55643
      ${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_HASH_ARGUMENTS}
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_check
    PROPERTIES
      DEPENDS "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_RUN_TESTS}"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
  )

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_restart_failures
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d_restart_failures.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 1
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml"
      --checkpoint
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_WORK_DIRECTORY}/selected_mpi_eb_restart.chk"
      --work-directory
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY}/corruptions"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_restart_failures
    PROPERTIES
      DEPENDS
        regression_selected_mpi_reactive_eb_patch_tree_2d_checkpoint_stop
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 300
  )

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_aliases
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d_aliases.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 1
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --reference-input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_selected_reference.nml"
      --checkpoint-input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_checkpoint_stop.nml"
      --restart-input
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/restart_restart.nml"
      --work-directory
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_RESTART_INVALID_WORK_DIRECTORY}/aliases"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_aliases
    PROPERTIES
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      RESOURCE_LOCK selected_mpi_reactive_eb_patch_tree_2d_restart
      TIMEOUT 180
  )

  add_test(
    NAME regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_np1
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 1 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_fixture_mpi_reactive_eb_patch_tree_2d>
      "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/fixture_selected.nml"
      "fixture_selected_np1.csv"
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_np1
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 1
      PASS_REGULAR_EXPRESSION "Chemistry integrator: explicit"
      TIMEOUT 300
  )
  add_test(
    NAME regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_valid
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_mpi_reactive_eb_patch_tree_2d.py"
      --reference
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_WORK_DIRECTORY}/fixture_selected_np1.csv"
      --species H2 H
      --expected-levels 0 1 2 3
      --final-time 1.0e-12
  )
  set_tests_properties(
    regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_valid
    PROPERTIES DEPENDS
      regression_selected_fixture_mpi_reactive_eb_patch_tree_2d_np1
  )

  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_model_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_fixed.nml"
      --expected "Selected sparse MPI EB requires chemistry_model='selected'"
      --forbidden-output
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}/selected_mpi_eb_tree_fixed.csv"
  )
  add_test(
    NAME regression_fixed_mpi_reactive_eb_patch_tree_2d_selected_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
      --expected "Unknown reactive 2D chemistry model"
      --forbidden-output
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}/selected_mpi_eb_tree_selected.csv"
  )
  add_test(
    NAME regression_selected_mpi_reactive_eb_patch_tree_2d_alias_rejected
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_mpi_reactive_eb_patch_tree_2d/h2o2_full_selected.nml"
      --expected "Selected sparse MPI EB input and output paths alias"
      --forbidden-output
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}/selected_mpi_eb_tree_selected.csv"
  )
  set_tests_properties(
    regression_selected_mpi_reactive_eb_patch_tree_2d_model_rejected
    regression_fixed_mpi_reactive_eb_patch_tree_2d_selected_rejected
    regression_selected_mpi_reactive_eb_patch_tree_2d_alias_rejected
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_REACTIVE_EB_PATCH_TREE_2D_INVALID_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      TIMEOUT 60
  )
endif()

