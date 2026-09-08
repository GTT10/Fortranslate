if(PELEF_ENABLE_MPI)
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_fixture"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_full"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_transport"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_restart"
  )
  set(
    SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/selected_mpi_amr_reactive_3d_invalid"
  )
  file(
    MAKE_DIRECTORY
    "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
    "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
  )
  configure_file(
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_path_alias.nml"
    "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_path_alias.nml"
    COPYONLY
  )
  configure_file(
    "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_checkpoint_input_alias.nml"
    "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_checkpoint_input_alias.nml"
    COPYONLY
  )

  set(SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_RUN_TESTS)
  foreach(selected_mpi_amr_3d_ranks IN ITEMS 1 2 4)
    set(
      fixture_test_name
      regression_selected_mpi_amr_reactive_3d_fixture_np${selected_mpi_amr_3d_ranks}
    )
    add_test(
      NAME ${fixture_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/fixture.nml"
        "selected_fixture_mpi_np${selected_mpi_amr_3d_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${fixture_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_3d_ranks}"
        TIMEOUT 180
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_RUN_TESTS
      ${fixture_test_name}
    )
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_fixture_np4
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_fixture_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
      --coarse
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np4_coarse.csv"
      --fine
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np4_fine.csv"
      --species H2 H
      --nx 4 --ny 4 --nz 4
      --i-lower 2 --i-upper 3
      --j-lower 2 --j-upper 3
      --k-lower 2 --k-upper 3
      --ratio 2
      --final-time 1.0e-7
      --require-uniform
      --activity-species H2
      --initial-mass-fraction 0.8888888888888889
      --minimum-change 1.0e-6
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_fixture_check
    PROPERTIES DEPENDS regression_selected_mpi_amr_reactive_3d_fixture_np4
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_fixture_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_amr_reactive_3d"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_WORK_DIRECTORY}/selected_fixture_mpi_np4"
      --final-time 1.0e-7
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_fixture_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_fixture_run_a;${SELECTED_MPI_AMR_REACTIVE_3D_FIXTURE_RUN_TESTS}"
  )

  set(SELECTED_MPI_AMR_REACTIVE_3D_FULL_RUN_TESTS)
  foreach(selected_mpi_amr_3d_ranks IN ITEMS 1 2 4 8)
    set(
      full_test_name
      regression_selected_mpi_amr_reactive_3d_full_np${selected_mpi_amr_3d_ranks}
    )
    add_test(
      NAME ${full_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_chemistry.nml"
        "selected_full_mpi_np${selected_mpi_amr_3d_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${full_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_3d_ranks}"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_3D_FULL_RUN_TESTS ${full_test_name}
    )
  endforeach()
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_full_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
      --coarse
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np8_coarse.csv"
      --fine
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np8_fine.csv"
      --species H2 H O O2 OH H2O HO2 H2O2 AR N2
      --nx 8 --ny 8 --nz 8
      --i-lower 3 --i-upper 6
      --j-lower 3 --j-upper 6
      --k-lower 3 --k-upper 6
      --ratio 2
      --final-time 1.0e-9
      --activity-species H2
      --initial-mass-fraction 0.028502723329253223
      --minimum-change 1.0e-12
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_full_check
    PROPERTIES DEPENDS regression_selected_mpi_amr_reactive_3d_full_np8
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_full_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_amr_hotspot_chemistry_3d"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_FULL_WORK_DIRECTORY}/selected_full_mpi_np8"
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --expected-coarse-sha256
        3313a9a6eff5d6175625e4bbe3613a782cdd794314fd9fbdd53ddaa309fc9516
      --expected-fine-sha256
        4ff0785c3462a49648fd694d090eb23bffd5af05a43a33360cd187417bd8a31d
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_full_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_full_selected_run;${SELECTED_MPI_AMR_REACTIVE_3D_FULL_RUN_TESTS}"
  )

  set(SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RUN_TESTS)
  foreach(selected_mpi_amr_3d_transport_ranks IN ITEMS 1 2 4 8)
    set(
      transport_test_name
      regression_selected_mpi_amr_reactive_3d_transport_np${selected_mpi_amr_3d_transport_ranks}
    )
    add_test(
      NAME ${transport_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_transport_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport.nml"
        "selected_transport_mpi_np${selected_mpi_amr_3d_transport_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${transport_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_transport_ranks}
        PASS_REGULAR_EXPRESSION
          "Completed transport fine substeps: 16"
        TIMEOUT 300
    )
    list(
      APPEND SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RUN_TESTS
      ${transport_test_name}
    )
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_check
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_selected_amr_reactive_3d.py"
      --coarse
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np8_coarse.csv"
      --fine
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np8_fine.csv"
      --species H2 H O O2 OH H2O HO2 H2O2 AR N2
      --nx 8 --ny 8 --nz 8
      --i-lower 3 --i-upper 6
      --j-lower 3 --j-upper 6
      --k-lower 3 --k-upper 6
      --ratio 2
      --final-time 1.0e-9
      --activity-species H2
      --initial-mass-fraction 0.028502723329253223
      --minimum-change 1.0e-12
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_check
    PROPERTIES DEPENDS regression_selected_mpi_amr_reactive_3d_transport_np8
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_transport_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_amr_hotspot_transport_3d"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_WORK_DIRECTORY}/selected_transport_mpi_np8"
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --expected-coarse-sha256
        $<IF:$<CONFIG:Debug>,5687cb97f834b759b3932bbd3f0b6e30fc5b48e6029e3f47cadfb3ec6a8d74fb,f0713f3a68f48b6537afaafd57209ea8ed4767f9c7d093fe2961bd3a1fa37fae>
      --expected-fine-sha256
        eaf132bbf89974d739c26a35f6b76137f9774438cfe9eebfc12436aa364cff58
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_transport_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_transport_selected_run;${SELECTED_MPI_AMR_REACTIVE_3D_TRANSPORT_RUN_TESTS}"
  )

  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_checkpoint_np2
    COMMAND
      "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
      ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
      "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart_checkpoint.nml"
      selected_full_mpi_checkpoint_np2
      ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_checkpoint_np2
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
      ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
      PROCESSORS 2
      PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
      TIMEOUT 300
  )
  set(SELECTED_MPI_AMR_REACTIVE_3D_RESTART_TESTS)
  foreach(selected_mpi_amr_3d_restart_ranks IN ITEMS 1 2 4 8)
    set(
      restart_test_name
      regression_selected_mpi_amr_reactive_3d_restart_np${selected_mpi_amr_3d_restart_ranks}
    )
    add_test(
      NAME ${restart_test_name}
      COMMAND
        "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
        "${selected_mpi_amr_3d_restart_ranks}" ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_restart.nml"
        "selected_full_mpi_restart_np${selected_mpi_amr_3d_restart_ranks}"
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      ${restart_test_name}
      PROPERTIES
        WORKING_DIRECTORY
          "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS ${selected_mpi_amr_3d_restart_ranks}
        PASS_REGULAR_EXPRESSION "MPI ranks: ${selected_mpi_amr_3d_restart_ranks}"
        DEPENDS regression_selected_mpi_amr_reactive_3d_checkpoint_np2
        TIMEOUT 300
    )
    list(APPEND SELECTED_MPI_AMR_REACTIVE_3D_RESTART_TESTS ${restart_test_name})
  endforeach()
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_restart_np8
    PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_restart_exact
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart_reference"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np1"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np2"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np4"
      --candidate-prefix
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_full_mpi_restart_np8"
      --checkpoint
        "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --expected-schema 3
      --expected-bundle-sha256
        f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
      --expected-integrator implicit
      --expected-checkpoint-sha256
        f44f69a7e1a1657b549b26373feb277840732e5453b0734321cc8c7cdf91c423
      --expected-coarse-sha256
        a03884837909bab719f0783d23ea9150a81f2afc628fee1a01c645c1b47e0ff3
      --expected-fine-sha256
        8d5a992d07b89832624e37af6ddacd510df63633707678693b9053057d4efcbe
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_restart_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_restart_reference;regression_selected_mpi_amr_reactive_3d_checkpoint_np2;${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_TESTS}"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_checkpoint_exact
    COMMAND "${CMAKE_COMMAND}" -E compare_files
      "${SELECTED_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
      "${SELECTED_MPI_AMR_REACTIVE_3D_RESTART_WORK_DIRECTORY}/selected_amr_full_restart.chk"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_checkpoint_exact
    PROPERTIES DEPENDS
      "regression_selected_amr_reactive_3d_checkpoint_stop;regression_selected_mpi_amr_reactive_3d_checkpoint_np2"
  )

  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_restart_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      --argument
        "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/invalid_restart.nml"
      --expected "Could not open 3D AMR checkpoint"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart_coarse.csv"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_restart_fine.csv"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_path_alias_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      --argument invalid_path_alias.nml
      --expected "input and output paths alias"
      --preserved-input invalid_path_alias.nml
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_alias_fine.csv"
  )
  add_test(
    NAME regression_selected_mpi_amr_reactive_3d_checkpoint_alias_rejected
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_mpi_expected_failure.py"
      --launcher "${MPIEXEC_EXECUTABLE}"
      "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
      --ranks 2
      --executable $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      --argument invalid_checkpoint_input_alias.nml
      --expected "checkpoint path aliases another file"
      --preserved-input invalid_checkpoint_input_alias.nml
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_alias_coarse.csv"
      --forbidden-output
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}/invalid_selected_amr_checkpoint_alias_fine.csv"
  )
  set_tests_properties(
    regression_selected_mpi_amr_reactive_3d_restart_rejected
    regression_selected_mpi_amr_reactive_3d_path_alias_rejected
    regression_selected_mpi_amr_reactive_3d_checkpoint_alias_rejected
    PROPERTIES
      WORKING_DIRECTORY
        "${SELECTED_MPI_AMR_REACTIVE_3D_INVALID_WORK_DIRECTORY}"
      PROCESSORS 2
      TIMEOUT 60
  )
endif()

