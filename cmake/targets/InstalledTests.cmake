if(PELEF_ENABLE_TESTS)
  # Exercise the installed 3D AMR transport/restart contract in a private
  # prefix.  The ordinary 0.245 tests above run build-tree executables; these
  # tests deliberately resolve the binaries through the install tree so that
  # install rules and restart schema selection are covered as well.
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT
    "${CMAKE_CURRENT_BINARY_DIR}/installed_amr_reactive_3d_transport_restart_0245"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/prefix"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/fixed"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/selected"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/fixed-mpi"
  )
  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}/selected-mpi"
  )
  file(
    MAKE_DIRECTORY
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}"
    "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}"
  )
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256
      "6a1b1ba7bece94e4f59e0da95fb6dc8ddc5bc5d356f3d88f56d37e184e0b573d"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256
      "0b0955e45d7d5e6296bdd69270cd78ed02eba4f814a3789d7218d861fd0eac6f"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256
      "45d51c690e8f3b1d20bb224deae24cbbb5ec4dc0c648a9b1367d6963acaeb37b"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256
      "6fa40282c6a12598bf7147efe2253540e49a66abb5e907aac8b5226443fdf342"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256
      "34d781bd9faef30a0fe5f859300fac459380403d67b00795f1c3da95e294cb39"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256
      "e4c3f6fe27ed9c8e62d33912c42adc6bdef4ee1c3e8cd0b6aebfdf34c42956ed"
    )
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256
      "15c2e5a3a087a38fe34ece1fa575bbc7857ec7c651ae9ee84f3fd8d52615afc5"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256
      "c81eb517ab4111f66eafa09b70d1f9543251e206d62143cfe02e50ca74e4f441"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256
      "45d51c690e8f3b1d20bb224deae24cbbb5ec4dc0c648a9b1367d6963acaeb37b"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256
      "52526217eb88e869167b8ef8b958f3cf6176b5346746a5d9ff1d1757ca526df3"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256
      "21470b7ac704274637266990f8bd77bbeacdba5c2808c2943fdc17fb106ee35a"
    )
    set(
      PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256
      "e4c3f6fe27ed9c8e62d33912c42adc6bdef4ee1c3e8cd0b6aebfdf34c42956ed"
    )
  else()
    message(
      FATAL_ERROR
      "installed 0.245 transport restart tests require "
      "CMAKE_BUILD_TYPE=Debug or Release"
    )
  endif()

  set(
    PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST
    regression_installed_amr_reactive_3d_transport_restart_0245_install
  )
  add_test(
    NAME ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
    COMMAND
      "${CMAKE_COMMAND}"
      "-DTEST_BUILD_DIR=${CMAKE_BINARY_DIR}"
      "-DTEST_INSTALL_PREFIX=${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}"
      "-DTEST_CONFIG=$<CONFIG>"
      -P "${PROJECT_SOURCE_DIR}/tests/cmake/install_clean_prefix.cmake"
  )
  set_tests_properties(
    ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
    PROPERTIES
      WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_ROOT}"
      TIMEOUT 600
  )

  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/reference.log"
      --
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d"
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_reference.nml"
  )
  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/checkpoint.log"
      --
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d"
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
  )
  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
      --output
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/restart.log"
      --
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d"
        "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run
    PROPERTIES
      WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}"
      TIMEOUT 300
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    PROPERTIES DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
    PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run
    PROPERTIES
      DEPENDS
        regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint
      PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
  )
  add_test(
    NAME regression_installed_amr_reactive_3d_transport_restart_0245_fixed_exact
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
      --reference-prefix
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_reference_0245"
      --candidate-prefix
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_restarted_0245"
      --checkpoint
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
      --expected-checkpoint-sha256
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256}"
      --expected-coarse-sha256
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256}"
      --expected-fine-sha256
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256}"
      --expected-schema 4
      --final-time 1.0e-9
      --reconstruction characteristic_plm
      --limiter mc
      --reference-log
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/reference.log"
      --candidate-log
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/restart.log"
  )
  set_tests_properties(
    regression_installed_amr_reactive_3d_transport_restart_0245_fixed_exact
    PROPERTIES DEPENDS
      "regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference;regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint;regression_installed_amr_reactive_3d_transport_restart_0245_fixed_run"
  )

  set(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED OFF)
  if(PELEF_MECHANISM_BUNDLE)
    file(SHA256 "${PELEF_MECHANISM_BUNDLE}"
      pelef_installed_transport_0245_bundle_sha256)
    if(pelef_installed_transport_0245_bundle_sha256 STREQUAL
        "f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3")
      set(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED ON)
    endif()
  endif()

  if(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED)
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/reference.log"
        --
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_reference.nml"
    )
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/checkpoint.log"
        --
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
    )
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_run
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/restart.log"
        --
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_run
      PROPERTIES
        WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}"
        TIMEOUT 600
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      PROPERTIES DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
      PROPERTIES PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_run
      PROPERTIES
        DEPENDS
          regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint
        PASS_REGULAR_EXPRESSION "Restarted from checkpoint"
    )
    add_test(
      NAME regression_installed_amr_reactive_3d_transport_restart_0245_selected_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_reference_0245"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_restarted_0245"
        --checkpoint
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_0245.chk"
        --expected-checkpoint-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256}"
        --expected-coarse-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256}"
        --expected-fine-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256}"
        --expected-schema 5
        --expected-bundle-sha256
          f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
        --expected-integrator implicit
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --reference-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/reference.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/restart.log"
    )
    set_tests_properties(
      regression_installed_amr_reactive_3d_transport_restart_0245_selected_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference;regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint;regression_installed_amr_reactive_3d_transport_restart_0245_selected_run"
    )
  endif()

  if(PELEF_ENABLE_MPI)
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/checkpoint-np2.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d"
          "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart_checkpoint.nml"
          installed_fixed_transport_checkpoint_np2
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
        TIMEOUT 300
    )
    set(PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_TESTS)
    foreach(installed_fixed_mpi_rank IN ITEMS 1 2 4 8)
      set(
        installed_fixed_mpi_restart_test
        regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_np${installed_fixed_mpi_rank}
      )
      add_test(
        NAME ${installed_fixed_mpi_restart_test}
        COMMAND
          "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
          --output
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np${installed_fixed_mpi_rank}.log"
          -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          ${installed_fixed_mpi_rank} ${MPIEXEC_PREFLAGS}
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d"
            "${PROJECT_SOURCE_DIR}/cases/amr_reactive_3d/amr_hotspot_transport_restart.nml"
            "installed_fixed_transport_restart_mpi_np${installed_fixed_mpi_rank}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        ${installed_fixed_mpi_restart_test}
        PROPERTIES
          WORKING_DIRECTORY "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${installed_fixed_mpi_rank}
          PASS_REGULAR_EXPRESSION "MPI ranks: ${installed_fixed_mpi_rank}"
          DEPENDS
            regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2
          TIMEOUT 300
      )
      list(APPEND PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_TESTS
        ${installed_fixed_mpi_restart_test})
    endforeach()
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_np8
      PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_reference_0245"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np1"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np2"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np4"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/installed_fixed_transport_restart_mpi_np8"
        --checkpoint
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
        --expected-checkpoint-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_CHECKPOINT_SHA256}"
        --expected-coarse-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_COARSE_SHA256}"
        --expected-fine-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_FINE_SHA256}"
        --expected-schema 4
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --reference-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/reference.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np1.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np2.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np4.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/restart-np8.log"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_fixed_reference;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2;${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_TESTS}"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_exact
      COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_FIXED_MPI_WORK}/amr_hotspot_transport_restart_0245_fixed.chk"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_fixed_checkpoint;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_fixed_checkpoint_np2"
    )
  endif()

  if(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_ENABLED AND PELEF_ENABLE_MPI)
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
        --output
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/checkpoint-np2.log"
        -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" 2
        ${MPIEXEC_PREFLAGS}
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d_selected"
          "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart_checkpoint.nml"
          installed_selected_transport_checkpoint_np2
        ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2
      PROPERTIES
        WORKING_DIRECTORY
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}"
        ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
        PROCESSORS 2
        PASS_REGULAR_EXPRESSION "Stopped after checkpoint"
        DEPENDS ${PELEF_INSTALLED_AMR_TRANSPORT_0245_INSTALL_TEST}
        TIMEOUT 600
    )
    set(PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_TESTS)
    foreach(installed_selected_mpi_rank IN ITEMS 1 2 4 8)
      set(
        installed_selected_mpi_restart_test
        regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_np${installed_selected_mpi_rank}
      )
      add_test(
        NAME ${installed_selected_mpi_restart_test}
        COMMAND
          "${Python3_EXECUTABLE}"
          "${PROJECT_SOURCE_DIR}/tools/run_and_capture.py"
          --output
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np${installed_selected_mpi_rank}.log"
          -- "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}"
          ${installed_selected_mpi_rank} ${MPIEXEC_PREFLAGS}
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_PREFIX}/bin/pelef_mpi_amr_reactive_3d_selected"
            "${PROJECT_SOURCE_DIR}/cases/selected_amr_reactive_3d/h2o2_full_transport_restart.nml"
            "installed_selected_transport_restart_mpi_np${installed_selected_mpi_rank}"
          ${MPIEXEC_POSTFLAGS}
      )
      set_tests_properties(
        ${installed_selected_mpi_restart_test}
        PROPERTIES
          WORKING_DIRECTORY
            "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}"
          ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1"
          PROCESSORS ${installed_selected_mpi_rank}
          PASS_REGULAR_EXPRESSION "MPI ranks: ${installed_selected_mpi_rank}"
          DEPENDS
            regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2
          TIMEOUT 600
      )
      list(APPEND PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_TESTS
        ${installed_selected_mpi_restart_test})
    endforeach()
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_np8
      PROPERTIES PASS_REGULAR_EXPRESSION "Fine x planes per rank: min=0"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_exact
      COMMAND
        "${Python3_EXECUTABLE}"
        "${PROJECT_SOURCE_DIR}/tools/compare_mpi_amr_reactive_3d.py"
        --reference-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_reference_0245"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np1"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np2"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np4"
        --candidate-prefix
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/installed_selected_transport_restart_mpi_np8"
        --checkpoint
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/selected_h2o2_full_transport_restart_0245.chk"
        --expected-checkpoint-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_CHECKPOINT_SHA256}"
        --expected-coarse-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_COARSE_SHA256}"
        --expected-fine-sha256
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_FINE_SHA256}"
        --expected-schema 5
        --expected-bundle-sha256
          f65e1c02e77618d188bc95f0868f3749d6345afb355fda924297521f69ce04c3
        --expected-integrator implicit
        --final-time 1.0e-9
        --reconstruction characteristic_plm
        --limiter mc
        --reference-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/reference.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np1.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np2.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np4.log"
        --candidate-log
          "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/restart-np8.log"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_selected_reference;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2;${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_TESTS}"
    )
    add_test(
      NAME regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_exact
      COMMAND "${CMAKE_COMMAND}" -E compare_files
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_WORK}/selected_h2o2_full_transport_restart_0245.chk"
        "${PELEF_INSTALLED_AMR_TRANSPORT_0245_SELECTED_MPI_WORK}/selected_h2o2_full_transport_restart_0245.chk"
    )
    set_tests_properties(
      regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_exact
      PROPERTIES DEPENDS
        "regression_installed_amr_reactive_3d_transport_restart_0245_selected_checkpoint;regression_installed_mpi_amr_reactive_3d_transport_restart_0245_selected_checkpoint_np2"
    )
  endif()

  add_test(
    NAME quality_project_contract
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_project_contract.py"
      --source "${PROJECT_SOURCE_DIR}"
  )
  add_test(
    NAME quality_mpi_f08_contract_accept
    COMMAND "${CMAKE_COMMAND}" -DTEST_MPI_F08=ON -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_mpi_f08_contract.cmake"
  )
  add_test(
    NAME quality_mpi_f08_contract_reject
    COMMAND "${CMAKE_COMMAND}" -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_mpi_f08_rejection.cmake"
  )
  add_test(
    NAME quality_mpi_launcher_contract_reject
    COMMAND "${CMAKE_COMMAND}" -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_mpi_launcher_rejection.cmake"
  )
  add_test(
    NAME quality_sundials_double_precision_accept
    COMMAND "${CMAKE_COMMAND}"
      -DTEST_PRECISION=double
      -DTEST_BINARY_DIR=${PROJECT_BINARY_DIR}/tests/cmake
      -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_sundials_precision_contract.cmake"
  )
  add_test(
    NAME quality_sundials_non_double_precision_reject
    COMMAND "${CMAKE_COMMAND}"
      -DTEST_BINARY_DIR=${PROJECT_BINARY_DIR}/tests/cmake
      -P
      "${PROJECT_SOURCE_DIR}/tests/cmake/check_sundials_precision_rejection.cmake"
  )

  set(
    PELEF_STACK_CHECK_TARGETS
    $<TARGET_FILE:pelef>
    $<TARGET_FILE:pelef2d>
    $<TARGET_FILE:pelef3d>
    $<TARGET_FILE:pelef_ms>
    $<TARGET_FILE:pelef0d>
    $<TARGET_FILE:pelef0d_h2o2>
    $<TARGET_FILE:pelef0d_h2o2_full>
    $<TARGET_FILE:pelef_reactive_1d>
    $<TARGET_FILE:pelef_amr_reactive_1d>
    $<TARGET_FILE:pelef_reactive_2d>
    $<TARGET_FILE:pelef_reactive_3d>
    $<TARGET_FILE:pelef_reactive_eb_3d>
    $<TARGET_FILE:pelef_amr_reactive_3d>
    $<TARGET_FILE:pelef_reactive_eb_2d>
    $<TARGET_FILE:pelef_reactive_eb_amr_2d>
    $<TARGET_FILE:pelef_reactive_eb_patch_tree_2d>
    $<TARGET_FILE:pelef_transport_probe>
    $<TARGET_FILE:test_amr_eb_multilevel_2d>
    $<TARGET_FILE:test_selected_reactor_config>
    $<TARGET_FILE:test_mechanism_bundle>
    $<TARGET_FILE:test_selected_fixture_reactor>
    $<TARGET_FILE:test_selected_fixture_reactive_1d>
    $<TARGET_FILE:test_selected_fixture_amr_reactive_1d>
    $<TARGET_FILE:test_selected_fixture_reactive_2d>
    $<TARGET_FILE:test_selected_fixture_reactive_eb_2d>
    $<TARGET_FILE:test_selected_fixture_reactive_3d>
    $<TARGET_FILE:test_selected_fixture_amr_reactive_3d>
    $<TARGET_FILE:test_selected_fixture_reactive_eb_3d>
    $<TARGET_FILE:test_selected_full_bundle>
    $<TARGET_FILE:test_selected_full_reactor>
    $<TARGET_FILE:test_selected_full_reactive_1d>
    $<TARGET_FILE:test_selected_full_amr_reactive_1d>
    $<TARGET_FILE:test_selected_full_reactive_2d>
    $<TARGET_FILE:test_selected_full_reactive_eb_2d>
    $<TARGET_FILE:test_selected_full_reactive_3d>
    $<TARGET_FILE:test_selected_full_amr_reactive_3d>
    $<TARGET_FILE:test_selected_full_reactive_eb_3d>
    $<TARGET_FILE:test_selected_elementary_bundle>
    $<TARGET_FILE:test_selected_elementary_reactive_eb_3d>
    $<TARGET_FILE:test_selected_unsupported_element_bundle>
    $<TARGET_FILE:test_selected_unsupported_element_reactive_eb_3d>
  )
  if(PELEF_MECHANISM_BUNDLE)
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:pelef_mechanism_probe>
      $<TARGET_FILE:pelef0d_selected>
      $<TARGET_FILE:pelef_reactive_1d_selected>
      $<TARGET_FILE:pelef_amr_reactive_1d_selected>
      $<TARGET_FILE:pelef_reactive_2d_selected>
      $<TARGET_FILE:pelef_reactive_eb_2d_selected>
      $<TARGET_FILE:pelef_reactive_eb_amr_2d_selected>
      $<TARGET_FILE:pelef_reactive_3d_selected>
      $<TARGET_FILE:pelef_amr_reactive_3d_selected>
      $<TARGET_FILE:pelef_reactive_eb_3d_selected>
    )
    if(PELEF_ENABLE_MPI)
      list(
        APPEND PELEF_STACK_CHECK_TARGETS
        $<TARGET_FILE:pelef_mpi_reactive_1d_selected>
        $<TARGET_FILE:pelef_mpi_amr_reactive_1d_selected>
        $<TARGET_FILE:pelef_mpi_amr_reactive_3d_selected>
        $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d_selected>
      )
    endif()
  endif()
  if(PELEF_ENABLE_SUNDIALS)
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:test_sundials_constant_volume_reactor>
      $<TARGET_FILE:test_sundials_constant_volume_multi_context>
    )
  endif()
  if(PELEF_ENABLE_MPI)
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:pelef_mpi_1d>
      $<TARGET_FILE:pelef_mpi_multispecies_1d>
      $<TARGET_FILE:pelef_mpi_h2o2_batch>
      $<TARGET_FILE:pelef_mpi_reactive_transport_1d>
      $<TARGET_FILE:pelef_mpi_reactive_1d>
      $<TARGET_FILE:pelef_mpi_amr_patch_1d>
      $<TARGET_FILE:pelef_mpi_eb_amr_patch_2d>
      $<TARGET_FILE:pelef_mpi_amr_eb_patch_tree_2d>
      $<TARGET_FILE:pelef_mpi_reactive_eb_patch_tree_2d>
      $<TARGET_FILE:pelef_mpi_amr_reactive_1d>
      $<TARGET_FILE:pelef_mpi_amr_reactive_3d>
      $<TARGET_FILE:test_mpi_reactive_integrator_policy>
      $<TARGET_FILE:test_mpi_amr_sparse_integrator_policy>
      $<TARGET_FILE:test_mpi_amr_eb_sparse_integrator_policy>
      $<TARGET_FILE:test_mpi_amr_selected_context_mismatch>
      $<TARGET_FILE:test_mpi_amr_selected_context_3d>
      $<TARGET_FILE:test_mpi_amr_eb_selected_context_mismatch>
      $<TARGET_FILE:test_mpi_amr_eb_selected_checkpoint_presence>
    )
    list(
      APPEND PELEF_STACK_CHECK_TARGETS
      $<TARGET_FILE:test_selected_fixture_mpi_reactive_1d>
      $<TARGET_FILE:test_selected_full_mpi_reactive_1d>
      $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_1d>
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_1d>
      $<TARGET_FILE:test_selected_fixture_mpi_amr_reactive_3d>
      $<TARGET_FILE:test_selected_full_mpi_amr_reactive_3d>
      $<TARGET_FILE:test_selected_fixture_mpi_reactive_eb_patch_tree_2d>
      $<TARGET_FILE:test_selected_full_mpi_reactive_eb_patch_tree_2d>
    )
  endif()
  add_test(
    NAME quality_no_executable_stack
    COMMAND "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/check_no_executable_stack.py"
      ${PELEF_STACK_CHECK_TARGETS}
  )
endif()
