set(pelef_selected_spray_option)
if(PELEF_ENABLE_SPRAY)
  set(pelef_selected_spray_option SPRAY_3D_NAME pelef_spray_3d_selected)
endif()
if(PELEF_MECHANISM_BUNDLE)
  set(pelef_selected_test_option)
  set(pelef_selected_mpi_option)
  if(PELEF_ENABLE_TESTS)
    list(APPEND pelef_selected_test_option ADD_TEST)
  endif()
  if(PELEF_ENABLE_MPI)
    list(
      APPEND pelef_selected_mpi_option
      MPI_REACTIVE_1D_NAME pelef_mpi_reactive_1d_selected
      MPI_AMR_REACTIVE_1D_NAME pelef_mpi_amr_reactive_1d_selected
      MPI_AMR_REACTIVE_3D_NAME pelef_mpi_amr_reactive_3d_selected
      MPI_REACTIVE_EB_PATCH_TREE_2D_NAME
        pelef_mpi_reactive_eb_patch_tree_2d_selected
    )
  endif()
  if(PELEF_MECHANISM_SOURCE)
    pelef_add_mechanism_bundle(
      NAME pelef_selected
      BUNDLE "${PELEF_MECHANISM_BUNDLE}"
      SOURCE "${PELEF_MECHANISM_SOURCE}"
      PROBE_NAME pelef_mechanism_probe
      REACTOR_NAME pelef0d_selected
      REACTIVE_1D_NAME pelef_reactive_1d_selected
      AMR_REACTIVE_1D_NAME pelef_amr_reactive_1d_selected
      REACTIVE_2D_NAME pelef_reactive_2d_selected
      REACTIVE_EB_2D_NAME pelef_reactive_eb_2d_selected
      REACTIVE_EB_AMR_2D_NAME pelef_reactive_eb_amr_2d_selected
      ${pelef_selected_spray_option}
      REACTIVE_3D_NAME pelef_reactive_3d_selected
      AMR_REACTIVE_3D_NAME pelef_amr_reactive_3d_selected
      REACTIVE_EB_3D_NAME pelef_reactive_eb_3d_selected
      ${pelef_selected_mpi_option}
      INSTALL
      ${pelef_selected_test_option}
    )
  else()
    pelef_add_mechanism_bundle(
      NAME pelef_selected
      BUNDLE "${PELEF_MECHANISM_BUNDLE}"
      PROBE_NAME pelef_mechanism_probe
      REACTOR_NAME pelef0d_selected
      REACTIVE_1D_NAME pelef_reactive_1d_selected
      AMR_REACTIVE_1D_NAME pelef_amr_reactive_1d_selected
      REACTIVE_2D_NAME pelef_reactive_2d_selected
      REACTIVE_EB_2D_NAME pelef_reactive_eb_2d_selected
      REACTIVE_EB_AMR_2D_NAME pelef_reactive_eb_amr_2d_selected
      ${pelef_selected_spray_option}
      REACTIVE_3D_NAME pelef_reactive_3d_selected
      AMR_REACTIVE_3D_NAME pelef_amr_reactive_3d_selected
      REACTIVE_EB_3D_NAME pelef_reactive_eb_3d_selected
      ${pelef_selected_mpi_option}
      INSTALL
      ${pelef_selected_test_option}
    )
  endif()
endif()

add_executable(pelef app/pelef.F90)
target_link_libraries(pelef PRIVATE pelef_core)

add_executable(pelef2d app/pelef2d.F90)
target_link_libraries(pelef2d PRIVATE pelef_core)

add_executable(pelef3d app/pelef3d.F90)
target_link_libraries(pelef3d PRIVATE pelef_core)

add_executable(pelef_ms app/pelef_ms.F90)
target_link_libraries(pelef_ms PRIVATE pelef_core)

add_executable(pelef0d app/pelef0d.F90)
target_link_libraries(pelef0d PRIVATE pelef_core)

add_executable(pelef0d_h2o2 app/pelef0d_h2o2.F90)
target_link_libraries(pelef0d_h2o2 PRIVATE pelef_core)

add_executable(pelef0d_h2o2_full app/pelef0d_h2o2_full.F90)
target_link_libraries(pelef0d_h2o2_full PRIVATE pelef_core)

add_executable(pelef_reactive_1d app/pelef_reactive_1d.F90)
target_link_libraries(pelef_reactive_1d PRIVATE pelef_core)

add_executable(pelef_amr_reactive_1d app/pelef_amr_reactive_1d.F90)
target_link_libraries(pelef_amr_reactive_1d PRIVATE pelef_core)

add_executable(pelef_reactive_2d app/pelef_reactive_2d.F90)
target_link_libraries(pelef_reactive_2d PRIVATE pelef_core)

add_executable(pelef_reactive_3d app/pelef_reactive_3d.F90)
target_link_libraries(pelef_reactive_3d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_3d app/pelef_reactive_eb_3d.F90)
target_link_libraries(pelef_reactive_eb_3d PRIVATE pelef_core)

add_executable(pelef_amr_reactive_3d app/pelef_amr_reactive_3d.F90)
target_link_libraries(pelef_amr_reactive_3d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_2d app/pelef_reactive_eb_2d.F90)
target_link_libraries(pelef_reactive_eb_2d PRIVATE pelef_core)

add_executable(pelef_reactive_eb_amr_2d app/pelef_reactive_eb_amr_2d.F90)
target_link_libraries(pelef_reactive_eb_amr_2d PRIVATE pelef_core)

add_executable(
  pelef_reactive_eb_patch_tree_2d
  app/pelef_reactive_eb_patch_tree_2d.F90
)
target_link_libraries(pelef_reactive_eb_patch_tree_2d PRIVATE pelef_core)

add_executable(pelef_transport_probe app/pelef_transport_probe.F90)
target_link_libraries(pelef_transport_probe PRIVATE pelef_core)

install(
  TARGETS pelef pelef2d pelef3d pelef_ms pelef0d pelef0d_h2o2
    pelef_reactive_1d
    pelef_amr_reactive_1d
    pelef0d_h2o2_full pelef_reactive_2d pelef_reactive_3d
    pelef_reactive_eb_3d
    pelef_amr_reactive_3d
    pelef_reactive_eb_2d
    pelef_reactive_eb_amr_2d
    pelef_reactive_eb_patch_tree_2d
    pelef_transport_probe
  RUNTIME DESTINATION bin
)

