include(CMakeParseArguments)

function(_pelef_read_mechanism_bundle_build_contract bundle_path variable_prefix)
  if(NOT EXISTS "${bundle_path}")
    message(FATAL_ERROR "mechanism bundle does not exist: ${bundle_path}")
  endif()

  file(READ "${bundle_path}" bundle_json)
  string(
    JSON bundle_schema_version ERROR_VARIABLE bundle_json_error
    GET "${bundle_json}" schema_version
  )
  if(NOT bundle_json_error STREQUAL "NOTFOUND" OR
      NOT "${bundle_schema_version}" STREQUAL "1")
    message(
      FATAL_ERROR
      "selected mechanism build contract requires schema_version=1; "
      "legacy kinetics-only bundles are generator-only"
    )
  endif()

  foreach(
      json_key
      module_name
      loader_name
      kernel_name
      jacobian_name
      thermo_loader_name
      transport_loader_name
      symbol_prefix
      chemistry_integrator
    )
    string(
      JSON json_value ERROR_VARIABLE bundle_json_error
      GET "${bundle_json}" "${json_key}"
    )
    if(NOT bundle_json_error STREQUAL "NOTFOUND")
      message(
        FATAL_ERROR
        "selected mechanism build contract requires explicit ${json_key}; "
        "legacy kinetics-only bundles are generator-only"
      )
    endif()
    if(NOT json_value MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
      message(FATAL_ERROR "mechanism bundle has invalid ${json_key}")
    endif()
    set("${variable_prefix}_${json_key}" "${json_value}")
    set("${variable_prefix}_${json_key}" "${json_value}" PARENT_SCOPE)
  endforeach()

  set(
    selected_chemistry_integrator
    "${${variable_prefix}_chemistry_integrator}"
  )
  if(NOT selected_chemistry_integrator STREQUAL "explicit" AND
      NOT selected_chemistry_integrator STREQUAL "implicit")
    message(
      FATAL_ERROR
      "mechanism bundle chemistry_integrator must be 'explicit' or 'implicit'"
    )
  endif()

  foreach(json_array species reactions thermo transport)
    string(
      JSON json_array_length ERROR_VARIABLE bundle_json_error
      LENGTH "${bundle_json}" "${json_array}"
    )
    if(NOT bundle_json_error STREQUAL "NOTFOUND")
      message(
        FATAL_ERROR
        "selected mechanism build contract requires a ${json_array} array; "
        "legacy kinetics-only bundles are generator-only"
      )
    endif()
    set("bundle_n${json_array}" "${json_array_length}")
  endforeach()

  if(bundle_nspecies LESS 1 OR bundle_nreactions LESS 1)
    message(FATAL_ERROR "mechanism bundle must contain species and reactions")
  endif()
  if(NOT bundle_nthermo EQUAL bundle_nspecies OR
      NOT bundle_ntransport EQUAL bundle_nspecies)
    message(
      FATAL_ERROR
      "selected mechanism build contract requires one thermo and transport "
      "record per species"
    )
  endif()

  set("${variable_prefix}_json" "${bundle_json}" PARENT_SCOPE)
  set("${variable_prefix}_nspecies" "${bundle_nspecies}" PARENT_SCOPE)
  set("${variable_prefix}_nreactions" "${bundle_nreactions}" PARENT_SCOPE)
endfunction()

function(pelef_add_mechanism_bundle)
  set(options ADD_TEST INSTALL REQUIRE_ACTIVITY)
  set(
    one_value_args
    NAME BUNDLE SOURCE PROBE_NAME REACTOR_NAME REACTIVE_1D_NAME
    AMR_REACTIVE_1D_NAME REACTIVE_2D_NAME REACTIVE_EB_2D_NAME
    REACTIVE_EB_AMR_2D_NAME REACTIVE_3D_NAME AMR_REACTIVE_3D_NAME
    REACTIVE_EB_3D_NAME
    MPI_REACTIVE_1D_NAME MPI_AMR_REACTIVE_1D_NAME
    MPI_AMR_REACTIVE_3D_NAME
    MPI_REACTIVE_EB_PATCH_TREE_2D_NAME
  )
  cmake_parse_arguments(
    PELEF_BUNDLE
    "${options}"
    "${one_value_args}"
    ""
    ${ARGN}
  )

  if(PELEF_BUNDLE_UNPARSED_ARGUMENTS)
    message(
      FATAL_ERROR
      "unparsed pelef_add_mechanism_bundle arguments: "
      "${PELEF_BUNDLE_UNPARSED_ARGUMENTS}"
    )
  endif()
  if(NOT PELEF_BUNDLE_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism bundle NAME must be target-safe")
  endif()
  if(PELEF_BUNDLE_PROBE_NAME AND
      NOT PELEF_BUNDLE_PROBE_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism probe target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTOR_NAME AND
      NOT PELEF_BUNDLE_REACTOR_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactor target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTIVE_1D_NAME AND
      NOT PELEF_BUNDLE_REACTIVE_1D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactive 1D target name is invalid")
  endif()
  if(PELEF_BUNDLE_AMR_REACTIVE_1D_NAME AND
      NOT PELEF_BUNDLE_AMR_REACTIVE_1D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism AMR reactive 1D target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTIVE_2D_NAME AND
      NOT PELEF_BUNDLE_REACTIVE_2D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactive 2D target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_2D_NAME AND
      NOT PELEF_BUNDLE_REACTIVE_EB_2D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactive EB 2D target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME AND
      NOT PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactive EB AMR 2D target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTIVE_3D_NAME AND
      NOT PELEF_BUNDLE_REACTIVE_3D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactive 3D target name is invalid")
  endif()
  if(PELEF_BUNDLE_AMR_REACTIVE_3D_NAME AND
      NOT PELEF_BUNDLE_AMR_REACTIVE_3D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism AMR reactive 3D target name is invalid")
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_3D_NAME AND
      NOT PELEF_BUNDLE_REACTIVE_EB_3D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism reactive EB 3D target name is invalid")
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_1D_NAME AND
      NOT PELEF_BUNDLE_MPI_REACTIVE_1D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism MPI reactive 1D target name is invalid")
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME AND
      NOT PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism MPI AMR reactive 1D target name is invalid")
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME AND
      NOT PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(FATAL_ERROR "mechanism MPI AMR reactive 3D target name is invalid")
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME AND
      NOT PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
    message(
      FATAL_ERROR "mechanism sparse MPI reactive EB 2D target name is invalid"
    )
  endif()
  if(NOT PELEF_BUNDLE_BUNDLE)
    message(FATAL_ERROR "mechanism bundle BUNDLE is required")
  endif()
  if(NOT TARGET pelef_selected_runtime_core)
    message(
      FATAL_ERROR
      "pelef_selected_runtime_core must exist before adding a mechanism bundle"
    )
  endif()
  if(NOT Python3_EXECUTABLE)
    message(FATAL_ERROR "Python3 is required to generate a mechanism bundle")
  endif()

  get_filename_component(
    bundle_path "${PELEF_BUNDLE_BUNDLE}" ABSOLUTE
    BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  set_property(
    DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${bundle_path}"
  )
  _pelef_read_mechanism_bundle_build_contract("${bundle_path}" bundle)
  file(SHA256 "${bundle_path}" bundle_sha256)
  if(PELEF_BUNDLE_REACTOR_NAME AND bundle_nspecies LESS 2)
    message(FATAL_ERROR "selected reactor requires at least two species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_1D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected reactive 1D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_1D_NAME AND
      NOT TARGET pelef_selected_reactive_1d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_reactive_1d_runtime must exist before adding a "
      "selected reactive 1D application"
    )
  endif()
  if(PELEF_BUNDLE_AMR_REACTIVE_1D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected AMR reactive 1D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_AMR_REACTIVE_1D_NAME AND
      NOT TARGET pelef_selected_amr_reactive_1d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_amr_reactive_1d_runtime must exist before adding a "
      "selected AMR reactive 1D application"
    )
  endif()
  if(PELEF_BUNDLE_REACTIVE_2D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected reactive 2D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_2D_NAME AND
      NOT TARGET pelef_selected_reactive_2d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_reactive_2d_runtime must exist before adding a "
      "selected reactive 2D application"
    )
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_2D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected reactive EB 2D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_2D_NAME AND
      NOT TARGET pelef_selected_reactive_eb_2d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_reactive_eb_2d_runtime must exist before adding a "
      "selected reactive EB 2D application"
    )
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected reactive EB AMR 2D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME AND
      NOT TARGET pelef_selected_reactive_eb_amr_2d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_reactive_eb_amr_2d_runtime must exist before adding a "
      "selected reactive EB AMR 2D application"
    )
  endif()
  if(PELEF_BUNDLE_REACTIVE_3D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected reactive 3D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_3D_NAME AND
      NOT TARGET pelef_selected_reactive_3d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_reactive_3d_runtime must exist before adding a "
      "selected reactive 3D application"
    )
  endif()
  if(PELEF_BUNDLE_AMR_REACTIVE_3D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected AMR reactive 3D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_AMR_REACTIVE_3D_NAME AND
      NOT TARGET pelef_selected_amr_reactive_3d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_amr_reactive_3d_runtime must exist before adding a "
      "selected AMR reactive 3D application"
    )
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_3D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected reactive EB 3D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_REACTIVE_EB_3D_NAME AND
      NOT TARGET pelef_selected_reactive_eb_3d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_reactive_eb_3d_runtime must exist before adding a "
      "selected reactive EB 3D application"
    )
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_1D_NAME AND
      NOT bundle_nspecies EQUAL 2 AND NOT bundle_nspecies EQUAL 10)
    message(
      FATAL_ERROR
      "selected MPI reactive 1D requires the qualified 2- or 10-species "
      "initialization profile"
    )
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_1D_NAME AND NOT PELEF_ENABLE_MPI)
    message(FATAL_ERROR "selected MPI reactive 1D requires PELEF_ENABLE_MPI")
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_1D_NAME AND
      NOT TARGET pelef_selected_mpi_reactive_1d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_mpi_reactive_1d_runtime must exist before adding a "
      "selected MPI reactive 1D application"
    )
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected MPI AMR reactive 1D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME AND NOT PELEF_ENABLE_MPI)
    message(
      FATAL_ERROR "selected MPI AMR reactive 1D requires PELEF_ENABLE_MPI"
    )
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME AND
      NOT TARGET pelef_selected_mpi_amr_reactive_1d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_mpi_amr_reactive_1d_runtime must exist before adding a "
      "selected MPI AMR reactive 1D application"
    )
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(FATAL_ERROR "selected MPI AMR reactive 3D requires 2--32 species")
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME AND NOT PELEF_ENABLE_MPI)
    message(
      FATAL_ERROR "selected MPI AMR reactive 3D requires PELEF_ENABLE_MPI"
    )
  endif()
  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME AND
      NOT TARGET pelef_selected_mpi_amr_reactive_3d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_mpi_amr_reactive_3d_runtime must exist before adding a "
      "selected MPI AMR reactive 3D application"
    )
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME AND
      (bundle_nspecies LESS 2 OR bundle_nspecies GREATER 32))
    message(
      FATAL_ERROR "selected sparse MPI reactive EB 2D requires 2--32 species"
    )
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME AND
      NOT PELEF_ENABLE_MPI)
    message(
      FATAL_ERROR
      "selected sparse MPI reactive EB 2D requires PELEF_ENABLE_MPI"
    )
  endif()
  if(PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME AND
      NOT TARGET pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime)
    message(
      FATAL_ERROR
      "pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime must exist before "
      "adding a selected sparse MPI reactive EB 2D application"
    )
  endif()

  set(source_check_command)
  set(source_dependency)
  if(PELEF_BUNDLE_SOURCE)
    get_filename_component(
      source_path "${PELEF_BUNDLE_SOURCE}" ABSOLUTE
      BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    if(NOT EXISTS "${source_path}")
      message(FATAL_ERROR "mechanism source does not exist: ${source_path}")
    endif()
    string(JSON source_sha256 GET "${bundle_json}" source sha256)
    string(JSON source_file GET "${bundle_json}" source file)
    file(SHA256 "${source_path}" actual_source_sha256)
    get_filename_component(actual_source_file "${source_path}" NAME)
    if(NOT actual_source_sha256 STREQUAL source_sha256)
      message(FATAL_ERROR "selected mechanism source SHA-256 does not match bundle")
    endif()
    if(NOT actual_source_file STREQUAL source_file)
      message(FATAL_ERROR "selected mechanism source basename does not match bundle")
    endif()
    set_property(
      DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${source_path}"
    )
    set(
      source_check_command
      COMMAND
        "${CMAKE_COMMAND}"
        "-DPELEF_MECHANISM_SOURCE=${source_path}"
        "-DPELEF_MECHANISM_SHA256=${source_sha256}"
        "-DPELEF_MECHANISM_SOURCE_FILE=${source_file}"
        -P "${PROJECT_SOURCE_DIR}/cmake/CheckMechanismSource.cmake"
    )
    set(source_dependency "${source_path}")
  endif()

  set(bundle_directory "${CMAKE_CURRENT_BINARY_DIR}/${PELEF_BUNDLE_NAME}_bundle")
  set(generated_source "${bundle_directory}/${PELEF_BUNDLE_NAME}_mechanism.F90")
  add_custom_command(
    OUTPUT "${generated_source}"
    ${source_check_command}
    COMMAND
      "${Python3_EXECUTABLE}"
      "${PROJECT_SOURCE_DIR}/tools/generate_elementary_mechanism.py"
      --input "${bundle_path}"
      --output "${generated_source}"
    DEPENDS
      "${bundle_path}"
      "${PROJECT_SOURCE_DIR}/tools/generate_elementary_mechanism.py"
      "${PROJECT_SOURCE_DIR}/cmake/CheckMechanismSource.cmake"
      ${source_dependency}
    COMMENT "Generating ${PELEF_BUNDLE_NAME} mechanism bundle"
    VERBATIM
  )
  set_source_files_properties("${generated_source}" PROPERTIES GENERATED TRUE)

  set(mechanism_target "${PELEF_BUNDLE_NAME}_mechanism")
  if(TARGET "${mechanism_target}")
    message(FATAL_ERROR "mechanism target already exists: ${mechanism_target}")
  endif()
  add_library("${mechanism_target}" STATIC "${generated_source}")
  target_link_libraries("${mechanism_target}" PUBLIC pelef_selected_runtime_core)
  set(bundle_module_directory "${bundle_directory}/modules")
  set_target_properties(
    "${mechanism_target}"
    PROPERTIES Fortran_MODULE_DIRECTORY "${bundle_module_directory}"
  )
  target_include_directories(
    "${mechanism_target}"
    PUBLIC "$<BUILD_INTERFACE:${bundle_module_directory}>"
  )

  set(PELEF_BUNDLE_MODULE_NAME "${bundle_module_name}")
  set(PELEF_BUNDLE_LOADER_NAME "${bundle_loader_name}")
  set(PELEF_BUNDLE_KERNEL_NAME "${bundle_kernel_name}")
  set(PELEF_BUNDLE_JACOBIAN_NAME "${bundle_jacobian_name}")
  set(PELEF_BUNDLE_THERMO_LOADER_NAME "${bundle_thermo_loader_name}")
  set(PELEF_BUNDLE_TRANSPORT_LOADER_NAME "${bundle_transport_loader_name}")
  set(PELEF_BUNDLE_NSPECIES_SYMBOL "${bundle_symbol_prefix}_nspecies")
  set(PELEF_BUNDLE_NREACTIONS_SYMBOL "${bundle_symbol_prefix}_nreactions")
  set(
    PELEF_BUNDLE_CHEMISTRY_INTEGRATOR_SYMBOL
    "${bundle_symbol_prefix}_chemistry_integrator"
  )
  set(PELEF_BUNDLE_BUNDLE_SHA256 "${bundle_sha256}")
  if(PELEF_BUNDLE_REQUIRE_ACTIVITY)
    set(PELEF_BUNDLE_REQUIRE_ACTIVITY_LITERAL ".true.")
  else()
    set(PELEF_BUNDLE_REQUIRE_ACTIVITY_LITERAL ".false.")
  endif()
  set(configured_probe "${bundle_directory}/${PELEF_BUNDLE_NAME}_probe.F90")
  configure_file(
    "${PROJECT_SOURCE_DIR}/app/pelef_mechanism_probe.F90.in"
    "${configured_probe}"
    @ONLY
  )

  if(PELEF_BUNDLE_PROBE_NAME)
    set(probe_target "${PELEF_BUNDLE_PROBE_NAME}")
  else()
    set(probe_target "${PELEF_BUNDLE_NAME}_probe")
  endif()
  if(TARGET "${probe_target}")
    message(FATAL_ERROR "mechanism probe target already exists: ${probe_target}")
  endif()
  add_executable("${probe_target}" "${configured_probe}")
  target_link_libraries("${probe_target}" PRIVATE "${mechanism_target}")

  set(install_targets "${probe_target}")
  if(PELEF_BUNDLE_REACTOR_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTOR_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactor target already exists: ${PELEF_BUNDLE_REACTOR_NAME}"
      )
    endif()
    set(
      configured_reactor
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactor.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef0d_selected.F90.in"
      "${configured_reactor}"
      @ONLY
    )
    add_executable("${PELEF_BUNDLE_REACTOR_NAME}" "${configured_reactor}")
    target_link_libraries(
      "${PELEF_BUNDLE_REACTOR_NAME}" PRIVATE "${mechanism_target}"
    )
    if(PELEF_ENABLE_SUNDIALS)
      target_compile_definitions(
        "${PELEF_BUNDLE_REACTOR_NAME}" PRIVATE PELEF_HAVE_SUNDIALS=1
      )
      target_link_libraries(
        "${PELEF_BUNDLE_REACTOR_NAME}" PRIVATE pelef_selected_cvode_runtime
      )
    endif()
    list(APPEND install_targets "${PELEF_BUNDLE_REACTOR_NAME}")
  endif()

  if(PELEF_BUNDLE_REACTIVE_1D_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTIVE_1D_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactive 1D target already exists: "
        "${PELEF_BUNDLE_REACTIVE_1D_NAME}"
      )
    endif()
    set(
      configured_reactive_1d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactive_1d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_reactive_1d_selected.F90.in"
      "${configured_reactive_1d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_REACTIVE_1D_NAME}" "${configured_reactive_1d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_REACTIVE_1D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_reactive_1d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_REACTIVE_1D_NAME}")
  endif()

  if(PELEF_BUNDLE_AMR_REACTIVE_1D_NAME)
    if(TARGET "${PELEF_BUNDLE_AMR_REACTIVE_1D_NAME}")
      message(
        FATAL_ERROR
        "mechanism AMR reactive 1D target already exists: "
        "${PELEF_BUNDLE_AMR_REACTIVE_1D_NAME}"
      )
    endif()
    set(
      configured_amr_reactive_1d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_amr_reactive_1d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_amr_reactive_1d_selected.F90.in"
      "${configured_amr_reactive_1d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_AMR_REACTIVE_1D_NAME}"
      "${configured_amr_reactive_1d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_AMR_REACTIVE_1D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_amr_reactive_1d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_AMR_REACTIVE_1D_NAME}")
  endif()

  if(PELEF_BUNDLE_MPI_REACTIVE_1D_NAME)
    if(TARGET "${PELEF_BUNDLE_MPI_REACTIVE_1D_NAME}")
      message(
        FATAL_ERROR
        "mechanism MPI reactive 1D target already exists: "
        "${PELEF_BUNDLE_MPI_REACTIVE_1D_NAME}"
      )
    endif()
    set(
      configured_mpi_reactive_1d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_mpi_reactive_1d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_mpi_reactive_1d_selected.F90.in"
      "${configured_mpi_reactive_1d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_MPI_REACTIVE_1D_NAME}"
      "${configured_mpi_reactive_1d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_MPI_REACTIVE_1D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_mpi_reactive_1d_runtime
    )
    pelef_require_coherent_mpi_runtime(
      "${PELEF_BUNDLE_MPI_REACTIVE_1D_NAME}"
    )
    list(APPEND install_targets "${PELEF_BUNDLE_MPI_REACTIVE_1D_NAME}")
  endif()

  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME)
    if(TARGET "${PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME}")
      message(
        FATAL_ERROR
        "mechanism MPI AMR reactive 1D target already exists: "
        "${PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME}"
      )
    endif()
    set(
      configured_mpi_amr_reactive_1d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_mpi_amr_reactive_1d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_mpi_amr_reactive_1d_selected.F90.in"
      "${configured_mpi_amr_reactive_1d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME}"
      "${configured_mpi_amr_reactive_1d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_mpi_amr_reactive_1d_runtime
    )
    pelef_require_coherent_mpi_runtime(
      "${PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME}"
    )
    list(
      APPEND install_targets "${PELEF_BUNDLE_MPI_AMR_REACTIVE_1D_NAME}"
    )
  endif()

  if(PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME)
    if(TARGET "${PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME}")
      message(
        FATAL_ERROR
        "mechanism MPI AMR reactive 3D target already exists: "
        "${PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME}"
      )
    endif()
    set(
      configured_mpi_amr_reactive_3d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_mpi_amr_reactive_3d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_mpi_amr_reactive_3d_selected.F90.in"
      "${configured_mpi_amr_reactive_3d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME}"
      "${configured_mpi_amr_reactive_3d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_mpi_amr_reactive_3d_runtime
    )
    pelef_require_coherent_mpi_runtime(
      "${PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME}"
    )
    list(
      APPEND install_targets "${PELEF_BUNDLE_MPI_AMR_REACTIVE_3D_NAME}"
    )
  endif()

  if(PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME)
    if(TARGET "${PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME}")
      message(
        FATAL_ERROR
        "mechanism sparse MPI reactive EB 2D target already exists: "
        "${PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME}"
      )
    endif()
    set(
      configured_mpi_reactive_eb_patch_tree_2d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_mpi_reactive_eb_patch_tree_2d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_mpi_reactive_eb_patch_tree_2d_selected.F90.in"
      "${configured_mpi_reactive_eb_patch_tree_2d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME}"
      "${configured_mpi_reactive_eb_patch_tree_2d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_mpi_reactive_eb_patch_tree_2d_runtime
    )
    pelef_require_coherent_mpi_runtime(
      "${PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME}"
    )
    list(
      APPEND install_targets
      "${PELEF_BUNDLE_MPI_REACTIVE_EB_PATCH_TREE_2D_NAME}"
    )
  endif()

  if(PELEF_BUNDLE_REACTIVE_2D_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTIVE_2D_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactive 2D target already exists: "
        "${PELEF_BUNDLE_REACTIVE_2D_NAME}"
      )
    endif()
    set(
      configured_reactive_2d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactive_2d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_reactive_2d_selected.F90.in"
      "${configured_reactive_2d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_REACTIVE_2D_NAME}" "${configured_reactive_2d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_REACTIVE_2D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_reactive_2d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_REACTIVE_2D_NAME}")
  endif()

  if(PELEF_BUNDLE_REACTIVE_EB_2D_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTIVE_EB_2D_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactive EB 2D target already exists: "
        "${PELEF_BUNDLE_REACTIVE_EB_2D_NAME}"
      )
    endif()
    set(
      configured_reactive_eb_2d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactive_eb_2d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_reactive_eb_2d_selected.F90.in"
      "${configured_reactive_eb_2d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_REACTIVE_EB_2D_NAME}"
      "${configured_reactive_eb_2d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_REACTIVE_EB_2D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_reactive_eb_2d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_REACTIVE_EB_2D_NAME}")
  endif()

  if(PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactive EB AMR 2D target already exists: "
        "${PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME}"
      )
    endif()
    set(
      configured_reactive_eb_amr_2d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactive_eb_amr_2d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_reactive_eb_amr_2d_selected.F90.in"
      "${configured_reactive_eb_amr_2d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME}"
      "${configured_reactive_eb_amr_2d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_reactive_eb_amr_2d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_REACTIVE_EB_AMR_2D_NAME}")
  endif()

  if(PELEF_BUNDLE_REACTIVE_3D_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTIVE_3D_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactive 3D target already exists: "
        "${PELEF_BUNDLE_REACTIVE_3D_NAME}"
      )
    endif()
    set(
      configured_reactive_3d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactive_3d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_reactive_3d_selected.F90.in"
      "${configured_reactive_3d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_REACTIVE_3D_NAME}" "${configured_reactive_3d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_REACTIVE_3D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_reactive_3d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_REACTIVE_3D_NAME}")
  endif()

  if(PELEF_BUNDLE_AMR_REACTIVE_3D_NAME)
    if(TARGET "${PELEF_BUNDLE_AMR_REACTIVE_3D_NAME}")
      message(
        FATAL_ERROR
        "mechanism AMR reactive 3D target already exists: "
        "${PELEF_BUNDLE_AMR_REACTIVE_3D_NAME}"
      )
    endif()
    set(
      configured_amr_reactive_3d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_amr_reactive_3d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_amr_reactive_3d_selected.F90.in"
      "${configured_amr_reactive_3d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_AMR_REACTIVE_3D_NAME}"
      "${configured_amr_reactive_3d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_AMR_REACTIVE_3D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_amr_reactive_3d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_AMR_REACTIVE_3D_NAME}")
  endif()

  if(PELEF_BUNDLE_REACTIVE_EB_3D_NAME)
    if(TARGET "${PELEF_BUNDLE_REACTIVE_EB_3D_NAME}")
      message(
        FATAL_ERROR
        "mechanism reactive EB 3D target already exists: "
        "${PELEF_BUNDLE_REACTIVE_EB_3D_NAME}"
      )
    endif()
    set(
      configured_reactive_eb_3d
      "${bundle_directory}/${PELEF_BUNDLE_NAME}_reactive_eb_3d.F90"
    )
    configure_file(
      "${PROJECT_SOURCE_DIR}/app/pelef_reactive_eb_3d_selected.F90.in"
      "${configured_reactive_eb_3d}"
      @ONLY
    )
    add_executable(
      "${PELEF_BUNDLE_REACTIVE_EB_3D_NAME}"
      "${configured_reactive_eb_3d}"
    )
    target_link_libraries(
      "${PELEF_BUNDLE_REACTIVE_EB_3D_NAME}"
      PRIVATE
        "${mechanism_target}"
        pelef_selected_reactive_eb_3d_runtime
    )
    list(APPEND install_targets "${PELEF_BUNDLE_REACTIVE_EB_3D_NAME}")
  endif()

  if(PELEF_BUNDLE_ADD_TEST)
    add_test(NAME "unit_${PELEF_BUNDLE_NAME}_bundle" COMMAND "${probe_target}")
  endif()
  if(PELEF_BUNDLE_INSTALL)
    install(TARGETS ${install_targets} RUNTIME DESTINATION bin)
  endif()
endfunction()
