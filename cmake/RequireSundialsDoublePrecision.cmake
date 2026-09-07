include_guard(GLOBAL)

function(pelef_require_sundials_double_precision config_header)
  if(NOT EXISTS "${config_header}")
    message(
      FATAL_ERROR
      "SUNDIALS precision header is unavailable: ${config_header}"
    )
  endif()

  file(
    STRINGS "${config_header}" precision_definitions
    REGEX
      "^[ \t]*#[ \t]*define[ \t]+SUNDIALS_(SINGLE|DOUBLE|EXTENDED)_PRECISION[ \t]+1[ \t]*$"
  )
  set(has_single_precision FALSE)
  set(has_double_precision FALSE)
  set(has_extended_precision FALSE)
  foreach(precision_definition IN LISTS precision_definitions)
    if(precision_definition MATCHES "SUNDIALS_SINGLE_PRECISION")
      set(has_single_precision TRUE)
    elseif(precision_definition MATCHES "SUNDIALS_DOUBLE_PRECISION")
      set(has_double_precision TRUE)
    elseif(precision_definition MATCHES "SUNDIALS_EXTENDED_PRECISION")
      set(has_extended_precision TRUE)
    endif()
  endforeach()

  if(
      NOT has_double_precision
      OR has_single_precision
      OR has_extended_precision
    )
    if(precision_definitions)
      string(JOIN ", " precision_summary ${precision_definitions})
    else()
      set(precision_summary "none")
    endif()
    message(
      FATAL_ERROR
      "PeleF requires a double-precision SUNDIALS ABI. "
      "Active precision definitions in '${config_header}': "
      "${precision_summary}"
    )
  endif()
endfunction()
