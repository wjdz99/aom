#
# Copyright (c) 2020, Alliance for Open Media. All rights reserved
#
# This source code is subject to the terms of the BSD 2 Clause License and the
# Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License was
# not distributed with this source code in the LICENSE file, you can obtain it
# at www.aomedia.org/license/software. If the Alliance for Open Media Patent
# License 1.0 was not distributed with this source code in the PATENTS file, you
# can obtain it at www.aomedia.org/license/patent.
#
if(AOM_BUILD_CMAKE_TENSORFLOW_LITE_CMAKE_)
  return()
endif() # AOM_BUILD_CMAKE_TENSORFLOW_LITE_CMAKE_
set(AOM_BUILD_CMAKE_TENSORFLOW_LITE_CMAKE_ 1)

include(ExternalProject)
include(FindGit)

# Checks if Tensorflow has been checked out -- if not, uses the git submodule
# command to fetch it.
function(checkout_submodule_)
  # As a quick sanity check, see if at least 1 expected file or directory is
  # present in each submodule. If so, assume they are all checked out (if they
  # are not, then the base directory will be empty).
  if(EXISTS "${AOM_ROOT}/third_party/tensorflow/tensorflow")
    return()
  endif()
  if(NOT GIT_FOUND)
    message(
      FATAL_ERROR
        "Tensorflow-Lite not present; " "git could not be found; "
        "please check out submodules with 'git submodule update --init'")
  endif()
  # Note that "git submodule update --init" must be run from inside the git
  # repository; the --git-dir flag does not work.
  message("Checking out Tensorflow-Lite submodule")
  execute_process(COMMAND "${GIT_EXECUTABLE}" submodule update --init
                  WORKING_DIRECTORY "${AOM_ROOT}"
                  OUTPUT_VARIABLE submodule_out
                  ERROR_VARIABLE submodule_err
                  RESULT_VARIABLE submodule_result)
  if(NOT ${submodule_result} EQUAL 0)
    message(FATAL_ERROR "Unable to run 'git submodule update --init': "
                        "Return code: " ${submodule_result} ", STDOUT: "
                        ${submodule_out} ", STDERR: " ${submodule_err})
  endif()
endfunction()

function(target_link_tf_lite_dep_ target subdir libname)
  if(NOT (("${subdir}" STREQUAL "") OR ("${subdir}" MATCHES "/$")))
    message(
      FATAL_ERROR "sub-directory must be empty or end with a slash: ${subdir}")
  endif()

  set(STATIC_LIBRARY_DIR "")
  if(MSVC)
    set(STATIC_LIBRARY_DIR "$<CONFIG>/")
  endif()
  target_link_libraries(
    ${target}
    PRIVATE
      "${CMAKE_BINARY_DIR}/tensorflow_lite/${subdir}${STATIC_LIBRARY_DIR}${CMAKE_STATIC_LIBRARY_PREFIX}${libname}${CMAKE_STATIC_LIBRARY_SUFFIX}"
    )
endfunction()

# Add TF-lite libraries onto the target at link time. For enabling TF-lite for
# experiments, prefer the "experiment_requires_tf_lite" function.
function(target_link_tf_lite_libraries target)
  target_link_libraries(${target} ${AOM_LIB_LINK_TYPE} Threads::Threads)
  target_link_libraries(${target} PRIVATE ${CMAKE_DL_LIBS})
  target_link_tf_lite_dep_(${target} "" tensorflow-lite)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                           absl_flags)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                           absl_flags_internal)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                           absl_flags_registry)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                           absl_flags_config)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                           absl_flags_program_name)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                           absl_flags_marshalling)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/hash/
                           absl_hash)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/hash/
                           absl_city)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/status/
                           absl_status)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/types/
                           absl_bad_optional_access)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/
                           absl_cord)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/
                           absl_str_format_internal)
  target_link_tf_lite_dep_(${target}
                           _deps/abseil-cpp-build/absl/synchronization/
                           absl_synchronization)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                           absl_stacktrace)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                           absl_symbolize)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                           absl_debugging_internal)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                           absl_demangle_internal)
  target_link_tf_lite_dep_(${target}
                           _deps/abseil-cpp-build/absl/synchronization/
                           absl_graphcycles_internal)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_malloc_internal)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/time/
                           absl_time)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/
                           absl_strings)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/
                           absl_strings_internal)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_throw_delegate)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_base)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_dynamic_annotations)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_spinlock_wait)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/numeric/
                           absl_int128)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/time/
                           absl_civil_time)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/time/
                           absl_time_zone)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/types/
                           absl_bad_variant_access)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_raw_logging_internal)
  target_link_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                           absl_log_severity)
  target_link_tf_lite_dep_(${target} _deps/farmhash-build/ farmhash)
  target_link_tf_lite_dep_(${target} _deps/fft2d-build/ fft2d_fftsg2d)
  target_link_tf_lite_dep_(${target} _deps/fft2d-build/ fft2d_fftsg)
  target_link_tf_lite_dep_(${target} _deps/flatbuffers-build/ flatbuffers)
  target_link_tf_lite_dep_(${target} _deps/ruy-build/ ruy)
endfunction()

# If Tensorflow-Lite should be enabled, adds appropriate build rules / targets.
function(setup_tensorflow_lite)
  if("${AOM_ROOT}" STREQUAL "")
    message(FATAL_ERROR "AOM_ROOT variable must not be empty.")
  endif()
  # Cross-compile is not currently implemented.
  if(CMAKE_TOOLCHAIN_FILE)
    message("TOOLCHAIN: ${CMAKE_TOOLCHAIN_FILE}")
    message(WARNING "No cross-compile support for TensorFlow Lite; disabling")
    set(CONFIG_TENSORFLOW_LITE 0 PARENT_SCOPE)
    return()
  endif()
  if(MSVC)
    add_compile_definitions(NOMINMAX=1)
  endif()
  checkout_submodule_()

  # Allow code to reference TF.
  include_directories("${AOM_ROOT}/third_party/tensorflow")

  if(NOT (EXISTS "${CMAKE_BINARY_DIR}/tensorflow"))
    execute_process(COMMAND "${GIT_EXECUTABLE}" apply
                            "${AOM_ROOT}/build/cmake/tflite.patch"
                    WORKING_DIRECTORY "${AOM_ROOT}/third_party/tensorflow"
                    OUTPUT_VARIABLE submodule_out
                    ERROR_VARIABLE submodule_err
                    RESULT_VARIABLE submodule_result)
    if(NOT ${submodule_result} EQUAL 0)
      message(FATAL_ERROR "Unable to apply CMake patch to TF-Lite: "
                          "Return code: " ${submodule_result} ", STDOUT: "
                          ${submodule_out} ", STDERR: " ${submodule_err})
    endif()
    file(COPY "${AOM_ROOT}/third_party/tensorflow/"
         DESTINATION "${CMAKE_BINARY_DIR}/tensorflow")
    execute_process(COMMAND "${GIT_EXECUTABLE}" reset HEAD --hard
                    WORKING_DIRECTORY "${AOM_ROOT}/third_party/tensorflow"
                    OUTPUT_VARIABLE submodule_out
                    ERROR_VARIABLE submodule_err
                    RESULT_VARIABLE submodule_result)
    if(NOT ${submodule_result} EQUAL 0)
      message(FATAL_ERROR "Unable to apply CMake patch to TF-Lite: "
                          "Return code: " ${submodule_result} ", STDOUT: "
                          ${submodule_out} ", STDERR: " ${submodule_err})
    endif()

  endif()
  externalproject_add(
    tensorflow_lite
    SOURCE_DIR "${CMAKE_BINARY_DIR}/tensorflow/tensorflow/lite"
    PREFIX "${CMAKE_BINARY_DIR}/tensorflow_lite"
    BINARY_DIR "${CMAKE_BINARY_DIR}/tensorflow_lite"
    DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/tensorflow_lite"
    LOG_BUILD 1)

  # TF-Lite depends on this, and downloads it during compilation.
  include_directories(
    "${CMAKE_CURRENT_BINARY_DIR}/tensorflow_lite/flatbuffers/include/")

  add_dependencies(aom_av1_common tensorflow_lite)
  foreach(aom_app ${AOM_APP_TARGETS})
    add_dependencies(${aom_app} tensorflow_lite)
    target_link_tf_lite_libraries(${aom_app})
  endforeach()
endfunction()

# Signal that the experiment needs TF-lite enabled.
function(experiment_requires_tf_lite experiment_name)
  # Experiment is not enabled, no need to include TF-Lite in the build.
  if(${${experiment_name}} EQUAL "0")
    return()
  endif()
  # Disable the experiment so Gerrit will not test this case.
  if(CMAKE_TOOLCHAIN_FILE)
    message(WARNING "--- Cross-compile support not implemented for TF-Lite. "
                    "Disabling ${experiment_name}.")
    set(${experiment_name} 0 PARENT_SCOPE)
    return()
  endif()
  if(NOT CONFIG_TENSORFLOW_LITE)
    set(CONFIG_TENSORFLOW_LITE 1 PARENT_SCOPE)
  endif()
endfunction()
