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

function(add_tf_lite_dep_ target subdir libname)
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

# Add Tensorflow-Lite as a dependency on the target. For experiments, prefer the
# "experiment_requires_tf_lite" function.
function(add_all_tf_lite_dependencies target)
  add_dependencies(${target} tensorflow_lite)
  target_link_libraries(${target} ${AOM_LIB_LINK_TYPE} Threads::Threads)
  target_link_libraries(${target} PRIVATE ${CMAKE_DL_LIBS})
  add_tf_lite_dep_(${target} "" tensorflow-lite)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/ absl_flags)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                   absl_flags_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                   absl_flags_registry)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                   absl_flags_config)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                   absl_flags_program_name)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/flags/
                   absl_flags_marshalling)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/hash/ absl_hash)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/hash/ absl_city)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/status/ absl_status)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/types/
                   absl_bad_optional_access)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/ absl_cord)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/
                   absl_str_format_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/synchronization/
                   absl_synchronization)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                   absl_stacktrace)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                   absl_symbolize)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                   absl_debugging_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/debugging/
                   absl_demangle_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/synchronization/
                   absl_graphcycles_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                   absl_malloc_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/time/ absl_time)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/ absl_strings)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/strings/
                   absl_strings_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                   absl_throw_delegate)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/ absl_base)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                   absl_dynamic_annotations)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                   absl_spinlock_wait)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/numeric/ absl_int128)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/time/ absl_civil_time)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/time/ absl_time_zone)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/types/
                   absl_bad_variant_access)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                   absl_raw_logging_internal)
  add_tf_lite_dep_(${target} _deps/abseil-cpp-build/absl/base/
                   absl_log_severity)
  add_tf_lite_dep_(${target} _deps/farmhash-build/ farmhash)
  add_tf_lite_dep_(${target} _deps/fft2d-build/ fft2d_fftsg2d)
  add_tf_lite_dep_(${target} _deps/fft2d-build/ fft2d_fftsg)
  add_tf_lite_dep_(${target} _deps/flatbuffers-build/ flatbuffers)
  add_tf_lite_dep_(${target} _deps/ruy-build/ ruy)
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

  externalproject_add(
    tensorflow_lite
    SOURCE_DIR "${AOM_ROOT}/third_party/tensorflow/tensorflow/lite"
    PREFIX "${CMAKE_BINARY_DIR}/tensorflow_lite"
    BINARY_DIR "${CMAKE_BINARY_DIR}/tensorflow_lite"
    DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/tensorflow_lite"
    LOG_BUILD 1)

  # TF-Lite depends on this, and downloads it during compilation.
  include_directories(
    "${CMAKE_CURRENT_BINARY_DIR}/tensorflow_lite/flatbuffers/include/")

  add_all_tf_lite_dependencies(aom)
  add_all_tf_lite_dependencies(aom_av1_common)
  if(ENABLE_EXAMPLES OR ENABLE_TESTS OR ENABLE_TOOLS)
    add_all_tf_lite_dependencies(aom_common_app_util)
  endif()
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
