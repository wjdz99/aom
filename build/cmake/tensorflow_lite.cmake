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

# Checks if Tensorflow has been checked out -- if not, uses  the git submodule
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

function(add_tf_lite_dependency target)
  add_dependencies(${target} tensorflow_lite)
  target_link_libraries(
    ${target}
    PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/tensorflow_lite/libtensorflow-lite.a"
            ${AOM_LIB_LINK_TYPE} Threads::Threads
    PRIVATE ${CMAKE_DL_LIBS})
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
  checkout_submodule_()

  # Allow code to reference TF.
  include_directories("${AOM_ROOT}/third_party/tensorflow")
  # TF-Lite depends on this, and downloads it during compilation.
  include_directories(
    "${AOM_ROOT}/third_party/tensorflow_dependencies/flatbuffers/include/")

  # Create a target for tf_lite. For windows, we need to define NOMINMAX=1 to
  # disable windows.h from including min/max macros, which conflict with the C++
  # std functions.
  add_compile_definitions(NOMINMAX=1)
  externalproject_add(
    tensorflow_lite
    SOURCE_DIR "${AOM_ROOT}/third_party/tensorflow/tensorflow/lite"
    PREFIX "${CMAKE_BINARY_DIR}/tensorflow_lite"
    BINARY_DIR "${CMAKE_BINARY_DIR}/tensorflow_lite"
    DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/tensorflow_lite"
    LOG_BUILD 1)

  # Add tensorflow-lite as a dependency on all AOM applications.
  foreach(aom_app ${AOM_APP_TARGETS})
    add_tf_lite_dependency(${aom_app})
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
