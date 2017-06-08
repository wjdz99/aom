##
## Copyright (c) 2017, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
cmake_minimum_required(VERSION 3.5)

set(REQUIRED_ARGS "AOM_ROOT" "AOM_CONFIG_DIR" "AOM_TARGET_SYSTEM" "AOM_SYM_FILE"
    "CONFIG_AV1_DECODER" "CONFIG_AV1_ENCODER")

foreach (arg ${REQUIRED_ARGS})
  if ("${${arg}}" STREQUAL "")
    message(FATAL_ERROR "${arg} must not be empty.")
  endif ()
endforeach ()

if ("${AOM_TARGET_SYSTEM}" STREQUAL "Darwin")
  set(symbol_prefix "_")
elseif ("${AOM_TARGET_SYSTEM}" MATCHES "Windows\|MSYS")
  message(FATAL_ERROR "Windows DLL builds not supported yet.")
else ()
  set(symbol_suffix ";")
endif ()

