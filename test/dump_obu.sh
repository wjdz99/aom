#!/bin/sh
## Copyright (c) 2018, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
## This file tests the libaom dump_obu tool. To add new tests to this
## file, do the following:
##   1. Write a shell function (this is your test).
##   2. Add the function to dump_obu_tests (on a new line).
##
. $(dirname $0)/tools_common.sh

dump_obu_verify_environment() {
  echo HI
  if [ ! -e "${YUV_RAW_INPUT}" ]; then
    elog "The file ${YUV_RAW_INPUT##*/} must exist in LIBAOM_TEST_DATA_PATH."
    return 1
  fi
  if [ "$(dump_obu_available)" = "yes" ]; then
    if [ -z "$(aom_tool_path dump_obu)" ]; then
      elog "dump_obu not found in LIBAOM_BIN_PATH, its parent, or child tools/."
    fi
  fi
}

dump_obu_available() {
  if [ "$(av1_decode_available)" = "yes" ] && \
     [ "$(av1_encode_available)" = "yes" ]; then
    echo yes
  fi
}

dump_obu() {
  return 0
}

dump_obu_tests="dump_obu"

run_tests dump_obu_verify_environment "${dump_obu_tests}"
