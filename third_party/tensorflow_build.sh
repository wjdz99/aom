#!/bin/bash
#
# Copyright (c) 2019, Alliance for Open Media. All rights reserved
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at www.aomedia.org/license/software. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at www.aomedia.org/license/patent.
#
###########################################################################
#
# Script to build the static TensorFlow lite library.
#
# Tensorflow's build process generates the static libraries in the same
# directory as the source code. AOM, however, generates binaries/objects/etc.
# in a different directory. This script:
#
# 1.) Copies the Tensorflow code to a temporary directory
# 2.) Downloads dependencies / compiles it there
# 3.) Copies over the TensorFlow lite static library
#
# This script has numerous limitations:
#
# 1.) It is only tested on Linux
# 2.) The included static library is ~6mb in size
# 3.) It may require a full recompile of Tensorflow on every make, even if
#     no changes to the TensorFlow files are made
#
# As such, it is meant as a starting point for integration.

set -ue

function cleanup() {
  local ret=$?
  rm -rf "${temp_dir}"
  if [[ $ret -ne 0 ]]; then
    echo "TensorFlow Lite: unable to delete temporary directory ${temp_dir}"
  fi
}

if [[ $# -ne 2 ]]; then
  echo "Usage: aom_build.sh /path/to/tensorflow /output/file"
  exit 1
fi

trap cleanup EXIT INT TERM QUIT

temp_dir=$(mktemp -d -t tensorflowlite.XXXXXXXX)

cp -r "$1" "${temp_dir}"
"${temp_dir}/tensorflow/tensorflow/lite/tools/make/download_dependencies.sh"
"${temp_dir}/tensorflow/tensorflow/lite/tools/make/build_lib.sh"
cp "${temp_dir}/tensorflow/tensorflow/lite/tools/make/gen/linux_x86_64/lib/libtensorflow-lite.a" "$2"
