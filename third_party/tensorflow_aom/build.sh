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
# 2.) Copies the necessary dependencies into the temporary directory
# 2.) Compiles it there
# 3.) Copies over the TensorFlow lite static library

set -uxe

readonly temp_dir=$(mktemp -d -t tensorflowlite.XXXXXXXX)

function cleanup() {
  echo $temp_dir
  return
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

cp -r "$1" "${temp_dir}"
readonly tf_aom_dir=$(dirname "$0")
readonly temp_make_dir="${temp_dir}/tensorflow/tensorflow/lite/tools/make"
# Note that fft2d is not in a git repository, so it has been downloaded and
# archived in the downloads/ folder.
cp -r "${tf_aom_dir}/downloads" "${temp_make_dir}"
cp -r "${tf_aom_dir}/abseil-cpp" "${temp_make_dir}/downloads/absl"
cp -r "${tf_aom_dir}/eigen" "${temp_make_dir}/downloads/"
cp -r "${tf_aom_dir}/farmhash" "${temp_make_dir}/downloads/"
cp -r "${tf_aom_dir}/flatbuffers" "${temp_make_dir}/downloads/"
cp -r "${tf_aom_dir}/gemmlowp" "${temp_make_dir}/downloads/"
"${temp_dir}/tensorflow/tensorflow/lite/tools/make/build_lib.sh"
cp "${temp_dir}/tensorflow/tensorflow/lite/tools/make/gen/linux_x86_64/lib/libtensorflow-lite.a" "$2"
