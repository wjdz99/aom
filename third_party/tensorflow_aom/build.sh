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
#
# Note that we do not use the download_dependencies.sh script directly, as
# it does not perform checksumming on the downloaded code, and it downloads
# directly into the source directory.
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

trap cleanup EXIT INT TERM QUIT

function download_fft2d() {
  curl --location --silent "https://storage.googleapis.com/mirror.tensorflow.org/www.kurims.kyoto-u.ac.jp/~ooura/fft2d.tgz" \
  --output "${temp_make_dir}/downloads/fft2d.tgz"
  checksum=$(sha256sum "${temp_make_dir}/downloads/fft2d.tgz" | cut -d ' ' -f 1)
  if [[ ${checksum} != "ada7e99087c4ed477bfdf11413f2ba8db8a840ba9bbf8ac94f4f3972e2a7cec9" ]]; then
    echo "Bad checksum for fft2d.tgz: ${checksum}"
    exit 1
  fi
  tar -C "${temp_make_dir}/downloads" -xzf \
    "${temp_make_dir}/downloads/fft2d.tgz"
}

function tweak_download_files() {
  # These tweaks are taken from
  # tensorflow/lite/tools/make/download_dependencies.sh
  local fname="${temp_make_dir}/downloads/eigen/Eigen/src/Core/arch/NEON/Complex.h"
  sed -i -e 's/static uint32x4_t p4ui_CONJ_XOR = vld1q_u32( conj_XOR_DATA );/static uint32x4_t p4ui_CONJ_XOR;/' \
    "${fname}"
  sed -i -e 's/static uint32x2_t p2ui_CONJ_XOR = vld1_u32( conj_XOR_DATA );/static uint32x2_t p2ui_CONJ_XOR;/' \
  "${fname}"
  sed -i -e 's/static uint64x2_t p2ul_CONJ_XOR = vld1q_u64( p2ul_conj_XOR_DATA );/static uint64x2_t p2ul_CONJ_XOR;/' \
  "${fname}"
}

if [[ $# -ne 2 ]]; then
  echo "Usage: aom_build.sh /path/to/tensorflow /output/file"
  exit 1
fi

cp -r "$1" "${temp_dir}"
readonly tf_aom_dir=$(dirname "$0")
readonly temp_make_dir="${temp_dir}/tensorflow/tensorflow/lite/tools/make"
mkdir "${temp_make_dir}/downloads"
cp -r "${tf_aom_dir}/ARM_NEON_2_x86_SSE" "${temp_make_dir}/downloads/neon_2_sse"
cp -r "${tf_aom_dir}/abseil-cpp" "${temp_make_dir}/downloads/absl"
cp -r "${tf_aom_dir}/eigen" "${temp_make_dir}/downloads/"
cp -r "${tf_aom_dir}/farmhash" "${temp_make_dir}/downloads/"
cp -r "${tf_aom_dir}/flatbuffers" "${temp_make_dir}/downloads/"
cp -r "${tf_aom_dir}/gemmlowp" "${temp_make_dir}/downloads/"
# Note that fft2d is not in a git repository, so it has been downloaded
# separately.
download_fft2d
# There are specific tweaks applied by TF-Lite: re-apply them here.
tweak_download_files

"${temp_dir}/tensorflow/tensorflow/lite/tools/make/build_lib.sh"
cp "${temp_dir}/tensorflow/tensorflow/lite/tools/make/gen/linux_x86_64/lib/libtensorflow-lite.a" "$2"
