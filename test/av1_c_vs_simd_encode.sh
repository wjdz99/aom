#!/bin/sh
## Copyright (c) 2023, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
##  This script checks the bit exactness between C and SIMD
##  implementations of AV1 encoder.

TEST_BITRATES="1600 6400"
PRESETS="good rt"
TEST_CLIPS="yuv_raw_input yuv_480p_raw_input y4m_720p_input"
OUT_FILE_SUFFIX=".ivf"
SCRIPT_DIR=$(dirname "$0")
LIBAOM_SOURCE_DIR=$(cd ${SCRIPT_DIR}/..; pwd)
devnull='> /dev/null 2>&1'

# Clips used in test.
YUV_RAW_INPUT="${LIBAOM_TEST_DATA_PATH}/hantro_collage_w352h288.yuv"
YUV_480P_RAW_INPUT="${LIBAOM_TEST_DATA_PATH}/niklas_640_480_30.yuv"
Y4M_720P_INPUT="${LIBAOM_TEST_DATA_PATH}/niklas_1280_720_30.y4m"

# Number of frames to test.
AV1_ENCODE_C_VS_SIMD_TEST_FRAME_LIMIT=20

# Create a temporary directory for output files.
if [ -n "${TMPDIR}" ]; then
  AOM_TEST_TEMP_ROOT="${TMPDIR}"
elif [ -n "${TEMPDIR}" ]; then
  AOM_TEST_TEMP_ROOT="${TEMPDIR}"
else
  AOM_TEST_TEMP_ROOT=/tmp
fi

AOM_TEST_OUTPUT_DIR="${AOM_TEST_TEMP_ROOT}/av1_test_$$"

if ! mkdir -p "${AOM_TEST_OUTPUT_DIR}" || \
   [ ! -d "${AOM_TEST_OUTPUT_DIR}" ]; then
  echo "${0##*/}: Cannot create output directory, giving up."
  echo "${0##*/}:   AOM_TEST_OUTPUT_DIR=${AOM_TEST_OUTPUT_DIR}"
  exit 1
fi

elog() {
  echo "$@" 1>&2
}

# Echoes path to $1 when it's executable and exists in ${AOM_TEST_OUTPUT_DIR},
# or an empty string. Caller is responsible for testing the string once the
# function returns.
av1_enc_tool_path() {
  local target="$1"
  local tool_path="${AOM_TEST_OUTPUT_DIR}/build_target_${target}/aomenc"

  if [ ! -x "${tool_path}" ]; then
    tool_path=""
  fi
  echo "${tool_path}"
}

# Environment check: Make sure input and source directories are available.
av1_c_vs_simd_enc_verify_environment () {
  if [ ! -e "${YUV_RAW_INPUT}" ]; then
    elog "libaom test data must exist in LIBAOM_TEST_DATA_PATH."
    return 1
  fi
  if [ ! -e "${YUV_480P_RAW_INPUT}" ]; then
    elog "libaom test data must exist in LIBAOM_TEST_DATA_PATH."
    return 1
  fi
  if [ ! -e "${Y4M_720P_INPUT}" ]; then
    elog "libaom test data must exist in LIBAOM_TEST_DATA_PATH."
    return 1
  fi
  if [ ! -d "$LIBAOM_SOURCE_DIR" ]; then
    elog "LIBAOM_SOURCE_DIR does not exist."
    return 1
  fi
}

cleanup() {
  rm -rf  ${AOM_TEST_OUTPUT_DIR}
}

# Echo AOM_SIMD_CAPS_MASK for different instruction set architecture.
avx512f() {
   echo "0x0FF"
}

avx2() {
   echo "0x0FF"
}

avx() {
   echo "0x07F"
}

sse4_1() {
   echo "0x03F"
}

ssse3() {
   echo "0x01F"
}

sse3() {
   echo "0x00F"
}

sse2() {
   echo "0x007"
}

# Echo clip details to be used as input to aomenc.
yuv_raw_input() {
  echo ""${YUV_RAW_INPUT}"
       --width=352
       --height=288
       --bit-depth=8"
}

yuv_480p_raw_input() {
  echo ""${YUV_480P_RAW_INPUT}"
       --width=640
       --height=480
       --bit-depth=8"
}

y4m_720p_input() {
  echo ""${Y4M_720P_INPUT}"
       --bit-depth=8"
}

y4m_360p_10bit_input() {
  echo ""${Y4M_360P_10BIT_INPUT}"
       --bit-depth=10"
}

has_x86_isa_extn() {
  instruction_set=$1
  grep -q "$instruction_set" /proc/cpuinfo
  if [ $? -eq 1 ]; then
    return 1
  fi
}

# Echo good encode params for use with AV1 encoder.
av1_encode_good_params() {
  echo "--good \
  --ivf \
  --profile=0 \
  --static-thresh=0 \
  --threads=1 \
  --tile-columns=0 \
  --tile-rows=0 \
  --verbose \
  --end-usage=vbr \
  --kf-max-dist=160 \
  --kf-min-dist=0 \
  --max-q=63 \
  --min-q=0 \
  --overshoot-pct=100 \
  --undershoot-pct=100 \
  --passes=2 \
  --arnr-maxframes=7 \
  --arnr-strength=5 \
  --auto-alt-ref=1 \
  --drop-frame=0 \
  --frame-parallel=0 \
  --lag-in-frames=35 \
  --maxsection-pct=2000 \
  --minsection-pct=0 \
  --sharpness=0"
}

# Echo realtime encode params for use with AV1 encoder.
av1_encode_rt_params() {
  echo "--rt \
  --ivf \
  --profile=0 \
  --static-thresh=0 \
  --threads=1 \
  --tile-columns=0 \
  --tile-rows=0 \
  --verbose \
  --end-usage=cbr \
  --kf-max-dist=90000 \
  --max-q=58 \
  --min-q=2 \
  --overshoot-pct=50 \
  --undershoot-pct=50 \
  --passes=1 \
  --aq-mode=3 \
  --buf-initial-sz=500 \
  --buf-optimal-sz=600 \
  --buf-sz=1000 \
  --coeff-cost-upd-freq=3 \
  --dv-cost-upd-freq=3 \
  --mode-cost-upd-freq=3 \
  --mv-cost-upd-freq=3 \
  --deltaq-mode=0 \
  --enable-global-motion=0 \
  --enable-obmc=0 \
  --enable-order-hint=0 \
  --enable-ref-frame-mvs=0 \
  --enable-tpl-model=0 \
  --enable-warped-motion=0 \
  --lag-in-frames=0 \
  --max-intra-rate=300 \
  --noise-sensitivity=0"
}

# Configures for the given target in AOM_TEST_OUTPUT_DIR/build_target_${target}
# directory.
av1_enc_build() {
  local target="$1"
  local cmake_command="$2"
  local tmp_build_dir=${AOM_TEST_OUTPUT_DIR}/build_target_${target}
  mkdir -p $tmp_build_dir
  cd $tmp_build_dir

  echo "Building target: ${target}"
  local cmake_extra_args="-DCONFIG_EXCLUDE_SIMD_MISMATCH=1 \
             -DCMAKE_BUILD_TYPE=Release \
             -DENABLE_CCACHE=1 \
             '-DCMAKE_C_FLAGS_RELEASE=-O3 -g' \
             '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -g' \
             -DCONFIG_AV1_HIGHBITDEPTH=1"

  eval "$cmake_command" "${cmake_extra_args}" ${devnull}
  eval make -j$(nproc) ${devnull}
  echo "Done building target: ${target}"
}

compare_enc_output() {
  local target=$1
  local cpu=$2
  local clip=$3
  local bitrate=$4
  local preset=$5
  diff ${AOM_TEST_OUTPUT_DIR}/Out-generic-"${clip}"-${preset}-${bitrate}kbps-cpu${cpu}${OUT_FILE_SUFFIX} \
       ${AOM_TEST_OUTPUT_DIR}/Out-${target}-"${clip}"-${preset}-${bitrate}kbps-cpu${cpu}${OUT_FILE_SUFFIX} > /dev/null
  if [ $? -eq 1 ]; then
    elog "C vs ${target} encode mismatches for ${clip}, at ${bitrate} kbps, speed ${cpu}, ${preset} preset"
    return 1
  fi
}

av1_enc_test() {
  local encoder="$1"
  local target="$2"
  if [ -z "$(av1_enc_tool_path "${target}")" ]; then
    elog "aomenc not found. It must exist in ${AOM_TEST_OUTPUT_DIR}/build_target_${target} path"
    return 1
  fi

  for preset in ${PRESETS}; do
    if [ "${preset}" = "good" ]; then
      local min_cpu_used=0
      local max_cpu_used=6
      local test_params=av1_encode_good_params
    elif [ "${preset}" = "rt" ]; then
      local min_cpu_used=5
      local max_cpu_used=10
      local test_params=av1_encode_rt_params
    else
      elog "Invalid preset"
      return 1
    fi

    for cpu in $(seq $min_cpu_used $max_cpu_used); do
      for clip in ${TEST_CLIPS}; do
        for bitrate in ${TEST_BITRATES}; do
          eval "${encoder}" $($clip) $($test_params) \
          "--limit=${AV1_ENCODE_C_VS_SIMD_TEST_FRAME_LIMIT}" \
          "--cpu-used=${cpu}" "--target-bitrate=${bitrate}" "-o" \
          ${AOM_TEST_OUTPUT_DIR}/Out-${target}-"${clip}"-${preset}-${bitrate}kbps-cpu${cpu}${OUT_FILE_SUFFIX} \
          ${devnull}

          if [ "${target}" != "generic" ]; then
            compare_enc_output ${target} $cpu ${clip} $bitrate ${preset}
            if [ $? -eq 1 ]; then
              return 1
            fi
          fi
        done
      done
    done
  done
}

av1_test_generic() {
  local target="generic"
  local cmake_command="cmake $LIBAOM_SOURCE_DIR -DAOM_TARGET_CPU=${target}"

  echo "Build for: Generic"
  av1_enc_build "${target}" "${cmake_command}"
  local encoder="$(av1_enc_tool_path "${target}")"
  av1_enc_test $encoder "${target}"
}

# This function encodes AV1 bitstream by enabling SSE2, SSE3, SSSE3, SSE4_1, AVX, AVX2 as there are
# no functions with MMX, SSE and AVX512 specialization.
# The value of environment variable 'AOM_SIMD_CAPS_MASK' controls enabling of different instruction
# set extension optimizations. The value of the flag 'AOM_SIMD_CAPS_MASK' and the corresponding
# instruction set extension optimization enabled are as follows:
# AVX512 AVX2 AVX SSE4_1 SSSE3 SSE3 SSE2 SSE MMX
#   1     1    1    1      1    1    1    1   1  -> 0x1FF -> Enable AVX512 and lower variants
#   0     1    1    1      1    1    1    1   1  -> 0x0FF -> Enable AVX2 and lower variants
#   0     0    1    1      1    1    1    1   1  -> 0x07F -> Enable AVX and lower variants
#   0     0    0    1      1    1    1    1   1  -> 0x03F  -> Enable SSE4_1 and lower variants
#   0     0    0    0      1    1    1    1   1  -> 0x01F  -> Enable SSSE3 and lower variants
#   0     0    0    0      0    1    1    1   1  -> 0x00F  -> Enable SSE3 and lower variants
#   0     0    0    0      0    0    1    1   1  -> 0x007  -> Enable SSE2 and lower variants
#   0     0    0    0      0    0    0    1   1  -> 0x003  -> Enable SSE and lower variants
#   0     0    0    0      0    0    0    0   1  -> 0x001  -> Enable MMX
## NOTE: In x86_64 platform, it is not possible to enable sse/mmx/c using "AOM_SIMD_CAPS_MASK" as
#  all x86_64 platforms implement sse2.
av1_test_x86() {
  local arch=$1

  uname -m | grep -q "x86"
  if [ $? -eq 1 ]; then
    elog "Machine architecture is not x86 or x86_64"
    return 0
  fi

  if [ $arch = "x86" ]; then
    local target="x86-linux"
    local cmake_command="cmake \
    $LIBAOM_SOURCE_DIR \
    -DCMAKE_TOOLCHAIN_FILE=${LIBAOM_SOURCE_DIR}/build/cmake/toolchains/${target}.cmake"
  elif [ $arch = "x86_64" ]; then
    local target="x86_64-linux"
    local cmake_command="cmake $LIBAOM_SOURCE_DIR"
  fi

  local x86_isa_variants="avx2 avx sse4_1 ssse3 sse3 sse2"

  echo "Build for x86: ${target}"
  av1_enc_build "${target}" "${cmake_command}"
  local encoder="$(av1_enc_tool_path "${target}")"
  for isa in $x86_isa_variants; do
    has_x86_isa_extn $isa
    if [ $? -eq 1 ]; then
      echo "${isa} is not supported in this machine"
      continue
    fi
    export AOM_SIMD_CAPS_MASK=$($isa)
    av1_enc_test $encoder "${target}"
    if [ $? -eq 1 ]; then
      return 1
    fi
    unset AOM_SIMD_CAPS_MASK
  done
}

av1_test_arm() {
  local target="arm64-linux-gcc"
  local cmake_command="cmake $LIBAOM_SOURCE_DIR \
        -DCMAKE_TOOLCHAIN_FILE=$LIBAOM_SOURCE_DIR/build/cmake/toolchains/${target}.cmake \
        -DCMAKE_C_FLAGS=-Wno-maybe-uninitialized"
  echo "Build for arm64: ${target}"
  av1_enc_build "${target}" "${cmake_command}"

  local encoder="$(av1_enc_tool_path "${target}")"
  av1_enc_test "qemu-aarch64 -L /usr/aarch64-linux-gnu ${encoder}" "${target}"
  if [ $? -eq 1 ]; then
    return 1
  fi
}

av1_c_vs_simd_enc_test () {
  # Test Generic
  av1_test_generic

  # Test x86 (32 bit)
  echo "av1 test for x86 (32 bit): Started."
  av1_test_x86 "x86"
  if [ $? -eq 1 ]; then
    echo "av1 test for x86 (32 bit): Done, test failed."
  else
    echo "av1 test for x86 (32 bit): Done, all tests passed."
  fi

  # Test x86_64 (64 bit)
  if [ "$(eval uname -m)" = "x86_64" ]; then
    echo "av1 test for x86_64 (64 bit): Started."
    av1_test_x86 "x86_64"
    if [ $? -eq 1 ]; then
      echo "av1 test for x86_64 (64 bit): Done, test failed."
    else
      echo "av1 test for x86_64 (64 bit): Done, all tests passed."
    fi
  fi

  # Test ARM
  echo "av1_test_arm: Started."
  av1_test_arm
  if [ $? -eq 1 ]; then
    echo "av1 test for arm: Done, test failed."
  else
    echo "av1 test for arm: Done, all tests passed."
  fi
}

# Setup a trap function to clean up build, and output files after tests complete.
trap cleanup EXIT

av1_c_vs_simd_enc_verify_environment
if [ $? -eq 1 ]; then
  echo "Environment check failed."
  exit 1
fi
av1_c_vs_simd_enc_test
