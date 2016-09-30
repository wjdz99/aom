/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include "aom_dsp/aom_simd.h"
#include "aom_dsp/simd/v128_intrinsics_c.h"
#include "./simd_test_defs.h"

namespace simd_test {

// Wrap templates around intrinsics using immediate values
template <int shift>
v64 imm_v64_shl_n_byte(v64 a) {
  return v64_shl_n_byte(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_byte(v64 a) {
  return v64_shr_n_byte(a, shift);
}
template <int shift>
v64 imm_v64_shl_n_8(v64 a) {
  return v64_shl_n_8(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_u8(v64 a) {
  return v64_shr_n_u8(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_s8(v64 a) {
  return v64_shr_n_s8(a, shift);
}
template <int shift>
v64 imm_v64_shl_n_16(v64 a) {
  return v64_shl_n_16(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_u16(v64 a) {
  return v64_shr_n_u16(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_s16(v64 a) {
  return v64_shr_n_s16(a, shift);
}
template <int shift>
v64 imm_v64_shl_n_32(v64 a) {
  return v64_shl_n_32(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_u32(v64 a) {
  return v64_shr_n_u32(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_s32(v64 a) {
  return v64_shr_n_s32(a, shift);
}
template <int shift>
v64 imm_v64_align(v64 a, v64 b) {
  return v64_align(a, b, shift);
}

template <int shift>
v128 imm_v128_shl_n_byte(v128 a) {
  return v128_shl_n_byte(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_byte(v128 a) {
  return v128_shr_n_byte(a, shift);
}
template <int shift>
v128 imm_v128_shl_n_8(v128 a) {
  return v128_shl_n_8(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_u8(v128 a) {
  return v128_shr_n_u8(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_s8(v128 a) {
  return v128_shr_n_s8(a, shift);
}
template <int shift>
v128 imm_v128_shl_n_16(v128 a) {
  return v128_shl_n_16(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_u16(v128 a) {
  return v128_shr_n_u16(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_s16(v128 a) {
  return v128_shr_n_s16(a, shift);
}
template <int shift>
v128 imm_v128_shl_n_32(v128 a) {
  return v128_shl_n_32(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_u32(v128 a) {
  return v128_shr_n_u32(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_s32(v128 a) {
  return v128_shr_n_s32(a, shift);
}
template <int shift>
v128 imm_v128_align(v128 a, v128 b) {
  return v128_align(a, b, shift);
}

template <int shift>
c_v64 c_imm_v64_shl_n_byte(c_v64 a) {
  return c_v64_shl_n_byte(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_byte(c_v64 a) {
  return c_v64_shr_n_byte(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shl_n_8(c_v64 a) {
  return c_v64_shl_n_8(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_u8(c_v64 a) {
  return c_v64_shr_n_u8(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_s8(c_v64 a) {
  return c_v64_shr_n_s8(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shl_n_16(c_v64 a) {
  return c_v64_shl_n_16(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_u16(c_v64 a) {
  return c_v64_shr_n_u16(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_s16(c_v64 a) {
  return c_v64_shr_n_s16(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shl_n_32(c_v64 a) {
  return c_v64_shl_n_32(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_u32(c_v64 a) {
  return c_v64_shr_n_u32(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_s32(c_v64 a) {
  return c_v64_shr_n_s32(a, shift);
}
template <int shift>
c_v64 c_imm_v64_align(c_v64 a, c_v64 b) {
  return c_v64_align(a, b, shift);
}

template <int shift>
c_v128 c_imm_v128_shl_n_byte(c_v128 a) {
  return c_v128_shl_n_byte(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_byte(c_v128 a) {
  return c_v128_shr_n_byte(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shl_n_8(c_v128 a) {
  return c_v128_shl_n_8(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_u8(c_v128 a) {
  return c_v128_shr_n_u8(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_s8(c_v128 a) {
  return c_v128_shr_n_s8(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shl_n_16(c_v128 a) {
  return c_v128_shl_n_16(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_u16(c_v128 a) {
  return c_v128_shr_n_u16(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_s16(c_v128 a) {
  return c_v128_shr_n_s16(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shl_n_32(c_v128 a) {
  return c_v128_shl_n_32(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_u32(c_v128 a) {
  return c_v128_shr_n_u32(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_s32(c_v128 a) {
  return c_v128_shr_n_s32(a, shift);
}
template <int shift>
c_v128 c_imm_v128_align(c_v128 a, c_v128 b) {
  return c_v128_align(a, b, shift);
}

// Wrappers around the the SAD and SSD functions
static uint32_t v64_sad_u8(v64 a, v64 b) {
  return v64_sad_u8_sum(::v64_sad_u8(v64_sad_u8_init(), a, b));
}
static uint32_t v64_ssd_u8(v64 a, v64 b) {
  return v64_ssd_u8_sum(::v64_ssd_u8(v64_ssd_u8_init(), a, b));
}
static uint32_t v128_sad_u8(v128 a, v128 b) {
  return v128_sad_u8_sum(::v128_sad_u8(v128_sad_u8_init(), a, b));
}
static uint32_t v128_ssd_u8(v128 a, v128 b) {
  return v128_ssd_u8_sum(::v128_ssd_u8(v128_ssd_u8_init(), a, b));
}

static uint32_t c_v64_sad_u8(c_v64 a, c_v64 b) {
  return c_v64_sad_u8_sum(::c_v64_sad_u8(c_v64_sad_u8_init(), a, b));
}
static uint32_t c_v64_ssd_u8(c_v64 a, c_v64 b) {
  return c_v64_ssd_u8_sum(::c_v64_ssd_u8(c_v64_ssd_u8_init(), a, b));
}
static uint32_t c_v128_sad_u8(c_v128 a, c_v128 b) {
  return c_v128_sad_u8_sum(::c_v128_sad_u8(c_v128_sad_u8_init(), a, b));
}
static uint32_t c_v128_ssd_u8(c_v128 a, c_v128 b) {
  return c_v128_ssd_u8_sum(::c_v128_ssd_u8(c_v128_ssd_u8_init(), a, b));
}


// Add a macro layer since INSTANTIATE_TEST_CASE_P will quote the name
// so we need to expand it first with the prefix
#define INSTANTIATE(name, type, ...) \
  INSTANTIATE_TEST_CASE_P(name, type, ::testing::Values(__VA_ARGS__))

#define SIMD_TUPLE(name, mask, maskwidth) \
  std::tr1::make_tuple(name, c_##name, mask, maskwidth, #name)

#define IMM_SIMD_TUPLE(name, imm) \
  std::tr1::make_tuple(imm_##name<imm>, c_imm_##name<imm>, 0U, 0U, #name)

INSTANTIATE(ARCH_PREFIX, U32_V64V64,
            (SIMD_TUPLE(v64_sad_u8, 0U, 0U), SIMD_TUPLE(v64_ssd_u8, 0U, 0U)));

INSTANTIATE(ARCH_PREFIX, U32_V128V128, SIMD_TUPLE(v128_sad_u8, 0U, 0U),
            SIMD_TUPLE(v128_ssd_u8, 0U, 0U));

INSTANTIATE(
    ARCH_PREFIX, V128_V128V128, SIMD_TUPLE(v128_add_8, 0U, 0U),
    SIMD_TUPLE(v128_add_16, 0U, 0U), SIMD_TUPLE(v128_sadd_s16, 0U, 0U),
    SIMD_TUPLE(v128_add_32, 0U, 0U), SIMD_TUPLE(v128_sub_8, 0U, 0U),
    SIMD_TUPLE(v128_ssub_u8, 0U, 0U), SIMD_TUPLE(v128_ssub_s8, 0U, 0U),
    SIMD_TUPLE(v128_sub_16, 0U, 0U), SIMD_TUPLE(v128_ssub_s16, 0U, 0U),
    SIMD_TUPLE(v128_sub_32, 0U, 0U), SIMD_TUPLE(v128_ziplo_8, 0U, 0U),
    SIMD_TUPLE(v128_ziphi_8, 0U, 0U), SIMD_TUPLE(v128_ziplo_16, 0U, 0U),
    SIMD_TUPLE(v128_ziphi_16, 0U, 0U), SIMD_TUPLE(v128_ziplo_32, 0U, 0U),
    SIMD_TUPLE(v128_ziphi_32, 0U, 0U), SIMD_TUPLE(v128_ziplo_64, 0U, 0U),
    SIMD_TUPLE(v128_ziphi_64, 0U, 0U), SIMD_TUPLE(v128_unziphi_8, 0U, 0U),
    SIMD_TUPLE(v128_unziplo_8, 0U, 0U), SIMD_TUPLE(v128_unziphi_16, 0U, 0U),
    SIMD_TUPLE(v128_unziplo_16, 0U, 0U), SIMD_TUPLE(v128_unziphi_32, 0U, 0U),
    SIMD_TUPLE(v128_unziplo_32, 0U, 0U), SIMD_TUPLE(v128_pack_s32_s16, 0U, 0U),
    SIMD_TUPLE(v128_pack_s16_u8, 0U, 0U), SIMD_TUPLE(v128_pack_s16_s8, 0U, 0U),
    SIMD_TUPLE(v128_or, 0U, 0U), SIMD_TUPLE(v128_xor, 0U, 0U),
    SIMD_TUPLE(v128_and, 0U, 0U), SIMD_TUPLE(v128_andn, 0U, 0U),
    SIMD_TUPLE(v128_mullo_s16, 0U, 0U), SIMD_TUPLE(v128_mulhi_s16, 0U, 0U),
    SIMD_TUPLE(v128_mullo_s32, 0U, 0U), SIMD_TUPLE(v128_madd_s16, 0U, 0U),
    SIMD_TUPLE(v128_madd_us8, 0U, 0U), SIMD_TUPLE(v128_avg_u8, 0U, 0U),
    SIMD_TUPLE(v128_rdavg_u8, 0U, 0U), SIMD_TUPLE(v128_avg_u16, 0U, 0U),
    SIMD_TUPLE(v128_min_u8, 0U, 0U), SIMD_TUPLE(v128_max_u8, 0U, 0U),
    SIMD_TUPLE(v128_min_s8, 0U, 0U), SIMD_TUPLE(v128_max_s8, 0U, 0U),
    SIMD_TUPLE(v128_min_s16, 0U, 0U), SIMD_TUPLE(v128_max_s16, 0U, 0U),
    SIMD_TUPLE(v128_cmpgt_s8, 0U, 0U), SIMD_TUPLE(v128_cmplt_s8, 0U, 0U),
    SIMD_TUPLE(v128_cmpeq_8, 0U, 0U), SIMD_TUPLE(v128_cmpgt_s16, 0U, 0U),
    SIMD_TUPLE(v128_cmpeq_16, 0U, 0U));

INSTANTIATE(ARCH_PREFIX, V128_V128V128_Part2,
            SIMD_TUPLE(v128_cmplt_s16, 0U, 0U),
            SIMD_TUPLE(v128_shuffle_8, 15U, 8U), IMM_SIMD_TUPLE(v128_align, 1),
            IMM_SIMD_TUPLE(v128_align, 2), IMM_SIMD_TUPLE(v128_align, 3),
            IMM_SIMD_TUPLE(v128_align, 4), IMM_SIMD_TUPLE(v128_align, 5),
            IMM_SIMD_TUPLE(v128_align, 6), IMM_SIMD_TUPLE(v128_align, 7),
            IMM_SIMD_TUPLE(v128_align, 8), IMM_SIMD_TUPLE(v128_align, 9),
            IMM_SIMD_TUPLE(v128_align, 10), IMM_SIMD_TUPLE(v128_align, 11),
            IMM_SIMD_TUPLE(v128_align, 12), IMM_SIMD_TUPLE(v128_align, 13),
            IMM_SIMD_TUPLE(v128_align, 14), IMM_SIMD_TUPLE(v128_align, 15));

INSTANTIATE(
    ARCH_PREFIX, V64_V64V64, SIMD_TUPLE(v64_add_8, 0U, 0U),
    SIMD_TUPLE(v64_add_16, 0U, 0U), SIMD_TUPLE(v64_sadd_s16, 0U, 0U),
    SIMD_TUPLE(v64_add_32, 0U, 0U), SIMD_TUPLE(v64_sub_8, 0U, 0U),
    SIMD_TUPLE(v64_ssub_u8, 0U, 0U), SIMD_TUPLE(v64_ssub_s8, 0U, 0U),
    SIMD_TUPLE(v64_sub_16, 0U, 0U), SIMD_TUPLE(v64_ssub_s16, 0U, 0U),
    SIMD_TUPLE(v64_sub_32, 0U, 0U), SIMD_TUPLE(v64_ziplo_8, 0U, 0U),
    SIMD_TUPLE(v64_ziphi_8, 0U, 0U), SIMD_TUPLE(v64_ziplo_16, 0U, 0U),
    SIMD_TUPLE(v64_ziphi_16, 0U, 0U), SIMD_TUPLE(v64_ziplo_32, 0U, 0U),
    SIMD_TUPLE(v64_ziphi_32, 0U, 0U), SIMD_TUPLE(v64_pack_s32_s16, 0U, 0U),
    SIMD_TUPLE(v64_pack_s16_u8, 0U, 0U), SIMD_TUPLE(v64_pack_s16_s8, 0U, 0U),
    SIMD_TUPLE(v64_unziphi_8, 0U, 0U), SIMD_TUPLE(v64_unziplo_8, 0U, 0U),
    SIMD_TUPLE(v64_unziphi_16, 0U, 0U), SIMD_TUPLE(v64_unziplo_16, 0U, 0U),
    SIMD_TUPLE(v64_or, 0U, 0U), SIMD_TUPLE(v64_xor, 0U, 0U),
    SIMD_TUPLE(v64_and, 0U, 0U), SIMD_TUPLE(v64_andn, 0U, 0U),
    SIMD_TUPLE(v64_mullo_s16, 0U, 0U), SIMD_TUPLE(v64_mulhi_s16, 0U, 0U),
    SIMD_TUPLE(v64_mullo_s32, 0U, 0U), SIMD_TUPLE(v64_madd_s16, 0U, 0U),
    SIMD_TUPLE(v64_madd_us8, 0U, 0U), SIMD_TUPLE(v64_avg_u8, 0U, 0U),
    SIMD_TUPLE(v64_rdavg_u8, 0U, 0U), SIMD_TUPLE(v64_avg_u16, 0U, 0U),
    SIMD_TUPLE(v64_min_u8, 0U, 0U), SIMD_TUPLE(v64_max_u8, 0U, 0U),
    SIMD_TUPLE(v64_min_s8, 0U, 0U), SIMD_TUPLE(v64_max_s8, 0U, 0U),
    SIMD_TUPLE(v64_min_s16, 0U, 0U), SIMD_TUPLE(v64_max_s16, 0U, 0U),
    SIMD_TUPLE(v64_cmpgt_s8, 0U, 0U), SIMD_TUPLE(v64_cmplt_s8, 0U, 0U),
    SIMD_TUPLE(v64_cmpeq_8, 0U, 0U), SIMD_TUPLE(v64_cmpgt_s16, 0U, 0U),
    SIMD_TUPLE(v64_cmplt_s16, 0U, 0U), SIMD_TUPLE(v64_cmpeq_16, 0U, 0U),
    SIMD_TUPLE(v64_shuffle_8, 7U, 8U));

INSTANTIATE(ARCH_PREFIX, V64_V64V64_Part2, IMM_SIMD_TUPLE(v64_align, 1),
            IMM_SIMD_TUPLE(v64_align, 2), IMM_SIMD_TUPLE(v64_align, 3),
            IMM_SIMD_TUPLE(v64_align, 4), IMM_SIMD_TUPLE(v64_align, 5),
            IMM_SIMD_TUPLE(v64_align, 6), IMM_SIMD_TUPLE(v64_align, 7));

INSTANTIATE(
    ARCH_PREFIX, V128_V128, SIMD_TUPLE(v128_abs_s16, 0U, 0U),
    SIMD_TUPLE(v128_padd_s16, 0U, 0U),
    SIMD_TUPLE(v128_unpacklo_u16_s32, 0U, 0U),
    SIMD_TUPLE(v128_unpacklo_s16_s32, 0U, 0U),
    SIMD_TUPLE(v128_unpackhi_u16_s32, 0U, 0U),
    SIMD_TUPLE(v128_unpackhi_s16_s32, 0U, 0U),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 1), IMM_SIMD_TUPLE(v128_shr_n_byte, 2),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 3), IMM_SIMD_TUPLE(v128_shr_n_byte, 4),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 5), IMM_SIMD_TUPLE(v128_shr_n_byte, 6),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 7), IMM_SIMD_TUPLE(v128_shr_n_byte, 8),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 9), IMM_SIMD_TUPLE(v128_shr_n_byte, 10),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 11), IMM_SIMD_TUPLE(v128_shr_n_byte, 12),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 13), IMM_SIMD_TUPLE(v128_shr_n_byte, 14),
    IMM_SIMD_TUPLE(v128_shr_n_byte, 15), IMM_SIMD_TUPLE(v128_shl_n_byte, 1),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 2), IMM_SIMD_TUPLE(v128_shl_n_byte, 3),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 4), IMM_SIMD_TUPLE(v128_shl_n_byte, 5),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 6), IMM_SIMD_TUPLE(v128_shl_n_byte, 7),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 8), IMM_SIMD_TUPLE(v128_shl_n_byte, 9),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 10), IMM_SIMD_TUPLE(v128_shl_n_byte, 11),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 12), IMM_SIMD_TUPLE(v128_shl_n_byte, 13),
    IMM_SIMD_TUPLE(v128_shl_n_byte, 14), IMM_SIMD_TUPLE(v128_shl_n_byte, 15),
    IMM_SIMD_TUPLE(v128_shl_n_8, 1), IMM_SIMD_TUPLE(v128_shl_n_8, 2),
    IMM_SIMD_TUPLE(v128_shl_n_8, 3), IMM_SIMD_TUPLE(v128_shl_n_8, 4),
    IMM_SIMD_TUPLE(v128_shl_n_8, 5), IMM_SIMD_TUPLE(v128_shl_n_8, 6),
    IMM_SIMD_TUPLE(v128_shl_n_8, 7), IMM_SIMD_TUPLE(v128_shr_n_u8, 1),
    IMM_SIMD_TUPLE(v128_shr_n_u8, 2), IMM_SIMD_TUPLE(v128_shr_n_u8, 3),
    IMM_SIMD_TUPLE(v128_shr_n_u8, 4), IMM_SIMD_TUPLE(v128_shr_n_u8, 5),
    IMM_SIMD_TUPLE(v128_shr_n_u8, 6));

INSTANTIATE(
    ARCH_PREFIX, V128_V128_Part2, IMM_SIMD_TUPLE(v128_shr_n_u8, 7),
    IMM_SIMD_TUPLE(v128_shr_n_s8, 1), IMM_SIMD_TUPLE(v128_shr_n_s8, 2),
    IMM_SIMD_TUPLE(v128_shr_n_s8, 3), IMM_SIMD_TUPLE(v128_shr_n_s8, 4),
    IMM_SIMD_TUPLE(v128_shr_n_s8, 5), IMM_SIMD_TUPLE(v128_shr_n_s8, 6),
    IMM_SIMD_TUPLE(v128_shr_n_s8, 7), IMM_SIMD_TUPLE(v128_shl_n_16, 1),
    IMM_SIMD_TUPLE(v128_shl_n_16, 2), IMM_SIMD_TUPLE(v128_shl_n_16, 4),
    IMM_SIMD_TUPLE(v128_shl_n_16, 6), IMM_SIMD_TUPLE(v128_shl_n_16, 8),
    IMM_SIMD_TUPLE(v128_shl_n_16, 10), IMM_SIMD_TUPLE(v128_shl_n_16, 12),
    IMM_SIMD_TUPLE(v128_shl_n_16, 14), IMM_SIMD_TUPLE(v128_shr_n_u16, 1),
    IMM_SIMD_TUPLE(v128_shr_n_u16, 2), IMM_SIMD_TUPLE(v128_shr_n_u16, 4),
    IMM_SIMD_TUPLE(v128_shr_n_u16, 6), IMM_SIMD_TUPLE(v128_shr_n_u16, 8),
    IMM_SIMD_TUPLE(v128_shr_n_u16, 10), IMM_SIMD_TUPLE(v128_shr_n_u16, 12),
    IMM_SIMD_TUPLE(v128_shr_n_u16, 14), IMM_SIMD_TUPLE(v128_shr_n_s16, 1),
    IMM_SIMD_TUPLE(v128_shr_n_s16, 2), IMM_SIMD_TUPLE(v128_shr_n_s16, 4),
    IMM_SIMD_TUPLE(v128_shr_n_s16, 6), IMM_SIMD_TUPLE(v128_shr_n_s16, 8),
    IMM_SIMD_TUPLE(v128_shr_n_s16, 10), IMM_SIMD_TUPLE(v128_shr_n_s16, 12),
    IMM_SIMD_TUPLE(v128_shr_n_s16, 14), IMM_SIMD_TUPLE(v128_shl_n_32, 1),
    IMM_SIMD_TUPLE(v128_shl_n_32, 4), IMM_SIMD_TUPLE(v128_shl_n_32, 8),
    IMM_SIMD_TUPLE(v128_shl_n_32, 12), IMM_SIMD_TUPLE(v128_shl_n_32, 16),
    IMM_SIMD_TUPLE(v128_shl_n_32, 20), IMM_SIMD_TUPLE(v128_shl_n_32, 24),
    IMM_SIMD_TUPLE(v128_shl_n_32, 28), IMM_SIMD_TUPLE(v128_shr_n_u32, 1),
    IMM_SIMD_TUPLE(v128_shr_n_u32, 4), IMM_SIMD_TUPLE(v128_shr_n_u32, 8),
    IMM_SIMD_TUPLE(v128_shr_n_u32, 12), IMM_SIMD_TUPLE(v128_shr_n_u32, 16),
    IMM_SIMD_TUPLE(v128_shr_n_u32, 20), IMM_SIMD_TUPLE(v128_shr_n_u32, 24));

INSTANTIATE(
    ARCH_PREFIX, V128_V128_Part3, IMM_SIMD_TUPLE(v128_shr_n_u32, 28),
    IMM_SIMD_TUPLE(v128_shr_n_s32, 1), IMM_SIMD_TUPLE(v128_shr_n_s32, 4),
    IMM_SIMD_TUPLE(v128_shr_n_s32, 8), IMM_SIMD_TUPLE(v128_shr_n_s32, 12),
    IMM_SIMD_TUPLE(v128_shr_n_s32, 16), IMM_SIMD_TUPLE(v128_shr_n_s32, 20),
    IMM_SIMD_TUPLE(v128_shr_n_s32, 24), IMM_SIMD_TUPLE(v128_shr_n_s32, 28));

INSTANTIATE(
    ARCH_PREFIX, V64_V64, SIMD_TUPLE(v64_abs_s16, 0U, 0U),
    SIMD_TUPLE(v64_unpacklo_u8_s16, 0U, 0U),
    SIMD_TUPLE(v64_unpackhi_u8_s16, 0U, 0U),
    SIMD_TUPLE(v64_unpacklo_u16_s32, 0U, 0U),
    SIMD_TUPLE(v64_unpacklo_s16_s32, 0U, 0U),
    SIMD_TUPLE(v64_unpackhi_u16_s32, 0U, 0U),
    SIMD_TUPLE(v64_unpackhi_s16_s32, 0U, 0U), IMM_SIMD_TUPLE(v64_shr_n_byte, 1),
    IMM_SIMD_TUPLE(v64_shr_n_byte, 2), IMM_SIMD_TUPLE(v64_shr_n_byte, 3),
    IMM_SIMD_TUPLE(v64_shr_n_byte, 4), IMM_SIMD_TUPLE(v64_shr_n_byte, 5),
    IMM_SIMD_TUPLE(v64_shr_n_byte, 6), IMM_SIMD_TUPLE(v64_shr_n_byte, 7),
    IMM_SIMD_TUPLE(v64_shl_n_byte, 1), IMM_SIMD_TUPLE(v64_shl_n_byte, 2),
    IMM_SIMD_TUPLE(v64_shl_n_byte, 3), IMM_SIMD_TUPLE(v64_shl_n_byte, 4),
    IMM_SIMD_TUPLE(v64_shl_n_byte, 5), IMM_SIMD_TUPLE(v64_shl_n_byte, 6),
    IMM_SIMD_TUPLE(v64_shl_n_byte, 7), IMM_SIMD_TUPLE(v64_shl_n_8, 1),
    IMM_SIMD_TUPLE(v64_shl_n_8, 2), IMM_SIMD_TUPLE(v64_shl_n_8, 3),
    IMM_SIMD_TUPLE(v64_shl_n_8, 4), IMM_SIMD_TUPLE(v64_shl_n_8, 5),
    IMM_SIMD_TUPLE(v64_shl_n_8, 6), IMM_SIMD_TUPLE(v64_shl_n_8, 7),
    IMM_SIMD_TUPLE(v64_shr_n_u8, 1), IMM_SIMD_TUPLE(v64_shr_n_u8, 2),
    IMM_SIMD_TUPLE(v64_shr_n_u8, 3), IMM_SIMD_TUPLE(v64_shr_n_u8, 4),
    IMM_SIMD_TUPLE(v64_shr_n_u8, 5), IMM_SIMD_TUPLE(v64_shr_n_u8, 6),
    IMM_SIMD_TUPLE(v64_shr_n_u8, 7), IMM_SIMD_TUPLE(v64_shr_n_s8, 1),
    IMM_SIMD_TUPLE(v64_shr_n_s8, 2), IMM_SIMD_TUPLE(v64_shr_n_s8, 3),
    IMM_SIMD_TUPLE(v64_shr_n_s8, 4), IMM_SIMD_TUPLE(v64_shr_n_s8, 5),
    IMM_SIMD_TUPLE(v64_shr_n_s8, 6), IMM_SIMD_TUPLE(v64_shr_n_s8, 7),
    IMM_SIMD_TUPLE(v64_shl_n_16, 1), IMM_SIMD_TUPLE(v64_shl_n_16, 2),
    IMM_SIMD_TUPLE(v64_shl_n_16, 4), IMM_SIMD_TUPLE(v64_shl_n_16, 6),
    IMM_SIMD_TUPLE(v64_shl_n_16, 8), IMM_SIMD_TUPLE(v64_shl_n_16, 10),
    IMM_SIMD_TUPLE(v64_shl_n_16, 12), IMM_SIMD_TUPLE(v64_shl_n_16, 14));

INSTANTIATE(
    ARCH_PREFIX, V64_V64_Part2, IMM_SIMD_TUPLE(v64_shr_n_u16, 1),
    IMM_SIMD_TUPLE(v64_shr_n_u16, 2), IMM_SIMD_TUPLE(v64_shr_n_u16, 4),
    IMM_SIMD_TUPLE(v64_shr_n_u16, 6), IMM_SIMD_TUPLE(v64_shr_n_u16, 8),
    IMM_SIMD_TUPLE(v64_shr_n_u16, 10), IMM_SIMD_TUPLE(v64_shr_n_u16, 12),
    IMM_SIMD_TUPLE(v64_shr_n_u16, 14), IMM_SIMD_TUPLE(v64_shr_n_s16, 1),
    IMM_SIMD_TUPLE(v64_shr_n_s16, 2), IMM_SIMD_TUPLE(v64_shr_n_s16, 4),
    IMM_SIMD_TUPLE(v64_shr_n_s16, 6), IMM_SIMD_TUPLE(v64_shr_n_s16, 8),
    IMM_SIMD_TUPLE(v64_shr_n_s16, 10), IMM_SIMD_TUPLE(v64_shr_n_s16, 12),
    IMM_SIMD_TUPLE(v64_shr_n_s16, 14), IMM_SIMD_TUPLE(v64_shl_n_32, 1),
    IMM_SIMD_TUPLE(v64_shl_n_32, 4), IMM_SIMD_TUPLE(v64_shl_n_32, 8),
    IMM_SIMD_TUPLE(v64_shl_n_32, 12), IMM_SIMD_TUPLE(v64_shl_n_32, 16),
    IMM_SIMD_TUPLE(v64_shl_n_32, 20), IMM_SIMD_TUPLE(v64_shl_n_32, 24),
    IMM_SIMD_TUPLE(v64_shl_n_32, 28), IMM_SIMD_TUPLE(v64_shr_n_u32, 1),
    IMM_SIMD_TUPLE(v64_shr_n_u32, 4), IMM_SIMD_TUPLE(v64_shr_n_u32, 8),
    IMM_SIMD_TUPLE(v64_shr_n_u32, 12), IMM_SIMD_TUPLE(v64_shr_n_u32, 16),
    IMM_SIMD_TUPLE(v64_shr_n_u32, 20), IMM_SIMD_TUPLE(v64_shr_n_u32, 24),
    IMM_SIMD_TUPLE(v64_shr_n_u32, 28), IMM_SIMD_TUPLE(v64_shr_n_s32, 1),
    IMM_SIMD_TUPLE(v64_shr_n_s32, 4), IMM_SIMD_TUPLE(v64_shr_n_s32, 8),
    IMM_SIMD_TUPLE(v64_shr_n_s32, 12), IMM_SIMD_TUPLE(v64_shr_n_s32, 16),
    IMM_SIMD_TUPLE(v64_shr_n_s32, 20), IMM_SIMD_TUPLE(v64_shr_n_s32, 24),
    IMM_SIMD_TUPLE(v64_shr_n_s32, 28));

INSTANTIATE(ARCH_PREFIX, V128_V64V64, SIMD_TUPLE(v128_from_v64, 0U, 0U),
            SIMD_TUPLE(v128_zip_8, 0U, 0U), SIMD_TUPLE(v128_zip_16, 0U, 0U),
            SIMD_TUPLE(v128_zip_32, 0U, 0U), SIMD_TUPLE(v128_mul_s16, 0U, 0U));

INSTANTIATE(ARCH_PREFIX, V128_V64, SIMD_TUPLE(v128_unpack_u8_s16, 0U, 0U),
            SIMD_TUPLE(v128_unpack_u16_s32, 0U, 0U),
            SIMD_TUPLE(v128_unpack_s16_s32, 0U, 0U));

INSTANTIATE(ARCH_PREFIX, V128_V128U32, SIMD_TUPLE(v128_shl_8, 7U, 32U),
            SIMD_TUPLE(v128_shr_u8, 7U, 32U), SIMD_TUPLE(v128_shr_s8, 7U, 32U),
            SIMD_TUPLE(v128_shl_16, 15U, 32U),
            SIMD_TUPLE(v128_shr_u16, 15U, 32U),
            SIMD_TUPLE(v128_shr_s16, 15U, 32U),
            SIMD_TUPLE(v128_shl_32, 31U, 32U),
            SIMD_TUPLE(v128_shr_u32, 31U, 32U),
            SIMD_TUPLE(v128_shr_s32, 31U, 32U));

INSTANTIATE(ARCH_PREFIX, V64_V64U32, SIMD_TUPLE(v64_shl_8, 7U, 32U),
            SIMD_TUPLE(v64_shr_u8, 7U, 32U), SIMD_TUPLE(v64_shr_s8, 7U, 32U),
            SIMD_TUPLE(v64_shl_16, 15U, 32U), SIMD_TUPLE(v64_shr_u16, 15U, 32U),
            SIMD_TUPLE(v64_shr_s16, 15U, 32U), SIMD_TUPLE(v64_shl_32, 31U, 32U),
            SIMD_TUPLE(v64_shr_u32, 31U, 32U),
            SIMD_TUPLE(v64_shr_s32, 31U, 32U));

INSTANTIATE(ARCH_PREFIX, U64_V64, SIMD_TUPLE(v64_hadd_u8, 0U, 0U));

INSTANTIATE(ARCH_PREFIX, S64_V64, SIMD_TUPLE(v64_hadd_s16, 0U, 0U));

INSTANTIATE(ARCH_PREFIX, U64_V128, SIMD_TUPLE(v128_hadd_u8, 0U, 0U));

}  // namespace simd_test
