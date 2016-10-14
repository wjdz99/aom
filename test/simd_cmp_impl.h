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

#include <string>
#include "./aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/register_state_check.h"
#include "aom_dsp/aom_simd.h"
#include "aom_dsp/simd/v256_intrinsics_c.h"

// Machine tuned code goes into this file. This file is included from
// simd_cmp_sse2.cc, simd_cmp_ssse3.cc etc which define the macros
// ARCH (=neon, sse2, ssse3, etc), SIMD_NAMESPACE and ARCH_POSTFIX().

using libaom_test::ACMRandom;

namespace SIMD_NAMESPACE {

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
v256 imm_v256_shl_n_byte(v256 a) {
  return v256_shl_n_byte(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_byte(v256 a) {
  return v256_shr_n_byte(a, shift);
}
template <int shift>
v256 imm_v256_shl_n_8(v256 a) {
  return v256_shl_n_8(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_u8(v256 a) {
  return v256_shr_n_u8(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_s8(v256 a) {
  return v256_shr_n_s8(a, shift);
}
template <int shift>
v256 imm_v256_shl_n_16(v256 a) {
  return v256_shl_n_16(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_u16(v256 a) {
  return v256_shr_n_u16(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_s16(v256 a) {
  return v256_shr_n_s16(a, shift);
}
template <int shift>
v256 imm_v256_shl_n_32(v256 a) {
  return v256_shl_n_32(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_u32(v256 a) {
  return v256_shr_n_u32(a, shift);
}
template <int shift>
v256 imm_v256_shr_n_s32(v256 a) {
  return v256_shr_n_s32(a, shift);
}
template <int shift>
v256 imm_v256_align(v256 a, v256 b) {
  return v256_align(a, b, shift);
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

template <int shift>
c_v256 c_imm_v256_shl_n_byte(c_v256 a) {
  return c_v256_shl_n_byte(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_byte(c_v256 a) {
  return c_v256_shr_n_byte(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shl_n_8(c_v256 a) {
  return c_v256_shl_n_8(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_u8(c_v256 a) {
  return c_v256_shr_n_u8(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_s8(c_v256 a) {
  return c_v256_shr_n_s8(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shl_n_16(c_v256 a) {
  return c_v256_shl_n_16(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_u16(c_v256 a) {
  return c_v256_shr_n_u16(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_s16(c_v256 a) {
  return c_v256_shr_n_s16(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shl_n_32(c_v256 a) {
  return c_v256_shl_n_32(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_u32(c_v256 a) {
  return c_v256_shr_n_u32(a, shift);
}
template <int shift>
c_v256 c_imm_v256_shr_n_s32(c_v256 a) {
  return c_v256_shr_n_s32(a, shift);
}
template <int shift>
c_v256 c_imm_v256_align(c_v256 a, c_v256 b) {
  return c_v256_align(a, b, shift);
}

// Wrappers around the the SAD and SSD functions
uint32_t v64_sad_u8(v64 a, v64 b) {
  return v64_sad_u8_sum(::v64_sad_u8(v64_sad_u8_init(), a, b));
}
uint32_t v64_ssd_u8(v64 a, v64 b) {
  return v64_ssd_u8_sum(::v64_ssd_u8(v64_ssd_u8_init(), a, b));
}
uint32_t v128_sad_u8(v128 a, v128 b) {
  return v128_sad_u8_sum(::v128_sad_u8(v128_sad_u8_init(), a, b));
}
uint32_t v128_ssd_u8(v128 a, v128 b) {
  return v128_ssd_u8_sum(::v128_ssd_u8(v128_ssd_u8_init(), a, b));
}
uint32_t v256_sad_u8(v256 a, v256 b) {
  return v256_sad_u8_sum(::v256_sad_u8(v256_sad_u8_init(), a, b));
}
uint32_t v256_ssd_u8(v256 a, v256 b) {
  return v256_ssd_u8_sum(::v256_ssd_u8(v256_ssd_u8_init(), a, b));
}

uint32_t c_v64_sad_u8(c_v64 a, c_v64 b) {
  return c_v64_sad_u8_sum(::c_v64_sad_u8(c_v64_sad_u8_init(), a, b));
}
uint32_t c_v64_ssd_u8(c_v64 a, c_v64 b) {
  return c_v64_ssd_u8_sum(::c_v64_ssd_u8(c_v64_ssd_u8_init(), a, b));
}
uint32_t c_v128_sad_u8(c_v128 a, c_v128 b) {
  return c_v128_sad_u8_sum(::c_v128_sad_u8(c_v128_sad_u8_init(), a, b));
}
uint32_t c_v128_ssd_u8(c_v128 a, c_v128 b) {
  return c_v128_ssd_u8_sum(::c_v128_ssd_u8(c_v128_ssd_u8_init(), a, b));
}
uint32_t c_v256_sad_u8(c_v256 a, c_v256 b) {
  return c_v256_sad_u8_sum(::c_v256_sad_u8(c_v256_sad_u8_init(), a, b));
}
uint32_t c_v256_ssd_u8(c_v256 a, c_v256 b) {
  return c_v256_ssd_u8_sum(::c_v256_ssd_u8(c_v256_ssd_u8_init(), a, b));
}

#define MAP(name) \
  { (void *)#name, (void *)c_##name, (void *)name }

// Map reference functions to machine tuned functions. Since the
// functions depend on machine tuned types, the non-machine tuned
// instantiations of the test can't refer to these functions directly,
// so we refer to them by name and do the mapping here.
void map(const char *name, void **ref, void **simd) {
  static void *m[][3] = { MAP(v64_sad_u8),
                          MAP(v64_ssd_u8),
                          MAP(v128_sad_u8),
                          MAP(v128_ssd_u8),
                          MAP(v256_sad_u8),
                          MAP(v256_ssd_u8),
                          MAP(v128_add_8),
                          MAP(v128_add_16),
                          MAP(v128_sadd_s16),
                          MAP(v128_add_32),
                          MAP(v128_sub_8),
                          MAP(v128_ssub_u8),
                          MAP(v128_ssub_s8),
                          MAP(v128_sub_16),
                          MAP(v128_ssub_s16),
                          MAP(v128_sub_32),
                          MAP(v128_ziplo_8),
                          MAP(v128_ziphi_8),
                          MAP(v128_ziplo_16),
                          MAP(v128_ziphi_16),
                          MAP(v128_ziplo_32),
                          MAP(v128_ziphi_32),
                          MAP(v128_ziplo_64),
                          MAP(v128_ziphi_64),
                          MAP(v128_unziphi_8),
                          MAP(v128_unziplo_8),
                          MAP(v128_unziphi_16),
                          MAP(v128_unziplo_16),
                          MAP(v128_unziphi_32),
                          MAP(v128_unziplo_32),
                          MAP(v128_pack_s32_s16),
                          MAP(v128_pack_s16_u8),
                          MAP(v128_pack_s16_s8),
                          MAP(v128_or),
                          MAP(v128_xor),
                          MAP(v128_and),
                          MAP(v128_andn),
                          MAP(v128_mullo_s16),
                          MAP(v128_mulhi_s16),
                          MAP(v128_mullo_s32),
                          MAP(v128_madd_s16),
                          MAP(v128_madd_us8),
                          MAP(v128_avg_u8),
                          MAP(v128_rdavg_u8),
                          MAP(v128_avg_u16),
                          MAP(v128_min_u8),
                          MAP(v128_max_u8),
                          MAP(v128_min_s8),
                          MAP(v128_max_s8),
                          MAP(v128_min_s16),
                          MAP(v128_max_s16),
                          MAP(v128_cmpgt_s8),
                          MAP(v128_cmplt_s8),
                          MAP(v128_cmpeq_8),
                          MAP(v128_cmpgt_s16),
                          MAP(v128_cmpeq_16),
                          MAP(v128_cmplt_s16),
                          MAP(v128_shuffle_8),
                          MAP(imm_v128_align<1>),
                          MAP(imm_v128_align<2>),
                          MAP(imm_v128_align<3>),
                          MAP(imm_v128_align<4>),
                          MAP(imm_v128_align<5>),
                          MAP(imm_v128_align<6>),
                          MAP(imm_v128_align<7>),
                          MAP(imm_v128_align<8>),
                          MAP(imm_v128_align<9>),
                          MAP(imm_v128_align<10>),
                          MAP(imm_v128_align<11>),
                          MAP(imm_v128_align<12>),
                          MAP(imm_v128_align<13>),
                          MAP(imm_v128_align<14>),
                          MAP(imm_v128_align<15>),
                          MAP(v64_add_8),
                          MAP(v64_add_16),
                          MAP(v64_sadd_s16),
                          MAP(v64_add_32),
                          MAP(v64_sub_8),
                          MAP(v64_ssub_u8),
                          MAP(v64_ssub_s8),
                          MAP(v64_sub_16),
                          MAP(v64_ssub_s16),
                          MAP(v64_sub_32),
                          MAP(v64_ziplo_8),
                          MAP(v64_ziphi_8),
                          MAP(v64_ziplo_16),
                          MAP(v64_ziphi_16),
                          MAP(v64_ziplo_32),
                          MAP(v64_ziphi_32),
                          MAP(v64_pack_s32_s16),
                          MAP(v64_pack_s16_u8),
                          MAP(v64_pack_s16_s8),
                          MAP(v64_unziphi_8),
                          MAP(v64_unziplo_8),
                          MAP(v64_unziphi_16),
                          MAP(v64_unziplo_16),
                          MAP(v64_or),
                          MAP(v64_xor),
                          MAP(v64_and),
                          MAP(v64_andn),
                          MAP(v64_mullo_s16),
                          MAP(v64_mulhi_s16),
                          MAP(v64_mullo_s32),
                          MAP(v64_madd_s16),
                          MAP(v64_madd_us8),
                          MAP(v64_avg_u8),
                          MAP(v64_rdavg_u8),
                          MAP(v64_avg_u16),
                          MAP(v64_min_u8),
                          MAP(v64_max_u8),
                          MAP(v64_min_s8),
                          MAP(v64_max_s8),
                          MAP(v64_min_s16),
                          MAP(v64_max_s16),
                          MAP(v64_cmpgt_s8),
                          MAP(v64_cmplt_s8),
                          MAP(v64_cmpeq_8),
                          MAP(v64_cmpgt_s16),
                          MAP(v64_cmplt_s16),
                          MAP(v64_cmpeq_16),
                          MAP(v64_shuffle_8),
                          MAP(imm_v64_align<1>),
                          MAP(imm_v64_align<2>),
                          MAP(imm_v64_align<3>),
                          MAP(imm_v64_align<4>),
                          MAP(imm_v64_align<5>),
                          MAP(imm_v64_align<6>),
                          MAP(imm_v64_align<7>),
                          MAP(v128_abs_s16),
                          MAP(v128_padd_s16),
                          MAP(v128_unpacklo_u16_s32),
                          MAP(v128_unpacklo_s16_s32),
                          MAP(v128_unpackhi_u16_s32),
                          MAP(v128_unpackhi_s16_s32),
                          MAP(imm_v128_shr_n_byte<1>),
                          MAP(imm_v128_shr_n_byte<2>),
                          MAP(imm_v128_shr_n_byte<3>),
                          MAP(imm_v128_shr_n_byte<4>),
                          MAP(imm_v128_shr_n_byte<5>),
                          MAP(imm_v128_shr_n_byte<6>),
                          MAP(imm_v128_shr_n_byte<7>),
                          MAP(imm_v128_shr_n_byte<8>),
                          MAP(imm_v128_shr_n_byte<9>),
                          MAP(imm_v128_shr_n_byte<10>),
                          MAP(imm_v128_shr_n_byte<11>),
                          MAP(imm_v128_shr_n_byte<12>),
                          MAP(imm_v128_shr_n_byte<13>),
                          MAP(imm_v128_shr_n_byte<14>),
                          MAP(imm_v128_shr_n_byte<15>),
                          MAP(imm_v128_shl_n_byte<1>),
                          MAP(imm_v128_shl_n_byte<2>),
                          MAP(imm_v128_shl_n_byte<3>),
                          MAP(imm_v128_shl_n_byte<4>),
                          MAP(imm_v128_shl_n_byte<5>),
                          MAP(imm_v128_shl_n_byte<6>),
                          MAP(imm_v128_shl_n_byte<7>),
                          MAP(imm_v128_shl_n_byte<8>),
                          MAP(imm_v128_shl_n_byte<9>),
                          MAP(imm_v128_shl_n_byte<10>),
                          MAP(imm_v128_shl_n_byte<11>),
                          MAP(imm_v128_shl_n_byte<12>),
                          MAP(imm_v128_shl_n_byte<13>),
                          MAP(imm_v128_shl_n_byte<14>),
                          MAP(imm_v128_shl_n_byte<15>),
                          MAP(imm_v128_shl_n_8<1>),
                          MAP(imm_v128_shl_n_8<2>),
                          MAP(imm_v128_shl_n_8<3>),
                          MAP(imm_v128_shl_n_8<4>),
                          MAP(imm_v128_shl_n_8<5>),
                          MAP(imm_v128_shl_n_8<6>),
                          MAP(imm_v128_shl_n_8<7>),
                          MAP(imm_v128_shr_n_u8<1>),
                          MAP(imm_v128_shr_n_u8<2>),
                          MAP(imm_v128_shr_n_u8<3>),
                          MAP(imm_v128_shr_n_u8<4>),
                          MAP(imm_v128_shr_n_u8<5>),
                          MAP(imm_v128_shr_n_u8<6>),
                          MAP(imm_v128_shr_n_u8<7>),
                          MAP(imm_v128_shr_n_s8<1>),
                          MAP(imm_v128_shr_n_s8<2>),
                          MAP(imm_v128_shr_n_s8<3>),
                          MAP(imm_v128_shr_n_s8<4>),
                          MAP(imm_v128_shr_n_s8<5>),
                          MAP(imm_v128_shr_n_s8<6>),
                          MAP(imm_v128_shr_n_s8<7>),
                          MAP(imm_v128_shl_n_16<1>),
                          MAP(imm_v128_shl_n_16<2>),
                          MAP(imm_v128_shl_n_16<4>),
                          MAP(imm_v128_shl_n_16<6>),
                          MAP(imm_v128_shl_n_16<8>),
                          MAP(imm_v128_shl_n_16<10>),
                          MAP(imm_v128_shl_n_16<12>),
                          MAP(imm_v128_shl_n_16<14>),
                          MAP(imm_v128_shr_n_u16<1>),
                          MAP(imm_v128_shr_n_u16<2>),
                          MAP(imm_v128_shr_n_u16<4>),
                          MAP(imm_v128_shr_n_u16<6>),
                          MAP(imm_v128_shr_n_u16<8>),
                          MAP(imm_v128_shr_n_u16<10>),
                          MAP(imm_v128_shr_n_u16<12>),
                          MAP(imm_v128_shr_n_u16<14>),
                          MAP(imm_v128_shr_n_s16<1>),
                          MAP(imm_v128_shr_n_s16<2>),
                          MAP(imm_v128_shr_n_s16<4>),
                          MAP(imm_v128_shr_n_s16<6>),
                          MAP(imm_v128_shr_n_s16<8>),
                          MAP(imm_v128_shr_n_s16<10>),
                          MAP(imm_v128_shr_n_s16<12>),
                          MAP(imm_v128_shr_n_s16<14>),
                          MAP(imm_v128_shl_n_32<1>),
                          MAP(imm_v128_shl_n_32<4>),
                          MAP(imm_v128_shl_n_32<8>),
                          MAP(imm_v128_shl_n_32<12>),
                          MAP(imm_v128_shl_n_32<16>),
                          MAP(imm_v128_shl_n_32<20>),
                          MAP(imm_v128_shl_n_32<24>),
                          MAP(imm_v128_shl_n_32<28>),
                          MAP(imm_v128_shr_n_u32<1>),
                          MAP(imm_v128_shr_n_u32<4>),
                          MAP(imm_v128_shr_n_u32<8>),
                          MAP(imm_v128_shr_n_u32<12>),
                          MAP(imm_v128_shr_n_u32<16>),
                          MAP(imm_v128_shr_n_u32<20>),
                          MAP(imm_v128_shr_n_u32<24>),
                          MAP(imm_v128_shr_n_u32<28>),
                          MAP(imm_v128_shr_n_s32<1>),
                          MAP(imm_v128_shr_n_s32<4>),
                          MAP(imm_v128_shr_n_s32<8>),
                          MAP(imm_v128_shr_n_s32<12>),
                          MAP(imm_v128_shr_n_s32<16>),
                          MAP(imm_v128_shr_n_s32<20>),
                          MAP(imm_v128_shr_n_s32<24>),
                          MAP(imm_v128_shr_n_s32<28>),
                          MAP(v64_abs_s16),
                          MAP(v64_unpacklo_u8_s16),
                          MAP(v64_unpackhi_u8_s16),
                          MAP(v64_unpacklo_u16_s32),
                          MAP(v64_unpacklo_s16_s32),
                          MAP(v64_unpackhi_u16_s32),
                          MAP(v64_unpackhi_s16_s32),
                          MAP(imm_v64_shr_n_byte<1>),
                          MAP(imm_v64_shr_n_byte<2>),
                          MAP(imm_v64_shr_n_byte<3>),
                          MAP(imm_v64_shr_n_byte<4>),
                          MAP(imm_v64_shr_n_byte<5>),
                          MAP(imm_v64_shr_n_byte<6>),
                          MAP(imm_v64_shr_n_byte<7>),
                          MAP(imm_v64_shl_n_byte<1>),
                          MAP(imm_v64_shl_n_byte<2>),
                          MAP(imm_v64_shl_n_byte<3>),
                          MAP(imm_v64_shl_n_byte<4>),
                          MAP(imm_v64_shl_n_byte<5>),
                          MAP(imm_v64_shl_n_byte<6>),
                          MAP(imm_v64_shl_n_byte<7>),
                          MAP(imm_v64_shl_n_8<1>),
                          MAP(imm_v64_shl_n_8<2>),
                          MAP(imm_v64_shl_n_8<3>),
                          MAP(imm_v64_shl_n_8<4>),
                          MAP(imm_v64_shl_n_8<5>),
                          MAP(imm_v64_shl_n_8<6>),
                          MAP(imm_v64_shl_n_8<7>),
                          MAP(imm_v64_shr_n_u8<1>),
                          MAP(imm_v64_shr_n_u8<2>),
                          MAP(imm_v64_shr_n_u8<3>),
                          MAP(imm_v64_shr_n_u8<4>),
                          MAP(imm_v64_shr_n_u8<5>),
                          MAP(imm_v64_shr_n_u8<6>),
                          MAP(imm_v64_shr_n_u8<7>),
                          MAP(imm_v64_shr_n_s8<1>),
                          MAP(imm_v64_shr_n_s8<2>),
                          MAP(imm_v64_shr_n_s8<3>),
                          MAP(imm_v64_shr_n_s8<4>),
                          MAP(imm_v64_shr_n_s8<5>),
                          MAP(imm_v64_shr_n_s8<6>),
                          MAP(imm_v64_shr_n_s8<7>),
                          MAP(imm_v64_shl_n_16<1>),
                          MAP(imm_v64_shl_n_16<2>),
                          MAP(imm_v64_shl_n_16<4>),
                          MAP(imm_v64_shl_n_16<6>),
                          MAP(imm_v64_shl_n_16<8>),
                          MAP(imm_v64_shl_n_16<10>),
                          MAP(imm_v64_shl_n_16<12>),
                          MAP(imm_v64_shl_n_16<14>),
                          MAP(imm_v64_shr_n_u16<1>),
                          MAP(imm_v64_shr_n_u16<2>),
                          MAP(imm_v64_shr_n_u16<4>),
                          MAP(imm_v64_shr_n_u16<6>),
                          MAP(imm_v64_shr_n_u16<8>),
                          MAP(imm_v64_shr_n_u16<10>),
                          MAP(imm_v64_shr_n_u16<12>),
                          MAP(imm_v64_shr_n_u16<14>),
                          MAP(imm_v64_shr_n_s16<1>),
                          MAP(imm_v64_shr_n_s16<2>),
                          MAP(imm_v64_shr_n_s16<4>),
                          MAP(imm_v64_shr_n_s16<6>),
                          MAP(imm_v64_shr_n_s16<8>),
                          MAP(imm_v64_shr_n_s16<10>),
                          MAP(imm_v64_shr_n_s16<12>),
                          MAP(imm_v64_shr_n_s16<14>),
                          MAP(imm_v64_shl_n_32<1>),
                          MAP(imm_v64_shl_n_32<4>),
                          MAP(imm_v64_shl_n_32<8>),
                          MAP(imm_v64_shl_n_32<12>),
                          MAP(imm_v64_shl_n_32<16>),
                          MAP(imm_v64_shl_n_32<20>),
                          MAP(imm_v64_shl_n_32<24>),
                          MAP(imm_v64_shl_n_32<28>),
                          MAP(imm_v64_shr_n_u32<1>),
                          MAP(imm_v64_shr_n_u32<4>),
                          MAP(imm_v64_shr_n_u32<8>),
                          MAP(imm_v64_shr_n_u32<12>),
                          MAP(imm_v64_shr_n_u32<16>),
                          MAP(imm_v64_shr_n_u32<20>),
                          MAP(imm_v64_shr_n_u32<24>),
                          MAP(imm_v64_shr_n_u32<28>),
                          MAP(imm_v64_shr_n_s32<1>),
                          MAP(imm_v64_shr_n_s32<4>),
                          MAP(imm_v64_shr_n_s32<8>),
                          MAP(imm_v64_shr_n_s32<12>),
                          MAP(imm_v64_shr_n_s32<16>),
                          MAP(imm_v64_shr_n_s32<20>),
                          MAP(imm_v64_shr_n_s32<24>),
                          MAP(imm_v64_shr_n_s32<28>),
                          MAP(v128_from_v64),
                          MAP(v128_zip_8),
                          MAP(v128_zip_16),
                          MAP(v128_zip_32),
                          MAP(v128_mul_s16),
                          MAP(v128_unpack_u8_s16),
                          MAP(v128_unpack_u16_s32),
                          MAP(v128_unpack_s16_s32),
                          MAP(v128_shl_8),
                          MAP(v128_shr_u8),
                          MAP(v128_shr_s8),
                          MAP(v128_shl_16),
                          MAP(v128_shr_u16),
                          MAP(v128_shr_s16),
                          MAP(v128_shl_32),
                          MAP(v128_shr_u32),
                          MAP(v128_shr_s32),
                          MAP(v64_shl_8),
                          MAP(v64_shr_u8),
                          MAP(v64_shr_s8),
                          MAP(v64_shl_16),
                          MAP(v64_shr_u16),
                          MAP(v64_shr_s16),
                          MAP(v64_shl_32),
                          MAP(v64_shr_u32),
                          MAP(v64_shr_s32),
                          MAP(v64_hadd_u8),
                          MAP(v64_hadd_s16),
                          MAP(v128_hadd_u8),
                          MAP(v64_dotp_s16),
                          MAP(v128_dotp_s16),
                          MAP(v256_hadd_u8),
                          MAP(v256_dotp_s16),
                          MAP(v256_add_8),
                          MAP(v256_add_16),
                          MAP(v256_sadd_s16),
                          MAP(v256_add_32),
                          MAP(v256_sub_8),
                          MAP(v256_ssub_u8),
                          MAP(v256_ssub_s8),
                          MAP(v256_sub_16),
                          MAP(v256_ssub_s16),
                          MAP(v256_sub_32),
                          MAP(v256_ziplo_8),
                          MAP(v256_ziphi_8),
                          MAP(v256_ziplo_16),
                          MAP(v256_ziphi_16),
                          MAP(v256_ziplo_32),
                          MAP(v256_ziphi_32),
                          MAP(v256_ziplo_64),
                          MAP(v256_ziphi_64),
                          MAP(v256_unziphi_8),
                          MAP(v256_unziplo_8),
                          MAP(v256_unziphi_16),
                          MAP(v256_unziplo_16),
                          MAP(v256_unziphi_32),
                          MAP(v256_unziplo_32),
                          MAP(v256_pack_s32_s16),
                          MAP(v256_pack_s16_u8),
                          MAP(v256_pack_s16_s8),
                          MAP(v256_or),
                          MAP(v256_xor),
                          MAP(v256_and),
                          MAP(v256_andn),
                          MAP(v256_mullo_s16),
                          MAP(v256_mulhi_s16),
                          MAP(v256_mullo_s32),
                          MAP(v256_madd_s16),
                          MAP(v256_madd_us8),
                          MAP(v256_avg_u8),
                          MAP(v256_rdavg_u8),
                          MAP(v256_avg_u16),
                          MAP(v256_min_u8),
                          MAP(v256_max_u8),
                          MAP(v256_min_s8),
                          MAP(v256_max_s8),
                          MAP(v256_min_s16),
                          MAP(v256_max_s16),
                          MAP(v256_cmpgt_s8),
                          MAP(v256_cmplt_s8),
                          MAP(v256_cmpeq_8),
                          MAP(v256_cmpgt_s16),
                          MAP(v256_cmplt_s16),
                          MAP(v256_cmpeq_16),
                          MAP(v256_shuffle_8),
                          MAP(v256_pshuffle_8),
                          MAP(imm_v256_align<1>),
                          MAP(imm_v256_align<2>),
                          MAP(imm_v256_align<3>),
                          MAP(imm_v256_align<4>),
                          MAP(imm_v256_align<5>),
                          MAP(imm_v256_align<6>),
                          MAP(imm_v256_align<7>),
                          MAP(imm_v256_align<8>),
                          MAP(imm_v256_align<9>),
                          MAP(imm_v256_align<10>),
                          MAP(imm_v256_align<11>),
                          MAP(imm_v256_align<12>),
                          MAP(imm_v256_align<13>),
                          MAP(imm_v256_align<14>),
                          MAP(imm_v256_align<15>),
                          MAP(imm_v256_align<16>),
                          MAP(imm_v256_align<17>),
                          MAP(imm_v256_align<18>),
                          MAP(imm_v256_align<19>),
                          MAP(imm_v256_align<20>),
                          MAP(imm_v256_align<21>),
                          MAP(imm_v256_align<22>),
                          MAP(imm_v256_align<23>),
                          MAP(imm_v256_align<24>),
                          MAP(imm_v256_align<25>),
                          MAP(imm_v256_align<26>),
                          MAP(imm_v256_align<27>),
                          MAP(imm_v256_align<28>),
                          MAP(imm_v256_align<29>),
                          MAP(imm_v256_align<30>),
                          MAP(imm_v256_align<31>),
                          MAP(v256_from_v128),
                          MAP(v256_zip_8),
                          MAP(v256_zip_16),
                          MAP(v256_zip_32),
                          MAP(v256_mul_s16),
                          MAP(v256_unpack_u8_s16),
                          MAP(v256_unpack_u16_s32),
                          MAP(v256_unpack_s16_s32),
                          MAP(v256_shl_8),
                          MAP(v256_shr_u8),
                          MAP(v256_shr_s8),
                          MAP(v256_shl_16),
                          MAP(v256_shr_u16),
                          MAP(v256_shr_s16),
                          MAP(v256_shl_32),
                          MAP(v256_shr_u32),
                          MAP(v256_shr_s32),
                          MAP(v256_abs_s16),
                          MAP(v256_padd_s16),
                          MAP(v256_unpacklo_u16_s32),
                          MAP(v256_unpacklo_s16_s32),
                          MAP(v256_unpackhi_u16_s32),
                          MAP(v256_unpackhi_s16_s32),
                          MAP(imm_v256_shr_n_byte<1>),
                          MAP(imm_v256_shr_n_byte<2>),
                          MAP(imm_v256_shr_n_byte<3>),
                          MAP(imm_v256_shr_n_byte<4>),
                          MAP(imm_v256_shr_n_byte<5>),
                          MAP(imm_v256_shr_n_byte<6>),
                          MAP(imm_v256_shr_n_byte<7>),
                          MAP(imm_v256_shr_n_byte<8>),
                          MAP(imm_v256_shr_n_byte<9>),
                          MAP(imm_v256_shr_n_byte<10>),
                          MAP(imm_v256_shr_n_byte<11>),
                          MAP(imm_v256_shr_n_byte<12>),
                          MAP(imm_v256_shr_n_byte<13>),
                          MAP(imm_v256_shr_n_byte<14>),
                          MAP(imm_v256_shr_n_byte<15>),
                          MAP(imm_v256_shr_n_byte<16>),
                          MAP(imm_v256_shr_n_byte<17>),
                          MAP(imm_v256_shr_n_byte<18>),
                          MAP(imm_v256_shr_n_byte<19>),
                          MAP(imm_v256_shr_n_byte<20>),
                          MAP(imm_v256_shr_n_byte<21>),
                          MAP(imm_v256_shr_n_byte<22>),
                          MAP(imm_v256_shr_n_byte<23>),
                          MAP(imm_v256_shr_n_byte<24>),
                          MAP(imm_v256_shr_n_byte<25>),
                          MAP(imm_v256_shr_n_byte<26>),
                          MAP(imm_v256_shr_n_byte<27>),
                          MAP(imm_v256_shr_n_byte<28>),
                          MAP(imm_v256_shr_n_byte<29>),
                          MAP(imm_v256_shr_n_byte<30>),
                          MAP(imm_v256_shr_n_byte<31>),
                          MAP(imm_v256_shl_n_byte<1>),
                          MAP(imm_v256_shl_n_byte<2>),
                          MAP(imm_v256_shl_n_byte<3>),
                          MAP(imm_v256_shl_n_byte<4>),
                          MAP(imm_v256_shl_n_byte<5>),
                          MAP(imm_v256_shl_n_byte<6>),
                          MAP(imm_v256_shl_n_byte<7>),
                          MAP(imm_v256_shl_n_byte<8>),
                          MAP(imm_v256_shl_n_byte<9>),
                          MAP(imm_v256_shl_n_byte<10>),
                          MAP(imm_v256_shl_n_byte<11>),
                          MAP(imm_v256_shl_n_byte<12>),
                          MAP(imm_v256_shl_n_byte<13>),
                          MAP(imm_v256_shl_n_byte<14>),
                          MAP(imm_v256_shl_n_byte<15>),
                          MAP(imm_v256_shl_n_byte<16>),
                          MAP(imm_v256_shl_n_byte<17>),
                          MAP(imm_v256_shl_n_byte<18>),
                          MAP(imm_v256_shl_n_byte<19>),
                          MAP(imm_v256_shl_n_byte<20>),
                          MAP(imm_v256_shl_n_byte<21>),
                          MAP(imm_v256_shl_n_byte<22>),
                          MAP(imm_v256_shl_n_byte<23>),
                          MAP(imm_v256_shl_n_byte<24>),
                          MAP(imm_v256_shl_n_byte<25>),
                          MAP(imm_v256_shl_n_byte<26>),
                          MAP(imm_v256_shl_n_byte<27>),
                          MAP(imm_v256_shl_n_byte<28>),
                          MAP(imm_v256_shl_n_byte<29>),
                          MAP(imm_v256_shl_n_byte<30>),
                          MAP(imm_v256_shl_n_byte<31>),
                          MAP(imm_v256_shl_n_8<1>),
                          MAP(imm_v256_shl_n_8<2>),
                          MAP(imm_v256_shl_n_8<3>),
                          MAP(imm_v256_shl_n_8<4>),
                          MAP(imm_v256_shl_n_8<5>),
                          MAP(imm_v256_shl_n_8<6>),
                          MAP(imm_v256_shl_n_8<7>),
                          MAP(imm_v256_shr_n_u8<1>),
                          MAP(imm_v256_shr_n_u8<2>),
                          MAP(imm_v256_shr_n_u8<3>),
                          MAP(imm_v256_shr_n_u8<4>),
                          MAP(imm_v256_shr_n_u8<5>),
                          MAP(imm_v256_shr_n_u8<6>),
                          MAP(imm_v256_shr_n_u8<7>),
                          MAP(imm_v256_shr_n_s8<1>),
                          MAP(imm_v256_shr_n_s8<2>),
                          MAP(imm_v256_shr_n_s8<3>),
                          MAP(imm_v256_shr_n_s8<4>),
                          MAP(imm_v256_shr_n_s8<5>),
                          MAP(imm_v256_shr_n_s8<6>),
                          MAP(imm_v256_shr_n_s8<7>),
                          MAP(imm_v256_shl_n_16<1>),
                          MAP(imm_v256_shl_n_16<2>),
                          MAP(imm_v256_shl_n_16<4>),
                          MAP(imm_v256_shl_n_16<6>),
                          MAP(imm_v256_shl_n_16<8>),
                          MAP(imm_v256_shl_n_16<10>),
                          MAP(imm_v256_shl_n_16<12>),
                          MAP(imm_v256_shl_n_16<14>),
                          MAP(imm_v256_shr_n_u16<1>),
                          MAP(imm_v256_shr_n_u16<2>),
                          MAP(imm_v256_shr_n_u16<4>),
                          MAP(imm_v256_shr_n_u16<6>),
                          MAP(imm_v256_shr_n_u16<8>),
                          MAP(imm_v256_shr_n_u16<10>),
                          MAP(imm_v256_shr_n_u16<12>),
                          MAP(imm_v256_shr_n_u16<14>),
                          MAP(imm_v256_shr_n_s16<1>),
                          MAP(imm_v256_shr_n_s16<2>),
                          MAP(imm_v256_shr_n_s16<4>),
                          MAP(imm_v256_shr_n_s16<6>),
                          MAP(imm_v256_shr_n_s16<8>),
                          MAP(imm_v256_shr_n_s16<10>),
                          MAP(imm_v256_shr_n_s16<12>),
                          MAP(imm_v256_shr_n_s16<14>),
                          MAP(imm_v256_shl_n_32<1>),
                          MAP(imm_v256_shl_n_32<4>),
                          MAP(imm_v256_shl_n_32<8>),
                          MAP(imm_v256_shl_n_32<12>),
                          MAP(imm_v256_shl_n_32<16>),
                          MAP(imm_v256_shl_n_32<20>),
                          MAP(imm_v256_shl_n_32<24>),
                          MAP(imm_v256_shl_n_32<28>),
                          MAP(imm_v256_shr_n_u32<1>),
                          MAP(imm_v256_shr_n_u32<4>),
                          MAP(imm_v256_shr_n_u32<8>),
                          MAP(imm_v256_shr_n_u32<12>),
                          MAP(imm_v256_shr_n_u32<16>),
                          MAP(imm_v256_shr_n_u32<20>),
                          MAP(imm_v256_shr_n_u32<24>),
                          MAP(imm_v256_shr_n_u32<28>),
                          MAP(imm_v256_shr_n_s32<1>),
                          MAP(imm_v256_shr_n_s32<4>),
                          MAP(imm_v256_shr_n_s32<8>),
                          MAP(imm_v256_shr_n_s32<12>),
                          MAP(imm_v256_shr_n_s32<16>),
                          MAP(imm_v256_shr_n_s32<20>),
                          MAP(imm_v256_shr_n_s32<24>),
                          MAP(imm_v256_shr_n_s32<28>),
                          { NULL, NULL, NULL } };

  unsigned int i;
  for (i = 0; m[i][0] && strcmp(name, (const char *)m[i][0]); i++) {
  }

  *ref = m[i][1];
  *simd = m[i][2];
}

// Used for printing errors in test_simd1 and test_simd2
std::string print(uint8_t *a, int size) {
  std::string text = "0x";
  for (int i = 0; i < size; i++) {
    char buf[3];
    snprintf(buf, sizeof(buf), "%02x",
             a[!CONFIG_BIG_ENDIAN ? size - 1 - i : i]);
    text += buf;
  }

  return text;
}

// Used in test_simd1 and test_simd2 to restrict argument ranges
void setmask(uint8_t *s, int size, uint32_t mask, uint32_t maskwidth) {
  switch (maskwidth) {
    case 0: break;
    case 8:
      for (int i = 0; i < size; i++) s[i] &= mask;
      break;
    case 16: {
      uint16_t *t = (uint16_t *)s;
      for (int i = 0; i < size / 2; i++) t[i] &= mask;
      break;
    }
    case 32: {
      uint32_t *t = (uint32_t *)s;
      for (int i = 0; i < size / 4; i++) t[i] &= mask;
      break;
    }
    case 64: {
      uint64_t *t = (uint64_t *)s;
      for (int i = 0; i < size / 8; i++) t[i] &= mask;
      break;
    }
    default: FAIL() << "Unsupported mask width" << std::endl; break;
  }
}

// Used to test for type equality
template <class T, class U>
struct is_same {
  enum { value = 0 };
};

template <class T>
struct is_same<T, T> {
  enum { value = 1 };
};

// We need a store function for uint64_t
void u64_store_aligned(void *p, uint64_t a) {
  v64_store_aligned(p, v64_from_64(a));
}

void c_u64_store_aligned(void *p, uint64_t a) {
  c_v64_store_aligned(p, c_v64_from_64(a));
}

// compare_simd1 and compare_simd2 compare intrinsics taking 1 or 2
// arguments respectively with their corresponding C reference.
// Ideally, the loads and stores should be function templates, but
// v64 and v128 could be typedef'ed to the same type (which is the
// case on x86) and then we can't instantiate both v64 and v128.  So
// this workaround takes the combination of all types as template
// arguments and the functions as void pointers so we get no
// matching issues, and then cast the pointers into the proper
// function pointers.
template <typename ret, typename arg, typename c_ret, typename c_arg>
int compare_simd1(void *store, void *load, void *simd, void *d, void *c_store,
                  void *c_load, void *c_simd, void *ref_d, const void *a) {
  void (*my_store)(void *, ret) = (void (*)(void *, ret))store;
  arg (*my_load)(const void *) = (arg(*)(const void *))load;
  ret (*my_simd)(arg) = (ret(*)(arg))simd;
  void (*my_c_store)(void *, c_ret) = (void (*)(void *, c_ret))c_store;
  c_arg (*my_c_load)(const void *) = (c_arg(*)(const void *))c_load;
  c_ret (*my_c_simd)(c_arg) = (c_ret(*)(c_arg))c_simd;

  // Call reference and intrinsic
  ASM_REGISTER_STATE_CHECK(my_c_store(ref_d, my_c_simd(my_c_load(a))));
  ASM_REGISTER_STATE_CHECK(my_store(d, my_simd(my_load(a))));

  // Compare results
  return memcmp(ref_d, d, sizeof(c_ret));
}

template <typename ret, typename arg1, typename arg2, typename c_ret,
          typename c_arg1, typename c_arg2>
int compare_simd2(void *store, void *load1, void *load2, void *simd, void *d,
                  void *c_store, void *c_load1, void *c_load2, void *c_simd,
                  void *ref_d, const void *a, const void *b) {
  void (*my_store)(void *, ret) = (void (*)(void *, ret))store;
  arg1 (*my_load1)(const void *) = (arg1(*)(const void *))load1;
  arg2 (*my_load2)(const void *) = (arg2(*)(const void *))load2;
  ret (*my_simd)(arg1, arg2) = (ret(*)(arg1, arg2))simd;
  void (*my_c_store)(void *, c_ret) = (void (*)(void *, c_ret))c_store;
  c_arg1 (*my_c_load1)(const void *) = (c_arg1(*)(const void *))c_load1;
  c_arg2 (*my_c_load2)(const void *) = (c_arg2(*)(const void *))c_load2;
  c_ret (*my_c_simd)(c_arg1, c_arg2) = (c_ret(*)(c_arg1, c_arg2))c_simd;

  // Call reference and intrinsic
  ASM_REGISTER_STATE_CHECK(
      my_c_store(ref_d, my_c_simd(my_c_load1(a), my_c_load2(b))));
  ASM_REGISTER_STATE_CHECK(my_store(d, my_simd(my_load1(a), my_load2(b))));

  // Compare results
  return memcmp(ref_d, d, sizeof(c_ret));
}

template <typename c_ret, typename c_arg>
void test_simd1(uint32_t iterations, uint32_t mask, uint32_t maskwidth,
                const char *name) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  void *ref_simd;
  void *simd;
  int error = 0;
  DECLARE_ALIGNED(sizeof(c_arg), uint16_t, s[sizeof(c_arg) / sizeof(uint16_t)]);
  DECLARE_ALIGNED(sizeof(c_ret), uint8_t, d[sizeof(c_ret)]);
  DECLARE_ALIGNED(sizeof(c_ret), uint8_t, ref_d[sizeof(c_ret)]);
  memset(ref_d, 0, sizeof(ref_d));
  memset(d, 0, sizeof(d));

  map(name, &ref_simd, &simd);
  if (simd == NULL || ref_simd == NULL) {
    FAIL() << "Internal errer: Unknown intrinsic function " << name;
  }

  for (unsigned int count = 0; count < iterations && !error; count++) {
    for (unsigned int c = 0; c < sizeof(c_arg) / sizeof(uint16_t); c++)
      s[c] = rnd.Rand16();

    if (maskwidth) setmask((uint8_t *)s, sizeof(c_arg), mask, maskwidth);

    if (is_same<c_ret, c_v64>::value && is_same<c_arg, c_v64>::value) {
      // V64_V64
      error = compare_simd1<v64, v64, c_ret, c_arg>(
          (void *)v64_store_aligned, (void *)v64_load_aligned, simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, uint64_t>::value &&
               is_same<c_arg, c_v64>::value) {
      // U64_V64
      error = compare_simd1<uint64_t, v64, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v64_load_aligned, simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, int64_t>::value && is_same<c_arg, c_v64>::value) {
      // S64_V64
      error = compare_simd1<int64_t, v64, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v64_load_aligned, simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, uint64_t>::value &&
               is_same<c_arg, c_v256>::value) {
      // U64_V256
      error = compare_simd1<uint64_t, v256, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v256_load_aligned, simd, d,
          (void *)c_u64_store_aligned, (void *)c_v256_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, uint64_t>::value &&
               is_same<c_arg, c_v128>::value) {
      // U64_V128
      error = compare_simd1<uint64_t, v128, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v128_load_aligned, simd, d,
          (void *)c_u64_store_aligned, (void *)c_v128_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, c_v128>::value && is_same<c_arg, c_v128>::value) {
      // V128_V128
      error = compare_simd1<v128, v128, c_ret, c_arg>(
          (void *)v128_store_aligned, (void *)v128_load_aligned, simd, d,
          (void *)c_v128_store_aligned, (void *)c_v128_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, c_v256>::value && is_same<c_arg, c_v256>::value) {
      // V256_V256
      error = compare_simd1<v256, v256, c_ret, c_arg>(
          (void *)v256_store_aligned, (void *)v256_load_aligned, simd, d,
          (void *)c_v256_store_aligned, (void *)c_v256_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, c_v256>::value && is_same<c_arg, c_v128>::value) {
      // V256_V128
      error = compare_simd1<v256, v128, c_ret, c_arg>(
          (void *)v256_store_aligned, (void *)v128_load_aligned, simd, d,
          (void *)c_v256_store_aligned, (void *)c_v128_load_aligned, ref_simd,
          ref_d, s);
    } else if (is_same<c_ret, c_v128>::value && is_same<c_arg, c_v64>::value) {
      // V128_V64
      error = compare_simd1<v128, v64, c_ret, c_arg>(
          (void *)v128_store_aligned, (void *)v64_load_aligned, simd, d,
          (void *)c_v128_store_aligned, (void *)c_v64_load_aligned, ref_simd,
          ref_d, s);
    } else {
      FAIL() << "Internal errer: Unknown intrinsic function "
             << typeid(c_ret).name() << " " << name << "("
             << typeid(c_arg).name() << ")";
    }
  }

  EXPECT_EQ(0, error) << "Error: mismatch for " << name << "("
                      << print((uint8_t *)s, sizeof(s)) << ") -> "
                      << print(d, sizeof(d)) << " (simd), "
                      << print(ref_d, sizeof(ref_d)) << " (ref)" << std::endl;
}

template <typename c_ret, typename c_arg1, typename c_arg2>
void test_simd2(uint32_t iterations, uint32_t mask, uint32_t maskwidth,
                const char *name) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  void *ref_simd;
  void *simd;
  int error = 0;
  DECLARE_ALIGNED(sizeof(c_arg1), uint16_t,
                  s1[sizeof(c_arg1) / sizeof(uint16_t)]);
  DECLARE_ALIGNED(sizeof(c_arg2), uint16_t,
                  s2[sizeof(c_arg2) / sizeof(uint16_t)]);
  DECLARE_ALIGNED(sizeof(c_ret), uint8_t, d[sizeof(c_ret)]);
  DECLARE_ALIGNED(sizeof(c_ret), uint8_t, ref_d[sizeof(c_ret)]);
  memset(ref_d, 0, sizeof(ref_d));
  memset(d, 0, sizeof(d));

  map(name, &ref_simd, &simd);
  if (simd == NULL || ref_simd == NULL) {
    FAIL() << "Internal errer: Unknown intrinsic function " << name;
  }

  for (unsigned int count = 0; count < iterations && !error; count++) {
    for (unsigned int c = 0; c < sizeof(c_arg1) / sizeof(uint16_t); c++)
      s1[c] = rnd.Rand16();

    for (unsigned int c = 0; c < sizeof(c_arg2) / sizeof(uint16_t); c++)
      s2[c] = rnd.Rand16();

    if (maskwidth) setmask((uint8_t *)s2, sizeof(c_arg2), mask, maskwidth);

    if (is_same<c_ret, c_v64>::value && is_same<c_arg1, c_v64>::value &&
        is_same<c_arg2, c_v64>::value) {
      // V64_V64V64
      error = compare_simd2<v64, v64, v64, c_ret, c_arg1, c_arg2>(
          (void *)v64_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, simd, d, (void *)c_v64_store_aligned,
          (void *)c_v64_load_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v128>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, c_v128>::value) {
      // V128_V128V128
      error = compare_simd2<v128, v128, v128, c_ret, c_arg1, c_arg2>(
          (void *)v128_store_aligned, (void *)v128_load_aligned,
          (void *)v128_load_aligned, simd, d, (void *)c_v128_store_aligned,
          (void *)c_v128_load_aligned, (void *)c_v128_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v256>::value &&
               is_same<c_arg1, c_v256>::value &&
               is_same<c_arg2, c_v256>::value) {
      // V256_V256V256
      error = compare_simd2<v256, v256, v256, c_ret, c_arg1, c_arg2>(
          (void *)v256_store_aligned, (void *)v256_load_aligned,
          (void *)v256_load_aligned, simd, d, (void *)c_v256_store_aligned,
          (void *)c_v256_load_aligned, (void *)c_v256_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, uint32_t>::value &&
               is_same<c_arg1, c_v64>::value && is_same<c_arg2, c_v64>::value) {
      // U32_V64V64
      error = compare_simd2<uint32_t, v64, v64, c_ret, c_arg1, c_arg2>(
          (void *)u32_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, simd, d, (void *)c_u32_store_aligned,
          (void *)c_v64_load_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, int64_t>::value &&
               is_same<c_arg1, c_v64>::value && is_same<c_arg2, c_v64>::value) {
      // S64_V64V64
      error = compare_simd2<int64_t, v64, v64, c_ret, c_arg1, c_arg2>(
          (void *)u64_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, simd, d, (void *)c_u64_store_aligned,
          (void *)c_v64_load_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, uint32_t>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, c_v128>::value) {
      // U32_V128V128
      error = compare_simd2<uint32_t, v128, v128, c_ret, c_arg1, c_arg2>(
          (void *)u32_store_aligned, (void *)v128_load_aligned,
          (void *)v128_load_aligned, simd, d, (void *)c_u32_store_aligned,
          (void *)c_v128_load_aligned, (void *)c_v128_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, int64_t>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, c_v128>::value) {
      // S64_V128V128
      error = compare_simd2<int64_t, v128, v128, c_ret, c_arg1, c_arg2>(
          (void *)u64_store_aligned, (void *)v128_load_aligned,
          (void *)v128_load_aligned, simd, d, (void *)c_u64_store_aligned,
          (void *)c_v128_load_aligned, (void *)c_v128_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, int64_t>::value &&
               is_same<c_arg1, c_v256>::value &&
               is_same<c_arg2, c_v256>::value) {
      // S64_V256V256
      error = compare_simd2<int64_t, v256, v256, c_ret, c_arg1, c_arg2>(
          (void *)u64_store_aligned, (void *)v256_load_aligned,
          (void *)v256_load_aligned, simd, d, (void *)c_u64_store_aligned,
          (void *)c_v256_load_aligned, (void *)c_v256_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, uint32_t>::value &&
               is_same<c_arg1, c_v256>::value &&
               is_same<c_arg2, c_v256>::value) {
      // U32_V256V256
      error = compare_simd2<uint32_t, v256, v256, c_ret, c_arg1, c_arg2>(
          (void *)u32_store_aligned, (void *)v256_load_aligned,
          (void *)v256_load_aligned, simd, d, (void *)c_u32_store_aligned,
          (void *)c_v256_load_aligned, (void *)c_v256_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v256>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, c_v128>::value) {
      // V256_V128V128
      error = compare_simd2<v256, v128, v128, c_ret, c_arg1, c_arg2>(
          (void *)v256_store_aligned, (void *)v128_load_aligned,
          (void *)v128_load_aligned, simd, d, (void *)c_v256_store_aligned,
          (void *)c_v128_load_aligned, (void *)c_v128_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v128>::value && is_same<c_arg1, c_v64>::value &&
               is_same<c_arg2, c_v64>::value) {
      // V128_V64V64
      error = compare_simd2<v128, v64, v64, c_ret, c_arg1, c_arg2>(
          (void *)v128_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, simd, d, (void *)c_v128_store_aligned,
          (void *)c_v64_load_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v256>::value &&
               is_same<c_arg1, c_v256>::value &&
               is_same<c_arg2, uint32_t>::value) {
      // V256_V256U32
      error = compare_simd2<v256, v256, uint32_t, c_ret, c_arg1, c_arg2>(
          (void *)v256_store_aligned, (void *)v256_load_aligned,
          (void *)u32_load_aligned, simd, d, (void *)c_v256_store_aligned,
          (void *)c_v256_load_aligned, (void *)c_u32_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v128>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, uint32_t>::value) {
      // V128_V128U32
      error = compare_simd2<v128, v128, uint32_t, c_ret, c_arg1, c_arg2>(
          (void *)v128_store_aligned, (void *)v128_load_aligned,
          (void *)u32_load_aligned, simd, d, (void *)c_v128_store_aligned,
          (void *)c_v128_load_aligned, (void *)c_u32_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else if (is_same<c_ret, c_v64>::value && is_same<c_arg1, c_v64>::value &&
               is_same<c_arg2, uint32_t>::value) {
      // V64_V64U32
      error = compare_simd2<v64, v64, uint32_t, c_ret, c_arg1, c_arg2>(
          (void *)v64_store_aligned, (void *)v64_load_aligned,
          (void *)u32_load_aligned, simd, d, (void *)c_v64_store_aligned,
          (void *)c_v64_load_aligned, (void *)c_u32_load_aligned,
          (void *)ref_simd, ref_d, s1, s2);
    } else {
      FAIL() << "Internal errer: Unknown intrinsic function "
             << typeid(c_ret).name() << " " << name << "("
             << typeid(c_arg1).name() << ", " << typeid(c_arg2).name() << ")";
    }
  }

  EXPECT_EQ(0, error) << "Error: mismatch for " << name << "("
                      << print((uint8_t *)s1, sizeof(s1)) << ", "
                      << print((uint8_t *)s2, sizeof(s2)) << ") -> "
                      << print(d, sizeof(d)) << " (simd), "
                      << print(ref_d, sizeof(ref_d)) << " (ref)" << std::endl;
}

// Instantiations to make the functions callable from another files
template void test_simd1<c_v128, c_v128>(uint32_t, uint32_t, uint32_t,
                                         const char *);
template void test_simd1<c_v128, c_v64>(uint32_t, uint32_t, uint32_t,
                                        const char *);
template void test_simd1<c_v256, c_v128>(uint32_t, uint32_t, uint32_t,
                                         const char *);
template void test_simd1<c_v256, c_v256>(uint32_t, uint32_t, uint32_t,
                                         const char *);
template void test_simd1<c_v64, c_v64>(uint32_t, uint32_t, uint32_t,
                                       const char *);
template void test_simd1<int64_t, c_v64>(uint32_t, uint32_t, uint32_t,
                                         const char *);
template void test_simd1<uint64_t, c_v128>(uint32_t, uint32_t, uint32_t,
                                           const char *);
template void test_simd1<uint64_t, c_v256>(uint32_t, uint32_t, uint32_t,
                                           const char *);
template void test_simd1<uint64_t, c_v64>(uint32_t, uint32_t, uint32_t,
                                          const char *);
template void test_simd2<c_v128, c_v128, c_v128>(uint32_t, uint32_t, uint32_t,
                                                 const char *);
template void test_simd2<c_v128, c_v128, uint32_t>(uint32_t, uint32_t, uint32_t,
                                                   const char *);
template void test_simd2<c_v128, c_v64, c_v64>(uint32_t, uint32_t, uint32_t,
                                               const char *);
template void test_simd2<c_v256, c_v128, c_v128>(uint32_t, uint32_t, uint32_t,
                                                 const char *);
template void test_simd2<c_v256, c_v256, c_v256>(uint32_t, uint32_t, uint32_t,
                                                 const char *);
template void test_simd2<c_v256, c_v256, uint32_t>(uint32_t, uint32_t, uint32_t,
                                                   const char *);
template void test_simd2<c_v64, c_v64, c_v64>(uint32_t, uint32_t, uint32_t,
                                              const char *);
template void test_simd2<c_v64, c_v64, uint32_t>(uint32_t, uint32_t, uint32_t,
                                                 const char *);
template void test_simd2<int64_t, c_v128, c_v128>(uint32_t, uint32_t, uint32_t,
                                                  const char *);
template void test_simd2<int64_t, c_v256, c_v256>(uint32_t, uint32_t, uint32_t,
                                                  const char *);
template void test_simd2<int64_t, c_v64, c_v64>(uint32_t, uint32_t, uint32_t,
                                                const char *);
template void test_simd2<uint32_t, c_v128, c_v128>(uint32_t, uint32_t, uint32_t,
                                                   const char *);
template void test_simd2<uint32_t, c_v256, c_v256>(uint32_t, uint32_t, uint32_t,
                                                   const char *);
template void test_simd2<uint32_t, c_v64, c_v64>(uint32_t, uint32_t, uint32_t,
                                                 const char *);

}  // namespace SIMD_NAMESPACE
