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

#ifndef _AOM_SIMD_H
#define _AOM_SIMD_H

#ifndef SIMD_INLINE
#ifdef __GNUC__
#define SIMD_INLINE static inline __attribute__((always_inline))
#elif __STDC_VERSION__ >= 199901L
#define SIMD_INLINE static inline
#else
#define SIMD_INLINE static
#endif
#endif

#include <stdint.h>

#if defined(_WIN32)
#include <intrin.h>
#endif

static const int simd_check = 1;

#if defined(__ARM_NEON__)
static const int simd_available = 1;
#include "simd/v128_intrinsics_arm.h"
#elif (defined(__SSE2__) || _M_IX86_FP == 2)
static const int simd_available = 1;
#include "simd/v128_intrinsics_x86.h"
#else
static const int simd_available = 0;
#include "simd/v128_intrinsics.h"
#endif

extern int aom_use_simd;
SIMD_INLINE void aom_init_use_simd() {
  /* SIMD optimisations supported only for little endian architectures */
  const uint16_t t = 0x100;
  aom_use_simd = simd_available && !*(const uint8_t *)&t;
}

#endif /* _AOM_SIMD_H */
