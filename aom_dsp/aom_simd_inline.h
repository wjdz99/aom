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

#ifndef AOM_DSP_AOM_SIMD_INLINE_H_
#define AOM_DSP_AOM_SIMD_INLINE_H_

#ifndef SIMD_INLINE
#ifdef __GNUC__
#define SIMD_INLINE static inline __attribute__((always_inline))
#elif __STDC_VERSION__ >= 199901L
#define SIMD_INLINE static inline
#elif defined(_MSC_VER)
#define SIMD_INLINE static __inline
#else
#define SIMD_INLINE static
#endif
#endif

#endif /* AOM_DSP_AOM_SIMD_INLINE_H_ */
