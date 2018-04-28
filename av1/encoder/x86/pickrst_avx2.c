/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <immintrin.h>  // avx2

#include "./av1_rtcd.h"
#include "./aom_dsp_rtcd.h"

#include "av1/common/restoration.h"

int64_t av1_get_pixel_proj_error_avx2(const uint8_t *src8, int width,
                                      int height, int src_stride,
                                      const uint8_t *dat8, int dat_stride,
                                      int use_highbitdepth, int32_t *flt0,
                                      int flt0_stride, int32_t *flt1,
                                      int flt1_stride, int *xqd,
                                      const sgr_params_type *params) {
  int64_t err = 0;
  int xq[2];
  decode_xq(xqd, xq, params);
  if (!use_highbitdepth) {
    const uint8_t *src = src8;
    const uint8_t *dat = dat8;
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        const int32_t u =
            (int32_t)(dat[i * dat_stride + j] << SGRPROJ_RST_BITS);
        int32_t v = u << SGRPROJ_PRJ_BITS;
        if (params->r0 > 0) v += xq[0] * (flt0[i * flt0_stride + j] - u);
        if (params->r1 > 0) v += xq[1] * (flt1[i * flt1_stride + j] - u);
        const int32_t e =
            ROUND_POWER_OF_TWO(v, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS) -
            src[i * src_stride + j];
        err += e * e;
      }
    }
  } else {
    const uint16_t *src = CONVERT_TO_SHORTPTR(src8);
    const uint16_t *dat = CONVERT_TO_SHORTPTR(dat8);
    const int32_t half = 1 << (SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS - 1);
    const __m256i vhalf = _mm256_set1_epi32(half);
    const int n16 = width / 16;
    __m256i verr = _mm256_setzero_si256();
    if (params->r0 > 0 && params->r1 > 0) {
      const int xq0 = xq[0];
      const int xq1 = xq[1];
      const __m256i vxq0 = _mm256_set1_epi32(xq0);
      const __m256i vxq1 = _mm256_set1_epi32(xq1);
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < n16; ++j) {
          __m256i d = _mm256_loadu_si256((const __m256i *)(&dat[j * 16]));
          __m256i s = _mm256_loadu_si256((const __m256i *)(&src[j * 16]));
          __m256i dl = _mm256_unpacklo_epi16(d, _mm256_setzero_si256());
          __m256i dh = _mm256_unpackhi_epi16(d, _mm256_setzero_si256());
          __m256i sl = _mm256_unpacklo_epi16(s, _mm256_setzero_si256());
          __m256i sh = _mm256_unpackhi_epi16(s, _mm256_setzero_si256());

          __m256i tmp;
          tmp = _mm256_permute2x128_si256(dl, dh, 0x20);
          dh = _mm256_permute2x128_si256(dl, dh, 0x31);
          dl = tmp;
          tmp = _mm256_permute2x128_si256(sl, sh, 0x20);
          sh = _mm256_permute2x128_si256(sl, sh, 0x31);
          sl = tmp;

          __m256i vflt0l =
              _mm256_loadu_si256((const __m256i *)&flt0[j * 16 + 0]);
          __m256i vflt0h =
              _mm256_loadu_si256((const __m256i *)&flt0[j * 16 + 8]);
          __m256i vflt1l =
              _mm256_loadu_si256((const __m256i *)&flt1[j * 16 + 0]);
          __m256i vflt1h =
              _mm256_loadu_si256((const __m256i *)&flt1[j * 16 + 8]);
          __m256i ul = _mm256_slli_epi32(dl, SGRPROJ_RST_BITS);
          __m256i uh = _mm256_slli_epi32(dh, SGRPROJ_RST_BITS);
          dl = _mm256_sub_epi32(dl, sl);
          dh = _mm256_sub_epi32(dh, sh);

          __m256i v0l = _mm256_sub_epi32(vflt0l, ul);
          __m256i v0h = _mm256_sub_epi32(vflt0h, uh);
          __m256i v1l = _mm256_sub_epi32(vflt1l, ul);
          __m256i v1h = _mm256_sub_epi32(vflt1h, uh);

          v0l = _mm256_mullo_epi32(vxq0, v0l);
          v0h = _mm256_mullo_epi32(vxq0, v0h);
          v1l = _mm256_mullo_epi32(vxq1, v1l);
          v1h = _mm256_mullo_epi32(vxq1, v1h);

          v0l = _mm256_add_epi32(v0l, v1l);
          v0h = _mm256_add_epi32(v0h, v1h);
          v0l = _mm256_add_epi32(v0l, vhalf);
          v0h = _mm256_add_epi32(v0h, vhalf);

          v0l = _mm256_srai_epi32(v0l, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS);
          v0h = _mm256_srai_epi32(v0h, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS);

          v0l = _mm256_add_epi32(v0l, dl);
          v0h = _mm256_add_epi32(v0h, dh);

          v0l = _mm256_mullo_epi32(v0l, v0l);
          v0h = _mm256_mullo_epi32(v0h, v0h);

          verr = _mm256_add_epi32(verr, v0l);
          verr = _mm256_add_epi32(verr, v0h);
        }
        for (int j = n16 * 16; j < width; ++j) {
          const int32_t d = dat[j];
          const int32_t s = src[j];
          const int32_t u = (int32_t)(d << SGRPROJ_RST_BITS);
          int32_t v0 = flt0[j] - u;
          int32_t v1 = flt1[j] - u;
          int32_t v = half;
          v += xq0 * v0;
          v += xq1 * v1;
          const int32_t e =
              (v >> (SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS)) + d - s;
          err += e * e;
        }
        dat += dat_stride;
        flt0 += flt0_stride;
        flt1 += flt1_stride;
        src += src_stride;
      }
    } else if (params->r0 > 0 || params->r1 > 0) {
      int exq;
      int32_t *flt;
      int flt_stride;
      if (params->r0 > 0) {
        exq = xq[0];
        flt = flt0;
        flt_stride = flt0_stride;
      } else {
        exq = xq[1];
        flt = flt1;
        flt_stride = flt1_stride;
      }
      const __m256i vxq = _mm256_set1_epi32(exq);
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < n16; ++j) {
          __m256i d = _mm256_loadu_si256((const __m256i *)(&dat[j * 16]));
          __m256i s = _mm256_loadu_si256((const __m256i *)(&src[j * 16]));
          __m256i dl = _mm256_unpacklo_epi16(d, _mm256_setzero_si256());
          __m256i dh = _mm256_unpackhi_epi16(d, _mm256_setzero_si256());
          __m256i sl = _mm256_unpacklo_epi16(s, _mm256_setzero_si256());
          __m256i sh = _mm256_unpackhi_epi16(s, _mm256_setzero_si256());

          __m256i tmp;
          tmp = _mm256_permute2x128_si256(dl, dh, 0x20);
          dh = _mm256_permute2x128_si256(dl, dh, 0x31);
          dl = tmp;
          tmp = _mm256_permute2x128_si256(sl, sh, 0x20);
          sh = _mm256_permute2x128_si256(sl, sh, 0x31);
          sl = tmp;

          __m256i vflt0l =
              _mm256_loadu_si256((const __m256i *)&flt[j * 16 + 0]);
          __m256i vflt0h =
              _mm256_loadu_si256((const __m256i *)&flt[j * 16 + 8]);
          __m256i ul = _mm256_slli_epi32(dl, SGRPROJ_RST_BITS);
          __m256i uh = _mm256_slli_epi32(dh, SGRPROJ_RST_BITS);
          dl = _mm256_sub_epi32(dl, sl);
          dh = _mm256_sub_epi32(dh, sh);

          __m256i v0l = _mm256_sub_epi32(vflt0l, ul);
          __m256i v0h = _mm256_sub_epi32(vflt0h, uh);

          v0l = _mm256_mullo_epi32(vxq, v0l);
          v0h = _mm256_mullo_epi32(vxq, v0h);

          v0l = _mm256_add_epi32(v0l, vhalf);
          v0h = _mm256_add_epi32(v0h, vhalf);

          v0l = _mm256_srai_epi32(v0l, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS);
          v0h = _mm256_srai_epi32(v0h, SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS);

          v0l = _mm256_add_epi32(v0l, dl);
          v0h = _mm256_add_epi32(v0h, dh);

          v0l = _mm256_mullo_epi32(v0l, v0l);
          v0h = _mm256_mullo_epi32(v0h, v0h);

          verr = _mm256_add_epi32(verr, v0l);
          verr = _mm256_add_epi32(verr, v0h);
        }
        for (int j = n16 * 16; j < width; ++j) {
          const int32_t d = dat[j];
          const int32_t s = src[j];
          const int32_t u = (int32_t)(d << SGRPROJ_RST_BITS);
          int32_t v = half;
          v += exq * (flt[j] - u);
          const int32_t e =
              (v >> (SGRPROJ_RST_BITS + SGRPROJ_PRJ_BITS)) + d - s;
          err += e * e;
        }
        dat += dat_stride;
        flt += flt_stride;
        src += src_stride;
      }
    } else {
      for (int i = 0; i < height; ++i) {
        for (int j = 0; j < n16; ++j) {
          __m256i d = _mm256_loadu_si256((const __m256i *)(&dat[j * 16]));
          __m256i s = _mm256_loadu_si256((const __m256i *)(&src[j * 16]));
          __m256i dl = _mm256_unpacklo_epi16(d, _mm256_setzero_si256());
          __m256i dh = _mm256_unpackhi_epi16(d, _mm256_setzero_si256());
          __m256i sl = _mm256_unpacklo_epi16(s, _mm256_setzero_si256());
          __m256i sh = _mm256_unpackhi_epi16(s, _mm256_setzero_si256());
          dl = _mm256_sub_epi32(dl, sl);
          dh = _mm256_sub_epi32(dh, sh);
          verr = _mm256_add_epi32(verr, dl);
          verr = _mm256_add_epi32(verr, dh);
        }
        for (int j = n16 * 16; j < width; ++j) {
          const int32_t d = dat[j];
          const int32_t s = src[j];
          const int32_t e = d - s;
          err += e * e;
        }
        dat += dat_stride;
        src += src_stride;
      }
    }
    err += _mm256_extract_epi32(verr, 0);
    err += _mm256_extract_epi32(verr, 1);
    err += _mm256_extract_epi32(verr, 2);
    err += _mm256_extract_epi32(verr, 3);
    err += _mm256_extract_epi32(verr, 4);
    err += _mm256_extract_epi32(verr, 5);
    err += _mm256_extract_epi32(verr, 6);
    err += _mm256_extract_epi32(verr, 7);
  }
  return err;
}
