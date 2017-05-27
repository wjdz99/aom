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

#include <math.h>

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"

#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"

#define DST(x, y) dst[(x) + (y)*stride]
#define AVG3(a, b, c) (((a) + 2 * (b) + (c) + 2) >> 2)
#define AVG2(a, b) (((a) + (b) + 1) >> 1)

static INLINE void d207_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  int r, c;
  (void)above;
  // first column
  for (r = 0; r < bs - 1; ++r) dst[r * stride] = AVG2(left[r], left[r + 1]);
  dst[(bs - 1) * stride] = left[bs - 1];
  dst++;

  // second column
  for (r = 0; r < bs - 2; ++r)
    dst[r * stride] = AVG3(left[r], left[r + 1], left[r + 2]);
  dst[(bs - 2) * stride] = AVG3(left[bs - 2], left[bs - 1], left[bs - 1]);
  dst[(bs - 1) * stride] = left[bs - 1];
  dst++;

  // rest of last row
  for (c = 0; c < bs - 2; ++c) dst[(bs - 1) * stride + c] = left[bs - 1];

  for (r = bs - 2; r >= 0; --r)
    for (c = 0; c < bs - 2; ++c)
      dst[r * stride + c] = dst[(r + 1) * stride + c - 2];
}

static INLINE void d207e_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                   const uint8_t *above, const uint8_t *left) {
  int r, c;
  (void)above;

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = c & 1 ? AVG3(left[(c >> 1) + r], left[(c >> 1) + r + 1],
                            left[(c >> 1) + r + 2])
                     : AVG2(left[(c >> 1) + r], left[(c >> 1) + r + 1]);
    }
    dst += stride;
  }
}

static INLINE void d63_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                 const uint8_t *above, const uint8_t *left) {
  int r, c;
  int size;
  (void)left;
  for (c = 0; c < bs; ++c) {
    dst[c] = AVG2(above[c], above[c + 1]);
    dst[stride + c] = AVG3(above[c], above[c + 1], above[c + 2]);
  }
  for (r = 2, size = bs - 2; r < bs; r += 2, --size) {
    memcpy(dst + (r + 0) * stride, dst + (r >> 1), size);
    memset(dst + (r + 0) * stride + size, above[bs - 1], bs - size);
    memcpy(dst + (r + 1) * stride, dst + stride + (r >> 1), size);
    memset(dst + (r + 1) * stride + size, above[bs - 1], bs - size);
  }
}

static INLINE void d63e_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  int r, c;
  (void)left;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = r & 1 ? AVG3(above[(r >> 1) + c], above[(r >> 1) + c + 1],
                            above[(r >> 1) + c + 2])
                     : AVG2(above[(r >> 1) + c], above[(r >> 1) + c + 1]);
    }
    dst += stride;
  }
}

static INLINE void d45_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                 const uint8_t *above, const uint8_t *left) {
  const uint8_t above_right = above[bs - 1];
  const uint8_t *const dst_row0 = dst;
  int x, size;
  (void)left;

  for (x = 0; x < bs - 1; ++x) {
    dst[x] = AVG3(above[x], above[x + 1], above[x + 2]);
  }
  dst[bs - 1] = above_right;
  dst += stride;
  for (x = 1, size = bs - 2; x < bs; ++x, --size) {
    memcpy(dst, dst_row0 + x, size);
    memset(dst + size, above_right, x + 1);
    dst += stride;
  }
}

static INLINE void d45e_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  int r, c;
  (void)left;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = AVG3(above[r + c], above[r + c + 1],
                    above[r + c + 1 + (r + c + 2 < bs * 2)]);
    }
    dst += stride;
  }
}

static INLINE void d117_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  int r, c;

  // first row
  for (c = 0; c < bs; c++) dst[c] = AVG2(above[c - 1], above[c]);
  dst += stride;

  // second row
  dst[0] = AVG3(left[0], above[-1], above[0]);
  for (c = 1; c < bs; c++) dst[c] = AVG3(above[c - 2], above[c - 1], above[c]);
  dst += stride;

  // the rest of first col
  dst[0] = AVG3(above[-1], left[0], left[1]);
  for (r = 3; r < bs; ++r)
    dst[(r - 2) * stride] = AVG3(left[r - 3], left[r - 2], left[r - 1]);

  // the rest of the block
  for (r = 2; r < bs; ++r) {
    for (c = 1; c < bs; c++) dst[c] = dst[-2 * stride + c - 1];
    dst += stride;
  }
}

static INLINE void d135_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  int i;
#if CONFIG_TX64X64
#if defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ > 7
  // silence a spurious -Warray-bounds warning, possibly related to:
  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56273
  uint8_t border[133];
#else
  uint8_t border[64 + 64 - 1];  // outer border from bottom-left to top-right
#endif
#else
#if defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ > 7
  // silence a spurious -Warray-bounds warning, possibly related to:
  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56273
  uint8_t border[69];
#else
  uint8_t border[32 + 32 - 1];  // outer border from bottom-left to top-right
#endif
#endif  // CONFIG_TX64X64

  // dst(bs, bs - 2)[0], i.e., border starting at bottom-left
  for (i = 0; i < bs - 2; ++i) {
    border[i] = AVG3(left[bs - 3 - i], left[bs - 2 - i], left[bs - 1 - i]);
  }
  border[bs - 2] = AVG3(above[-1], left[0], left[1]);
  border[bs - 1] = AVG3(left[0], above[-1], above[0]);
  border[bs - 0] = AVG3(above[-1], above[0], above[1]);
  // dst[0][2, size), i.e., remaining top border ascending
  for (i = 0; i < bs - 2; ++i) {
    border[bs + 1 + i] = AVG3(above[i], above[i + 1], above[i + 2]);
  }

  for (i = 0; i < bs; ++i) {
    memcpy(dst + i * stride, border + bs - 1 - i, bs);
  }
}

static INLINE void d153_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                  const uint8_t *above, const uint8_t *left) {
  int r, c;
  dst[0] = AVG2(above[-1], left[0]);
  for (r = 1; r < bs; r++) dst[r * stride] = AVG2(left[r - 1], left[r]);
  dst++;

  dst[0] = AVG3(left[0], above[-1], above[0]);
  dst[stride] = AVG3(above[-1], left[0], left[1]);
  for (r = 2; r < bs; r++)
    dst[r * stride] = AVG3(left[r - 2], left[r - 1], left[r]);
  dst++;

  for (c = 0; c < bs - 2; c++)
    dst[c] = AVG3(above[c - 1], above[c], above[c + 1]);
  dst += stride;

  for (r = 1; r < bs; ++r) {
    for (c = 0; c < bs - 2; c++) dst[c] = dst[-stride + c - 2];
    dst += stride;
  }
}

static INLINE void v_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                               const uint8_t *above, const uint8_t *left) {
  int r;
  (void)left;

  for (r = 0; r < bs; r++) {
    memcpy(dst, above, bs);
    dst += stride;
  }
}

static INLINE void h_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                               const uint8_t *above, const uint8_t *left) {
  int r;
  (void)above;

  for (r = 0; r < bs; r++) {
    memset(dst, left[r], bs);
    dst += stride;
  }
}

#if CONFIG_ALT_INTRA
static INLINE int abs_diff(int a, int b) { return (a > b) ? a - b : b - a; }

static INLINE uint16_t paeth_predictor_single(uint16_t left, uint16_t top,
                                              uint16_t top_left) {
  const int base = top + left - top_left;
  const int p_left = abs_diff(base, left);
  const int p_top = abs_diff(base, top);
  const int p_top_left = abs_diff(base, top_left);

  // Return nearest to base of left, top and top_left.
  return (p_left <= p_top && p_left <= p_top_left)
             ? left
             : (p_top <= p_top_left) ? top : top_left;
}

static INLINE void paeth_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                   const uint8_t *above, const uint8_t *left) {
  int r, c;
  const uint8_t ytop_left = above[-1];

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++)
      dst[c] = (uint8_t)paeth_predictor_single(left[r], above[c], ytop_left);
    dst += stride;
  }
}

// Weights are quadratic from 'bs' to '1'.
// Scale is same as 'bs'.
// TODO(urvang): Integerize the weights at a suitable precision.
#if CONFIG_TX64X64
static const double sm_weights_fwd[6][64] = {
#else
static const double sm_weights_fwd[5][32] = {
#endif  // CONFIG_TX64X64
  // bs = 2
  { 2, 1 },
  // bs = 4
  { 4, 2.33333, 1.33333, 1 },
  // bs = 8
  { 8, 6.14286, 4.57143, 3.28571, 2.28571, 1.57143, 1.14286, 1 },
  // bs = 16
  { 16, 14.0667, 12.2667, 10.6, 9.06667, 7.66667, 6.4, 5.26667, 4.26667, 3.4,
    2.66667, 2.06667, 1.6, 1.26667, 1.06667, 1 },
  // bs = 32
  { 32,      30.0323, 28.129,  26.2903, 24.5161, 22.8065, 21.1613, 19.5806,
    18.0645, 16.6129, 15.2258, 13.9032, 12.6452, 11.4516, 10.3226, 9.25806,
    8.25806, 7.32258, 6.45161, 5.64516, 4.90323, 4.22581, 3.6129,  3.06452,
    2.58065, 2.16129, 1.80645, 1.51613, 1.29032, 1.12903, 1.03226, 1 },
#if CONFIG_TX64X64
  // bs = 64
  { 64,      62.0159, 60.0635, 58.1429, 56.254,  54.3968, 52.5714, 50.7778,
    49.0159, 47.2857, 45.5873, 43.9206, 42.2857, 40.6825, 39.1111, 37.5714,
    36.0635, 34.5873, 33.1429, 31.7302, 30.3492, 29,      27.6825, 26.3968,
    25.1429, 23.9206, 22.7302, 21.5714, 20.4444, 19.3492, 18.2857, 17.254,
    16.254,  15.2857, 14.3492, 13.4444, 12.5714, 11.7302, 10.9206, 10.1429,
    9.39683, 8.68254, 8,       7.34921, 6.73016, 6.14286, 5.5873,  5.06349,
    4.57143, 4.11111, 3.68254, 3.28571, 2.92063, 2.5873,  2.28571, 2.01587,
    1.77778, 1.57143, 1.39683, 1.25397, 1.14286, 1.06349, 1.01587, 1 },
#endif  // CONFIG_TX64X64
};

#if CONFIG_TX64X64
static const double sm_weights_rev[6][64] = {
#else
static const double sm_weights_rev[5][32] = {
#endif  // CONFIG_TX64X64
  // bs = 2
  { 0, 1 },
  // bs = 4
  { 0, 1.66667, 2.66667, 3 },
  // bs = 8
  { 0, 1.85714, 3.42857, 4.71429, 5.71429, 6.42857, 6.85714, 7 },
  // bs = 16
  { 0, 1.93333, 3.73333, 5.4, 6.93333, 8.33333, 9.6, 10.7333, 11.7333, 12.6,
    13.3333, 13.9333, 14.4, 14.7333, 14.9333, 15 },
  // bs = 32
  { 0,       1.96774, 3.87097, 5.70968, 7.48387, 9.19355, 10.8387, 12.4194,
    13.9355, 15.3871, 16.7742, 18.0968, 19.3548, 20.5484, 21.6774, 22.7419,
    23.7419, 24.6774, 25.5484, 26.3548, 27.0968, 27.7742, 28.3871, 28.9355,
    29.4194, 29.8387, 30.1935, 30.4839, 30.7097, 30.871,  30.9677, 31 },
#if CONFIG_TX64X64
  // bs = 64
  { 0,       1.98413, 3.93651, 5.85714, 7.74603, 9.60317, 11.4286, 13.2222,
    14.9841, 16.7143, 18.4127, 20.0794, 21.7143, 23.3175, 24.8889, 26.4286,
    27.9365, 29.4127, 30.8571, 32.2698, 33.6508, 35,      36.3175, 37.6032,
    38.8571, 40.0794, 41.2698, 42.4286, 43.5556, 44.6508, 45.7143, 46.746,
    47.746,  48.7143, 49.6508, 50.5556, 51.4286, 52.2698, 53.0794, 53.8571,
    54.6032, 55.3175, 56,      56.6508, 57.2698, 57.8571, 58.4127, 58.9365,
    59.4286, 59.8889, 60.3175, 60.7143, 61.0794, 61.4127, 61.7143, 61.9841,
    62.2222, 62.4286, 62.6032, 62.746,  62.8571, 62.9365, 62.9841, 63 },
#endif  // CONFIG_TX64X64
};

static INLINE void smooth_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  const uint8_t below_pred = left[bs - 1];   // estimated by bottom-left pixel
  const uint8_t right_pred = above[bs - 1];  // estimated by top-right pixel
  const int arr_index = (int)lround(log2(bs)) - 1;
  const double *const fwd_weights = sm_weights_fwd[arr_index];
  const double *const rev_weights = sm_weights_rev[arr_index];
  const double scale = 2.0 * bs;
  int r;
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c) {
      const int pixels[] = { above[c], below_pred, left[r], right_pred };
      const double weights[] = { fwd_weights[r], rev_weights[r], fwd_weights[c],
                                 rev_weights[c] };
      double this_pred = 0;
      int i;
      for (i = 0; i < 4; ++i) {
        this_pred += weights[i] * pixels[i];
      }
      dst[c] = clip_pixel(lround(this_pred / scale));
    }
    dst += stride;
  }
}

static INLINE void least_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  uint8_t *temp;
  int r;
  temp = (uint8_t*) malloc(bs);
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c)
    	temp[c] = (above[c]*((r+1)^3)+left[r]*((c+1)^3))/(((r+1)^3)+((c+1)^3));
    memcpy(dst, temp, bs);
    dst += stride;
  }//for
}

static INLINE void wcalic_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  uint8_t *temp;
  int r;
  temp = (uint8_t*) malloc(bs);
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c)
    	temp[c] = (above[c]*((r+1)^3)+left[r]*((c+1)^3))/(((r+1)^3)+((c+1)^3));
    memcpy(dst, temp, bs);
    dst += stride;
  }//for
}

static INLINE void isle_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  uint8_t *temp;
  int r;
  temp = (uint8_t*) malloc(bs);
  switch (bs){
	  case 4:
		temp[0] = (int)(0.149806*above[-1]+0.176843*above[0]+0.188342*above[1]+0.0459267*above[2]+0.0188894*above[3]+0.176007*left[0]+0.187323*left[1]+0.0447558*left[2]+0.018554*left[3]);
		temp[1] = (int)(0.0365993*above[-1]+0.200869*above[0]+0.23798*above[1]+0.223417*above[2]+0.0591478*above[3]+0.0710809*left[0]+0.0902342*left[1]+0.0673923*left[2]+0.0329107*left[3]);
		temp[2] = (int)(0.0227616*above[-1]+0.0683466*above[0]+0.243396*above[1]+0.277401*above[2]+0.231816*above[3]+0.0451048*left[0]+0.0656805*left[1]+0.059505*left[2]+0.0371618*left[3]);
		temp[3] = (int)(0.0418371*above[-1]+0.0775787*above[0]+0.149306*above[1]+0.378025*above[2]+0.342284*above[3]+0.0650567*left[0]+0.0866076*left[1]+0.0640531*left[2]+0.0408335*left[3]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0353905*above[-1]+0.0741366*above[0]+0.092907*above[1]+0.0710178*above[2]+0.0322717*above[3]+0.197696*left[0]+0.233239*left[1]+0.218764*left[2]+0.0564578*left[3]);
		temp[1] = (int)(0.046192*above[-1]+0.102838*above[0]+0.157393*above[1]+0.146637*above[2]+0.089991*above[3]+0.101043*left[0]+0.15502*left[1]+0.143672*left[2]+0.0888212*left[3]);
		temp[2] = (int)(0.0345212*above[-1]+0.0967754*above[0]+0.172136*above[1]+0.222433*above[2]+0.160178*above[3]+0.0739959*left[0]+0.115602*left[1]+0.120132*left[2]+0.0806576*left[3]); 
		temp[3] = (int)(0.0329124*above[-1]+0.0809815*above[0]+0.173462*above[1]+0.25446*above[2]+0.206391*above[3]+0.0668332*left[0]+0.104711*left[1]+0.108575*left[2]+0.0746541*left[3]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0181858*above[-1]+0.0408267*above[0]+0.0609626*above[1]+0.0590855*above[2]+0.0364446*above[3]+0.0618371*left[0]+0.235456*left[1]+0.272629*left[2]+0.228978*left[3]);
		temp[1] = (int)(0.0293539*above[-1]+0.0695521*above[0]+0.110938*above[1]+0.120645*above[2]+0.0804468*above[3]+0.088949*left[0]+0.162465*left[1]+0.216047*left[2]+0.156452*left[3]);
		temp[2] = (int)(0.0327864*above[-1]+0.0804687*above[0]+0.140936*above[1]+0.171405*above[2]+0.123723*above[3]+0.0788817*left[0]+0.138539*left[1]+0.168051*left[2]+0.121955*left[3]);
		temp[3] = (int)(0.0324419*above[-1]+0.0823681*above[0]+0.151597*above[1]+0.196851*above[2]+0.146925*above[3]+0.0747049*left[0]+0.127471*left[1]+0.149932*left[2]+0.107669*left[3]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0133633*above[-1]+0.0312002*above[0]+0.0491957*above[1]+0.052381*above[2]+0.0345441*above[3]+0.0414506*left[0]+0.106618*left[1]+0.357037*left[2]+0.32895*left[3]);
		temp[1] = (int)(0.0242045*above[-1]+0.0576405*above[0]+0.0948127*above[1]+0.107002*above[2]+0.0735662*above[3]+0.0689807*left[0]+0.158961*left[1]+0.246071*left[2]+0.201295*left[3]);
		temp[2] = (int)(0.0299398*above[-1]+0.0732499*above[0]+0.126379*above[1]+0.151905*above[2]+0.108595*above[3]+0.0777886*left[0]+0.145704*left[1]+0.192022*left[2]+0.144173*left[3]);
		temp[3] = (int)(0.0317227*above[-1]+0.0786956*above[0]+0.139638*above[1]+0.173387*above[2]+0.126414*above[3]+0.077125*left[0]+0.137238*left[1]+0.170002*left[2]+0.124599*left[3]);
		memcpy(dst, temp, bs);
		dst += stride;
		break;
	  case 8:
	    	temp[0] = (int)(0.14833*above[-1]+0.173909*above[0]+0.183211*above[1]+0.0388265*above[2]+0.0152115*above[3]+0.00701522*above[4]+0.0037912*above[5]+0.00240563*above[6]+0.0012997*above[7]+0.173909*left[0]+0.183211*left[1]+0.0388265*left[2]+0.0152115*left[3]+0.00701522*left[4]+0.0037912*left[5]+0.0012997*left[7]+0.00240563*left[6]);
		temp[1] = (int)(0.0324404*above[-1]+0.192703*above[0]+0.22334*above[1]+0.202718*above[2]+0.0478104*above[3]+0.0200378*above[4]+0.0100258*above[5]+0.00610205*above[6]+0.00323684*above[7]+0.068135*left[0]+0.0825451*left[1]+0.0571483*left[2]+0.0252024*left[3]+0.0130013*left[4]+0.00743934*left[5]+0.00268731*left[7]+0.00489635*left[6]);
		temp[2] = (int)(0.0138637*above[-1]+0.0524606*above[0]+0.215812*above[1]+0.23432*above[2]+0.20864*above[3]+0.0514628*above[4]+0.0227482*above[5]+0.0126481*above[6]+0.0064743*above[7]+0.0317157*left[0]+0.0456611*left[1]+0.0398651*left[2]+0.0268342*left[3]+0.0159375*left[4]+0.0100095*left[5]+0.00391516*left[7]+0.00696385*left[6]);
		temp[3] = (int)(0.00707833*above[-1]+0.0246328*above[0]+0.0653219*above[1]+0.22292*above[2]+0.23871*above[3]+0.211817*above[4]+0.0543876*above[5]+0.0263705*above[6]+0.0125742*above[7]+0.0175328*left[0]+0.0273918*left[1]+0.0277564*left[2]+0.0223789*left[3]+0.0160289*left[4]+0.0111902*left[5]+0.00481798*left[7]+0.0083269*left[6]);
		temp[4] = (int)(0.004214*above[-1]+0.0137676*above[0]+0.0329025*above[1]+0.070531*above[2]+0.226671*above[3]+0.242021*above[4]+0.215706*above[5]+0.060428*above[6]+0.0259426*above[7]+0.0108895*left[0]+0.0179754*left[1]+0.0199078*left[2]+0.0180316*left[3]+0.0145758*left[4]+0.0112867*left[5]+0.00534808*left[7]+0.00897815*left[6]);
		temp[5] = (int)(0.00281734*above[-1]+0.00888846*above[0]+0.0198025*above[1]+0.0373051*above[2]+0.0743666*above[3]+0.230933*above[4]+0.248369*above[5]+0.227819*above[6]+0.0603386*above[7]+0.0074889*left[0]+0.0128197*left[1]+0.0150649*left[2]+0.0147268*left[3]+0.0129388*left[4]+0.0108144*left[5]+0.00555796*left[7]+0.00910077*left[6]);
		temp[6] = (int)(0.00210921*above[-1]+0.00651279*above[0]+0.0139624*above[1]+0.024243*above[2]+0.0421435*above[3]+0.0811453*above[4]+0.243584*above[5]+0.27762*above[6]+0.231169*above[7]+0.00570671*left[0]+0.0100006*left[1]+0.0122043*left[2]+0.0125227*left[3]+0.0116027*left[4]+0.0101963*left[5]+0.00554138*left[7]+0.00891529*left[6]);
		temp[7] = (int)(0.00178045*above[-1]+0.00544615*above[0]+0.0114824*above[1]+0.0192911*above[2]+0.0314956*above[3]+0.0547233*above[4]+0.111792*above[5]+0.361636*above[6]+0.332372*above[7]+0.00485575*left[0]+0.00860028*left[1]+0.0106755*left[2]+0.0111946*left[3]+0.010622*left[4]+0.00954774*left[5]+0.00532326*left[7]+0.0084952*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0324404*above[-1]+0.068135*above[0]+0.0825451*above[1]+0.0571483*above[2]+0.0252024*above[3]+0.0130013*above[4]+0.00743934*above[5]+0.00489635*above[6]+0.00268731*above[7]+0.192703*left[0]+0.22334*left[1]+0.202718*left[2]+0.0478104*left[3]+0.0200378*left[4]+0.0100258*left[5]+0.00323684*left[7]+0.00610205*left[6]);
		temp[1] = (int)(0.039125*above[-1]+0.0888304*above[0]+0.133119*above[1]+0.113016*above[2]+0.0730757*above[3]+0.034334*above[4]+0.018945*above[5]+0.0121049*above[6]+0.006558*above[7]+0.0888304*left[0]+0.133119*left[1]+0.113016*left[2]+0.0730757*left[3]+0.034334*left[4]+0.018945*left[5]+0.006558*left[7]+0.0121049*left[6]);
		temp[2] = (int)(0.0198519*above[-1]+0.0692403*above[0]+0.124529*above[1]+0.152273*above[2]+0.12408*above[3]+0.0802103*above[4]+0.0397624*above[5]+0.0242077*above[6]+0.0127889*above[7]+0.0496902*left[0]+0.0771072*left[1]+0.0769953*left[2]+0.0589075*left[3]+0.0388662*left[4]+0.0242461*left[5]+0.00934474*left[7]+0.0167204*left[6]);
		temp[3] = (int)(0.0117099*above[-1]+0.0381781*above[0]+0.0915304*above[1]+0.137714*above[2]+0.160797*above[3]+0.130413*above[4]+0.0860637*above[5]+0.0468496*above[6]+0.0238807*above[7]+0.0300926*left[0]+0.0492325*left[1]+0.0533991*left[2]+0.0466457*left[3]+0.0356584*left[4]+0.0258523*left[5]+0.0112292*left[7]+0.0193819*left[6]);
		temp[4] = (int)(0.00753287*above[-1]+0.0239139*above[0]+0.053479*above[1]+0.10161*above[2]+0.145171*above[3]+0.167418*above[4]+0.138014*above[5]+0.0973115*above[6]+0.0454999*above[7]+0.0199016*left[0]+0.0337392*left[1]+0.0389439*left[2]+0.0370218*left[3]+0.0313349*left[4]+0.0250944*left[5]+0.0121837*left[7]+0.0203345*left[6]);
		temp[5] = (int)(0.00530997*above[-1]+0.0164755*above[0]+0.0355813*above[1]+0.0622246*above[2]+0.109283*above[3]+0.153529*above[4]+0.179221*above[5]+0.158522*above[6]+0.09517*above[7]+0.0142998*left[0]+0.0248898*left[1]+0.0300147*left[2]+0.0302755*left[3]+0.027455*left[4]+0.0235716*left[5]+0.0124297*left[7]+0.0201992*left[6]);
		temp[6] = (int)(0.00412606*above[-1]+0.0126157*above[0]+0.0265669*above[1]+0.0444868*above[2]+0.0717503*above[3]+0.121965*above[4]+0.174871*above[5]+0.222967*above[6]+0.158802*above[7]+0.0112553*left[0]+0.0199374*left[1]+0.024747*left[2]+0.025936*left[3]+0.0245781*left[4]+0.0220498*left[5]+0.0122555*left[7]+0.0195796*left[6]);
		temp[7] = (int)(0.00361763*above[-1]+0.0109922*above[0]+0.0229077*above[1]+0.0376167*above[2]+0.0589794*above[3]+0.0949259*above[4]+0.167501*above[5]+0.250914*above[6]+0.203052*above[7]+0.00992278*left[0]+0.0177104*left[1]+0.0222563*left[2]+0.0237086*left[3]+0.022889*left[4]+0.0209182*left[5]+0.0118866*left[7]+0.0188547*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0138637*above[-1]+0.0317157*above[0]+0.0456611*above[1]+0.0398651*above[2]+0.0268342*above[3]+0.0159375*above[4]+0.0100095*above[5]+0.00696385*above[6]+0.00391516*above[7]+0.0524606*left[0]+0.215812*left[1]+0.23432*left[2]+0.20864*left[3]+0.0514628*left[4]+0.0227482*left[5]+0.0064743*left[7]+0.0126481*left[6]);
		temp[1] = (int)(0.0198519*above[-1]+0.0496902*above[0]+0.0771072*above[1]+0.0769953*above[2]+0.0589075*above[3]+0.0388662*above[4]+0.0242461*above[5]+0.0167204*above[6]+0.00934474*above[7]+0.0692403*left[0]+0.124529*left[1]+0.152273*left[2]+0.12408*left[3]+0.0802103*left[4]+0.0397624*left[5]+0.0127889*left[7]+0.0242077*left[6]);
		temp[2] = (int)(0.0178816*above[-1]+0.0489296*above[0]+0.0856684*above[1]+0.0996215*above[2]+0.091367*above[3]+0.0687655*above[4]+0.0466089*above[5]+0.0317267*above[6]+0.0175895*above[7]+0.0489296*left[0]+0.0856684*left[1]+0.0996215*left[2]+0.091367*left[3]+0.0687655*left[4]+0.0466089*left[5]+0.0175895*left[7]+0.0317267*left[6]);
		temp[3] = (int)(0.0126896*above[-1]+0.0386718*above[0]+0.0749722*above[1]+0.102659*above[2]+0.111318*above[3]+0.100357*above[4]+0.0770646*above[5]+0.0562627*above[6]+0.0309235*above[7]+0.0344499*left[0]+0.059922*left[1]+0.0713459*left[2]+0.0691087*left[3]+0.0582008*left[4]+0.0450512*left[5]+0.0202976*left[7]+0.0348618*left[6]);
		temp[4] = (int)(0.00919032*above[-1]+0.0280397*above[0]+0.0583276*above[1]+0.0887475*above[2]+0.113211*above[3]+0.120703*above[4]+0.110702*above[5]+0.0911764*above[6]+0.053396*above[7]+0.0250085*left[0]+0.0440046*left[1]+0.053734*left[2]+0.054607*left[3]+0.0493315*left[4]+0.0416703*left[5]+0.0211999*left[7]+0.0349204*left[6]);
		temp[5] = (int)(0.00699176*above[-1]+0.0212213*above[0]+0.043966*above[1]+0.0706853*above[2]+0.0996234*above[3]+0.124625*above[4]+0.135479*above[5]+0.133242*above[6]+0.085496*above[7]+0.0191712*left[0]+0.0341529*left[1]+0.0426975*left[2]+0.0450014*left[3]+0.0427045*left[4]+0.0381682*left[5]+0.0209851*left[7]+0.0336699*left[6]);
		temp[6] = (int)(0.0057245*above[-1]+0.017267*above[0]+0.0355041*above[1]+0.0566357*above[2]+0.0837809*above[3]+0.115701*above[4]+0.14794*above[5]+0.170061*above[6]+0.119904*above[7]+0.0157962*left[0]+0.0284126*left[1]+0.0361293*left[2]+0.03903*left[3]+0.0382131*left[4]+0.0353455*left[5]+0.0203289*left[7]+0.0321283*left[6]);
		temp[7] = (int)(0.00518237*above[-1]+0.0155903*above[0]+0.0319378*above[1]+0.0507298*above[2]+0.0750059*above[3]+0.10902*above[4]+0.151161*above[5]+0.190086*above[6]+0.139725*above[7]+0.0143384*left[0]+0.0258933*left[1]+0.0331554*left[2]+0.036175*left[3]+0.0358599*left[4]+0.0336166*left[5]+0.0196723*left[7]+0.0309121*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00707833*above[-1]+0.0175328*above[0]+0.0273918*above[1]+0.0277564*above[2]+0.0223789*above[3]+0.0160289*above[4]+0.0111902*above[5]+0.0083269*above[6]+0.00481798*above[7]+0.0246328*left[0]+0.0653219*left[1]+0.22292*left[2]+0.23871*left[3]+0.211817*left[4]+0.0543876*left[5]+0.0125742*left[7]+0.0263705*left[6]);
		temp[1] = (int)(0.0117099*above[-1]+0.0300926*above[0]+0.0492325*above[1]+0.0533991*above[2]+0.0466457*above[3]+0.0356584*above[4]+0.0258523*above[5]+0.0193819*above[6]+0.0112292*above[7]+0.0381781*left[0]+0.0915304*left[1]+0.137714*left[2]+0.160797*left[3]+0.130413*left[4]+0.0860637*left[5]+0.0238807*left[7]+0.0468496*left[6]);
		temp[2] = (int)(0.0126896*above[-1]+0.0344499*above[0]+0.059922*above[1]+0.0713459*above[2]+0.0691087*above[3]+0.0582008*above[4]+0.0450512*above[5]+0.0348618*above[6]+0.0202976*above[7]+0.0386718*left[0]+0.0749722*left[1]+0.102659*left[2]+0.111318*left[3]+0.100357*left[4]+0.0770646*left[5]+0.0309235*left[7]+0.0562627*left[6]);
		temp[3] = (int)(0.0113304*above[-1]+0.0322829*above[0]+0.0599166*above[1]+0.0783775*above[2]+0.0849751*above[3]+0.0799697*above[4]+0.0681784*above[5]+0.056082*above[6]+0.0334537*above[7]+0.0322829*left[0]+0.0599166*left[1]+0.0783775*left[2]+0.0849751*left[3]+0.0799697*left[4]+0.0681784*left[5]+0.0334537*left[7]+0.056082*left[6]);
		temp[4] = (int)(0.00935239*above[-1]+0.0275423*above[0]+0.0535881*above[1]+0.0759109*above[2]+0.0910895*above[3]+0.0962589*above[4]+0.0917877*above[5]+0.082812*above[6]+0.0515967*above[7]+0.026244*left[0]+0.0480411*left[1]+0.0623507*left[2]+0.068202*left[3]+0.0665294*left[4]+0.060207*left[5]+0.0329615*left[7]+0.053086*left[6]);
		temp[5] = (int)(0.00777657*above[-1]+0.0230769*above[0]+0.0460151*above[1]+0.0684525*above[2]+0.088983*above[3]+0.104171*above[4]+0.111538*above[5]+0.111534*above[6]+0.073408*above[7]+0.021738*left[0]+0.0397344*left[1]+0.0517134*left[2]+0.0572968*left[3]+0.0573582*left[4]+0.0538939*left[5]+0.0313424*left[7]+0.0494002*left[6]);
		temp[6] = (int)(0.00675256*above[-1]+0.02007*above[0]+0.0402461*above[1]+0.0612318*above[2]+0.0834567*above[3]+0.105625*above[4]+0.124406*above[5]+0.136016*above[6]+0.0936722*above[7]+0.018876*left[0]+0.0345543*left[1]+0.0451835*left[2]+0.0505697*left[3]+0.051469*left[4]+0.0494078*left[5]+0.0296582*left[7]+0.0462358*left[6]);
		temp[7] = (int)(0.00630253*above[-1]+0.0187411*above[0]+0.0376507*above[1]+0.0576339*above[2]+0.0801525*above[3]+0.105254*above[4]+0.130255*above[5]+0.148799*above[6]+0.104788*above[7]+0.0176184*left[0]+0.0322692*left[1]+0.0422647*left[2]+0.0474693*left[3]+0.0485913*left[4]+0.0469903*left[5]+0.0285108*left[7]+0.0442836*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.004214*above[-1]+0.0108895*above[0]+0.0179754*above[1]+0.0199078*above[2]+0.0180316*above[3]+0.0145758*above[4]+0.0112867*above[5]+0.00897815*above[6]+0.00534808*above[7]+0.0137676*left[0]+0.0329025*left[1]+0.070531*left[2]+0.226671*left[3]+0.242021*left[4]+0.215706*left[5]+0.0259426*left[7]+0.060428*left[6]);
		temp[1] = (int)(0.00753287*above[-1]+0.0199016*above[0]+0.0337392*above[1]+0.0389439*above[2]+0.0370218*above[3]+0.0313349*above[4]+0.0250944*above[5]+0.0203345*above[6]+0.0121837*above[7]+0.0239139*left[0]+0.053479*left[1]+0.10161*left[2]+0.145171*left[3]+0.167418*left[4]+0.138014*left[5]+0.0454999*left[7]+0.0973115*left[6]);
		temp[2] = (int)(0.00919032*above[-1]+0.0250085*above[0]+0.0440046*above[1]+0.053734*above[2]+0.054607*above[3]+0.0493315*above[4]+0.0416703*above[5]+0.0349204*above[6]+0.0211999*above[7]+0.0280397*left[0]+0.0583276*left[1]+0.0887475*left[2]+0.113211*left[3]+0.120703*left[4]+0.110702*left[5]+0.053396*left[7]+0.0911764*left[6]);
		temp[3] = (int)(0.00935239*above[-1]+0.026244*above[0]+0.0480411*above[1]+0.0623507*above[2]+0.068202*above[3]+0.0665294*above[4]+0.060207*above[5]+0.053086*above[6]+0.0329615*above[7]+0.0275423*left[0]+0.0535881*left[1]+0.0759109*left[2]+0.0910895*left[3]+0.0962589*left[4]+0.0917877*left[5]+0.0515967*left[7]+0.082812*left[6]);
		temp[4] = (int)(0.00871133*above[-1]+0.0250295*above[0]+0.0474089*above[1]+0.0649814*above[2]+0.0762769*above[3]+0.0805064*above[4]+0.0787616*above[5]+0.073962*above[6]+0.0473531*above[7]+0.0250295*left[0]+0.0474089*left[1]+0.0649814*left[2]+0.0762769*left[3]+0.0805064*left[4]+0.0787616*left[5]+0.0473531*left[7]+0.073962*left[6]);
		temp[5] = (int)(0.00789288*above[-1]+0.0230012*above[0]+0.0445815*above[1]+0.063667*above[2]+0.0792225*above[3]+0.0898359*above[4]+0.0948226*above[5]+0.0950102*above[6]+0.0628165*above[7]+0.0224426*left[0]+0.0420122*left[1]+0.0568149*left[2]+0.0661339*left[3]+0.0699664*left[4]+0.0694026*left[5]+0.0430403*left[7]+0.0664424*left[6]);
		temp[6] = (int)(0.0072576*above[-1]+0.0212807*above[0]+0.0417593*above[1]+0.0611318*above[2]+0.0791992*above[3]+0.0948036*above[4]+0.106303*above[5]+0.11218*above[6]+0.076097*above[7]+0.0205419*left[0]+0.0382636*left[1]+0.051453*left[2]+0.0597011*left[3]+0.0632608*left[4]+0.0631627*left[5]+0.0396999*left[7]+0.0609813*left[6]);
		temp[7] = (int)(0.00696276*above[-1]+0.020467*above[0]+0.0403546*above[1]+0.0597063*above[2]+0.0787947*above[3]+0.0968315*above[4]+0.111833*above[5]+0.121017*above[6]+0.0830924*above[7]+0.0196701*left[0]+0.0365597*left[1]+0.0490274*left[2]+0.0567628*left[3]+0.0601042*left[4]+0.0600662*left[5]+0.0378631*left[7]+0.0580941*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00281734*above[-1]+0.0074889*above[0]+0.0128197*above[1]+0.0150649*above[2]+0.0147268*above[3]+0.0129388*above[4]+0.0108144*above[5]+0.00910077*above[6]+0.00555796*above[7]+0.00888846*left[0]+0.0198025*left[1]+0.0373051*left[2]+0.0743666*left[3]+0.230933*left[4]+0.248369*left[5]+0.0603386*left[7]+0.227819*left[6]);
		temp[1] = (int)(0.00530997*above[-1]+0.0142998*above[0]+0.0248898*above[1]+0.0300147*above[2]+0.0302755*above[3]+0.027455*above[4]+0.0235716*above[5]+0.0201992*above[6]+0.0124297*above[7]+0.0164755*left[0]+0.0355813*left[1]+0.0622246*left[2]+0.109283*left[3]+0.153529*left[4]+0.179221*left[5]+0.09517*left[7]+0.158522*left[6]);
		temp[2] = (int)(0.00699176*above[-1]+0.0191712*above[0]+0.0341529*above[1]+0.0426975*above[2]+0.0450014*above[3]+0.0427045*above[4]+0.0381682*above[5]+0.0336699*above[6]+0.0209851*above[7]+0.0212213*left[0]+0.043966*left[1]+0.0706853*left[2]+0.0996234*left[3]+0.124625*left[4]+0.135479*left[5]+0.085496*left[7]+0.133242*left[6]);
		temp[3] = (int)(0.00777657*above[-1]+0.021738*above[0]+0.0397344*above[1]+0.0517134*above[2]+0.0572968*above[3]+0.0573582*above[4]+0.0538939*above[5]+0.0494002*above[6]+0.0313424*above[7]+0.0230769*left[0]+0.0460151*left[1]+0.0684525*left[2]+0.088983*left[3]+0.104171*left[4]+0.111538*left[5]+0.073408*left[7]+0.111534*left[6]);
		temp[4] = (int)(0.00789288*above[-1]+0.0224426*above[0]+0.0420122*above[1]+0.0568149*above[2]+0.0661339*above[3]+0.0699664*above[4]+0.0694026*above[5]+0.0664424*above[6]+0.0430403*above[7]+0.0230012*left[0]+0.0445815*left[1]+0.063667*left[2]+0.0792225*left[3]+0.0898359*left[4]+0.0948226*left[5]+0.0628165*left[7]+0.0950102*left[6]);
		temp[5] = (int)(0.00767584*above[-1]+0.0220993*above[0]+0.0421529*above[1]+0.0588286*above[2]+0.0714676*above[3]+0.0794965*above[4]+0.0829862*above[5]+0.0828336*above[6]+0.0547511*above[7]+0.0220993*left[0]+0.0421529*left[1]+0.0588286*left[2]+0.0714676*left[3]+0.0794965*left[4]+0.0829862*left[5]+0.0547511*left[7]+0.0828336*left[6]);
		temp[6] = (int)(0.00740601*above[-1]+0.021485*above[0]+0.0414837*above[1]+0.0591575*above[2]+0.0741187*above[3]+0.0856136*above[4]+0.092938*above[5]+0.095789*above[6]+0.0642954*above[7]+0.02118*left[0]+0.0400462*left[1]+0.0551844*left[2]+0.0661262*left[3]+0.0726835*left[4]+0.0752417*left[5]+0.0493398*left[7]+0.0747632*left[6]);
		temp[7] = (int)(0.00726438*above[-1]+0.0211432*above[0]+0.0410472*above[1]+0.0591184*above[2]+0.0751571*above[3]+0.0883947*above[4]+0.0977736*above[5]+0.102315*above[6]+0.0691739*above[7]+0.0207148*left[0]+0.0390161*left[1]+0.0534576*left[2]+0.0636385*left[3]+0.0695138*left[4]+0.0715984*left[5]+0.0467303*left[7]+0.0709069*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00210921*above[-1]+0.00570671*above[0]+0.0100006*above[1]+0.0122043*above[2]+0.0125227*above[3]+0.0116027*above[4]+0.0101963*above[5]+0.00891529*above[6]+0.00554138*above[7]+0.00651279*left[0]+0.0139624*left[1]+0.024243*left[2]+0.0421435*left[3]+0.0811453*left[4]+0.243584*left[5]+0.231169*left[7]+0.27762*left[6]);
		temp[1] = (int)(0.00412606*above[-1]+0.0112553*above[0]+0.0199374*above[1]+0.024747*above[2]+0.025936*above[3]+0.0245781*above[4]+0.0220498*above[5]+0.0195796*above[6]+0.0122555*above[7]+0.0126157*left[0]+0.0265669*left[1]+0.0444868*left[2]+0.0717503*left[3]+0.121965*left[4]+0.174871*left[5]+0.158802*left[7]+0.222967*left[6]);
		temp[2] = (int)(0.0057245*above[-1]+0.0157962*above[0]+0.0284126*above[1]+0.0361293*above[2]+0.03903*above[3]+0.0382131*above[4]+0.0353455*above[5]+0.0321283*above[6]+0.0203289*above[7]+0.017267*left[0]+0.0355041*left[1]+0.0566357*left[2]+0.0837809*left[3]+0.115701*left[4]+0.14794*left[5]+0.119904*left[7]+0.170061*left[6]);
		temp[3] = (int)(0.00675256*above[-1]+0.018876*above[0]+0.0345543*above[1]+0.0451835*above[2]+0.0505697*above[3]+0.051469*above[4]+0.0494078*above[5]+0.0462358*above[6]+0.0296582*above[7]+0.02007*left[0]+0.0402461*left[1]+0.0612318*left[2]+0.0834567*left[3]+0.105625*left[4]+0.124406*left[5]+0.0936722*left[7]+0.136016*left[6]);
		temp[4] = (int)(0.0072576*above[-1]+0.0205419*above[0]+0.0382636*above[1]+0.051453*above[2]+0.0597011*above[3]+0.0632608*above[4]+0.0631627*above[5]+0.0609813*above[6]+0.0396999*above[7]+0.0212807*left[0]+0.0417593*left[1]+0.0611318*left[2]+0.0791992*left[3]+0.0948036*left[4]+0.106303*left[5]+0.076097*left[7]+0.11218*left[6]);
		temp[5] = (int)(0.00740601*above[-1]+0.02118*above[0]+0.0400462*above[1]+0.0551844*above[2]+0.0661262*above[3]+0.0726835*above[4]+0.0752417*above[5]+0.0747632*above[6]+0.0493398*above[7]+0.021485*left[0]+0.0414837*left[1]+0.0591575*left[2]+0.0741187*left[3]+0.0856136*left[4]+0.092938*left[5]+0.0642954*left[7]+0.095789*left[6]);
		temp[6] = (int)(0.00738871*above[-1]+0.0212829*above[0]+0.0406736*above[1]+0.0570576*above[2]+0.0700236*above[3]+0.0791112*above[4]+0.084151*above[5]+0.0854328*above[6]+0.0569545*above[7]+0.0212829*left[0]+0.0406736*left[1]+0.0570576*left[2]+0.0700236*left[3]+0.0791112*left[4]+0.084151*left[5]+0.0569545*left[7]+0.0854328*left[6]);
		temp[7] = (int)(0.00735585*above[-1]+0.0212585*above[0]+0.0408307*above[1]+0.0577619*above[2]+0.0716962*above[3]+0.0820653*above[4]+0.0884211*above[5]+0.0906761*above[6]+0.0607349*above[7]+0.0211206*left[0]+0.0401785*left[1]+0.0559504*left[2]+0.0680275*left[3]+0.0760837*left[4]+0.0801718*left[5]+0.0537025*left[7]+0.0808167*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00178045*above[-1]+0.00485575*above[0]+0.00860028*above[1]+0.0106755*above[2]+0.0111946*above[3]+0.010622*above[4]+0.00954774*above[5]+0.0084952*above[6]+0.00532326*above[7]+0.00544615*left[0]+0.0114824*left[1]+0.0192911*left[2]+0.0314956*left[3]+0.0547233*left[4]+0.111792*left[5]+0.332372*left[7]+0.361636*left[6]);
		temp[1] = (int)(0.00361763*above[-1]+0.00992278*above[0]+0.0177104*above[1]+0.0222563*above[2]+0.0237086*above[3]+0.022889*above[4]+0.0209182*above[5]+0.0188547*above[6]+0.0118866*above[7]+0.0109922*left[0]+0.0229077*left[1]+0.0376167*left[2]+0.0589794*left[3]+0.0949259*left[4]+0.167501*left[5]+0.203052*left[7]+0.250914*left[6]);
		temp[2] = (int)(0.00518237*above[-1]+0.0143384*above[0]+0.0258933*above[1]+0.0331554*above[2]+0.036175*above[3]+0.0358599*above[4]+0.0336166*above[5]+0.0309121*above[6]+0.0196723*above[7]+0.0155903*left[0]+0.0319378*left[1]+0.0507298*left[2]+0.0750059*left[3]+0.10902*left[4]+0.151161*left[5]+0.139725*left[7]+0.190086*left[6]);
		temp[3] = (int)(0.00630253*above[-1]+0.0176184*above[0]+0.0322692*above[1]+0.0422647*above[2]+0.0474693*above[3]+0.0485913*above[4]+0.0469903*above[5]+0.0442836*above[6]+0.0285108*above[7]+0.0187411*left[0]+0.0376507*left[1]+0.0576339*left[2]+0.0801525*left[3]+0.105254*left[4]+0.130255*left[5]+0.104788*left[7]+0.148799*left[6]);
		temp[4] = (int)(0.00696276*above[-1]+0.0196701*above[0]+0.0365597*above[1]+0.0490274*above[2]+0.0567628*above[3]+0.0601042*above[4]+0.0600662*above[5]+0.0580941*above[6]+0.0378631*above[7]+0.020467*left[0]+0.0403546*left[1]+0.0597063*left[2]+0.0787947*left[3]+0.0968315*left[4]+0.111833*left[5]+0.0830924*left[7]+0.121017*left[6]);
		temp[5] = (int)(0.00726438*above[-1]+0.0207148*above[0]+0.0390161*above[1]+0.0534576*above[2]+0.0636385*above[3]+0.0695138*above[4]+0.0715984*above[5]+0.0709069*above[6]+0.0467303*above[7]+0.0211432*left[0]+0.0410472*left[1]+0.0591184*left[2]+0.0751571*left[3]+0.0883947*left[4]+0.0977736*left[5]+0.0691739*left[7]+0.102315*left[6]);
		temp[6] = (int)(0.00735585*above[-1]+0.0211206*above[0]+0.0401785*above[1]+0.0559504*above[2]+0.0680275*above[3]+0.0760837*above[4]+0.0801718*above[5]+0.0808167*above[6]+0.0537025*above[7]+0.0212585*left[0]+0.0408307*left[1]+0.0577619*left[2]+0.0716962*left[3]+0.0820653*left[4]+0.0884211*left[5]+0.0607349*left[7]+0.0906761*left[6]);
		temp[7] = (int)(0.0073668*above[-1]+0.0212206*above[0]+0.0405609*above[1]+0.0569233*above[2]+0.0699158*above[3]+0.0790867*above[4]+0.084248*above[5]+0.0856419*above[6]+0.0571306*above[7]+0.0212206*left[0]+0.0405609*left[1]+0.0569233*left[2]+0.0699158*left[3]+0.0790867*left[4]+0.084248*left[5]+0.0571306*left[7]+0.0856419*left[6]);
		memcpy(dst, temp, bs);
		dst += stride;
		break;	
	case 16:
		temp[0] = (int)(0.148309*above[-1]+0.173851*above[0]+0.173851*left[0]+0.183099*above[1]+0.183099*left[1]+0.0386645*above[2]+0.0386645*left[2]+0.0149999*above[3]+0.0149999*left[3]+0.00674458*above[4]+0.00674458*left[4]+0.00343045*above[5]+0.00343045*left[5]+0.00188662*above[6]+0.00188662*left[6]+0.00109963*above[7]+0.00109963*left[7]+0.000670179*above[8]+0.000670179*left[8]+0.000423584*above[9]+0.000423584*left[9]+0.000276361*above[10]+0.000276361*left[10]+0.000185898*above[11]+0.000185898*left[11]+0.000129328*above[12]+0.000129328*left[12]+0.0000941577*above[13]+0.0000941577*left[13]+0.0000437494*above[15]+0.0000437494*left[15]+0.0000736128*above[14]+0.0000736128*left[14]);
		temp[1] = (int)(0.0323946*above[-1]+0.19257*above[0]+0.0680031*left[0]+0.223084*above[1]+0.0822919*left[1]+0.202346*above[2]+0.056786*left[2]+0.0473198*above[3]+0.0247336*left[3]+0.0194023*above[4]+0.0124091*left[4]+0.00916511*above[5]+0.00666212*left[5]+0.00484079*above[6]+0.00379783*left[6]+0.00274582*above[7]+0.00226735*left[7]+0.0016424*above[8]+0.00140565*left[8]+0.00102404*above[9]+0.000899686*left[9]+0.000661318*above[10]+0.000592646*left[10]+0.000441333*above[11]+0.000401653*left[11]+0.000305124*above[12]+0.000281093*left[12]+0.000221072*above[13]+0.000205604*left[13]+0.000102212*above[15]+0.0000959861*left[15]+0.000172246*above[14]+0.000161268*left[14]);
		temp[2] = (int)(0.0137858*above[-1]+0.0522344*above[0]+0.0314926*left[0]+0.215371*above[1]+0.0452355*left[1]+0.233674*above[2]+0.039263*left[2]+0.207774*above[3]+0.0260676*left[3]+0.050316*above[4]+0.014989*left[4]+0.0211519*above[5]+0.00879662*left[5]+0.0102323*above[6]+0.00529943*left[6]+0.0055154*above[7]+0.00328629*left[7]+0.0031857*above[8]+0.00209342*left[8]+0.00193777*above[9]+0.00136725*left[9]+0.00122864*above[10]+0.000914717*left[10]+0.000808498*above[11]+0.000627534*left[11]+0.000552867*above[12]+0.000443454*left[12]+0.000397181*above[13]+0.000326834*left[13]+0.000182055*above[15]+0.000153771*left[15]+0.000307621*above[14]+0.000257732*left[14]);
		temp[3] = (int)(0.00695955*above[-1]+0.0242854*above[0]+0.0171945*left[0]+0.0646388*above[1]+0.0267522*left[1]+0.221902*above[2]+0.0268659*left[2]+0.237312*above[3]+0.0212704*left[3]+0.209909*above[4]+0.0146968*left[4]+0.0516223*above[5]+0.0095474*left[5]+0.0219791*above[6]+0.00616285*left[6]+0.0107722*above[7]+0.0040097*left[7]+0.00587801*above[8]+0.0026454*left[8]+0.00343664*above[9]+0.00177408*left[9]+0.00211772*above[10]+0.00121152*left[10]+0.00136394*above[11]+0.000844796*left[11]+0.00091737*above[12]+0.000604816*left[12]+0.000650702*above[13]+0.000450344*left[13]+0.000294432*above[15]+0.000214108*left[15]+0.000499507*above[14]+0.000357699*left[14]);
		temp[4] = (int)(0.00404048*above[-1]+0.013256*above[0]+0.0103991*left[0]+0.0318839*above[1]+0.0170594*left[1]+0.0689792*above[2]+0.018658*left[2]+0.224473*above[3]+0.0165205*left[3]+0.238894*above[4]+0.0128267*left[4]+0.210914*above[5]+0.00922672*left[5]+0.0522794*above[6]+0.00640032*left[6]+0.0224208*above[7]+0.0043913*left[7]+0.0110779*above[8]+0.00301519*left[8]+0.00609679*above[9]+0.00208537*left[9]+0.0036003*above[10]+0.00145916*left[10]+0.00224784*above[11]+0.00103753*left[11]+0.00147668*above[12]+0.000754591*left[12]+0.00102885*above[13]+0.000568887*left[13]+0.000457204*above[15]+0.000273924*left[15]+0.000780012*above[14]+0.00045583*left[14]);
		temp[5] = (int)(0.00256612*above[-1]+0.00814027*above[0]+0.00678583*left[0]+0.0182888*above[1]+0.0115253*left[1]+0.0349352*above[2]+0.0133423*left[2]+0.0708733*above[3]+0.0127167*left[3]+0.225683*above[4]+0.0107162*left[4]+0.239687*above[5]+0.0083377*left[5]+0.211448*above[6]+0.00618347*left[6]+0.0526496*above[7]+0.00447179*left[7]+0.0226856*above[8]+0.00320064*left[8]+0.0112753*above[9]+0.00228803*left[9]+0.0062526*above[10]+0.00164424*left[10]+0.00373308*above[11]+0.00119484*left[11]+0.00237316*above[12]+0.000884539*left[12]+0.00161337*above[13]+0.000676291*left[13]+0.000699541*above[15]+0.000330357*left[15]+0.00120259*above[14]+0.000547293*left[14]);
		temp[6] = (int)(0.00174105*above[-1]+0.0054025*above[0]+0.00468857*left[0]+0.0116713*above[1]+0.00815837*left[1]+0.0205292*above[2]+0.00982467*left[2]+0.0363791*above[3]+0.00986169*left[3]+0.0718257*above[4]+0.00881624*left[4]+0.226326*above[5]+0.0072859*left[5]+0.240134*above[6]+0.00571382*left[6]+0.211768*above[7]+0.00433653*left[7]+0.0528873*above[8]+0.00323119*left[8]+0.0228723*above[9]+0.00238802*left[9]+0.0114328*above[10]+0.0017641*left[10]+0.00639817*above[11]+0.00131166*left[11]+0.00388347*above[12]+0.000989559*left[12]+0.00255284*above[13]+0.000768087*left[13]+0.00107062*above[15]+0.000381043*left[15]+0.00185959*above[14]+0.000628246*left[14]);
		temp[7] = (int)(0.00124197*above[-1]+0.00379662*above[0]+0.00338791*left[0]+0.00799194*above[1]+0.00599897*left[1]+0.0133754*above[2]+0.00743362*left[2]+0.0216631*above[3]+0.00774936*left[3]+0.0371491*above[4]+0.00724038*left[4]+0.0723613*above[5]+0.00626976*left[5]+0.22671*above[6]+0.00514688*left[6]+0.240419*above[7]+0.00407381*left[7]+0.211991*above[8]+0.0031499*left[8]+0.0530745*above[9]+0.0024034*left[9]+0.0230432*above[10]+0.00182449*left[10]+0.0116053*above[11]+0.00138827*left[11]+0.00659256*above[12]+0.0010678*left[12]+0.00413464*above[13]+0.000841778*left[13]+0.00165671*above[15]+0.000424306*left[15]+0.00291839*above[14]+0.000696133*left[14]);
		temp[8] = (int)(0.000921862*above[-1]+0.00278816*above[0]+0.00253851*left[0]+0.00576475*above[1]+0.00455373*left[1]+0.00932795*above[2]+0.00576458*left[2]+0.0142903*above[3]+0.00618339*left[3]+0.0223026*above[4]+0.00597601*left[4]+0.0376084*above[5]+0.0053687*left[5]+0.0727033*above[6]+0.0045747*left[6]+0.226977*above[7]+0.00375293*left[7]+0.240642*above[8]+0.00299904*left[8]+0.212193*above[9]+0.00235683*left[9]+0.0532752*above[10]+0.00183614*left[10]+0.023264*above[11]+0.00142888*left[11]+0.0118739*above[12]+0.00112016*left[12]+0.0069629*above[13]+0.000896763*left[13]+0.0026167*above[15]+0.000459232*left[15]+0.00470101*above[14]+0.000749752*left[14]);
		temp[9] = (int)(0.000707116*above[-1]+0.00212187*above[0]+0.00196106*left[0]+0.00433048*above[1]+0.00355299*left[1]+0.00684161*above[2]+0.00457209*left[2]+0.0100868*above[3]+0.00501376*left[3]+0.0148377*above[4]+0.0049758*left[4]+0.0227109*above[5]+0.00460343*left[5]+0.0379274*above[6]+0.00404449*left[6]+0.0729686*above[7]+0.00342002*left[7]+0.227216*above[8]+0.00281294*left[8]+0.240877*above[9]+0.00227018*left[9]+0.212448*above[10]+0.00181158*left[10]+0.0535783*above[11]+0.00143994*left[11]+0.0236574*above[12]+0.00114953*left[12]+0.0124466*above[13]+0.000933992*left[13]+0.00426178*above[15]+0.00048563*left[15]+0.00787864*above[14]+0.000789127*left[14]);
		temp[10] = (int)(0.000558143*above[-1]+0.00166482*above[0]+0.00155644*left[0]+0.00336477*above[1]+0.00284188*left[1]+0.0052228*above[2]+0.00370435*left[2]+0.0074902*above[3]+0.00413363*left[3]+0.0105724*above[4]+0.00419001*left[4]+0.0152176*above[5]+0.00396973*left[5]+0.0230265*above[6]+0.00357685*left[6]+0.0382103*above[7]+0.00310286*left[7]+0.0732455*above[8]+0.00261644*left[8]+0.227513*above[9]+0.00216187*left[9]+0.241226*above[10]+0.0017629*left[10]+0.21289*above[11]+0.00142872*left[11]+0.0541847*above[12]+0.00115992*left[12]+0.0245849*above[13]+0.000955493*left[13]+0.00725305*above[15]+0.000503896*left[15]+0.0139985*above[14]+0.00081523*left[14]);
		temp[11] = (int)(0.000452457*above[-1]+0.00134331*above[0]+0.00126719*left[0]+0.00269471*above[1]+0.00232802*left[1]+0.00412693*above[2]+0.00306574*left[2]+0.00579754*above[3]+0.00346903*left[3]+0.00794119*above[4]+0.00357676*left[4]+0.0109471*above[5]+0.00345497*left[5]+0.0155527*above[6]+0.00317852*left[6]+0.023353*above[7]+0.00281699*left[7]+0.0385584*above[8]+0.00242635*left[8]+0.0736495*above[9]+0.0020461*left[9]+0.228019*above[10]+0.00170057*left[10]+0.241905*above[11]+0.0014022*left[11]+0.213871*above[12]+0.00115568*left[12]+0.0557626*above[13]+0.000963848*left[13]+0.0131451*above[15]+0.00051481*left[15]+0.0273687*above[14]+0.000829621*left[14]);
		temp[12] = (int)(0.000376794*above[-1]+0.0011146*above[0]+0.00105889*left[0]+0.002223*above[1]+0.00195486*left[1]+0.00336947*above[2]+0.00259535*left[2]+0.00465964*above[3]+0.00296962*left[3]+0.00624092*above[4]+0.00310398*left[4]+0.00833729*above[5]+0.00304554*left[5]+0.0113319*above[6]+0.00284977*left[6]+0.0159608*above[7]+0.00257054*left[7]+0.0238236*above[8]+0.00225352*left[8]+0.0391429*above[9]+0.00193321*left[9]+0.074425*above[10]+0.00163291*left[10]+0.229116*above[11]+0.00136642*left[11]+0.243572*above[12]+0.00114085*left[12]+0.216717*above[13]+0.000961647*left[13]+0.0263407*above[15]+0.000519273*left[15]+0.0611383*above[14]+0.000833999*left[14]);
		temp[13] = (int)(0.00032318*above[-1]+0.000953354*above[0]+0.000910597*left[0]+0.00189306*above[1]+0.00168742*left[1]+0.00284708*above[2]+0.00225439*left[2]+0.00389112*above[3]+0.00260175*left[3]+0.00512565*above[4]+0.0027484*left[4]+0.00669243*above[5]+0.00272962*left[5]+0.00881433*above[6]+0.0025882*left[6]+0.0118793*above[7]+0.0023671*left[7]+0.0166369*above[8]+0.00210428*left[8]+0.0247143*above[9]+0.00182991*left[9]+0.0403893*above[10]+0.00156569*left[10]+0.0762843*above[11]+0.00132572*left[11]+0.232115*above[12]+0.0011184*left[12]+0.249103*above[13]+0.000950785*left[13]+0.0606003*above[15]+0.000517907*left[15]+0.228304*above[14]+0.000829569*left[14]);
		temp[14] = (int)(0.0002867*above[-1]+0.000844128*above[0]+0.000809266*left[0]+0.00167114*above[1]+0.00150355*left[1]+0.00249996*above[2]+0.00201747*left[2]+0.00338943*above[3]+0.00234219*left[3]+0.00441481*above[4]+0.00249239*left[4]+0.00567639*above[5]+0.00249628*left[5]+0.00732166*above[6]+0.00238877*left[6]+0.00958849*above[7]+0.00220579*left[7]+0.0128954*above[8]+0.00197999*left[8]+0.0180521*above[9]+0.00173823*left[9]+0.026808*above[10]+0.00150068*left[10]+0.0437131*above[11]+0.0012812*left[11]+0.0820786*above[12]+0.00108869*left[12]+0.244128*above[13]+0.000931043*left[13]+0.231335*above[15]+0.000510242*left[15]+0.27795*above[14]+0.000815766*left[14]);
		temp[15] = (int)(0.000261947*above[-1]+0.000770563*above[0]+0.000740021*left[0]+0.00152336*above[1]+0.00137657*left[1]+0.00227324*above[2]+0.00185085*left[2]+0.00307059*above[3]+0.00215474*left[3]+0.00397868*above[4]+0.0023008*left[4]+0.00507944*above[5]+0.00231349*left[5]+0.00648914*above[6]+0.0022234*left[6]+0.00838785*above[7]+0.00206236*left[7]+0.0110772*above[8]+0.00185968*left[8]+0.0151041*above[9]+0.0016399*left[9]+0.0215568*above[10]+0.00142179*left[10]+0.0328669*above[11]+0.00121855*left[11]+0.0555194*above[12]+0.00103899*left[12]+0.112236*above[13]+0.000891019*left[13]+0.332491*above[15]+0.000489699*left[15]+0.361886*above[14]+0.000782237*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0323946*above[-1]+0.0680031*above[0]+0.19257*left[0]+0.0822919*above[1]+0.223084*left[1]+0.056786*above[2]+0.202346*left[2]+0.0247336*above[3]+0.0473198*left[3]+0.0124091*above[4]+0.0194023*left[4]+0.00666212*above[5]+0.00916511*left[5]+0.00379783*above[6]+0.00484079*left[6]+0.00226735*above[7]+0.00274582*left[7]+0.00140565*above[8]+0.0016424*left[8]+0.000899686*above[9]+0.00102404*left[9]+0.000592646*above[10]+0.000661318*left[10]+0.000401653*above[11]+0.000441333*left[11]+0.000281093*above[12]+0.000305124*left[12]+0.000205604*above[13]+0.000221072*left[13]+0.0000959861*above[15]+0.000102212*left[15]+0.000161268*above[14]+0.000172246*left[14]);
		temp[1] = (int)(0.0390225*above[-1]+0.0885344*above[0]+0.0885344*left[0]+0.132549*above[1]+0.132549*left[1]+0.112195*above[2]+0.112195*left[2]+0.072003*above[3]+0.072003*left[3]+0.0329626*above[4]+0.0329626*left[4]+0.0171183*above[5]+0.0171183*left[5]+0.00948022*above[6]+0.00948022*left[6]+0.0055475*above[7]+0.0055475*left[7]+0.00338825*above[8]+0.00338825*left[8]+0.0021443*above[9]+0.0021443*left[9]+0.00140016*above[10]+0.00140016*left[10]+0.000942349*above[11]+0.000942349*left[11]+0.000655828*above[12]+0.000655828*left[12]+0.0004776*above[13]+0.0004776*left[13]+0.000221964*above[15]+0.000221964*left[15]+0.000373451*above[14]+0.000373451*left[14]);
		temp[2] = (int)(0.0196795*above[-1]+0.0687404*above[0]+0.0491946*left[0]+0.123559*above[1]+0.0761582*left[1]+0.15086*above[2]+0.0756443*left[2]+0.122206*above[3]+0.0571722*left[3]+0.0777651*above[4]+0.0366953*left[4]+0.0364226*above[5]+0.0214343*left[5]+0.0192711*above[6]+0.0128103*left[6]+0.0108613*above[7]+0.00786347*left[7]+0.00645823*above[8]+0.00496411*left[8]+0.00400526*above[9]+0.00321726*left[9]+0.00257484*above[10]+0.00213856*left[10]+0.00171177*above[11]+0.00145925*left[11]+0.00117969*above[12]+0.00102659*left[12]+0.000852513*above[13]+0.000753884*left[13]+0.000393085*above[15]+0.000353355*left[15]+0.000662989*above[14]+0.00059295*left[14]);
		temp[3] = (int)(0.0114494*above[-1]+0.0374185*above[0]+0.0293487*left[0]+0.0900431*above[1]+0.0478211*left[1]+0.135514*above[2]+0.0514215*left[2]+0.157811*above[3]+0.0441632*left[3]+0.126404*above[4]+0.0326455*left[4]+0.0803863*above[5]+0.0220969*left[5]+0.0381078*above[6]+0.0143877*left[6]+0.020384*above[7]+0.00936813*left[7]+0.0116156*above[8]+0.00616491*left[8]+0.0069842*above[9]+0.00411983*left[9]+0.00438486*above[10]+0.00280315*left[10]+0.00286177*above[11]+0.00194789*left[11]+0.00194363*above[12]+0.00139023*left[12]+0.00138863*above[13]+0.00103244*left[13]+0.000632822*above[15]+0.000489465*left[15]+0.00107124*above[14]+0.000818455*left[14]);
		temp[4] = (int)(0.00715666*above[-1]+0.0228081*above[0]+0.0188354*left[0]+0.0512877*above[1]+0.0317397*left[1]+0.0983017*above[2]+0.036198*left[2]+0.14055*above[3]+0.0336736*left[3]+0.160974*above[4]+0.0274233*left[4]+0.128444*above[5]+0.0204473*left[5]+0.0817366*above[6]+0.0144846*left[6]+0.0390242*above[7]+0.0100255*left[7]+0.021023*above[8]+0.00690626*left[8]+0.012076*above[9]+0.0047788*left[9]+0.00733052*above[10]+0.00334098*left[10]+0.00466141*above[11]+0.00237221*left[11]+0.00310208*above[12]+0.0017225*left[12]+0.0021815*above[13]+0.00129659*left[13]+0.000978168*above[15]+0.0006232*left[15]+0.00166422*above[14]+0.00103765*left[14]);
		temp[5] = (int)(0.00477171*above[-1]+0.0148783*above[0]+0.0127885*left[0]+0.0323691*above[1]+0.0220947*left[1]+0.0572505*above[2]+0.0262678*left[2]+0.102081*above[3]+0.0258629*left[3]+0.143002*above[4]+0.022529*left[4]+0.162602*above[5]+0.0180398*left[5]+0.129551*above[6]+0.013662*left[6]+0.0825093*above[7]+0.0100102*left[7]+0.0395806*above[8]+0.00721624*left[8]+0.0214402*above[9]+0.00517762*left[9]+0.0124067*above[10]+0.00372682*left[10]+0.00761322*above[11]+0.00270943*left[11]+0.00492851*above[12]+0.00200542*left[12]+0.00339306*above[13]+0.00153256*left[13]+0.00148874*above[15]+0.000748077*left[15]+0.0025501*above[14]+0.00123963*left[14]);
		temp[6] = (int)(0.00334569*above[-1]+0.0102734*above[0]+0.0090881*left[0]+0.0217722*above[1]+0.0159935*left[1]+0.0368328*above[2]+0.0196062*left[2]+0.0601725*above[3]+0.0201218*left[3]+0.104032*above[4]+0.0184197*left[4]+0.144333*above[5]+0.015563*left[5]+0.163532*above[6]+0.0124301*left[6]+0.130221*above[7]+0.0095618*left[7]+0.0830111*above[8]+0.00718944*left[8]+0.0399764*above[9]+0.00534387*left[9]+0.0217751*above[10]+0.00396137*left[10]+0.0127167*above[11]+0.00295128*left[11]+0.00793311*above[12]+0.00222892*left[12]+0.00530932*above[13]+0.00173093*left[13]+0.00226346*above[15]+0.000858901*left[15]+0.00391207*above[14]+0.00141606*left[14]);
		temp[7] = (int)(0.00244329*above[-1]+0.00741815*above[0]+0.00670381*left[0]+0.0154298*above[1]+0.0119631*left[1]+0.0252174*above[2]+0.0150085*left[2]+0.0391528*above[3]+0.0158935*left[3]+0.0617633*above[4]+0.0151086*left[4]+0.105147*above[5]+0.0133077*left[5]+0.145137*above[6]+0.0110911*left[6]+0.164134*above[7]+0.00888787*left[7]+0.130695*above[8]+0.00693671*left[8]+0.083409*above[9]+0.00532834*left[9]+0.0403401*above[10]+0.00406363*left[10]+0.0221419*above[11]+0.0031017*left[11]+0.0131285*above[12]+0.00239057*left[12]+0.00846198*above[13]+0.00188702*left[13]+0.00347222*above[15]+0.000952168*left[15]+0.00607322*above[14]+0.0015617*left[14]);
		temp[8] = (int)(0.00184535*above[-1]+0.0055554*above[0]+0.00510215*left[0]+0.0113957*above[1]+0.00920319*left[1]+0.0181617*above[2]+0.011754*left[2]+0.0271064*above[3]+0.0127527*left[3]+0.0404837*above[4]+0.0124849*left[4]+0.0627256*above[5]+0.0113648*left[5]+0.105868*above[6]+0.00980438*left[6]+0.145703*above[7]+0.00813026*left[7]+0.164608*above[8]+0.00655462*left[8]+0.131124*above[9]+0.00518656*left[9]+0.0838356*above[10]+0.00406167*left[10]+0.0408077*above[11]+0.00317283*left[11]+0.0227073*above[12]+0.00249411*left[12]+0.0139013*above[13]+0.0020005*left[13]+0.00541912*above[15]+0.00102618*left[15]+0.00963191*above[14]+0.00167454*left[14]);
		temp[9] = (int)(0.00143446*above[-1]+0.00429039*above[0]+0.00398984*left[0]+0.00870852*above[1]+0.00725781*left[1]+0.013618*above[2]+0.00940058*left[2]+0.01974*above[3]+0.0103967*left[3]+0.0282526*above[4]+0.0104195*left[4]+0.0413437*above[5]+0.00973943*left[5]+0.0634008*above[6]+0.00864312*left[6]+0.106431*above[7]+0.00737599*left[7]+0.146211*above[8]+0.00611513*left[8]+0.165108*above[9]+0.00496786*left[9]+0.131664*above[10]+0.00398531*left[10]+0.0844734*above[11]+0.00318084*left[11]+0.0416284*above[12]+0.0025473*left[12]+0.0238884*above[13]+0.00207443*left[13]+0.00867656*above[15]+0.0010809*left[15]+0.0157636*above[14]+0.0017553*left[14]);
		temp[10] = (int)(0.00114425*above[-1]+0.00340497*above[0]+0.00319773*left[0]+0.00685503*above[1]+0.00585623*left[1]+0.0105641*above[2]+0.00767093*left[2]+0.0149758*above[3]+0.00861533*left[3]+0.0207627*above[4]+0.00879906*left[4]+0.0290565*above[5]+0.00840443*left[5]+0.0420137*above[6]+0.00763465*left[6]+0.0640023*above[7]+0.00667437*left[7]+0.10702*above[8]+0.00566757*left[8]+0.14684*above[9]+0.00471147*left[9]+0.165841*above[10]+0.00386168*left[10]+0.132586*above[11]+0.00314275*left[11]+0.0857233*above[12]+0.00255996*left[12]+0.0435106*above[13]+0.00211408*left[13]+0.0143884*above[15]+0.00111758*left[15]+0.0269761*above[14]+0.00180677*left[14]);
		temp[11] = (int)(0.000935606*above[-1]+0.00277285*above[0]+0.00262458*left[0]+0.00554657*above[1]+0.00483276*left[1]+0.00845041*above[2]+0.00638802*left[2]+0.0117742*above[3]+0.00726445*left[3]+0.01593*above[4]+0.00753445*left[4]+0.0215582*above[5]+0.00732512*left[5]+0.0297692*above[6]+0.00678398*left[6]+0.0427078*above[7]+0.00605157*left[7]+0.0647398*above[8]+0.00524415*left[8]+0.10787*above[9]+0.0044466*left[9]+0.147896*above[10]+0.00371342*left[10]+0.167241*above[11]+0.00307434*left[11]+0.134575*above[12]+0.00254228*left[12]+0.0888548*above[13]+0.00212577*left[13]+0.0250382*above[15]+0.00113829*left[15]+0.0488747*above[14]+0.00183296*left[14]);
		temp[12] = (int)(0.000784865*above[-1]+0.00231863*above[0]+0.0022084*left[0]+0.00461449*above[1]+0.0040842*left[1]+0.00696743*above[2]+0.00543809*left[2]+0.0095778*above[3]+0.00624657*left[3]+0.0127161*above[4]+0.0065598*left[4]+0.0167733*above[5]+0.00646985*left[5]+0.022377*above[6]+0.0060871*left[6]+0.030635*above[7]+0.00552072*left[7]+0.0437*above[8]+0.00486532*left[8]+0.0659608*above[9]+0.00419418*left[9]+0.10947*above[10]+0.00355825*left[10]+0.150123*above[11]+0.00298898*left[11]+0.170551*above[12]+0.00250362*left[12]+0.140042*above[13]+0.00211576*left[13]+0.0462868*above[15]+0.00114538*left[15]+0.0987234*above[14]+0.00183816*left[14]);
		temp[13] = (int)(0.000677653*above[-1]+0.00199697*above[0]+0.0019112*left[0]+0.0039589*above[1]+0.00354651*left[1]+0.00593667*above[2]+0.00474894*left[2]+0.00807758*above[3]+0.00549761*left[3]+0.0105724*above[4]+0.00582929*left[4]+0.0136802*above[5]+0.00581399*left[5]+0.0177891*above[6]+0.00553772*left[6]+0.0235359*above[7]+0.00508807*left[7]+0.0320532*above[8]+0.0045437*left[8]+0.0455451*above[9]+0.00396836*left[9]+0.0684997*above[10]+0.00340894*left[10]+0.113172*above[11]+0.00289677*left[11]+0.155901*above[12]+0.00245126*left[12]+0.180666*above[13]+0.00208903*left[13]+0.0956629*above[15]+0.00114076*left[15]+0.159453*above[14]+0.00182585*left[14]);
		temp[14] = (int)(0.000605415*above[-1]+0.00178108*above[0]+0.00171022*left[0]+0.00352152*above[1]+0.00318095*left[1]+0.00525609*above[2]+0.00427605*left[2]+0.00710174*above[3]+0.00497678*left[3]+0.00920543*above[4]+0.00531224*left[4]+0.0117573*above[5]+0.00533923*left[5]+0.0150271*above[6]+0.00512877*left[6]+0.0194313*above[7]+0.00475467*left[7]+0.0256631*above[8]+0.00428489*left[8]+0.034963*above[9]+0.00377619*left[9]+0.049746*above[10]+0.00327195*left[10]+0.0749414*above[11]+0.00280262*left[11]+0.123826*above[12]+0.00238837*left[12]+0.175919*above[13]+0.00204731*left[13]+0.159091*above[15]+0.00112466*left[15]+0.223567*above[14]+0.00179678*left[14]);
		temp[15] = (int)(0.000564426*above[-1]+0.00165924*above[0]+0.00159556*left[0]+0.00327678*above[1]+0.00297077*left[1]+0.00488069*above[2]+0.00400044*left[2]+0.0065742*above[3]+0.0046671*left[3]+0.00848531*above[4]+0.00499643*left[4]+0.0107758*above[5]+0.00503904*left[5]+0.0136692*above[6]+0.0048587*left[6]+0.0175004*above[7]+0.00452228*left[7]+0.0228085*above[8]+0.00409208*left[8]+0.0305232*above[9]+0.00362081*left[9]+0.0423661*above[10]+0.00314944*left[10]+0.0618234*above[11]+0.00270732*left[11]+0.0965391*above[12]+0.00231444*left[12]+0.16836*above[13]+0.0019891*left[13]+0.203248*above[15]+0.00109561*left[15]+0.25136*above[14]+0.00174892*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0137858*above[-1]+0.0314926*above[0]+0.0522344*left[0]+0.0452355*above[1]+0.215371*left[1]+0.039263*above[2]+0.233674*left[2]+0.0260676*above[3]+0.207774*left[3]+0.014989*above[4]+0.050316*left[4]+0.00879662*above[5]+0.0211519*left[5]+0.00529943*above[6]+0.0102323*left[6]+0.00328629*above[7]+0.0055154*left[7]+0.00209342*above[8]+0.0031857*left[8]+0.00136725*above[9]+0.00193777*left[9]+0.000914717*above[10]+0.00122864*left[10]+0.000627534*above[11]+0.000808498*left[11]+0.000443454*above[12]+0.000552867*left[12]+0.000326834*above[13]+0.000397181*left[13]+0.000153771*above[15]+0.000182055*left[15]+0.000257732*above[14]+0.000307621*left[14]);
		temp[1] = (int)(0.0196795*above[-1]+0.0491946*above[0]+0.0687404*left[0]+0.0761582*above[1]+0.123559*left[1]+0.0756443*above[2]+0.15086*left[2]+0.0571722*above[3]+0.122206*left[3]+0.0366953*above[4]+0.0777651*left[4]+0.0214343*above[5]+0.0364226*left[5]+0.0128103*above[6]+0.0192711*left[6]+0.00786347*above[7]+0.0108613*left[7]+0.00496411*above[8]+0.00645823*left[8]+0.00321726*above[9]+0.00400526*left[9]+0.00213856*above[10]+0.00257484*left[10]+0.00145925*above[11]+0.00171177*left[11]+0.00102659*above[12]+0.00117969*left[12]+0.000753884*above[13]+0.000852513*left[13]+0.000353355*above[15]+0.000393085*left[15]+0.00059295*above[14]+0.000662989*left[14]);
		temp[2] = (int)(0.0175951*above[-1]+0.0481027*above[0]+0.0481027*left[0]+0.0840749*above[1]+0.0840749*left[1]+0.0973284*above[2]+0.0973284*left[2]+0.0883771*above[3]+0.0883771*left[3]+0.0649553*above[4]+0.0649553*left[4]+0.0415661*above[5]+0.0415661*left[5]+0.024556*above[6]+0.024556*left[6]+0.0148576*above[7]+0.0148576*left[7]+0.00923662*above[8]+0.00923662*left[8]+0.00590699*above[9]+0.00590699*left[9]+0.00388235*above[10]+0.00388235*left[10]+0.00262416*above[11]+0.00262416*left[11]+0.00183162*above[12]+0.00183162*left[12]+0.00133654*above[13]+0.00133654*left[13]+0.000622301*above[15]+0.000622301*left[15]+0.00104643*above[14]+0.00104643*left[14]);
		temp[3] = (int)(0.0122635*above[-1]+0.0374348*above[0]+0.0332276*left[0]+0.0725672*above[1]+0.0575882*left[1]+0.0991457*above[2]+0.0680412*left[2]+0.106642*above[3]+0.0649013*left[3]+0.0942438*above[4]+0.0530088*left[4]+0.0687312*above[5]+0.0384651*left[5]+0.0440494*above[6]+0.0259663*left[6]+0.0262247*above[7]+0.0169969*left[7]+0.0160044*above[8]+0.0111565*left[8]+0.0100453*above[9]+0.00741428*left[9]+0.00649592*above[10]+0.00501347*left[10]+0.0043307*above[11]+0.0034629*left[11]+0.00298822*above[12]+0.00245803*left[12]+0.00216033*above[13]+0.00181695*left[13]+0.000996094*above[15]+0.000857059*left[15]+0.00168012*above[14]+0.00143541*left[14]);
		temp[4] = (int)(0.00858612*above[-1]+0.0262726*above[0]+0.0232882*left[0]+0.0548536*above[1]+0.0407568*left[1]+0.0835774*above[2]+0.0492257*left[2]+0.10615*above[3]+0.0490327*left[3]+0.111177*above[4]+0.0427185*left[4]+0.0972379*above[5]+0.0337013*left[5]+0.0707482*above[6]+0.0247949*left[6]+0.0454376*above[7]+0.0175105*left[7]+0.0272039*above[8]+0.0121522*left[8]+0.0167165*above[9]+0.00842234*left[9]+0.0105849*above[10]+0.00588114*left[10]+0.0069291*above[11]+0.00416553*left[11]+0.00470799*above[12]+0.00301594*left[12]+0.00336124*above[13]+0.00226388*left[13]+0.00152945*above[15]+0.00108459*left[15]+0.00259047*above[14]+0.00180778*left[14]);
		temp[5] = (int)(0.00614548*above[-1]+0.0187248*above[0]+0.0167823*left[0]+0.0389934*above[1]+0.0297015*left[1]+0.0631216*above[2]+0.0366582*left[2]+0.0889807*above[3]+0.037781*left[3]+0.109737*above[4]+0.0345164*left[4]+0.113602*above[5]+0.0288554*left[5]+0.098911*above[6]+0.0226066*left[6]+0.0719295*above[7]+0.0169482*left[7]+0.0462963*above[8]+0.0123832*left[8]+0.0278527*above[9]+0.00894586*left[9]+0.0172336*above[10]+0.00645828*left[10]+0.0110281*above[11]+0.00469875*left[11]+0.00734749*above[12]+0.00347632*left[12]+0.0051616*above[13]+0.00265409*left[13]+0.00230862*above[15]+0.00129362*left[15]+0.00393142*above[14]+0.00214472*left[14]);
		temp[6] = (int)(0.00452352*above[-1]+0.0136894*above[0]+0.0124384*left[0]+0.0282725*above[1]+0.0222453*left[1]+0.0453629*above[2]+0.0279715*left[2]+0.0673882*above[3]+0.0296332*left[3]+0.0918807*above[4]+0.0280731*left[4]+0.111745*above[5]+0.0245123*left[5]+0.115022*above[6]+0.020139*left[6]+0.0999436*above[7]+0.0158338*left[7]+0.0727086*above[8]+0.0120886*left[8]+0.0469143*above[9]+0.00907419*left[9]+0.0283772*above[10]+0.00676637*left[10]+0.0177188*above[11]+0.00505767*left[11]+0.0115265*above[12]+0.00382604*left[12]+0.00793588*above[13]+0.0029732*left[13]+0.00347288*above[15]+0.00147556*left[15]+0.00595507*above[14]+0.00243274*left[14]);
		temp[7] = (int)(0.00342044*above[-1]+0.0102856*above[0]+0.0094632*left[0]+0.0210419*above[1]+0.0170774*left[1]+0.0332966*above[2]+0.0218103*left[2]+0.0488072*above[3]+0.0236271*left[3]+0.0697842*above[4]+0.0230344*left[4]+0.0935808*above[5]+0.0208073*left[5]+0.112983*above[6]+0.0177456*left[6]+0.115956*above[7]+0.0144978*left[7]+0.100682*above[8]+0.0114863*left[8]+0.0733312*above[9]+0.00892007*left[9]+0.0474833*above[10]+0.00685455*left[10]+0.0289484*above[11]+0.0052584*left[11]+0.0183548*above[12]+0.00406594*left[12]+0.0123338*above[13]+0.00321585*left[13]+0.00525303*above[15]+0.00162516*left[15]+0.00908611*above[14]+0.00266439*left[14]);
		temp[8] = (int)(0.00265119*above[-1]+0.0079287*above[0]+0.00737327*left[0]+0.0160835*above[1]+0.0134069*left[1]+0.0250941*above[2]+0.0173449*left[2]+0.0361388*above[3]+0.0191351*left[3]+0.0508335*above[4]+0.0190915*left[4]+0.0712637*above[5]+0.0177218*left[5]+0.0946974*above[6]+0.0155762*left[6]+0.113865*above[7]+0.0131321*left[7]+0.116696*above[8]+0.0107342*left[8]+0.101353*above[9]+0.00858701*left[9]+0.0739954*above[10]+0.00678042*left[10]+0.0482053*above[11]+0.00532876*left[11]+0.0298111*above[12]+0.00420692*left[12]+0.019516*above[13]+0.00338427*left[13]+0.00804256*above[15]+0.00174048*left[15]+0.0140611*above[14]+0.00283801*left[14]);
		temp[9] = (int)(0.00210237*above[-1]+0.00625848*above[0]+0.00587248*left[0]+0.0126044*above[1]+0.0107456*left[1]+0.0194251*above[2]+0.0140515*left[2]+0.0274954*above[3]+0.0157359*left[3]+0.0378998*above[4]+0.0160002*left[4]+0.0521652*above[5]+0.0151867*left[5]+0.0723154*above[6]+0.0136821*left[6]+0.095578*above[7]+0.0118406*left[7]+0.114658*above[8]+0.00993771*left[8]+0.117475*above[9]+0.00815672*left[9]+0.102187*above[10]+0.00659797*left[10]+0.0749683*above[11]+0.00530036*left[11]+0.0494377*above[12]+0.00426556*left[12]+0.0315498*above[13]+0.00348611*left[13]+0.0125323*above[15]+0.00182246*left[15]+0.0221947*above[14]+0.00295664*left[14]);
		temp[10] = (int)(0.00170372*above[-1]+0.00505233*above[0]+0.00477628*left[0]+0.0101144*above[1]+0.00878588*left[1]+0.0154263*above[2]+0.0115917*left[2]+0.0215102*above[3]+0.0131436*left[3]+0.0290783*above[4]+0.0135752*left[4]+0.0391516*above[5]+0.0131244*left[5]+0.0532122*above[6]+0.0120694*left[6]+0.0732555*above[7]+0.0106759*left[7]+0.0964941*above[8]+0.00916307*left[8]+0.11563*above[9]+0.0076891*left[9]+0.118593*above[10]+0.00635247*left[10]+0.103571*above[11]+0.00520358*left[11]+0.0768071*above[12]+0.00426056*left[12]+0.0521352*above[13]+0.0035322*left[13]+0.0199699*above[15]+0.00187421*left[15]+0.0358354*above[14]+0.0030266*left[14]);
		temp[11] = (int)(0.00141108*above[-1]+0.00417123*above[0]+0.00396786*left[0]+0.00830907*above[1]+0.00733089*left[1]+0.0125637*above[2]+0.00974388*left[2]+0.0172989*above[3]+0.011163*left[3]+0.0229958*above[4]+0.0116802*left[4]+0.0303215*above[5]+0.0114662*left[5]+0.0402657*above[6]+0.010726*left[6]+0.0542932*above[7]+0.00966259*left[7]+0.0743949*above[8]+0.00845134*left[8]+0.0977924*above[9]+0.00722656*left[9]+0.117216*above[10]+0.00607972*left[10]+0.120653*above[11]+0.00506513*left[11]+0.10642*above[12]+0.0042102*left[12]+0.0811335*above[13]+0.00353455*left[13]+0.0326043*above[15]+0.00190008*left[15]+0.0592086*above[14]+0.00305601*left[14]);
		temp[12] = (int)(0.00119654*above[-1]+0.00352787*above[0]+0.00337294*left[0]+0.00699903*above[1]+0.00625422*left[1]+0.0105082*above[2]+0.00836346*left[2]+0.0143196*above[3]+0.00966285*left[3]+0.0187719*above[4]+0.0102186*left[4]+0.0243154*above[5]+0.0101575*left[5]+0.0315986*above[6]+0.00963563*left[6]+0.0416055*above[7]+0.00881174*left[7]+0.0558101*above[8]+0.00782804*left[8]+0.0762313*above[9]+0.00679884*left[9]+0.100149*above[10]+0.00580713*left[10]+0.120407*above[11]+0.00490708*left[11]+0.125226*above[12]+0.00413086*left[12]+0.113588*above[13]+0.00350479*left[13]+0.0544807*above[15]+0.00190485*left[15]+0.0931468*above[14]+0.0030533*left[14]);
		temp[13] = (int)(0.00104284*above[-1]+0.00306844*above[0]+0.00294537*left[0]+0.00606831*above[1]+0.00547684*left[1]+0.00906064*above[2]+0.00735886*left[2]+0.0122473*above[3]+0.00855861*left[3]+0.0158801*above[4]+0.00912639*left[4]+0.0202805*above[5]+0.00916083*left[5]+0.0258935*above[6]+0.00878554*left[6]+0.0333775*above[7]+0.00812913*left[7]+0.0437467*above[8]+0.00731005*left[8]+0.0585362*above[9]+0.00642703*left[9]+0.079879*above[10]+0.00555516*left[10]+0.105271*above[11]+0.00474674*left[11]+0.127993*above[12]+0.00403591*left[12]+0.137449*above[13]+0.00345277*left[13]+0.0861027*above[15]+0.00189278*left[15]+0.134439*above[14]+0.00302591*left[14]);
		temp[14] = (int)(0.00094038*above[-1]+0.00276312*above[0]+0.00265954*left[0]+0.00545265*above[1]+0.00495497*left[1]+0.00811058*above[2]+0.00667947*left[2]+0.0109022*above[3]+0.00780382*left[3]+0.0140297*above[4]+0.00836906*left[4]+0.0177426*above[5]+0.00845697*left[5]+0.0223748*above[6]+0.0081713*left[6]+0.0284054*above[7]+0.00762164*left[7]+0.0365602*above[8]+0.00691091*left[8]+0.0479762*above[9]+0.00612703*left[9]+0.0644111*above[10]+0.005339*left[10]+0.0884147*above[11]+0.00459687*left[11]+0.118299*above[12]+0.00393515*left[12]+0.149288*above[13]+0.00338568*left[13]+0.120179*above[15]+0.00186688*left[15]+0.170726*above[14]+0.00297912*left[14]);
		temp[15] = (int)(0.000886281*above[-1]+0.00260262*above[0]+0.00250796*left[0]+0.00513122*above[1]+0.00467642*left[1]+0.00762026*above[2]+0.00631274*left[2]+0.0102191*above[3]+0.00738963*left[3]+0.0131088*above[4]+0.0079441*left[4]+0.01651*above[5]+0.00805046*left[5]+0.0207128*above[6]+0.00780335*left[6]+0.0261282*above[7]+0.00730343*left[7]+0.0333733*above[8]+0.006646*left[8]+0.0434175*above[9]+0.00591332*left[9]+0.0578345*above[10]+0.00517079*left[10]+0.079166*above[11]+0.00446665*left[11]+0.111264*above[12]+0.00383491*left[12]+0.15223*above[13]+0.0033075*left[13]+0.139855*above[15]+0.00182837*left[15]+0.190516*above[14]+0.0029154*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00695955*above[-1]+0.0171945*above[0]+0.0242854*left[0]+0.0267522*above[1]+0.0646388*left[1]+0.0268659*above[2]+0.221902*left[2]+0.0212704*above[3]+0.237312*left[3]+0.0146968*above[4]+0.209909*left[4]+0.0095474*above[5]+0.0516223*left[5]+0.00616285*above[6]+0.0219791*left[6]+0.0040097*above[7]+0.0107722*left[7]+0.0026454*above[8]+0.00587801*left[8]+0.00177408*above[9]+0.00343664*left[9]+0.00121152*above[10]+0.00211772*left[10]+0.000844796*above[11]+0.00136394*left[11]+0.000604816*above[12]+0.00091737*left[12]+0.000450344*above[13]+0.000650702*left[13]+0.000214108*above[15]+0.000294432*left[15]+0.000357699*above[14]+0.000499507*left[14]);
		temp[1] = (int)(0.0114494*above[-1]+0.0293487*above[0]+0.0374185*left[0]+0.0478211*above[1]+0.0900431*left[1]+0.0514215*above[2]+0.135514*left[2]+0.0441632*above[3]+0.157811*left[3]+0.0326455*above[4]+0.126404*left[4]+0.0220969*above[5]+0.0803863*left[5]+0.0143877*above[6]+0.0381078*left[6]+0.00936813*above[7]+0.020384*left[7]+0.00616491*above[8]+0.0116156*left[8]+0.00411983*above[9]+0.0069842*left[9]+0.00280315*above[10]+0.00438486*left[10]+0.00194789*above[11]+0.00286177*left[11]+0.00139023*above[12]+0.00194363*left[12]+0.00103244*above[13]+0.00138863*left[13]+0.000489465*above[15]+0.000632822*left[15]+0.000818455*above[14]+0.00107124*left[14]);
		temp[2] = (int)(0.0122635*above[-1]+0.0332276*above[0]+0.0374348*left[0]+0.0575882*above[1]+0.0725672*left[1]+0.0680412*above[2]+0.0991457*left[2]+0.0649013*above[3]+0.106642*left[3]+0.0530088*above[4]+0.0942438*left[4]+0.0384651*above[5]+0.0687312*left[5]+0.0259663*above[6]+0.0440494*left[6]+0.0169969*above[7]+0.0262247*left[7]+0.0111565*above[8]+0.0160044*left[8]+0.00741428*above[9]+0.0100453*left[9]+0.00501347*above[10]+0.00649592*left[10]+0.0034629*above[11]+0.0043307*left[11]+0.00245803*above[12]+0.00298822*left[12]+0.00181695*above[13]+0.00216033*left[13]+0.000857059*above[15]+0.000996094*left[15]+0.00143541*above[14]+0.00168012*left[14]);
		temp[3] = (int)(0.0107097*above[-1]+0.0304918*above[0]+0.0304918*left[0]+0.0564672*above[1]+0.0564672*left[1]+0.0734235*above[2]+0.0734235*left[2]+0.0785477*above[3]+0.0785477*left[3]+0.0718625*above[4]+0.0718625*left[4]+0.0576593*above[5]+0.0576593*left[5]+0.0416114*above[6]+0.0416114*left[6]+0.0281277*above[7]+0.0281277*left[7]+0.0185087*above[8]+0.0185087*left[8]+0.0122379*above[9]+0.0122379*left[9]+0.00821095*above[10]+0.00821095*left[10]+0.00562527*above[11]+0.00562527*left[11]+0.0039624*above[12]+0.0039624*left[12]+0.00290958*above[13]+0.00290958*left[13]+0.00136251*above[15]+0.00136251*left[15]+0.0022872*above[14]+0.0022872*left[14]);
		temp[4] = (int)(0.00849466*above[-1]+0.02505*above[0]+0.023787*left[0]+0.0487387*above[1]+0.043362*left[1]+0.068828*above[2]+0.0557638*left[2]+0.0816942*above[3]+0.0599101*left[3]+0.0841074*above[4]+0.0565006*left[4]+0.0756396*above[5]+0.0479159*left[5]+0.0602614*above[6]+0.0373317*left[6]+0.0434347*above[7]+0.0273529*left[7]+0.0294326*above[8]+0.0193369*left[8]+0.0194689*above[9]+0.0134806*left[9]+0.0129721*above[10]+0.00941414*left[10]+0.0088038*above[11]+0.00665045*left[11]+0.00614256*above[12]+0.00479736*left[12]+0.00447218*above[13]+0.00358748*left[13]+0.0020744*above[15]+0.0017109*left[15]+0.00349282*above[14]+0.00285592*left[14]);
		temp[5] = (int)(0.00661179*above[-1]+0.0196662*above[0]+0.0184274*left[0]+0.039304*above[1]+0.0335061*left[1]+0.0584707*above[2]+0.043134*left[2]+0.0754303*above[3]+0.0468444*left[3]+0.0862052*above[4]+0.0452762*left[4]+0.0872268*above[5]+0.0399532*left[5]+0.0778307*above[6]+0.032782*left[6]+0.0618315*above[7]+0.0254497*left[7]+0.0445897*above[8]+0.019024*left[8]+0.0303135*above[9]+0.013921*left[9]+0.0201753*above[10]+0.010111*left[10]+0.0135789*above[11]+0.007372*left[11]+0.00937502*above[12]+0.00545389*left[12]+0.00675695*above[13]+0.00415951*left[13]+0.00309768*above[15]+0.00202355*left[15]+0.00523551*above[14]+0.00335712*left[14]);
		temp[6] = (int)(0.00515946*above[-1]+0.0153688*above[0]+0.0143838*left[0]+0.0308913*above[1]+0.0262061*left[1]+0.0470704*above[2]+0.0339342*left[2]+0.0638237*above[3]+0.0373124*left[3]+0.0791518*above[4]+0.0368319*left[4]+0.0888283*above[5]+0.0335043*left[5]+0.0891101*above[6]+0.0285665*left[6]+0.0792169*above[7]+0.0231568*left[7]+0.0628872*above[8]+0.0180888*left[8]+0.0454326*above[9]+0.0137916*left[9]+0.0310305*above[10]+0.0103847*left[10]+0.0208372*above[11]+0.0078059*left[11]+0.0142537*above[12]+0.00592228*left[12]+0.0101617*above[13]+0.00460812*left[13]+0.00459552*above[15]+0.00228806*left[15]+0.00780184*above[14]+0.00377203*left[14]);
		temp[7] = (int)(0.00407126*above[-1]+0.0121103*above[0]+0.0113727*left[0]+0.024327*above[1]+0.0207952*left[1]+0.0371854*above[2]+0.027124*left[2]+0.0514831*above[3]+0.0301909*left[3]+0.0669496*above[4]+0.0303493*left[4]+0.0814031*above[5]+0.0282908*left[5]+0.0904878*above[6]+0.0248556*left[6]+0.0903737*above[7]+0.0208414*left[7]+0.080223*above[8]+0.0168649*left[8]+0.0637377*above[9]+0.0133089*left[9]+0.0462082*above[10]+0.010345*left[10]+0.0318036*above[11]+0.00799795*left[11]+0.0216872*above[12]+0.00621544*left[12]+0.0153146*above[13]+0.00493123*left[13]+0.00682525*above[15]+0.00249803*left[15]+0.0116451*above[14]+0.00409271*left[14]);
		temp[8] = (int)(0.00325892*above[-1]+0.00967144*above[0]+0.00912642*left[0]+0.0193703*above[1]+0.0167546*left[1]+0.0295132*above[2]+0.0220144*left[2]+0.0408891*above[3]+0.0247819*left[3]+0.0541632*above[4]+0.025307*left[4]+0.0689305*above[5]+0.024073*left[5]+0.0829125*above[6]+0.02167*left[6]+0.091688*above[7]+0.0186731*left[7]+0.0913842*above[8]+0.0155532*left[8]+0.0811371*above[9]+0.0126341*left[9]+0.0646364*above[10]+0.0100953*left[10]+0.047173*above[11]+0.0080045*left[11]+0.0329368*above[12]+0.0063596*left[12]+0.0231804*above[13]+0.00513836*left[13]+0.0101966*above[15]+0.00265265*left[15]+0.0174839*above[14]+0.00432059*left[14]);
		temp[9] = (int)(0.00264983*above[-1]+0.00784413*above[0]+0.00743961*left[0]+0.0156539*above[1]+0.0137108*left[1]+0.0237217*above[2]+0.0181388*left[2]+0.0326866*above[3]+0.0206269*left[3]+0.0432448*above[4]+0.0213526*left[4]+0.055962*above[5]+0.0206607*left[5]+0.0703608*above[6]+0.0189767*left[6]+0.0841138*above[7]+0.0167253*left[7]+0.092769*above[8]+0.0142695*left[8]+0.0924371*above[9]+0.0118781*left[9]+0.0822513*above[10]+0.00971993*left[10]+0.0659139*above[11]+0.00788026*left[11]+0.0487562*above[12]+0.00638642*left[12]+0.0351115*above[13]+0.00524616*left[13]+0.0153619*above[15]+0.00275555*left[15]+0.0264303*above[14]+0.00446414*left[14]);
		temp[10] = (int)(0.00219073*above[-1]+0.00646949*above[0]+0.00616524*left[0]+0.0128644*above[1]+0.0114026*left[1]+0.0193812*above[2]+0.0151783*left[2]+0.0265093*above[3]+0.0174151*left[3]+0.0348235*above[4]+0.0182407*left[4]+0.0449464*above[5]+0.017907*left[5]+0.0573898*above[6]+0.016728*left[6]+0.0716411*above[7]+0.0150241*left[7]+0.0853532*above[8]+0.0130792*left[8]+0.0940679*above[9]+0.0111148*left[9]+0.0939052*above[10]+0.0092828*left[10]+0.0840285*above[11]+0.00767272*left[11]+0.0682109*above[12]+0.00632762*left[12]+0.0520114*above[13]+0.00527461*left[13]+0.0233275*above[15]+0.00281338*left[15]+0.0400749*above[14]+0.00453612*left[14]);
		temp[11] = (int)(0.0018443*above[-1]+0.00543462*above[0]+0.00520132*left[0]+0.010771*above[1]+0.00965012*left[1]+0.016138*above[2]+0.0129152*left[2]+0.0219084*above[3]+0.0149336*left[3]+0.0285284*above[4]+0.0157997*left[4]+0.0365189*above[5]+0.0157019*left[5]+0.046464*above[6]+0.014878*left[6]+0.0588527*above[7]+0.0135754*left[7]+0.0731643*above[8]+0.0120193*left[8]+0.0870586*above[9]+0.0103933*left[9]+0.096104*above[10]+0.00883149*left[10]+0.0964778*above[11]+0.00742094*left[11]+0.0874621*above[12]+0.00621221*left[12]+0.0731873*above[13]+0.00524416*left[13]+0.0354996*above[15]+0.0028344*left[15]+0.0596841*above[14]+0.00455126*left[14]);
		temp[12] = (int)(0.0015853*above[-1]+0.00466268*above[0]+0.00447905*left[0]+0.00921464*above[1]+0.00833245*left[1]+0.0137391*above[2]+0.0112031*left[2]+0.0185249*above[3]+0.0130389*left[3]+0.0239181*above[4]+0.0139117*left[4]+0.0303272*above[5]+0.0139671*left[5]+0.0382489*above[6]+0.0133905*left[6]+0.0482572*above[7]+0.012378*left[7]+0.0608476*above[8]+0.0111123*left[8]+0.0755248*above[9]+0.00974744*left[9]+0.0900019*above[10]+0.00840135*left[10]+0.099947*above[11]+0.00715625*left[11]+0.101729*above[12]+0.00606527*left[12]+0.0951685*above[13]+0.00517392*left[13]+0.0527833*above[15]+0.00282722*left[15]+0.0850275*above[14]+0.00452433*left[14]);
		temp[13] = (int)(0.00139768*above[-1]+0.00410464*above[0]+0.00395474*left[0]+0.00809305*above[1]+0.00737297*left[1]+0.012019*above[2]+0.00994961*left[2]+0.0161145*above[3]+0.0116403*left[3]+0.0206558*above[4]+0.0125025*left[4]+0.0259662*above[5]+0.0126532*left[5]+0.0324449*above[6]+0.0122429*left[6]+0.0405946*above[7]+0.0114325*left[7]+0.0510162*above[8]+0.0103748*left[8]+0.0642591*above[9]+0.00920223*left[9]+0.0799232*above[10]+0.00801958*left[10]+0.0958854*above[11]+0.00690363*left[11]+0.108099*above[12]+0.00590767*left[12]+0.113638*above[13]+0.00508048*left[13]+0.0738925*above[15]+0.00279977*left[15]+0.11263*above[14]+0.00446865*left[14]);
		temp[14] = (int)(0.00127291*above[-1]+0.00373434*above[0]+0.00360535*left[0]+0.00735115*above[1]+0.0067316*left[1]+0.0108872*above[2]+0.00910711*left[2]+0.0145397*above[3]+0.0106927*left[3]+0.0185418*above[4]+0.011537*left[4]+0.0231638*above[5]+0.01174*left[5]+0.0287377*above[6]+0.0114304*left[6]+0.0356919*above[7]+0.0107473*left[7]+0.0445871*above[8]+0.00982443*left[8]+0.0561229*above[9]+0.00877963*left[9]+0.071006*above[10]+0.00770849*left[10]+0.089075*above[11]+0.00668322*left[11]+0.108522*above[12]+0.00575614*left[12]+0.125632*above[13]+0.00497718*left[13]+0.0936672*above[15]+0.00275844*left[15]+0.136338*above[14]+0.00439496*left[14]);
		temp[15] = (int)(0.00120963*above[-1]+0.00354722*above[0]+0.00342751*left[0]+0.00697833*above[1]+0.00640339*left[1]+0.0103237*above[2]+0.00867191*left[2]+0.0137655*above[3]+0.0101964*left[3]+0.017519*above[4]+0.0110218*left[4]+0.0218323*above[5]+0.0112406*left[5]+0.0270102*above[6]+0.0109721*left[6]+0.0334497*above[7]+0.0103453*left[7]+0.0416854*above[8]+0.00948538*left[8]+0.0524407*above[9]+0.0085029*left[9]+0.0666375*above[10]+0.00748859*left[10]+0.0851984*above[11]+0.00651179*left[11]+0.107702*above[12]+0.0056236*left[12]+0.131107*above[13]+0.00487362*left[13]+0.104576*above[15]+0.00270745*left[15]+0.148792*above[14]+0.00431055*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00404048*above[-1]+0.0103991*above[0]+0.013256*left[0]+0.0170594*above[1]+0.0318839*left[1]+0.018658*above[2]+0.0689792*left[2]+0.0165205*above[3]+0.224473*left[3]+0.0128267*above[4]+0.238894*left[4]+0.00922672*above[5]+0.210914*left[5]+0.00640032*above[6]+0.0522794*left[6]+0.0043913*above[7]+0.0224208*left[7]+0.00301519*above[8]+0.0110779*left[8]+0.00208537*above[9]+0.00609679*left[9]+0.00145916*above[10]+0.0036003*left[10]+0.00103753*above[11]+0.00224784*left[11]+0.000754591*above[12]+0.00147668*left[12]+0.000568887*above[13]+0.00102885*left[13]+0.000273924*above[15]+0.000457204*left[15]+0.00045583*above[14]+0.000780012*left[14]);
		temp[1] = (int)(0.00715666*above[-1]+0.0188354*above[0]+0.0228081*left[0]+0.0317397*above[1]+0.0512877*left[1]+0.036198*above[2]+0.0983017*left[2]+0.0336736*above[3]+0.14055*left[3]+0.0274233*above[4]+0.160974*left[4]+0.0204473*above[5]+0.128444*left[5]+0.0144846*above[6]+0.0817366*left[6]+0.0100255*above[7]+0.0390242*left[7]+0.00690626*above[8]+0.021023*left[8]+0.0047788*above[9]+0.012076*left[9]+0.00334098*above[10]+0.00733052*left[10]+0.00237221*above[11]+0.00466141*left[11]+0.0017225*above[12]+0.00310208*left[12]+0.00129659*above[13]+0.0021815*left[13]+0.0006232*above[15]+0.000978168*left[15]+0.00103765*above[14]+0.00166422*left[14]);
		temp[2] = (int)(0.00858612*above[-1]+0.0232882*above[0]+0.0262726*left[0]+0.0407568*above[1]+0.0548536*left[1]+0.0492257*above[2]+0.0835774*left[2]+0.0490327*above[3]+0.10615*left[3]+0.0427185*above[4]+0.111177*left[4]+0.0337013*above[5]+0.0972379*left[5]+0.0247949*above[6]+0.0707482*left[6]+0.0175105*above[7]+0.0454376*left[7]+0.0121522*above[8]+0.0272039*left[8]+0.00842234*above[9]+0.0167165*left[9]+0.00588114*above[10]+0.0105849*left[10]+0.00416553*above[11]+0.0069291*left[11]+0.00301594*above[12]+0.00470799*left[12]+0.00226388*above[13]+0.00336124*left[13]+0.00108459*above[15]+0.00152945*left[15]+0.00180778*above[14]+0.00259047*left[14]);
		temp[3] = (int)(0.00849466*above[-1]+0.023787*above[0]+0.02505*left[0]+0.043362*above[1]+0.0487387*left[1]+0.0557638*above[2]+0.068828*left[2]+0.0599101*above[3]+0.0816942*left[3]+0.0565006*above[4]+0.0841074*left[4]+0.0479159*above[5]+0.0756396*left[5]+0.0373317*above[6]+0.0602614*left[6]+0.0273529*above[7]+0.0434347*left[7]+0.0193369*above[8]+0.0294326*left[8]+0.0134806*above[9]+0.0194689*left[9]+0.00941414*above[10]+0.0129721*left[10]+0.00665045*above[11]+0.0088038*left[11]+0.00479736*above[12]+0.00614256*left[12]+0.00358748*above[13]+0.00447218*left[13]+0.0017109*above[15]+0.0020744*left[15]+0.00285592*above[14]+0.00349282*left[14]);
		temp[4] = (int)(0.00756475*above[-1]+0.0217221*above[0]+0.0217221*left[0]+0.0410481*above[1]+0.0410481*left[1]+0.0558855*above[2]+0.0558855*left[2]+0.0645998*above[3]+0.0645998*left[3]+0.0660979*above[4]+0.0660979*left[4]+0.0608414*above[5]+0.0608414*left[5]+0.0509839*above[6]+0.0509839*left[6]+0.0395265*above[7]+0.0395265*left[7]+0.0289505*above[8]+0.0289505*left[8]+0.0205286*above[9]+0.0205286*left[9]+0.0144013*above[10]+0.0144013*left[10]+0.0101623*above[11]+0.0101623*left[11]+0.00730406*above[12]+0.00730406*left[12]+0.00543833*above[13]+0.00543833*left[13]+0.00257919*above[15]+0.00257919*left[15]+0.00431323*above[14]+0.00431323*left[14]);
		temp[5] = (int)(0.00639724*above[-1]+0.018657*above[0]+0.0181599*left[0]+0.0361445*above[1]+0.0338711*left[1]+0.051417*above[2]+0.0454175*left[2]+0.0632114*above[3]+0.0519751*left[3]+0.0697683*above[4]+0.0532956*left[4]+0.0697649*above[5]+0.0499554*left[5]+0.0634715*above[6]+0.0433519*left[6]+0.052901*above[7]+0.0352686*left[7]+0.0409564*above[8]+0.0272935*left[8]+0.0300526*above[9]+0.0204201*left[9]+0.0214182*above[10]+0.0150117*left[10]+0.0151668*above[11]+0.0110032*left[11]+0.0108799*above[12]+0.00815215*left[12]+0.0080678*above[13]+0.00621424*left[13]+0.00380315*above[15]+0.00301767*left[15]+0.00637317*above[14]+0.00500984*left[14]);
		temp[6] = (int)(0.005308*above[-1]+0.0155841*above[0]+0.0149949*left[0]+0.0306054*above[1]+0.0278268*left[1]+0.044745*above[2]+0.0371155*left[2]+0.0575364*above[3]+0.0424081*left[3]+0.0675762*above[4]+0.043741*left[4]+0.0729098*above[5]+0.0416504*left[5]+0.0720596*above[6]+0.037102*left[6]+0.0651842*above[7]+0.0312606*left[7]+0.0542194*above[8]+0.0251934*left[8]+0.0420163*above[9]+0.0196552*left[9]+0.0309562*above[10]+0.0150297*left[10]+0.0222492*above[11]+0.0114044*left[11]+0.0160053*above[12]+0.00869794*left[12]+0.0118419*above[13]+0.00678555*left[13]+0.00555141*above[15]+0.00337389*left[15]+0.00932169*above[14]+0.00556043*left[14]);
		temp[7] = (int)(0.00438708*above[-1]+0.0129125*above[0]+0.0123742*left[0]+0.0254981*above[1]+0.0229349*left[1]+0.0377903*above[2]+0.0305782*left[2]+0.0499119*above[3]+0.0350272*left[3]+0.0612738*above[4]+0.0364021*left[4]+0.0703146*above[5]+0.0351513*left[5]+0.0749569*above[6]+0.0319737*left[6]+0.0736352*above[7]+0.0276776*left[7]+0.0664477*above[8]+0.0230165*left[8]+0.0552901*above[9]+0.0185645*left[9]+0.0429894*above[10]+0.0146682*left[10]+0.0319165*above[11]+0.0114718*left[11]+0.0232883*above[12]+0.00898397*left[12]+0.0172754*above[13]+0.00716246*left[13]+0.00807385*above[15]+0.00364229*left[15]+0.0135774*above[14]+0.00596101*left[14]);
		temp[8] = (int)(0.00363993*above[-1]+0.0107179*above[0]+0.0102673*left[0]+0.0211966*above[1]+0.0190413*left[1]+0.0315644*above[2]+0.0254361*left[2]+0.0422131*above[3]+0.0292663*left[3]+0.0531677*above[4]+0.0306621*left[4]+0.0637142*above[5]+0.029984*left[5]+0.0721946*above[6]+0.0277518*left[6]+0.0764626*above[7]+0.0245517*left[7]+0.0749062*above[8]+0.0209366*left[8]+0.067594*above[9]+0.0173494*left[9]+0.0564062*above[10]+0.014087*left[10]+0.0441686*above[11]+0.0113075*left[11]+0.0332721*above[12]+0.00906517*left[12]+0.0250287*above[13]+0.00737046*left[13]+0.0117356*above[15]+0.00382605*left[15]+0.0197339*above[14]+0.00622174*left[14]);
		temp[9] = (int)(0.00304489*above[-1]+0.00896065*above[0]+0.0085956*left[0]+0.0177138*above[1]+0.0159643*left[1]+0.0263889*above[2]+0.0213893*left[2]+0.0354159*above[3]+0.0247374*left[3]+0.0451125*above[4]+0.0261255*left[4]+0.055406*above[5]+0.0258387*left[5]+0.0655072*above[6]+0.024272*left[6]+0.0737047*above[7]+0.0218645*left[7]+0.0778175*above[8]+0.0190343*left[8]+0.0762136*above[9]+0.0161288*left[9]+0.0689556*above[10]+0.0133983*left[10]+0.0579338*above[11]+0.0109958*left[11]+0.0460115*above[12]+0.00899645*left[12]+0.0357237*above[13]+0.00744213*left[13]+0.0170481*above[15]+0.00393453*left[15]+0.0285658*above[14]+0.00636179*left[14]);
		temp[10] = (int)(0.00257581*above[-1]+0.00757255*above[0]+0.00727946*left[0]+0.0149501*above[1]+0.0135441*left[1]+0.0222351*above[2]+0.0182074*left[2]+0.0298176*above[3]+0.0211684*left[3]+0.0380732*above[4]+0.0225257*left[4]+0.0472443*above[5]+0.022505*left[5]+0.0572001*above[6]+0.0214118*left[6]+0.0671115*above[7]+0.0195836*left[7]+0.0752434*above[8]+0.017345*left[8]+0.0794048*above[9]+0.0149733*left[9]+0.0779688*above[10]+0.0126789*left[10]+0.0710227*above[11]+0.0106026*left[11]+0.0605179*above[12]+0.00882689*left[12]+0.0495289*above[13]+0.0074114*left[13]+0.0246547*above[15]+0.00398073*left[15]+0.0408463*above[14]+0.00640481*left[14]);
		temp[11] = (int)(0.00220983*above[-1]+0.00648914*above[0]+0.00625264*left[0]+0.0127902*above[1]+0.0116552*left[1]+0.0189753*above[2]+0.0157199*left[2]+0.0253766*above[3]+0.0183679*left[3]+0.0323463*above[4]+0.0196805*left[4]+0.0402031*above[5]+0.0198383*left[5]+0.0491459*above[6]+0.0190822*left[6]+0.0590166*above[7]+0.0176783*left[7]+0.0689732*above[8]+0.0158845*left[8]+0.0772825*above[9]+0.0139264*left[9]+0.0817724*above[10]+0.011982*left[10]+0.0808597*above[11]+0.0101782*left[11]+0.0747224*above[12]+0.00859816*left[12]+0.0656059*above[13]+0.00731034*left[13]+0.0351062*above[15]+0.00397903*left[15]+0.0569015*above[14]+0.00637557*left[14]);
		temp[12] = (int)(0.0019296*above[-1]+0.00565986*above[0]+0.00546591*left[0]+0.0111373*above[1]+0.0102062*left[1]+0.0164792*above[2]+0.0138073*left[2]+0.0219644*above[3]+0.0162051*left[3]+0.0279013*above[4]+0.0174676*left[4]+0.0346008*above[5]+0.0177418*left[5]+0.0423519*above[6]+0.0172228*left[6]+0.0513386*above[7]+0.0161258*left[7]+0.0614026*above[8]+0.0146615*left[8]+0.0717181*above[9]+0.0130174*left[9]+0.0805881*above[10]+0.0113462*left[10]+0.0859067*above[11]+0.0097619*left[11]+0.0862092*above[12]+0.00834488*left[12]+0.0820369*above[13]+0.00716758*left[13]+0.0483125*above[15]+0.00394365*left[15]+0.0758971*above[14]+0.00629768*left[14]);
		temp[13] = (int)(0.00172356*above[-1]+0.00505064*above[0]+0.00488695*left[0]+0.00992422*above[1]+0.00913832*left[1]+0.0146495*above[2]+0.0123937*left[2]+0.0194642*above[3]+0.0145993*left[3]+0.0246361*above[4]+0.0158131*left[4]+0.0304469*above[5]+0.0161588*left[5]+0.0371919*above[6]+0.0157997*left[6]+0.0451589*above[7]+0.0149162*left[7]+0.0545472*above[8]+0.0136862*left[8]+0.0652322*above[9]+0.0122702*left[9]+0.0764434*above[10]+0.010802*left[10]+0.0865647*above[11]+0.00938493*left[11]+0.0936*above[12]+0.00809571*left[12]+0.0963879*above[13]+0.00700777*left[13]+0.0627814*above[15]+0.00388748*left[15]+0.0953884*above[14]+0.00619216*left[14]);
		temp[14] = (int)(0.00158612*above[-1]+0.00464478*above[0]+0.00450028*left[0]+0.00911751*above[1]+0.00842372*left[1]+0.013436*above[2]+0.0114444*left[2]+0.017811*above[3]+0.0135149*left[3]+0.0224822*above[4]+0.014687*left[4]+0.027706*above[5]+0.0150698*left[5]+0.0337637*above[6]+0.0148067*left[6]+0.0409665*above[7]+0.0140564*left[7]+0.0496359*above[8]+0.0129762*left[8]+0.0600221*above[9]+0.0117092*left[9]+0.0720632*above[10]+0.0103765*left[10]+0.0850351*above[11]+0.00907386*left[11]+0.0972664*above[12]+0.00787449*left[12]+0.106694*above[13]+0.00685133*left[13]+0.075359*above[15]+0.0038213*left[15]+0.111468*above[14]+0.00607652*left[14]);
		temp[15] = (int)(0.00151815*above[-1]+0.00444464*above[0]+0.00430847*left[0]+0.00872145*above[1]+0.00806765*left[1]+0.0128446*above[2]+0.0109677*left[2]+0.0170133*above[3]+0.012964*left[3]+0.0214552*above[4]+0.0141058*left[4]+0.0264158*above[5]+0.014496*left[5]+0.0321697*above[6]+0.0142694*left[6]+0.0390338*above[7]+0.0135754*left[7]+0.0473691*above[8]+0.012562*left[8]+0.0575489*above[9]+0.0113644*left[9]+0.0698317*above[10]+0.0100974*left[10]+0.0839893*above[11]+0.00885264*left[11]+0.0987507*above[12]+0.00770102*left[12]+0.111738*above[13]+0.00671423*left[13]+0.0820644*above[15]+0.00375301*left[15]+0.119855*above[14]+0.00596388*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00256612*above[-1]+0.00678583*above[0]+0.00814027*left[0]+0.0115253*above[1]+0.0182888*left[1]+0.0133423*above[2]+0.0349352*left[2]+0.0127167*above[3]+0.0708733*left[3]+0.0107162*above[4]+0.225683*left[4]+0.0083377*above[5]+0.239687*left[5]+0.00618347*above[6]+0.211448*left[6]+0.00447179*above[7]+0.0526496*left[7]+0.00320064*above[8]+0.0226856*left[8]+0.00228803*above[9]+0.0112753*left[9]+0.00164424*above[10]+0.0062526*left[10]+0.00119484*above[11]+0.00373308*left[11]+0.000884539*above[12]+0.00237316*left[12]+0.000676291*above[13]+0.00161337*left[13]+0.000330357*above[15]+0.000699541*left[15]+0.000547293*above[14]+0.00120259*left[14]);
		temp[1] = (int)(0.00477171*above[-1]+0.0127885*above[0]+0.0148783*left[0]+0.0220947*above[1]+0.0323691*left[1]+0.0262678*above[2]+0.0572505*left[2]+0.0258629*above[3]+0.102081*left[3]+0.022529*above[4]+0.143002*left[4]+0.0180398*above[5]+0.162602*left[5]+0.013662*above[6]+0.129551*left[6]+0.0100102*above[7]+0.0825093*left[7]+0.00721624*above[8]+0.0395806*left[8]+0.00517762*above[9]+0.0214402*left[9]+0.00372682*above[10]+0.0124067*left[10]+0.00270943*above[11]+0.00761322*left[11]+0.00200542*above[12]+0.00492851*left[12]+0.00153256*above[13]+0.00339306*left[13]+0.000748077*above[15]+0.00148874*left[15]+0.00123963*above[14]+0.0025501*left[14]);
		temp[2] = (int)(0.00614548*above[-1]+0.0167823*above[0]+0.0187248*left[0]+0.0297015*above[1]+0.0389934*left[1]+0.0366582*above[2]+0.0631216*left[2]+0.037781*above[3]+0.0889807*left[3]+0.0345164*above[4]+0.109737*left[4]+0.0288554*above[5]+0.113602*left[5]+0.0226066*above[6]+0.098911*left[6]+0.0169482*above[7]+0.0719295*left[7]+0.0123832*above[8]+0.0462963*left[8]+0.00894586*above[9]+0.0278527*left[9]+0.00645828*above[10]+0.0172336*left[10]+0.00469875*above[11]+0.0110281*left[11]+0.00347632*above[12]+0.00734749*left[12]+0.00265409*above[13]+0.0051616*left[13]+0.00129362*above[15]+0.00230862*left[15]+0.00214472*above[14]+0.00393142*left[14]);
		temp[3] = (int)(0.00661179*above[-1]+0.0184274*above[0]+0.0196662*left[0]+0.0335061*above[1]+0.039304*left[1]+0.043134*above[2]+0.0584707*left[2]+0.0468444*above[3]+0.0754303*left[3]+0.0452762*above[4]+0.0862052*left[4]+0.0399532*above[5]+0.0872268*left[5]+0.032782*above[6]+0.0778307*left[6]+0.0254497*above[7]+0.0618315*left[7]+0.019024*above[8]+0.0445897*left[8]+0.013921*above[9]+0.0303135*left[9]+0.010111*above[10]+0.0201753*left[10]+0.007372*above[11]+0.0135789*left[11]+0.00545389*above[12]+0.00937502*left[12]+0.00415951*above[13]+0.00675695*left[13]+0.00202355*above[15]+0.00309768*left[15]+0.00335712*above[14]+0.00523551*left[14]);
		temp[4] = (int)(0.00639724*above[-1]+0.0181599*above[0]+0.018657*left[0]+0.0338711*above[1]+0.0361445*left[1]+0.0454175*above[2]+0.051417*left[2]+0.0519751*above[3]+0.0632114*left[3]+0.0532956*above[4]+0.0697683*left[4]+0.0499554*above[5]+0.0697649*left[5]+0.0433519*above[6]+0.0634715*left[6]+0.0352686*above[7]+0.052901*left[7]+0.0272935*above[8]+0.0409564*left[8]+0.0204201*above[9]+0.0300526*left[9]+0.0150117*above[10]+0.0214182*left[10]+0.0110032*above[11]+0.0151668*left[11]+0.00815215*above[12]+0.0108799*left[12]+0.00621424*above[13]+0.0080678*left[13]+0.00301767*above[15]+0.00380315*left[15]+0.00500984*above[14]+0.00637317*left[14]);
		temp[5] = (int)(0.00581681*above[-1]+0.0167395*above[0]+0.0167395*left[0]+0.0318684*above[1]+0.0318684*left[1]+0.0442247*above[2]+0.0442247*left[2]+0.0530361*above[3]+0.0530361*left[3]+0.0575316*above[4]+0.0575316*left[4]+0.057347*above[5]+0.057347*left[5]+0.0529275*above[6]+0.0529275*left[6]+0.0455589*above[7]+0.0455589*left[7]+0.0369398*above[8]+0.0369398*left[8]+0.0285963*above[9]+0.0285963*left[9]+0.0214792*above[10]+0.0214792*left[10]+0.0159243*above[11]+0.0159243*left[11]+0.0118541*above[12]+0.0118541*left[12]+0.00904617*above[13]+0.00904617*left[13]+0.00438881*above[15]+0.00438881*left[15]+0.0072897*above[14]+0.0072897*left[14]);
		temp[6] = (int)(0.00511434*above[-1]+0.0148459*above[0]+0.0146088*left[0]+0.0286638*above[1]+0.0275499*left[1]+0.0407981*above[2]+0.0377376*left[2]+0.0507908*above[3]+0.0446934*left[3]+0.0578485*above[4]+0.0481179*left[4]+0.0610742*above[5]+0.0480203*left[5]+0.0599834*above[6]+0.0448632*left[6]+0.0549253*above[7]+0.0395345*left[7]+0.0471147*above[8]+0.0331318*left[8]+0.0381996*above[9]+0.0266695*left[9]+0.0296722*above[10]+0.0208636*left[10]+0.0224637*above[11]+0.016072*left[11]+0.0169056*above[12]+0.0123699*left[12]+0.0129591*above[13]+0.00969835*left[13]+0.00629542*above[15]+0.00483759*left[15]+0.0104556*above[14]+0.00796635*left[14]);
		temp[7] = (int)(0.00442633*above[-1]+0.0129097*above[0]+0.0125941*left[0]+0.0251352*above[1]+0.0236367*left[1]+0.0363634*above[2]+0.0321718*left[2]+0.0464862*above[3]+0.0378905*left[3]+0.0549985*above[4]+0.0407084*left[4]+0.0609895*above[5]+0.040775*left[5]+0.0634583*above[6]+0.0385105*left[6]+0.0618399*above[7]+0.0345669*left[7]+0.0564252*above[8]+0.0297043*left[8]+0.0483883*above[9]+0.0246343*left[9]+0.0393517*above[10]+0.0198981*left[10]+0.0307958*above[11]+0.0158191*left[11]+0.0236569*above[12]+0.0125308*left[12]+0.0183292*above[13]+0.010065*left[13]+0.00895582*above[15]+0.00514965*left[15]+0.0148549*above[14]+0.0084134*left[14]);
		temp[8] = (int)(0.0038106*above[-1]+0.0111401*above[0]+0.0108221*left[0]+0.0217847*above[1]+0.0202676*left[1]+0.0318069*above[2]+0.0275149*left[2]+0.0413369*above[3]+0.0323537*left[3]+0.050216*above[4]+0.0347937*left[4]+0.057837*above[5]+0.0350245*left[5]+0.0632021*above[6]+0.0334095*left[6]+0.0652439*above[7]+0.0304464*left[7]+0.0633506*above[8]+0.0266904*left[8]+0.0577817*above[9]+0.0226635*left[9]+0.0496937*above[10]+0.0187822*left[10]+0.0407055*above[11]+0.0153241*left[11]+0.0323138*above[12]+0.0124377*left[12]+0.0255494*above[13]+0.0102019*left[13]+0.0126581*above[15]+0.00533775*left[15]+0.0209221*above[14]+0.00865991*left[14]);
		temp[9] = (int)(0.00328462*above[-1]+0.00961168*above[0]+0.00932212*left[0]+0.0188327*above[1]+0.0174473*left[1]+0.0276192*above[2]+0.0236755*left[2]+0.0362152*above[3]+0.0278559*left[3]+0.0447062*above[4]+0.030036*left[4]+0.0528478*above[5]+0.0304043*left[5]+0.0599614*above[6]+0.0292681*left[6]+0.0649955*above[7]+0.027018*left[7]+0.0668464*above[8]+0.0240759*left[8]+0.0648793*above[9]+0.0208393*left[9]+0.0593444*above[10]+0.0176357*left[10]+0.0514036*above[11]+0.0146995*left[11]+0.0427057*above[12]+0.0121763*left[12]+0.0348823*above[13]+0.0101665*left[13]+0.0177487*above[15]+0.00542208*left[15]+0.0291204*above[14]+0.0087441*left[14]);
		temp[10] = (int)(0.00284798*above[-1]+0.00833544*above[0]+0.00808278*left[0]+0.0163413*above[1]+0.0151308*left[1]+0.0240049*above[2]+0.0205469*left[2]+0.0316008*above[3]+0.0242182*left[3]+0.0393374*above[4]+0.0262041*left[4]+0.0472301*above[5]+0.0266769*left[5]+0.0549773*above[6]+0.0258946*left[6]+0.0618582*above[7]+0.0241702*left[7]+0.0667944*above[8]+0.021835*left[8]+0.0686679*above[9]+0.019201*left[9]+0.0668433*above[10]+0.0165312*left[10]+0.0615859*above[11]+0.0140241*left[11]+0.0541027*above[12]+0.0118148*left[12]+0.0462248*above[13]+0.0100121*left[13]+0.024556*above[15]+0.00542596*left[15]+0.039777*above[14]+0.00870648*left[14]);
		temp[11] = (int)(0.00249397*above[-1]+0.0072976*above[0]+0.00708045*left[0]+0.014304*above[1]+0.0132627*left[1]+0.0210139*above[2]+0.0180336*left[2]+0.0276941*above[3]+0.0213052*left[3]+0.0345933*above[4]+0.023138*left[4]+0.0418643*above[5]+0.0236845*left[5]+0.0494777*above[6]+0.0231619*left[6]+0.0571005*above[7]+0.0218262*left[7]+0.0639943*above[8]+0.0199446*left[8]+0.0690752*above[9]+0.0177692*left[9]+0.0712332*above[10]+0.0155159*left[10]+0.0698577*above[11]+0.0133542*left[11]+0.0652709*above[12]+0.0114072*left[12]+0.0589057*above[13]+0.00978479*left[13]+0.0331989*above[15]+0.00537277*left[15]+0.0527904*above[14]+0.00858579*left[14]);
		temp[12] = (int)(0.00221527*above[-1]+0.00647935*above[0]+0.0062922*left[0]+0.0126933*above[1]+0.0117954*left[1]+0.018635*above[2]+0.0160623*left[2]+0.02455*above[3]+0.0190219*left[3]+0.0306886*above[4]+0.0207317*left[4]+0.0372573*above[5]+0.0213258*left[5]+0.0443756*above[6]+0.0209895*left[6]+0.0519931*above[7]+0.0199371*left[7]+0.0597684*above[8]+0.0183904*left[8]+0.0669661*above[9]+0.0165586*left[9]+0.0725172*above[10]+0.0146237*left[10]+0.0753398*above[11]+0.0127323*left[11]+0.0748754*above[12]+0.0109967*left[12]+0.0716798*above[13]+0.00952416*left[13]+0.0432742*above[15]+0.00528426*left[15]+0.0673047*above[14]+0.00841695*left[14]);
		temp[13] = (int)(0.00200665*above[-1]+0.00586656*above[0]+0.00570229*left[0]+0.0114858*above[1]+0.0106975*left[1]+0.0168472*above[2]+0.0145871*left[2]+0.0221743*above[3]+0.0173112*left[3]+0.0277059*above[4]+0.018924*left[4]+0.0336625*above[5]+0.0195446*left[5]+0.0402284*above[6]+0.0193351*left[6]+0.0475122*above[7]+0.0184804*left[7]+0.0554635*above[8]+0.0171705*left[8]+0.0637494*above[9]+0.0155855*left[9]+0.0716478*above[10]+0.0138832*left[10]+0.0780992*above[11]+0.0121933*left[11]+0.0820187*above[12]+0.0106188*left[12]+0.0830001*above[13]+0.00926353*left[13]+0.053574*above[15]+0.00517961*left[15]+0.0815323*above[14]+0.00823033*left[14]);
		temp[14] = (int)(0.00186655*above[-1]+0.00545523*above[0]+0.00530597*left[0]+0.0106756*above[1]+0.00995927*left[1]+0.015648*above[2]+0.0135934*left[2]+0.0205797*above[3]+0.0161555*left[3]+0.0256977*above[4]+0.0176967*left[4]+0.0312229*above[5]+0.0183266*left[5]+0.0373644*above[6]+0.0181921*left[6]+0.0443044*above[7]+0.0174596*left[7]+0.0521559*above[8]+0.0162996*left[8]+0.0608715*above[9]+0.0148734*left[9]+0.070103*above[10]+0.0133234*left[10]+0.0790607*above[11]+0.0117681*left[11]+0.0865271*above[12]+0.0103038*left[12]+0.0913767*above[13]+0.00903101*left[13]+0.062098*above[15]+0.00507495*left[15]+0.092942*above[14]+0.00805142*left[14]);
		temp[15] = (int)(0.00179841*above[-1]+0.00525562*above[0]+0.00511273*left[0]+0.0102839*above[1]+0.00959801*left[1]+0.0150714*above[2]+0.013104*left[2]+0.0198188*above[3]+0.0155808*left[3]+0.0247478*above[4]+0.0170783*left[4]+0.0300789*above[5]+0.0177021*left[5]+0.0360307*above[6]+0.0175928*left[6]+0.0428136*above[7]+0.016909*left[7]+0.0506044*above[8]+0.0158128*left[8]+0.0594815*above[9]+0.0144573*left[9]+0.0692868*above[10]+0.0129779*left[10]+0.0794182*above[11]+0.0114877*left[11]+0.0886659*above[12]+0.010079*left[12]+0.0955796*above[13]+0.00884995*left[13]+0.0665598*above[15]+0.00498285*left[15]+0.098841*above[14]+0.00790053*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00174105*above[-1]+0.00468857*above[0]+0.0054025*left[0]+0.00815837*above[1]+0.0116713*left[1]+0.00982467*above[2]+0.0205292*left[2]+0.00986169*above[3]+0.0363791*left[3]+0.00881624*above[4]+0.0718257*left[4]+0.0072859*above[5]+0.226326*left[5]+0.00571382*above[6]+0.240134*left[6]+0.00433653*above[7]+0.211768*left[7]+0.00323119*above[8]+0.0528873*left[8]+0.00238802*above[9]+0.0228723*left[9]+0.0017641*above[10]+0.0114328*left[10]+0.00131166*above[11]+0.00639817*left[11]+0.000989559*above[12]+0.00388347*left[12]+0.000768087*above[13]+0.00255284*left[13]+0.000381043*above[15]+0.00107062*left[15]+0.000628246*above[14]+0.00185959*left[14]);
		temp[1] = (int)(0.00334569*above[-1]+0.0090881*above[0]+0.0102734*left[0]+0.0159935*above[1]+0.0217722*left[1]+0.0196062*above[2]+0.0368328*left[2]+0.0201218*above[3]+0.0601725*left[3]+0.0184197*above[4]+0.104032*left[4]+0.015563*above[5]+0.144333*left[5]+0.0124301*above[6]+0.163532*left[6]+0.0095618*above[7]+0.130221*left[7]+0.00718944*above[8]+0.0830111*left[8]+0.00534387*above[9]+0.0399764*left[9]+0.00396137*above[10]+0.0217751*left[10]+0.00295128*above[11]+0.0127167*left[11]+0.00222892*above[12]+0.00793311*left[12]+0.00173093*above[13]+0.00530932*left[13]+0.000858901*above[15]+0.00226346*left[15]+0.00141606*above[14]+0.00391207*left[14]);
		temp[2] = (int)(0.00452352*above[-1]+0.0124384*above[0]+0.0136894*left[0]+0.0222453*above[1]+0.0282725*left[1]+0.0279715*above[2]+0.0453629*left[2]+0.0296332*above[3]+0.0673882*left[3]+0.0280731*above[4]+0.0918807*left[4]+0.0245123*above[5]+0.111745*left[5]+0.020139*above[6]+0.115022*left[6]+0.0158338*above[7]+0.0999436*left[7]+0.0120886*above[8]+0.0727086*left[8]+0.00907419*above[9]+0.0469143*left[9]+0.00676637*above[10]+0.0283772*left[10]+0.00505767*above[11]+0.0177188*left[11]+0.00382604*above[12]+0.0115265*left[12]+0.0029732*above[13]+0.00793588*left[13]+0.00147556*above[15]+0.00347288*left[15]+0.00243274*above[14]+0.00595507*left[14]);
		temp[3] = (int)(0.00515946*above[-1]+0.0143838*above[0]+0.0153688*left[0]+0.0262061*above[1]+0.0308913*left[1]+0.0339342*above[2]+0.0470704*left[2]+0.0373124*above[3]+0.0638237*left[3]+0.0368319*above[4]+0.0791518*left[4]+0.0335043*above[5]+0.0888283*left[5]+0.0285665*above[6]+0.0891101*left[6]+0.0231568*above[7]+0.0792169*left[7]+0.0180888*above[8]+0.0628872*left[8]+0.0137916*above[9]+0.0454326*left[9]+0.0103847*above[10]+0.0310305*left[10]+0.0078059*above[11]+0.0208372*left[11]+0.00592228*above[12]+0.0142537*left[12]+0.00460812*above[13]+0.0101617*left[13]+0.00228806*above[15]+0.00459552*left[15]+0.00377203*above[14]+0.00780184*left[14]);
		temp[4] = (int)(0.005308*above[-1]+0.0149949*above[0]+0.0155841*left[0]+0.0278268*above[1]+0.0306054*left[1]+0.0371155*above[2]+0.044745*left[2]+0.0424081*above[3]+0.0575364*left[3]+0.043741*above[4]+0.0675762*left[4]+0.0416504*above[5]+0.0729098*left[5]+0.037102*above[6]+0.0720596*left[6]+0.0312606*above[7]+0.0651842*left[7]+0.0251934*above[8]+0.0542194*left[8]+0.0196552*above[9]+0.0420163*left[9]+0.0150297*above[10]+0.0309562*left[10]+0.0114044*above[11]+0.0222492*left[11]+0.00869794*above[12]+0.0160053*left[12]+0.00678555*above[13]+0.0118419*left[13]+0.00337389*above[15]+0.00555141*left[15]+0.00556043*above[14]+0.00932169*left[14]);
		temp[5] = (int)(0.00511434*above[-1]+0.0146088*above[0]+0.0148459*left[0]+0.0275499*above[1]+0.0286638*left[1]+0.0377376*above[2]+0.0407981*left[2]+0.0446934*above[3]+0.0507908*left[3]+0.0481179*above[4]+0.0578485*left[4]+0.0480203*above[5]+0.0610742*left[5]+0.0448632*above[6]+0.0599834*left[6]+0.0395345*above[7]+0.0549253*left[7]+0.0331318*above[8]+0.0471147*left[8]+0.0266695*above[9]+0.0381996*left[9]+0.0208636*above[10]+0.0296722*left[10]+0.016072*above[11]+0.0224637*left[11]+0.0123699*above[12]+0.0169056*left[12]+0.00969835*above[13]+0.0129591*left[13]+0.00483759*above[15]+0.00629542*left[15]+0.00796635*above[14]+0.0104556*left[14]);
		temp[6] = (int)(0.00473139*above[-1]+0.0136262*above[0]+0.0136262*left[0]+0.0260211*above[1]+0.0260211*left[1]+0.0364214*above[2]+0.0364214*left[2]+0.044474*above[3]+0.044474*left[3]+0.0497674*above[4]+0.0497674*left[4]+0.0519404*above[5]+0.0519404*left[5]+0.0509207*above[6]+0.0509207*left[6]+0.0470962*above[7]+0.0470962*left[7]+0.0412944*above[8]+0.0412944*left[8]+0.0345672*above[9]+0.0345672*left[9]+0.027897*above[10]+0.027897*left[10]+0.0219802*above[11]+0.0219802*left[11]+0.0171699*above[12]+0.0171699*left[12]+0.0135812*above[13]+0.0135812*left[13]+0.00681733*above[15]+0.00681733*left[15]+0.0112076*above[14]+0.0112076*left[14]);
		temp[7] = (int)(0.00427437*above[-1]+0.0123779*above[0]+0.0122484*left[0]+0.0238467*above[1]+0.0232321*left[1]+0.0339128*above[2]+0.0321911*left[2]+0.0424109*above[3]+0.0388668*left[3]+0.0490077*above[4]+0.0430683*left[4]+0.0532178*above[5]+0.0446991*left[5]+0.0546012*above[6]+0.0438622*left[6]+0.0530178*above[7]+0.0409253*left[7]+0.048803*above[8]+0.0364929*left[8]+0.042746*above[9]+0.0312863*left[9]+0.0358732*above[10]+0.0259886*left[10]+0.0291539*above[11]+0.0211221*left[11]+0.0232877*above[12]+0.0170056*left[12]+0.0186896*above[13]+0.0138119*left[13]+0.00948863*above[15]+0.00713368*left[15]+0.0155505*above[14]+0.0116233*left[14]);
		temp[8] = (int)(0.00381399*above[-1]+0.0110825*above[0]+0.0108962*left[0]+0.0214731*above[1]+0.0205849*left[1]+0.0308664*above[2]+0.0283573*left[2]+0.039266*above[3]+0.0340237*left[3]+0.0465034*above[4]+0.0375138*left[4]+0.0521719*above[5]+0.0388584*left[5]+0.0557141*above[6]+0.0382284*left[6]+0.056631*above[7]+0.0359546*left[7]+0.0547382*above[8]+0.0324999*left[8]+0.0503398*above[9]+0.0283885*left[9]+0.0442059*above[10]+0.0241181*left[10]+0.0373563*above[11]+0.0200883*left[11]+0.0307725*above[12]+0.0165717*left[12]+0.0252433*above[13]+0.0137562*left[13]+0.013058*above[15]+0.00727687*left[15]+0.0212861*above[14]+0.0117676*left[14]);
		temp[9] = (int)(0.00338752*above[-1]+0.00986293*above[0]+0.00966119*left[0]+0.0191757*above[1]+0.0182114*left[1]+0.027749*above[2]+0.0250091*left[2]+0.0356994*above[3]+0.0299123*left[3]+0.0430183*above[4]+0.032914*left[4]+0.04947*above[5]+0.0341021*left[5]+0.0545847*above[6]+0.0336664*left[6]+0.0577551*above[7]+0.0318973*left[7]+0.0584456*above[8]+0.0291624*left[8]+0.0564472*above[9]+0.0258603*left[9]+0.0520514*above[10]+0.0223684*left[10]+0.0460279*above[11]+0.0189996*left[11]+0.0394182*above[12]+0.015984*left[12]+0.0333243*above[13]+0.0135055*left[13]+0.0177321*above[15]+0.00728561*left[15]+0.0286645*above[14]+0.0117091*left[14]);
		temp[10] = (int)(0.00301209*above[-1]+0.00877928*above[0]+0.0085828*left[0]+0.0171015*above[1]+0.0161608*left[1]+0.0248433*above[2]+0.0221609*left[2]+0.0321804*above[3]+0.0264736*left[3]+0.039216*above[4]+0.0291237*left[4]+0.0458829*above[5]+0.030221*left[5]+0.0518918*above[6]+0.0299519*left[6]+0.0567315*above[7]+0.0285697*left[7]+0.0597655*above[8]+0.0263738*left[8]+0.06044*above[9]+0.0236789*left[9]+0.0585389*above[10]+0.0207813*left[10]+0.0543589*above[11]+0.0179318*left[11]+0.0486988*above[12]+0.0153256*left[12]+0.0427526*above[13]+0.0131353*left[13]+0.0236416*above[15]+0.00719857*left[15]+0.0377692*above[14]+0.0115117*left[14]);
		temp[11] = (int)(0.00269428*above[-1]+0.00785692*above[0]+0.00767424*left[0]+0.0153193*above[1]+0.0144439*left[1]+0.0222997*above[2]+0.0197977*left[2]+0.0289964*above[3]+0.0236488*left[3]+0.0355747*above[4]+0.0260368*left[4]+0.0420882*above[5]+0.0270763*left[5]+0.0484254*above[6]+0.0269416*left[6]+0.0542636*above[7]+0.0258539*left[7]+0.0590693*above[8]+0.0240632*left[8]+0.0621932*above[9]+0.0218255*left[9]+0.0630781*above[10]+0.0193811*left[10]+0.0615149*above[11]+0.0169363*left[11]+0.0578294*above[12]+0.014658*left[12]+0.0529953*above[13]+0.0127065*left[13]+0.0307098*above[15]+0.00705109*left[15]+0.0483448*above[14]+0.0112315*left[14]);
		temp[12] = (int)(0.00243614*above[-1]+0.00710539*above[0]+0.00693829*left[0]+0.0138594*above[1]+0.0130581*left[1]+0.0201932*above[2]+0.0179001*left[2]+0.0263086*above[3]+0.0213929*left[3]+0.0323982*above[4]+0.0235824*left[4]+0.0385866*above[5]+0.0245806*left[5]+0.044892*above[6]+0.024548*left[6]+0.0511767*above[7]+0.0236796*left[7]+0.0571012*above[8]+0.0221889*left[8]+0.0621219*above[9]+0.0202912*left[9]+0.0655848*above[10]+0.0181872*left[10]+0.0669316*above[11]+0.0160514*left[11]+0.0659641*above[12]+0.0140291*left[12]+0.0631776*above[13]+0.0122687*left[13]+0.0384927*above[15]+0.00687407*left[15]+0.059636*above[14]+0.0109164*left[14]);
		temp[13] = (int)(0.00223887*above[-1]+0.00653014*above[0]+0.00637666*left[0]+0.0127386*above[1]+0.0120024*left[1]+0.0185667*above[2]+0.0164579*left[2]+0.0242115*above[3]+0.0196825*left[3]+0.0298743*above[4]+0.021724*left[4]+0.035716*above[5]+0.0226896*left[5]+0.0418336*above[6]+0.0227278*left[6]+0.0482245*above[7]+0.0220137*left[7]+0.0547363*above[8]+0.0207353*left[8]+0.061016*above[9]+0.0190798*left[9]+0.0665018*above[10]+0.0172209*left[10]+0.0705084*above[11]+0.0153107*left[11]+0.0724348*above[12]+0.0134786*left[12]+0.0722064*above[13]+0.0118635*left[13]+0.0460768*above[15]+0.00669384*left[15]+0.0703365*above[14]+0.0106068*left[14]);
		temp[14] = (int)(0.00210516*above[-1]+0.00614006*above[0]+0.00599605*left[0]+0.011978*above[1]+0.0112871*left[1]+0.0174608*above[2]+0.0154807*left[2]+0.0227799*above[3]+0.0185227*left[3]+0.028138*above[4]+0.0204613*left[4]+0.0337129*above[5]+0.0213998*left[5]+0.0396431*above[6]+0.021478*left[6]+0.0460064*above[7]+0.0208584*left[7]+0.05278*above[8]+0.0197128*left[8]+0.0597798*above[9]+0.0182111*left[9]+0.0665957*above[10]+0.01651*left[10]+0.0725652*above[11]+0.0147474*left[11]+0.0768642*above[12]+0.0130424*left[12]+0.0789466*above[13]+0.0115265*left[13]+0.0521492*above[15]+0.00653262*left[15]+0.0787322*above[14]+0.0103369*left[14]);
		temp[15] = (int)(0.00204088*above[-1]+0.00595289*above[0]+0.00581272*left[0]+0.0116141*above[1]+0.0109415*left[1]+0.016934*above[2]+0.0150061*left[2]+0.0221019*above[3]+0.0179549*left[3]+0.027321*above[4]+0.0198363*left[4]+0.032776*above[5]+0.020752*left[5]+0.0386228*above[6]+0.0208383*left[6]+0.044973*above[7]+0.0202526*left[7]+0.0518613*above[8]+0.0191605*left[8]+0.0591859*above[9]+0.0177241*left[9]+0.0666236*above[10]+0.016093*left[10]+0.0735496*above[11]+0.0143988*left[11]+0.0790544*above[12]+0.0127551*left[12]+0.0823577*above[13]+0.0112894*left[13]+0.0552925*above[15]+0.00640866*left[15]+0.0830485*above[14]+0.0101356*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00124197*above[-1]+0.00338791*above[0]+0.00379662*left[0]+0.00599897*above[1]+0.00799194*left[1]+0.00743362*above[2]+0.0133754*left[2]+0.00774936*above[3]+0.0216631*left[3]+0.00724038*above[4]+0.0371491*left[4]+0.00626976*above[5]+0.0723613*left[5]+0.00514688*above[6]+0.22671*left[6]+0.00407381*above[7]+0.240419*left[7]+0.0031499*above[8]+0.211991*left[8]+0.0024034*above[9]+0.0530745*left[9]+0.00182449*above[10]+0.0230432*left[10]+0.00138827*above[11]+0.0116053*left[11]+0.0010678*above[12]+0.00659256*left[12]+0.000841778*above[13]+0.00413464*left[13]+0.000424306*above[15]+0.00165671*left[15]+0.000696133*above[14]+0.00291839*left[14]);
		temp[1] = (int)(0.00244329*above[-1]+0.00670381*above[0]+0.00741815*left[0]+0.0119631*above[1]+0.0154298*left[1]+0.0150085*above[2]+0.0252174*left[2]+0.0158935*above[3]+0.0391528*left[3]+0.0151086*above[4]+0.0617633*left[4]+0.0133077*above[5]+0.105147*left[5]+0.0110911*above[6]+0.145137*left[6]+0.00888787*above[7]+0.164134*left[7]+0.00693671*above[8]+0.130695*left[8]+0.00532834*above[9]+0.083409*left[9]+0.00406363*above[10]+0.0403401*left[10]+0.0031017*above[11]+0.0221419*left[11]+0.00239057*above[12]+0.0131285*left[12]+0.00188702*above[13]+0.00846198*left[13]+0.000952168*above[15]+0.00347222*left[15]+0.0015617*above[14]+0.00607322*left[14]);
		temp[2] = (int)(0.00342044*above[-1]+0.0094632*above[0]+0.0102856*left[0]+0.0170774*above[1]+0.0210419*left[1]+0.0218103*above[2]+0.0332966*left[2]+0.0236271*above[3]+0.0488072*left[3]+0.0230344*above[4]+0.0697842*left[4]+0.0208073*above[5]+0.0935808*left[5]+0.0177456*above[6]+0.112983*left[6]+0.0144978*above[7]+0.115956*left[7]+0.0114863*above[8]+0.100682*left[8]+0.00892007*above[9]+0.0733312*left[9]+0.00685455*above[10]+0.0474833*left[10]+0.0052584*above[11]+0.0289484*left[11]+0.00406594*above[12]+0.0183548*left[12]+0.00321585*above[13]+0.0123338*left[13]+0.00162516*above[15]+0.00525303*left[15]+0.00266439*above[14]+0.00908611*left[14]);
		temp[3] = (int)(0.00407126*above[-1]+0.0113727*above[0]+0.0121103*left[0]+0.0207952*above[1]+0.024327*left[1]+0.027124*above[2]+0.0371854*left[2]+0.0301909*above[3]+0.0514831*left[3]+0.0303493*above[4]+0.0669496*left[4]+0.0282908*above[5]+0.0814031*left[5]+0.0248556*above[6]+0.0904878*left[6]+0.0208414*above[7]+0.0903737*left[7]+0.0168649*above[8]+0.080223*left[8]+0.0133089*above[9]+0.0637377*left[9]+0.010345*above[10]+0.0462082*left[10]+0.00799795*above[11]+0.0318036*left[11]+0.00621544*above[12]+0.0216872*left[12]+0.00493123*above[13]+0.0153146*left[13]+0.00249803*above[15]+0.00682525*left[15]+0.00409271*above[14]+0.0116451*left[14]);
		temp[4] = (int)(0.00438708*above[-1]+0.0123742*above[0]+0.0129125*left[0]+0.0229349*above[1]+0.0254981*left[1]+0.0305782*above[2]+0.0377903*left[2]+0.0350272*above[3]+0.0499119*left[3]+0.0364021*above[4]+0.0612738*left[4]+0.0351513*above[5]+0.0703146*left[5]+0.0319737*above[6]+0.0749569*left[6]+0.0276776*above[7]+0.0736352*left[7]+0.0230165*above[8]+0.0664477*left[8]+0.0185645*above[9]+0.0552901*left[9]+0.0146682*above[10]+0.0429894*left[10]+0.0114718*above[11]+0.0319165*left[11]+0.00898397*above[12]+0.0232883*left[12]+0.00716246*above[13]+0.0172754*left[13]+0.00364229*above[15]+0.00807385*left[15]+0.00596101*above[14]+0.0135774*left[14]);
		temp[5] = (int)(0.00442633*above[-1]+0.0125941*above[0]+0.0129097*left[0]+0.0236367*above[1]+0.0251352*left[1]+0.0321718*above[2]+0.0363634*left[2]+0.0378905*above[3]+0.0464862*left[3]+0.0407084*above[4]+0.0549985*left[4]+0.040775*above[5]+0.0609895*left[5]+0.0385105*above[6]+0.0634583*left[6]+0.0345669*above[7]+0.0618399*left[7]+0.0297043*above[8]+0.0564252*left[8]+0.0246343*above[9]+0.0483883*left[9]+0.0198981*above[10]+0.0393517*left[10]+0.0158191*above[11]+0.0307958*left[11]+0.0125308*above[12]+0.0236569*left[12]+0.010065*above[13]+0.0183292*left[13]+0.00514965*above[15]+0.00895582*left[15]+0.0084134*above[14]+0.0148549*left[14]);
		temp[6] = (int)(0.00427437*above[-1]+0.0122484*above[0]+0.0123779*left[0]+0.0232321*above[1]+0.0238467*left[1]+0.0321911*above[2]+0.0339128*left[2]+0.0388668*above[3]+0.0424109*left[3]+0.0430683*above[4]+0.0490077*left[4]+0.0446991*above[5]+0.0532178*left[5]+0.0438622*above[6]+0.0546012*left[6]+0.0409253*above[7]+0.0530178*left[7]+0.0364929*above[8]+0.048803*left[8]+0.0312863*above[9]+0.042746*left[9]+0.0259886*above[10]+0.0358732*left[10]+0.0211221*above[11]+0.0291539*left[11]+0.0170056*above[12]+0.0232877*left[12]+0.0138119*above[13]+0.0186896*left[13]+0.00713368*above[15]+0.00948863*left[15]+0.0116233*above[14]+0.0155505*left[14]);
		temp[7] = (int)(0.00401205*above[-1]+0.0115584*above[0]+0.0115584*left[0]+0.0221042*above[1]+0.0221042*left[1]+0.0310699*above[2]+0.0310699*left[2]+0.038297*above[3]+0.038297*left[3]+0.0435943*above[4]+0.0435943*left[4]+0.0467364*above[5]+0.0467364*left[5]+0.0475729*above[6]+0.0475729*left[6]+0.0461544*above[7]+0.0461544*left[7]+0.0428045*above[8]+0.0428045*left[8]+0.0380931*above[9]+0.0380931*left[9]+0.0327174*above[10]+0.0327174*left[10]+0.0273465*above[11]+0.0273465*left[11]+0.0225045*above[12]+0.0225045*left[12]+0.0185692*above[13]+0.0185692*left[13]+0.00972654*above[15]+0.00972654*left[15]+0.0157831*above[14]+0.0157831*left[14]);
		temp[8] = (int)(0.00370121*above[-1]+0.0107033*above[0]+0.0106251*left[0]+0.0205917*above[1]+0.020219*left[1]+0.029256*above[2]+0.0282015*left[2]+0.0366471*above[3]+0.0344378*left[3]+0.0426432*above[4]+0.0388344*left[4]+0.0470082*above[5]+0.0413166*left[5]+0.0494621*above[6]+0.041884*left[6]+0.0498056*above[7]+0.0406698*left[7]+0.0480494*above[8]+0.0379697*left[8]+0.0444875*above[9]+0.0342174*left[9]+0.03967*above[10]+0.0299148*left[10]+0.034285*above[11]+0.0255459*left[11]+0.0290098*above[12]+0.0215104*left[12]+0.0244505*above[13]+0.0181373*left[13]+0.0130649*above[15]+0.00973773*left[15]+0.0210749*above[14]+0.0156774*left[14]);
		temp[9] = (int)(0.00338324*above[-1]+0.00980877*above[0]+0.00968957*left[0]+0.018949*above[1]+0.0183792*left[1]+0.0271271*above[2]+0.0255086*left[2]+0.0343847*above[3]+0.0309667*left[3]+0.040692*above[4]+0.0347215*left[4]+0.0458795*above[5]+0.0367828*left[5]+0.0496598*above[6]+0.0372251*left[6]+0.0517086*above[7]+0.0362138*left[7]+0.0517917*above[8]+0.0340133*left[8]+0.0498948*above[9]+0.0309668*left[9]+0.0462962*above[10]+0.0274538*left[10]+0.0415409*above[11]+0.0238384*left[11]+0.0363309*above[12]+0.0204329*left[12]+0.0314502*above[13]+0.0175208*left[13]+0.0172567*above[15]+0.00958924*left[15]+0.0276133*above[14]+0.015344*left[14]);
		temp[10] = (int)(0.00308357*above[-1]+0.00895482*above[0]+0.00881805*left[0]+0.0173466*above[1]+0.0166921*left[1]+0.0249609*above[2]+0.0230956*left[2]+0.0319009*above[3]+0.0279377*left[3]+0.0382185*above[4]+0.0312229*left[4]+0.0438379*above[5]+0.0330059*left[5]+0.0485433*above[6]+0.0333966*left[6]+0.0520082*above[7]+0.032568*left[7]+0.0538783*above[8]+0.0307549*left[8]+0.0538984*above[9]+0.0282379*left[9]+0.0520413*above[10]+0.025315*left[10]+0.0485821*above[11]+0.0222709*left[11]+0.0440824*above[12]+0.0193555*left[12]+0.0393737*above[13]+0.0168145*left[13]+0.0223219*above[15]+0.00933997*left[15]+0.035359*above[14]+0.0148748*left[14]);
		temp[11] = (int)(0.00281716*above[-1]+0.00818971*above[0]+0.00804871*left[0]+0.0158923*above[1]+0.0152168*left[1]+0.0229445*above[2]+0.0210157*left[2]+0.0294859*above[3]+0.0253698*left[3]+0.0356262*above[4]+0.0283044*left[4]+0.0413756*above[5]+0.0298968*left[5]+0.0466174*above[6]+0.030271*left[6]+0.0511018*above[7]+0.0295964*left[7]+0.0544763*above[8]+0.0280823*left[8]+0.0563679*above[9]+0.0259642*left[9]+0.0565085*above[10]+0.0234847*left[10]+0.0548653*above[11]+0.0208737*left[11]+0.0517248*above[12]+0.0183369*left[12]+0.0477932*above[13]+0.0160899*left[13]+0.0281102*above[15]+0.00903906*left[15]+0.0440171*above[14]+0.0143439*left[14]);
		temp[12] = (int)(0.00259308*above[-1]+0.0075431*above[0]+0.00740438*left[0]+0.0146534*above[1]+0.0139884*left[1]+0.0212004*above[2]+0.0192988*left[2]+0.0273417*above[3]+0.0232715*left[3]+0.0332226*above[4]+0.0259431*left[4]+0.0389166*above[5]+0.0274009*left[5]+0.044398*above[6]+0.0277726*left[6]+0.04952*above[7]+0.02722*left[7]+0.0540091*above[8]+0.0259313*left[8]+0.0574926*above[9]+0.0241097*left[9]+0.0595792*above[10]+0.0219593*left[10]+0.0599829*above[11]+0.0196723*left[11]+0.0586643*above[12]+0.0174229*left[12]+0.056055*above[13]+0.0154035*left[13]+0.0342216*above[15]+0.0087268*left[15]+0.0529592*above[14]+0.0138115*left[14]);
		temp[13] = (int)(0.0024178*above[-1]+0.00703594*above[0]+0.0069016*left[0]+0.0136773*above[1]+0.013033*left[1]+0.019814*above[2]+0.01797*left[2]+0.0256113*above[3]+0.0216565*left[3]+0.031234*above[4]+0.024135*left[4]+0.0367959*above[5]+0.0254966*left[5]+0.0423387*above[6]+0.0258682*left[6]+0.0478083*above[7]+0.0254038*left[7]+0.053032*above[8]+0.0242756*left[8]+0.0577078*above[9]+0.0226644*left[9]+0.0614283*above[10]+0.0207482*left[10]+0.0637568*above[11]+0.0186937*left[11]+0.0643649*above[12]+0.0166533*left[12]+0.0633459*above[13]+0.0148025*left[13]+0.0399752*above[15]+0.00843655*left[15]+0.0612174*above[14]+0.0133269*left[14]);
		temp[14] = (int)(0.00229761*above[-1]+0.0066878*above[0]+0.00655715*left[0]+0.0130059*above[1]+0.0123793*left[1]+0.0188569*above[2]+0.0170621*left[2]+0.0244085*above[3]+0.0205547*left[3]+0.0298353*above[4]+0.0229021*left[4]+0.0352743*above[5]+0.0241968*left[5]+0.0408087*above[6]+0.0245641*left[6]+0.0464483*above[7]+0.0241522*left[7]+0.0521025*above[8]+0.0231232*left[8]+0.0575499*above[9]+0.0216438*left[9]+0.0624195*above[10]+0.0198761*left[10]+0.0662126*above[11]+0.017971*left[11]+0.0684033*above[12]+0.0160672*left[12]+0.0687914*above[13]+0.0143288*left[13]+0.0444754*above[15]+0.00819657*left[15]+0.0675897*above[14]+0.0129328*left[14]);
		temp[15] = (int)(0.00224033*above[-1]+0.00652213*above[0]+0.00639273*left[0]+0.0126872*above[1]+0.0120664*left[1]+0.0184041*above[2]+0.0166257*left[2]+0.0238423*above[3]+0.0200215*left[3]+0.0291804*above[4]+0.0223*left[4]+0.0345656*above[5]+0.0235542*left[5]+0.0400993*above[6]+0.023909*left[6]+0.0458202*above[7]+0.0235108*left[7]+0.0516767*above[8]+0.0225177*left[8]+0.0574896*above[9]+0.0210908*left[9]+0.0629145*above[10]+0.019386*left[10]+0.0674337*above[11]+0.0175471*left[11]+0.0704309*above[12]+0.0157066*left[12]+0.0715559*above[13]+0.0140226*left[13]+0.0467887*above[15]+0.00803147*left[15]+0.0708525*above[14]+0.0126673*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.000921862*above[-1]+0.00253851*above[0]+0.00278816*left[0]+0.00455373*above[1]+0.00576475*left[1]+0.00576458*above[2]+0.00932795*left[2]+0.00618339*above[3]+0.0142903*left[3]+0.00597601*above[4]+0.0223026*left[4]+0.0053687*above[5]+0.0376084*left[5]+0.0045747*above[6]+0.0727033*left[6]+0.00375293*above[7]+0.226977*left[7]+0.00299904*above[8]+0.240642*left[8]+0.00235683*above[9]+0.212193*left[9]+0.00183614*above[10]+0.0532752*left[10]+0.00142888*above[11]+0.023264*left[11]+0.00112016*above[12]+0.0118739*left[12]+0.000896763*above[13]+0.0069629*left[13]+0.000459232*above[15]+0.0026167*left[15]+0.000749752*above[14]+0.00470101*left[14]);
		temp[1] = (int)(0.00184535*above[-1]+0.00510215*above[0]+0.0055554*left[0]+0.00920319*above[1]+0.0113957*left[1]+0.011754*above[2]+0.0181617*left[2]+0.0127527*above[3]+0.0271064*left[3]+0.0124849*above[4]+0.0404837*left[4]+0.0113648*above[5]+0.0627256*left[5]+0.00980438*above[6]+0.105868*left[6]+0.00813026*above[7]+0.145703*left[7]+0.00655462*above[8]+0.164608*left[8]+0.00518656*above[9]+0.131124*left[9]+0.00406167*above[10]+0.0838356*left[10]+0.00317283*above[11]+0.0408077*left[11]+0.00249411*above[12]+0.0227073*left[12]+0.0020005*above[13]+0.0139013*left[13]+0.00102618*above[15]+0.00541912*left[15]+0.00167454*above[14]+0.00963191*left[14]);
		temp[2] = (int)(0.00265119*above[-1]+0.00737327*above[0]+0.0079287*left[0]+0.0134069*above[1]+0.0160835*left[1]+0.0173449*above[2]+0.0250941*left[2]+0.0191351*above[3]+0.0361388*left[3]+0.0190915*above[4]+0.0508335*left[4]+0.0177218*above[5]+0.0712637*left[5]+0.0155762*above[6]+0.0946974*left[6]+0.0131321*above[7]+0.113865*left[7]+0.0107342*above[8]+0.116696*left[8]+0.00858701*above[9]+0.101353*left[9]+0.00678042*above[10]+0.0739954*left[10]+0.00532876*above[11]+0.0482053*left[11]+0.00420692*above[12]+0.0298111*left[12]+0.00338427*above[13]+0.019516*left[13]+0.00174048*above[15]+0.00804256*left[15]+0.00283801*above[14]+0.0140611*left[14]);
		temp[3] = (int)(0.00325892*above[-1]+0.00912642*above[0]+0.00967144*left[0]+0.0167546*above[1]+0.0193703*left[1]+0.0220144*above[2]+0.0295132*left[2]+0.0247819*above[3]+0.0408891*left[3]+0.025307*above[4]+0.0541632*left[4]+0.024073*above[5]+0.0689305*left[5]+0.02167*above[6]+0.0829125*left[6]+0.0186731*above[7]+0.091688*left[7]+0.0155532*above[8]+0.0913842*left[8]+0.0126341*above[9]+0.0811371*left[9]+0.0100953*above[10]+0.0646364*left[10]+0.0080045*above[11]+0.047173*left[11]+0.0063596*above[12]+0.0329368*left[12]+0.00513836*above[13]+0.0231804*left[13]+0.00265265*above[15]+0.0101966*left[15]+0.00432059*above[14]+0.0174839*left[14]);
		temp[4] = (int)(0.00363993*above[-1]+0.0102673*above[0]+0.0107179*left[0]+0.0190413*above[1]+0.0211966*left[1]+0.0254361*above[2]+0.0315644*left[2]+0.0292663*above[3]+0.0422131*left[3]+0.0306621*above[4]+0.0531677*left[4]+0.029984*above[5]+0.0637142*left[5]+0.0277518*above[6]+0.0721946*left[6]+0.0245517*above[7]+0.0764626*left[7]+0.0209366*above[8]+0.0749062*left[8]+0.0173494*above[9]+0.067594*left[9]+0.014087*above[10]+0.0564062*left[10]+0.0113075*above[11]+0.0441686*left[11]+0.00906517*above[12]+0.0332721*left[12]+0.00737046*above[13]+0.0250287*left[13]+0.00382605*above[15]+0.0117356*left[15]+0.00622174*above[14]+0.0197339*left[14]);
		temp[5] = (int)(0.0038106*above[-1]+0.0108221*above[0]+0.0111401*left[0]+0.0202676*above[1]+0.0217847*left[1]+0.0275149*above[2]+0.0318069*left[2]+0.0323537*above[3]+0.0413369*left[3]+0.0347937*above[4]+0.050216*left[4]+0.0350245*above[5]+0.057837*left[5]+0.0334095*above[6]+0.0632021*left[6]+0.0304464*above[7]+0.0652439*left[7]+0.0266904*above[8]+0.0633506*left[8]+0.0226635*above[9]+0.0577817*left[9]+0.0187822*above[10]+0.0496937*left[10]+0.0153241*above[11]+0.0407055*left[11]+0.0124377*above[12]+0.0323138*left[12]+0.0102019*above[13]+0.0255494*left[13]+0.00533775*above[15]+0.0126581*left[15]+0.00865991*above[14]+0.0209221*left[14]);
		temp[6] = (int)(0.00381399*above[-1]+0.0108962*above[0]+0.0110825*left[0]+0.0205849*above[1]+0.0214731*left[1]+0.0283573*above[2]+0.0308664*left[2]+0.0340237*above[3]+0.039266*left[3]+0.0375138*above[4]+0.0465034*left[4]+0.0388584*above[5]+0.0521719*left[5]+0.0382284*above[6]+0.0557141*left[6]+0.0359546*above[7]+0.056631*left[7]+0.0324999*above[8]+0.0547382*left[8]+0.0283885*above[9]+0.0503398*left[9]+0.0241181*above[10]+0.0442059*left[10]+0.0200883*above[11]+0.0373563*left[11]+0.0165717*above[12]+0.0307725*left[12]+0.0137562*above[13]+0.0252433*left[13]+0.00727687*above[15]+0.013058*left[15]+0.0117676*above[14]+0.0212861*left[14]);
		temp[7] = (int)(0.00370121*above[-1]+0.0106251*above[0]+0.0107033*left[0]+0.020219*above[1]+0.0205917*left[1]+0.0282015*above[2]+0.029256*left[2]+0.0344378*above[3]+0.0366471*left[3]+0.0388344*above[4]+0.0426432*left[4]+0.0413166*above[5]+0.0470082*left[5]+0.041884*above[6]+0.0494621*left[6]+0.0406698*above[7]+0.0498056*left[7]+0.0379697*above[8]+0.0480494*left[8]+0.0342174*above[9]+0.0444875*left[9]+0.0299148*above[10]+0.03967*left[10]+0.0255459*above[11]+0.034285*left[11]+0.0215104*above[12]+0.0290098*left[12]+0.0181373*above[13]+0.0244505*left[13]+0.00973773*above[15]+0.0130649*left[15]+0.0156774*above[14]+0.0210749*left[14]);
		temp[8] = (int)(0.0035189*above[-1]+0.0101394*above[0]+0.0101394*left[0]+0.0194053*above[1]+0.0194053*left[1]+0.0273383*above[2]+0.0273383*left[2]+0.0338732*above[3]+0.0338732*left[3]+0.03894*above[4]+0.03894*left[4]+0.0424253*above[5]+0.0424253*left[5]+0.044217*above[6]+0.044217*left[6]+0.0442761*above[7]+0.0442761*left[7]+0.0427023*above[8]+0.0427023*left[8]+0.0397639*above[9]+0.0397639*left[9]+0.035875*above[10]+0.035875*left[10]+0.0315271*above[11]+0.0315271*left[11]+0.0272087*above[12]+0.0272087*left[12]+0.0233949*above[13]+0.0233949*left[13]+0.012804*above[15]+0.012804*left[15]+0.0204944*above[14]+0.0204944*left[14]);
		temp[9] = (int)(0.0033039*above[-1]+0.00954612*above[0]+0.00949495*left[0]+0.0183487*above[1]+0.0181041*left[1]+0.0260488*above[2]+0.0253531*left[2]+0.0326479*above[3]+0.0311759*left[3]+0.0381247*above[4]+0.0355441*left[4]+0.042381*above[5]+0.0384253*left[5]+0.0452642*above[6]+0.0398069*left[6]+0.0466248*above[7]+0.0397351*left[7]+0.0463923*above[8]+0.0383479*left[8]+0.0446415*above[9]+0.0358856*left[9]+0.0416236*above[10]+0.0326732*left[10]+0.0377439*above[11]+0.0290785*left[11]+0.0335018*above[12]+0.0254681*left[12]+0.0294813*above[13]+0.0222241*left[13]+0.0165202*above[15]+0.0123792*left[15]+0.0262519*above[14]+0.0197021*left[14]);
		temp[10] = (int)(0.00308295*above[-1]+0.00892536*above[0]+0.00884351*left[0]+0.0172095*above[1]+0.0168178*left[1]+0.0245705*above[2]+0.023454*left[2]+0.0310637*above[3]+0.0286907*left[3]+0.0367214*above[4]+0.0325289*left[4]+0.0414904*above[5]+0.0349871*left[5]+0.0452342*above[6]+0.0361098*left[6]+0.0477649*above[7]+0.0359874*left[7]+0.0489039*above[8]+0.0347716*left[8]+0.0485576*above[9]+0.0326769*left[9]+0.0467842*above[10]+0.0299669*left[10]+0.0438261*above[11]+0.0269285*left[11]+0.0400968*above[12]+0.0238473*left[12]+0.0362099*above[13]+0.0210388*left[13]+0.0208494*above[15]+0.0118713*left[15]+0.0328533*above[14]+0.0188153*left[14]);
		temp[11] = (int)(0.00287471*above[-1]+0.00833406*above[0]+0.00823547*left[0]+0.0161053*above[1]+0.015633*left[1]+0.023088*above[2]+0.0217397*left[2]+0.0293759*above[3]+0.0264997*left[3]+0.0350467*above[4]+0.0299336*left[4]+0.0401007*above[5]+0.03209*left[5]+0.0444484*above[6]+0.0330454*left[6]+0.0479201*above[7]+0.032912*left[7]+0.0503002*above[8]+0.0318431*left[8]+0.0513873*above[9]+0.0300312*left[9]+0.0510689*above[10]+0.0276968*left[10]+0.0493903*above[11]+0.025072*left[11]+0.0465969*above[12]+0.0223876*left[12]+0.0432306*above[13]+0.0199113*left[13]+0.0256243*above[15]+0.0113414*left[15]+0.0400129*above[14]+0.0179212*left[14]);
		temp[12] = (int)(0.00269231*above[-1]+0.00781281*above[0]+0.00770606*left[0]+0.0151216*above[1]+0.01461*left[1]+0.0217401*above[2]+0.0202776*left[2]+0.0277867*above[3]+0.0246587*left[3]+0.0333729*above[4]+0.0277858*left[4]+0.0385459*above[5]+0.0297248*left[5]+0.0432713*above[6]+0.0305692*left[6]+0.0474278*above[7]+0.0304407*left[7]+0.0508175*above[8]+0.0294897*left[8]+0.0531985*above[9]+0.0278906*left[9]+0.0543436*above[10]+0.0258337*left[10]+0.0541166*above[11]+0.0235137*left[11]+0.0525534*above[12]+0.0211237*left[12]+0.0500364*above[13]+0.0188976*left[13]+0.0305077*above[15]+0.010837*left[15]+0.047218*above[14]+0.0170873*left[14]);
		temp[13] = (int)(0.00254578*above[-1]+0.00739245*above[0]+0.00728224*left[0]+0.0143232*above[1]+0.0137949*left[1]+0.0206331*above[2]+0.0191212*left[2]+0.0264549*above[3]+0.0232153*left[3]+0.0319225*above[4]+0.0261166*left[4]+0.0371181*above[5]+0.0279005*left[5]+0.0420552*above[6]+0.0286692*left[6]+0.0466663*above[7]+0.0285481*left[7]+0.0507974*above[8]+0.0276831*left[8]+0.0542152*above[9]+0.0262354*left[9]+0.0566379*above[10]+0.0243745*left[10]+0.0577951*above[11]+0.0222703*left[11]+0.0575226*above[12]+0.0200907*left[12]+0.0560016*above[13]+0.0180458*left[13]+0.0349885*above[15]+0.0103963*left[15]+0.0537393*above[14]+0.0163684*left[14]);
		temp[14] = (int)(0.0024439*above[-1]+0.00709969*above[0]+0.00698803*left[0]+0.0137657*above[1]+0.0132302*left[1]+0.0198556*above[2]+0.0183226*left[2]+0.0255105*above[3]+0.0222219*left[3]+0.0308774*above[4]+0.0249711*left[4]+0.0360609*above[5]+0.0266509*left[5]+0.0411075*above[6]+0.0273674*left[6]+0.0459918*above[7]+0.0272476*left[7]+0.0506026*above[8]+0.0264339*left[8]+0.0547331*above[9]+0.0250793*left[9]+0.0580864*above[10]+0.0233407*left[10]+0.0603126*above[11]+0.0213728*left[11]+0.0610975*above[12]+0.0193282*left[12]+0.0604447*above[13]+0.0174019*left[13]+0.0384334*above[15]+0.0100525*left[15]+0.0587061*above[14]+0.0158134*left[14]);
		temp[15] = (int)(0.00239568*above[-1]+0.0069613*above[0]+0.00684856*left[0]+0.0135026*above[1]+0.0129619*left[1]+0.01949*above[2]+0.0179416*left[2]+0.0250685*above[3]+0.0217453*left[3]+0.0303911*above[4]+0.0244172*left[4]+0.0355725*above[5]+0.0260403*left[5]+0.0406744*above[6]+0.0267229*left[6]+0.0456911*above[7]+0.0265929*left[7]+0.050533*above[8]+0.0257921*left[8]+0.0550088*above[9]+0.0244706*left[9]+0.0588193*above[10]+0.0227803*left[10]+0.0615795*above[11]+0.02087*left[11]+0.0629032*above[12]+0.0188855*left[12]+0.0627016*above[13]+0.0170145*left[13]+0.0401957*above[15]+0.00983654*left[15]+0.0612411*above[14]+0.0154697*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.000707116*above[-1]+0.00196106*above[0]+0.00212187*left[0]+0.00355299*above[1]+0.00433048*left[1]+0.00457209*above[2]+0.00684161*left[2]+0.00501376*above[3]+0.0100868*left[3]+0.0049758*above[4]+0.0148377*left[4]+0.00460343*above[5]+0.0227109*left[5]+0.00404449*above[6]+0.0379274*left[6]+0.00342002*above[7]+0.0729686*left[7]+0.00281294*above[8]+0.227216*left[8]+0.00227018*above[9]+0.240877*left[9]+0.00181158*above[10]+0.212448*left[10]+0.00143994*above[11]+0.0535783*left[11]+0.00114953*above[12]+0.0236574*left[12]+0.000933992*above[13]+0.0124466*left[13]+0.00048563*above[15]+0.00426178*left[15]+0.000789127*above[14]+0.00787864*left[14]);
		temp[1] = (int)(0.00143446*above[-1]+0.00398984*above[0]+0.00429039*left[0]+0.00725781*above[1]+0.00870852*left[1]+0.00940058*above[2]+0.013618*left[2]+0.0103967*above[3]+0.01974*left[3]+0.0104195*above[4]+0.0282526*left[4]+0.00973943*above[5]+0.0413437*left[5]+0.00864312*above[6]+0.0634008*left[6]+0.00737599*above[7]+0.106431*left[7]+0.00611513*above[8]+0.146211*left[8]+0.00496786*above[9]+0.165108*left[9]+0.00398531*above[10]+0.131664*left[10]+0.00318084*above[11]+0.0844734*left[11]+0.0025473*above[12]+0.0416284*left[12]+0.00207443*above[13]+0.0238884*left[13]+0.0010809*above[15]+0.00867656*left[15]+0.0017553*above[14]+0.0157636*left[14]);
		temp[2] = (int)(0.00210237*above[-1]+0.00587248*above[0]+0.00625848*left[0]+0.0107456*above[1]+0.0126044*left[1]+0.0140515*above[2]+0.0194251*left[2]+0.0157359*above[3]+0.0274954*left[3]+0.0160002*above[4]+0.0378998*left[4]+0.0151867*above[5]+0.0521652*left[5]+0.0136821*above[6]+0.0723154*left[6]+0.0118406*above[7]+0.095578*left[7]+0.00993771*above[8]+0.114658*left[8]+0.00815672*above[9]+0.117475*left[9]+0.00659797*above[10]+0.102187*left[10]+0.00530036*above[11]+0.0749683*left[11]+0.00426556*above[12]+0.0494377*left[12]+0.00348611*above[13]+0.0315498*left[13]+0.00182246*above[15]+0.0125323*left[15]+0.00295664*above[14]+0.0221947*left[14]);
		temp[3] = (int)(0.00264983*above[-1]+0.00743961*above[0]+0.00784413*left[0]+0.0137108*above[1]+0.0156539*left[1]+0.0181388*above[2]+0.0237217*left[2]+0.0206269*above[3]+0.0326866*left[3]+0.0213526*above[4]+0.0432448*left[4]+0.0206607*above[5]+0.055962*left[5]+0.0189767*above[6]+0.0703608*left[6]+0.0167253*above[7]+0.0841138*left[7]+0.0142695*above[8]+0.092769*left[8]+0.0118781*above[9]+0.0924371*left[9]+0.00971993*above[10]+0.0822513*left[10]+0.00788026*above[11]+0.0659139*left[11]+0.00638642*above[12]+0.0487562*left[12]+0.00524616*above[13]+0.0351115*left[13]+0.00275555*above[15]+0.0153619*left[15]+0.00446414*above[14]+0.0264303*left[14]);
		temp[4] = (int)(0.00304489*above[-1]+0.0085956*above[0]+0.00896065*left[0]+0.0159643*above[1]+0.0177138*left[1]+0.0213893*above[2]+0.0263889*left[2]+0.0247374*above[3]+0.0354159*left[3]+0.0261255*above[4]+0.0451125*left[4]+0.0258387*above[5]+0.055406*left[5]+0.024272*above[6]+0.0655072*left[6]+0.0218645*above[7]+0.0737047*left[7]+0.0190343*above[8]+0.0778175*left[8]+0.0161288*above[9]+0.0762136*left[9]+0.0133983*above[10]+0.0689556*left[10]+0.0109958*above[11]+0.0579338*left[11]+0.00899645*above[12]+0.0460115*left[12]+0.00744213*above[13]+0.0357237*left[13]+0.00393453*above[15]+0.0170481*left[15]+0.00636179*above[14]+0.0285658*left[14]);
		temp[5] = (int)(0.00328462*above[-1]+0.00932212*above[0]+0.00961168*left[0]+0.0174473*above[1]+0.0188327*left[1]+0.0236755*above[2]+0.0276192*left[2]+0.0278559*above[3]+0.0362152*left[3]+0.030036*above[4]+0.0447062*left[4]+0.0304043*above[5]+0.0528478*left[5]+0.0292681*above[6]+0.0599614*left[6]+0.027018*above[7]+0.0649955*left[7]+0.0240759*above[8]+0.0668464*left[8]+0.0208393*above[9]+0.0648793*left[9]+0.0176357*above[10]+0.0593444*left[10]+0.0146995*above[11]+0.0514036*left[11]+0.0121763*above[12]+0.0427057*left[12]+0.0101665*above[13]+0.0348823*left[13]+0.00542208*above[15]+0.0177487*left[15]+0.0087441*above[14]+0.0291204*left[14]);
		temp[6] = (int)(0.00338752*above[-1]+0.00966119*above[0]+0.00986293*left[0]+0.0182114*above[1]+0.0191757*left[1]+0.0250091*above[2]+0.027749*left[2]+0.0299123*above[3]+0.0356994*left[3]+0.032914*above[4]+0.0430183*left[4]+0.0341021*above[5]+0.04947*left[5]+0.0336664*above[6]+0.0545847*left[6]+0.0318973*above[7]+0.0577551*left[7]+0.0291624*above[8]+0.0584456*left[8]+0.0258603*above[9]+0.0564472*left[9]+0.0223684*above[10]+0.0520514*left[10]+0.0189996*above[11]+0.0460279*left[11]+0.015984*above[12]+0.0394182*left[12]+0.0135055*above[13]+0.0333243*left[13]+0.00728561*above[15]+0.0177321*left[15]+0.0117091*above[14]+0.0286645*left[14]);
		temp[7] = (int)(0.00338324*above[-1]+0.00968957*above[0]+0.00980877*left[0]+0.0183792*above[1]+0.018949*left[1]+0.0255086*above[2]+0.0271271*left[2]+0.0309667*above[3]+0.0343847*left[3]+0.0347215*above[4]+0.040692*left[4]+0.0367828*above[5]+0.0458795*left[5]+0.0372251*above[6]+0.0496598*left[6]+0.0362138*above[7]+0.0517086*left[7]+0.0340133*above[8]+0.0517917*left[8]+0.0309668*above[9]+0.0498948*left[9]+0.0274538*above[10]+0.0462962*left[10]+0.0238384*above[11]+0.0415409*left[11]+0.0204329*above[12]+0.0363309*left[12]+0.0175208*above[13]+0.0314502*left[13]+0.00958924*above[15]+0.0172567*left[15]+0.015344*above[14]+0.0276133*left[14]);
		temp[8] = (int)(0.0033039*above[-1]+0.00949495*above[0]+0.00954612*left[0]+0.0181041*above[1]+0.0183487*left[1]+0.0253531*above[2]+0.0260488*left[2]+0.0311759*above[3]+0.0326479*left[3]+0.0355441*above[4]+0.0381247*left[4]+0.0384253*above[5]+0.042381*left[5]+0.0398069*above[6]+0.0452642*left[6]+0.0397351*above[7]+0.0466248*left[7]+0.0383479*above[8]+0.0463923*left[8]+0.0358856*above[9]+0.0446415*left[9]+0.0326732*above[10]+0.0416236*left[10]+0.0290785*above[11]+0.0377439*left[11]+0.0254681*above[12]+0.0335018*left[12]+0.0222241*above[13]+0.0294813*left[13]+0.0123792*above[15]+0.0165202*left[15]+0.0197021*above[14]+0.0262519*left[14]);
		temp[9] = (int)(0.00317868*above[-1]+0.00915997*above[0]+0.00915997*left[0]+0.0175384*above[1]+0.0175384*left[1]+0.0247408*above[2]+0.0247408*left[2]+0.0307489*above[3]+0.0307489*left[3]+0.0355584*above[4]+0.0355584*left[4]+0.0391294*above[5]+0.0391294*left[5]+0.041402*above[6]+0.041402*left[6]+0.0423335*above[7]+0.0423335*left[7]+0.0419425*above[8]+0.0419425*left[8]+0.0403436*above[9]+0.0403436*left[9]+0.0377592*above[10]+0.0377592*left[10]+0.0345034*above[11]+0.0345034*left[11]+0.0309485*above[12]+0.0309485*left[12]+0.0275484*above[13]+0.0275484*left[13]+0.0156617*above[15]+0.0156617*left[15]+0.0247699*above[14]+0.0247699*left[14]);
		temp[10] = (int)(0.00303164*above[-1]+0.00875456*above[0]+0.00871865*left[0]+0.0168168*above[1]+0.0166449*left[1]+0.0238595*above[2]+0.0233692*left[2]+0.029908*above[3]+0.0288644*left[3]+0.0349908*above[4]+0.0331427*left[4]+0.0390853*above[5]+0.0362069*left[5]+0.0421216*above[6]+0.038058*left[6]+0.0440092*above[7]+0.0387175*left[7]+0.0446779*above[8]+0.0382525*left[8]+0.044123*above[9]+0.0367934*left[9]+0.0424405*above[10]+0.0345371*left[10]+0.039841*above[11]+0.0317359*left[11]+0.0366428*above[12]+0.0286792*left[12]+0.033322*above[13]+0.0257338*left[13]+0.0193737*above[15]+0.0147764*left[15]+0.0304283*above[14]+0.0232939*left[14]);
		temp[11] = (int)(0.00288167*above[-1]+0.00833462*above[0]+0.00827486*left[0]+0.0160497*above[1]+0.0157635*left[1]+0.0228718*above[2]+0.0220544*left[2]+0.0288611*above[3]+0.0271167*left[3]+0.0340792*above[4]+0.0309756*left[4]+0.0385326*above[5]+0.033663*left[5]+0.0421674*above[6]+0.035219*left[6]+0.044883*above[7]+0.035703*left[7]+0.046561*above[8]+0.035208*left[8]+0.0471064*above[9]+0.0338686*left[9]+0.0464933*above[10]+0.0318609*left[10]+0.0448015*above[11]+0.0293943*left[11]+0.0422429*above[12]+0.0267031*left[12]+0.0392545*above[13]+0.0240946*left[13]+0.0233545*above[15]+0.0139301*left[15]+0.0364184*above[14]+0.021911*left[14]);
		temp[12] = (int)(0.00274343*above[-1]+0.00794406*above[0]+0.00786916*left[0]+0.0153258*above[1]+0.0149668*left[1]+0.0219124*above[2]+0.0208863*left[2]+0.0277899*above[3]+0.0255957*left[3]+0.0330476*above[4]+0.029129*left[4]+0.0377227*above[5]+0.0315375*left[5]+0.04179*above[6]+0.0328849*left[6]+0.0451639*above[7]+0.0332532*left[7]+0.0477138*above[8]+0.0327486*left[8]+0.0492926*above[9]+0.031505*left[9]+0.0497777*above[10]+0.0296822*left[10]+0.0491193*above[11]+0.0274602*left[11]+0.0473922*above[12]+0.0250361*left[12]+0.0449416*above[13]+0.0226763*left[13]+0.0273262*above[15]+0.0131708*left[15]+0.0423232*above[14]+0.0206856*left[14]);
		temp[13] = (int)(0.00262869*above[-1]+0.00761821*above[0]+0.00753404*left[0]+0.0147167*above[1]+0.0143132*left[1]+0.0210919*above[2]+0.0199377*left[2]+0.0268474*above[3]+0.0243753*left[3]+0.032093*above[4]+0.0276661*left[4]+0.0368923*above[5]+0.0298729*left[5]+0.0412492*above[6]+0.0310741*left[6]+0.0451037*above[7]+0.0313643*left[7]+0.0483354*above[8]+0.0308564*left[8]+0.050777*above[9]+0.0296822*left[9]+0.0522434*above[10]+0.0279897*left[10]+0.0525773*above[11]+0.0259393*left[11]+0.0517173*above[12]+0.0237037*left[12]+0.0498901*above[13]+0.0215215*left[13]+0.0308992*above[15]+0.0125369*left[15]+0.0475826*above[14]+0.0196713*left[14]);
		temp[14] = (int)(0.00254754*above[-1]+0.00738719*above[0]+0.00729756*left[0]+0.0142831*above[1]+0.0138534*left[1]+0.0205033*above[2]+0.0192733*left[2]+0.0261622*above[3]+0.0235252*left[3]+0.0313831*above[4]+0.0266525*left[4]+0.0362474*above[5]+0.0287247*left[5]+0.0407817*above[6]+0.0298284*left[6]+0.044949*above[7]+0.0300652*left[7]+0.048645*above[8]+0.0295516*left[8]+0.0517001*above[9]+0.0284178*left[9]+0.0538957*above[10]+0.0268048*left[10]+0.0550016*above[11]+0.0248609*left[11]+0.0548497*above[12]+0.0227445*left[12]+0.0535612*above[13]+0.0206768*left[13]+0.0336101*above[15]+0.0120638*left[15]+0.0515462*above[14]+0.0189192*left[14]);
		temp[15] = (int)(0.00250934*above[-1]+0.00727859*above[0]+0.00718607*left[0]+0.0140798*above[1]+0.0136361*left[1]+0.0202282*above[2]+0.0189582*left[2]+0.0258438*above[3]+0.02312*left[3]+0.031056*above[4]+0.0261659*left[4]+0.0359547*above[5]+0.0281686*left[5]+0.040577*above[6]+0.0292183*left[6]+0.0448972*above[7]+0.0294204*left[7]+0.0488197*above[8]+0.0288934*left[8]+0.0521757*above[9]+0.0277677*left[9]+0.0547323*above[10]+0.0261819*left[10]+0.0562253*above[11]+0.02428*left[11]+0.0564336*above[12]+0.0222144*left[12]+0.0554233*above[13]+0.0201984*left[13]+0.0349911*above[15]+0.0117881*left[15]+0.0535624*above[14]+0.018485*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.000558143*above[-1]+0.00155644*above[0]+0.00166482*left[0]+0.00284188*above[1]+0.00336477*left[1]+0.00370435*above[2]+0.0052228*left[2]+0.00413363*above[3]+0.0074902*left[3]+0.00419001*above[4]+0.0105724*left[4]+0.00396973*above[5]+0.0152176*left[5]+0.00357685*above[6]+0.0230265*left[6]+0.00310286*above[7]+0.0382103*left[7]+0.00261644*above[8]+0.0732455*left[8]+0.00216187*above[9]+0.227513*left[9]+0.0017629*above[10]+0.241226*left[10]+0.00142872*above[11]+0.21289*left[11]+0.00115992*above[12]+0.0541847*left[12]+0.000955493*above[13]+0.0245849*left[13]+0.000503896*above[15]+0.00725305*left[15]+0.00081523*above[14]+0.0139985*left[14]);
		temp[1] = (int)(0.00114425*above[-1]+0.00319773*above[0]+0.00340497*left[0]+0.00585623*above[1]+0.00685503*left[1]+0.00767093*above[2]+0.0105641*left[2]+0.00861533*above[3]+0.0149758*left[3]+0.00879906*above[4]+0.0207627*left[4]+0.00840443*above[5]+0.0290565*left[5]+0.00763465*above[6]+0.0420137*left[6]+0.00667437*above[7]+0.0640023*left[7]+0.00566757*above[8]+0.10702*left[8]+0.00471147*above[9]+0.14684*left[9]+0.00386168*above[10]+0.165841*left[10]+0.00314275*above[11]+0.132586*left[11]+0.00255996*above[12]+0.0857233*left[12]+0.00211408*above[13]+0.0435106*left[13]+0.00111758*above[15]+0.0143884*left[15]+0.00180677*above[14]+0.0269761*left[14]);
		temp[2] = (int)(0.00170372*above[-1]+0.00477628*above[0]+0.00505233*left[0]+0.00878588*above[1]+0.0101144*left[1]+0.0115917*above[2]+0.0154263*left[2]+0.0131436*above[3]+0.0215102*left[3]+0.0135752*above[4]+0.0290783*left[4]+0.0131244*above[5]+0.0391516*left[5]+0.0120694*above[6]+0.0532122*left[6]+0.0106759*above[7]+0.0732555*left[7]+0.00916307*above[8]+0.0964941*left[8]+0.0076891*above[9]+0.11563*left[9]+0.00635247*above[10]+0.118593*left[10]+0.00520358*above[11]+0.103571*left[11]+0.00426056*above[12]+0.0768071*left[12]+0.0035322*above[13]+0.0521352*left[13]+0.00187421*above[15]+0.0199699*left[15]+0.0030266*above[14]+0.0358354*left[14]);
		temp[3] = (int)(0.00219073*above[-1]+0.00616524*above[0]+0.00646949*left[0]+0.0114026*above[1]+0.0128644*left[1]+0.0151783*above[2]+0.0193812*left[2]+0.0174151*above[3]+0.0265093*left[3]+0.0182407*above[4]+0.0348235*left[4]+0.017907*above[5]+0.0449464*left[5]+0.016728*above[6]+0.0573898*left[6]+0.0150241*above[7]+0.0716411*left[7]+0.0130792*above[8]+0.0853532*left[8]+0.0111148*above[9]+0.0940679*left[9]+0.0092828*above[10]+0.0939052*left[10]+0.00767272*above[11]+0.0840285*left[11]+0.00632762*above[12]+0.0682109*left[12]+0.00527461*above[13]+0.0520114*left[13]+0.00281338*above[15]+0.0233275*left[15]+0.00453612*above[14]+0.0400749*left[14]);
		temp[4] = (int)(0.00257581*above[-1]+0.00727946*above[0]+0.00757255*left[0]+0.0135441*above[1]+0.0149501*left[1]+0.0182074*above[2]+0.0222351*left[2]+0.0211684*above[3]+0.0298176*left[3]+0.0225257*above[4]+0.0380732*left[4]+0.022505*above[5]+0.0472443*left[5]+0.0214118*above[6]+0.0572001*left[6]+0.0195836*above[7]+0.0671115*left[7]+0.017345*above[8]+0.0752434*left[8]+0.0149733*above[9]+0.0794048*left[9]+0.0126789*above[10]+0.0779688*left[10]+0.0106026*above[11]+0.0710227*left[11]+0.00882689*above[12]+0.0605179*left[12]+0.0074114*above[13]+0.0495289*left[13]+0.00398073*above[15]+0.0246547*left[15]+0.00640481*above[14]+0.0408463*left[14]);
		temp[5] = (int)(0.00284798*above[-1]+0.00808278*above[0]+0.00833544*left[0]+0.0151308*above[1]+0.0163413*left[1]+0.0205469*above[2]+0.0240049*left[2]+0.0242182*above[3]+0.0316008*left[3]+0.0262041*above[4]+0.0393374*left[4]+0.0266769*above[5]+0.0472301*left[5]+0.0258946*above[6]+0.0549773*left[6]+0.0241702*above[7]+0.0618582*left[7]+0.021835*above[8]+0.0667944*left[8]+0.019201*above[9]+0.0686679*left[9]+0.0165312*above[10]+0.0668433*left[10]+0.0140241*above[11]+0.0615859*left[11]+0.0118148*above[12]+0.0541027*left[12]+0.0100121*above[13]+0.0462248*left[13]+0.00542596*above[15]+0.024556*left[15]+0.00870648*above[14]+0.039777*left[14]);
		temp[6] = (int)(0.00301209*above[-1]+0.0085828*above[0]+0.00877928*left[0]+0.0161608*above[1]+0.0171015*left[1]+0.0221609*above[2]+0.0248433*left[2]+0.0264736*above[3]+0.0321804*left[3]+0.0291237*above[4]+0.039216*left[4]+0.030221*above[5]+0.0458829*left[5]+0.0299519*above[6]+0.0518918*left[6]+0.0285697*above[7]+0.0567315*left[7]+0.0263738*above[8]+0.0597655*left[8]+0.0236789*above[9]+0.06044*left[9]+0.0207813*above[10]+0.0585389*left[10]+0.0179318*above[11]+0.0543589*left[11]+0.0153256*above[12]+0.0486988*left[12]+0.0131353*above[13]+0.0427526*left[13]+0.00719857*above[15]+0.0236416*left[15]+0.0115117*above[14]+0.0377692*left[14]);
		temp[7] = (int)(0.00308357*above[-1]+0.00881805*above[0]+0.00895482*left[0]+0.0166921*above[1]+0.0173466*left[1]+0.0230956*above[2]+0.0249609*left[2]+0.0279377*above[3]+0.0319009*left[3]+0.0312229*above[4]+0.0382185*left[4]+0.0330059*above[5]+0.0438379*left[5]+0.0333966*above[6]+0.0485433*left[6]+0.032568*above[7]+0.0520082*left[7]+0.0307549*above[8]+0.0538783*left[8]+0.0282379*above[9]+0.0538984*left[9]+0.025315*above[10]+0.0520413*left[10]+0.0222709*above[11]+0.0485821*left[11]+0.0193555*above[12]+0.0440824*left[12]+0.0168145*above[13]+0.0393737*left[13]+0.00933997*above[15]+0.0223219*left[15]+0.0148748*above[14]+0.035359*left[14]);
		temp[8] = (int)(0.00308295*above[-1]+0.00884351*above[0]+0.00892536*left[0]+0.0168178*above[1]+0.0172095*left[1]+0.023454*above[2]+0.0245705*left[2]+0.0286907*above[3]+0.0310637*left[3]+0.0325289*above[4]+0.0367214*left[4]+0.0349871*above[5]+0.0414904*left[5]+0.0361098*above[6]+0.0452342*left[6]+0.0359874*above[7]+0.0477649*left[7]+0.0347716*above[8]+0.0489039*left[8]+0.0326769*above[9]+0.0485576*left[9]+0.0299669*above[10]+0.0467842*left[10]+0.0269285*above[11]+0.0438261*left[11]+0.0238473*above[12]+0.0400968*left[12]+0.0210388*above[13]+0.0362099*left[13]+0.0118713*above[15]+0.0208494*left[15]+0.0188153*above[14]+0.0328533*left[14]);
		temp[9] = (int)(0.00303164*above[-1]+0.00871865*above[0]+0.00875456*left[0]+0.0166449*above[1]+0.0168168*left[1]+0.0233692*above[2]+0.0238595*left[2]+0.0288644*above[3]+0.029908*left[3]+0.0331427*above[4]+0.0349908*left[4]+0.0362069*above[5]+0.0390853*left[5]+0.038058*above[6]+0.0421216*left[6]+0.0387175*above[7]+0.0440092*left[7]+0.0382525*above[8]+0.0446779*left[8]+0.0367934*above[9]+0.044123*left[9]+0.0345371*above[10]+0.0424405*left[10]+0.0317359*above[11]+0.039841*left[11]+0.0286792*above[12]+0.0366428*left[12]+0.0257338*above[13]+0.033322*left[13]+0.0147764*above[15]+0.0193737*left[15]+0.0232939*above[14]+0.0304283*left[14]);
		temp[10] = (int)(0.0029494*above[-1]+0.00849975*above[0]+0.00849975*left[0]+0.0162786*above[1]+0.0162786*left[1]+0.0229819*above[2]+0.0229819*left[2]+0.0286163*above[3]+0.0286163*left[3]+0.0332137*above[4]+0.0332137*left[4]+0.0367791*above[5]+0.0367791*left[5]+0.0392934*above[6]+0.0392934*left[6]+0.040732*above[7]+0.040732*left[7]+0.0410915*above[8]+0.0410915*left[8]+0.0404166*above[9]+0.0404166*left[9]+0.0388187*above[10]+0.0388187*left[10]+0.0364821*above[11]+0.0364821*left[11]+0.0336617*above[12]+0.0336617*left[12]+0.0307442*above[13]+0.0307442*left[13]+0.0179831*above[15]+0.0179831*left[15]+0.0281869*above[14]+0.0281869*left[14]);
		temp[11] = (int)(0.0028533*above[-1]+0.00823634*above[0]+0.00820967*left[0]+0.0158142*above[1]+0.0156864*left[1]+0.0224259*above[2]+0.0220607*left[2]+0.028108*above[3]+0.0273282*left[3]+0.0329149*above[4]+0.0315252*left[4]+0.036865*above[5]+0.0346791*left[5]+0.0399384*above[6]+0.0368072*left[6]+0.0420883*above[7]+0.0379285*left[7]+0.0432632*above[8]+0.0380807*left[8]+0.0434354*above[9]+0.037336*left[9]+0.0426276*above[10]+0.0358109*left[10]+0.0409356*above[11]+0.0336692*left[11]+0.0385455*above[12]+0.0311235*left[12]+0.0358255*above[13]+0.0285001*left[13]+0.0213452*above[15]+0.0167317*left[15]+0.0332659*above[14]+0.0261931*left[14]);
		temp[12] = (int)(0.00275771*above[-1]+0.00797055*above[0]+0.00792487*left[0]+0.0153341*above[1]+0.0151151*left[1]+0.0218209*above[2]+0.0211949*left[2]+0.0274916*above[3]+0.0261526*left[3]+0.0324215*above[4]+0.0300292*left[4]+0.0366474*above[5]+0.0328683*left[5]+0.04016*above[6]+0.0347119*left[6]+0.042909*above[7]+0.0356076*left[7]+0.0448181*above[8]+0.0356184*left[8]+0.0458077*above[9]+0.0348317*left[9]+0.0458233*above[10]+0.0333648*left[10]+0.0448678*above[11]+0.0313666*left[11]+0.0430371*above[12]+0.029021*left[12]+0.0406482*above[13]+0.0266132*left[13]+0.0246346*above[15]+0.015657*left[15]+0.0381896*above[14]+0.0244932*left[14]);
		temp[13] = (int)(0.0026747*above[-1]+0.00773795*above[0]+0.00767929*left[0]+0.0149085*above[1]+0.0146273*left[1]+0.0212705*above[2]+0.0204663*left[2]+0.0269025*above[3]+0.0251802*left[3]+0.0318969*above[4]+0.0288134*left[4]+0.0363078*above[5]+0.0314204*left[5]+0.0401406*above[6]+0.0330596*left[6]+0.0433525*above[7]+0.0337967*left[7]+0.0458591*above[8]+0.03371*left[8]+0.0475496*above[9]+0.0328957*left[9]+0.0483119*above[10]+0.0314705*left[10]+0.0480668*above[11]+0.029573*left[11]+0.0468183*above[12]+0.0273682*left[12]+0.044815*above[13]+0.0251139*left[13]+0.0275481*above[15]+0.0147908*left[15]+0.0425182*above[14]+0.0231295*left[14]);
		temp[14] = (int)(0.00261459*above[-1]+0.00756893*above[0]+0.00750207*left[0]+0.0145974*above[1]+0.0142769*left[1]+0.0208636*above[2]+0.0199463*left[2]+0.0264576*above[3]+0.0244919*left[3]+0.0314836*above[4]+0.0279597*left[4]+0.0360087*above[5]+0.030411*left[5]+0.0400512*above[6]+0.0319142*left[6]+0.0435775*above[7]+0.0325457*left[7]+0.0465034*above[8]+0.0323929*left[8]+0.0487028*above[9]+0.0315571*left[9]+0.0500266*above[10]+0.0301548*left[10]+0.0503363*above[11]+0.0283182*left[11]+0.0495624*above[12]+0.0262016*left[12]+0.0478918*above[13]+0.0240456*left[13]+0.0297352*above[15]+0.0141663*left[15]+0.0457515*above[14]+0.0221502*left[14]);
		temp[15] = (int)(0.00258641*above[-1]+0.00748981*above[0]+0.00741884*left[0]+0.0144522*above[1]+0.0141119*left[1]+0.0206745*above[2]+0.0197007*left[2]+0.0262527*above[3]+0.024165*left[3]+0.0312964*above[4]+0.0275517*left[4]+0.035879*above[5]+0.029925*left[5]+0.0400251*above[6]+0.0313578*left[6]+0.0437057*above[7]+0.0319317*left[7]+0.0468377*above[8]+0.0317386*left[8]+0.0492882*above[9]+0.0308829*left[9]+0.0508912*above[10]+0.0294816*left[10]+0.0514795*above[11]+0.0276654*left[11]+0.0509465*above[12]+0.0255842*left[12]+0.0494471*above[13]+0.0234713*left[13]+0.0308438*above[15]+0.0138246*left[15]+0.047389*above[14]+0.0216173*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.000452457*above[-1]+0.00126719*above[0]+0.00134331*left[0]+0.00232802*above[1]+0.00269471*left[1]+0.00306574*above[2]+0.00412693*left[2]+0.00346903*above[3]+0.00579754*left[3]+0.00357676*above[4]+0.00794119*left[4]+0.00345497*above[5]+0.0109471*left[5]+0.00317852*above[6]+0.0155527*left[6]+0.00281699*above[7]+0.023353*left[7]+0.00242635*above[8]+0.0385584*left[8]+0.0020461*above[9]+0.0736495*left[9]+0.00170057*above[10]+0.228019*left[10]+0.0014022*above[11]+0.241905*left[11]+0.00115568*above[12]+0.213871*left[12]+0.000963848*above[13]+0.0557626*left[13]+0.00051481*above[15]+0.0131451*left[15]+0.000829621*above[14]+0.0273687*left[14]);
		temp[1] = (int)(0.000935606*above[-1]+0.00262458*above[0]+0.00277285*left[0]+0.00483276*above[1]+0.00554657*left[1]+0.00638802*above[2]+0.00845041*left[2]+0.00726445*above[3]+0.0117742*left[3]+0.00753445*above[4]+0.01593*left[4]+0.00732512*above[5]+0.0215582*left[5]+0.00678398*above[6]+0.0297692*left[6]+0.00605157*above[7]+0.0427078*left[7]+0.00524415*above[8]+0.0647398*left[8]+0.0044466*above[9]+0.10787*left[9]+0.00371342*above[10]+0.147896*left[10]+0.00307434*above[11]+0.167241*left[11]+0.00254228*above[12]+0.134575*left[12]+0.00212577*above[13]+0.0888548*left[13]+0.00113829*above[15]+0.0250382*left[15]+0.00183296*above[14]+0.0488747*left[14]);
		temp[2] = (int)(0.00141108*above[-1]+0.00396786*above[0]+0.00417123*left[0]+0.00733089*above[1]+0.00830907*left[1]+0.00974388*above[2]+0.0125637*left[2]+0.011163*above[3]+0.0172989*left[3]+0.0116802*above[4]+0.0229958*left[4]+0.0114662*above[5]+0.0303215*left[5]+0.010726*above[6]+0.0402657*left[6]+0.00966259*above[7]+0.0542932*left[7]+0.00845134*above[8]+0.0743949*left[8]+0.00722656*above[9]+0.0977924*left[9]+0.00607972*above[10]+0.117216*left[10]+0.00506513*above[11]+0.120653*left[11]+0.0042102*above[12]+0.10642*left[12]+0.00353455*above[13]+0.0811335*left[13]+0.00190008*above[15]+0.0326043*left[15]+0.00305601*above[14]+0.0592086*left[14]);
		temp[3] = (int)(0.0018443*above[-1]+0.00520132*above[0]+0.00543462*left[0]+0.00965012*above[1]+0.010771*left[1]+0.0129152*above[2]+0.016138*left[2]+0.0149336*above[3]+0.0219084*left[3]+0.0157997*above[4]+0.0285284*left[4]+0.0157019*above[5]+0.0365189*left[5]+0.014878*above[6]+0.046464*left[6]+0.0135754*above[7]+0.0588527*left[7]+0.0120193*above[8]+0.0731643*left[8]+0.0103933*above[9]+0.0870586*left[9]+0.00883149*above[10]+0.096104*left[10]+0.00742094*above[11]+0.0964778*left[11]+0.00621221*above[12]+0.0874621*left[12]+0.00524416*above[13]+0.0731873*left[13]+0.0028344*above[15]+0.0354996*left[15]+0.00455126*above[14]+0.0596841*left[14]);
		temp[4] = (int)(0.00220983*above[-1]+0.00625264*above[0]+0.00648914*left[0]+0.0116552*above[1]+0.0127902*left[1]+0.0157199*above[2]+0.0189753*left[2]+0.0183679*above[3]+0.0253766*left[3]+0.0196805*above[4]+0.0323463*left[4]+0.0198383*above[5]+0.0402031*left[5]+0.0190822*above[6]+0.0491459*left[6]+0.0176783*above[7]+0.0590166*left[7]+0.0158845*above[8]+0.0689732*left[8]+0.0139264*above[9]+0.0772825*left[9]+0.011982*above[10]+0.0817724*left[10]+0.0101782*above[11]+0.0808597*left[11]+0.00859816*above[12]+0.0747224*left[12]+0.00731034*above[13]+0.0656059*left[13]+0.00397903*above[15]+0.0351062*left[15]+0.00637557*above[14]+0.0569015*left[14]);
		temp[5] = (int)(0.00249397*above[-1]+0.00708045*above[0]+0.0072976*left[0]+0.0132627*above[1]+0.014304*left[1]+0.0180336*above[2]+0.0210139*left[2]+0.0213052*above[3]+0.0276941*left[3]+0.023138*above[4]+0.0345933*left[4]+0.0236845*above[5]+0.0418643*left[5]+0.0231619*above[6]+0.0494777*left[6]+0.0218262*above[7]+0.0571005*left[7]+0.0199446*above[8]+0.0639943*left[8]+0.0177692*above[9]+0.0690752*left[9]+0.0155159*above[10]+0.0712332*left[10]+0.0133542*above[11]+0.0698577*left[11]+0.0114072*above[12]+0.0652709*left[12]+0.00978479*above[13]+0.0589057*left[13]+0.00537277*above[15]+0.0331989*left[15]+0.00858579*above[14]+0.0527904*left[14]);
		temp[6] = (int)(0.00269428*above[-1]+0.00767424*above[0]+0.00785692*left[0]+0.0144439*above[1]+0.0153193*left[1]+0.0197977*above[2]+0.0222997*left[2]+0.0236488*above[3]+0.0289964*left[3]+0.0260368*above[4]+0.0355747*left[4]+0.0270763*above[5]+0.0420882*left[5]+0.0269416*above[6]+0.0484254*left[6]+0.0258539*above[7]+0.0542636*left[7]+0.0240632*above[8]+0.0590693*left[8]+0.0218255*above[9]+0.0621932*left[9]+0.0193811*above[10]+0.0630781*left[10]+0.0169363*above[11]+0.0615149*left[11]+0.014658*above[12]+0.0578294*left[12]+0.0127065*above[13]+0.0529953*left[13]+0.00705109*above[15]+0.0307098*left[15]+0.0112315*above[14]+0.0483448*left[14]);
		temp[7] = (int)(0.00281716*above[-1]+0.00804871*above[0]+0.00818971*left[0]+0.0152168*above[1]+0.0158923*left[1]+0.0210157*above[2]+0.0229445*left[2]+0.0253698*above[3]+0.0294859*left[3]+0.0283044*above[4]+0.0356262*left[4]+0.0298968*above[5]+0.0413756*left[5]+0.030271*above[6]+0.0466174*left[6]+0.0295964*above[7]+0.0511018*left[7]+0.0280823*above[8]+0.0544763*left[8]+0.0259642*above[9]+0.0563679*left[9]+0.0234847*above[10]+0.0565085*left[10]+0.0208737*above[11]+0.0548653*left[11]+0.0183369*above[12]+0.0517248*left[12]+0.0160899*above[13]+0.0477932*left[13]+0.00903906*above[15]+0.0281102*left[15]+0.0143439*above[14]+0.0440171*left[14]);
		temp[8] = (int)(0.00287471*above[-1]+0.00823547*above[0]+0.00833406*left[0]+0.015633*above[1]+0.0161053*left[1]+0.0217397*above[2]+0.023088*left[2]+0.0264997*above[3]+0.0293759*left[3]+0.0299336*above[4]+0.0350467*left[4]+0.03209*above[5]+0.0401007*left[5]+0.0330454*above[6]+0.0444484*left[6]+0.032912*above[7]+0.0479201*left[7]+0.0318431*above[8]+0.0503002*left[8]+0.0300312*above[9]+0.0513873*left[9]+0.0276968*above[10]+0.0510689*left[10]+0.025072*above[11]+0.0493903*left[11]+0.0223876*above[12]+0.0465969*left[12]+0.0199113*above[13]+0.0432306*left[13]+0.0113414*above[15]+0.0256243*left[15]+0.0179212*above[14]+0.0400129*left[14]);
		temp[9] = (int)(0.00288167*above[-1]+0.00827486*above[0]+0.00833462*left[0]+0.0157635*above[1]+0.0160497*left[1]+0.0220544*above[2]+0.0228718*left[2]+0.0271167*above[3]+0.0288611*left[3]+0.0309756*above[4]+0.0340792*left[4]+0.033663*above[5]+0.0385326*left[5]+0.035219*above[6]+0.0421674*left[6]+0.035703*above[7]+0.044883*left[7]+0.035208*above[8]+0.046561*left[8]+0.0338686*above[9]+0.0471064*left[9]+0.0318609*above[10]+0.0464933*left[10]+0.0293943*above[11]+0.0448015*left[11]+0.0267031*above[12]+0.0422429*left[12]+0.0240946*above[13]+0.0392545*left[13]+0.0139301*above[15]+0.0233545*left[15]+0.021911*above[14]+0.0364184*left[14]);
		temp[10] = (int)(0.0028533*above[-1]+0.00820967*above[0]+0.00823634*left[0]+0.0156864*above[1]+0.0158142*left[1]+0.0220607*above[2]+0.0224259*left[2]+0.0273282*above[3]+0.028108*left[3]+0.0315252*above[4]+0.0329149*left[4]+0.0346791*above[5]+0.036865*left[5]+0.0368072*above[6]+0.0399384*left[6]+0.0379285*above[7]+0.0420883*left[7]+0.0380807*above[8]+0.0432632*left[8]+0.037336*above[9]+0.0434354*left[9]+0.0358109*above[10]+0.0426276*left[10]+0.0336692*above[11]+0.0409356*left[11]+0.0311235*above[12]+0.0385455*left[12]+0.0285001*above[13]+0.0358255*left[13]+0.0167317*above[15]+0.0213452*left[15]+0.0261931*above[14]+0.0332659*left[14]);
		temp[11] = (int)(0.00280402*above[-1]+0.00808107*above[0]+0.00808107*left[0]+0.0154792*above[1]+0.0154792*left[1]+0.0218638*above[2]+0.0218638*left[2]+0.027255*above[3]+0.027255*left[3]+0.0317042*above[4]+0.0317042*left[4]+0.035243*above[5]+0.035243*left[5]+0.0378789*above[6]+0.0378789*left[6]+0.0396051*above[7]+0.0396051*left[7]+0.0404162*above[8]+0.0404162*left[8]+0.0403267*above[9]+0.0403267*left[9]+0.0393882*above[10]+0.0393882*left[10]+0.0377033*above[11]+0.0377033*left[11]+0.0354399*above[12]+0.0354399*left[12]+0.0329188*above[13]+0.0329188*left[13]+0.0196163*above[15]+0.0196163*left[15]+0.0305674*above[14]+0.0305674*left[14]);
		temp[12] = (int)(0.00274682*above[-1]+0.00792658*above[0]+0.00790607*left[0]+0.0152137*above[1]+0.0151154*left[1]+0.021564*above[2]+0.0212829*left[2]+0.0270187*above[3]+0.0264172*left[3]+0.0316447*above[4]+0.0305693*left[4]+0.0354825*above[5]+0.0337813*left[5]+0.0385397*above[6]+0.0360823*left[6]+0.0407966*above[7]+0.0374939*left[7]+0.0422188*above[8]+0.0380417*left[8]+0.0427743*above[9]+0.0377674*left[9]+0.0424544*above[10]+0.03674*left[10]+0.0412948*above[11]+0.0350653*left[11]+0.0394038*above[12]+0.032898*left[12]+0.0370761*above[13]+0.0305264*left[13]+0.0223946*above[15]+0.0181806*left[15]+0.0347509*above[14]+0.0283336*left[14]);
		temp[13] = (int)(0.00269309*above[-1]+0.0077793*above[0]+0.00774389*left[0]+0.0149539*above[1]+0.0147841*left[1]+0.0212525*above[2]+0.0207669*left[2]+0.0267329*above[3]+0.025693*left[3]+0.0314754*above[4]+0.0296135*left[4]+0.0355306*above[5]+0.0325789*left[5]+0.0389111*above[6]+0.0346333*left[6]+0.0415935*above[7]+0.0358179*left[7]+0.0435265*above[8]+0.0361787*left[8]+0.0446455*above[9]+0.035775*left[9]+0.0448923*above[10]+0.0346866*left[10]+0.0442417*above[11]+0.0330207*left[11]+0.0427401*above[12]+0.0309236*left[12]+0.0406445*above[13]+0.028661*left[13]+0.0248248*above[15]+0.0170546*left[15]+0.0383894*above[14]+0.0265851*left[14]);
		temp[14] = (int)(0.00265267*above[-1]+0.00766776*above[0]+0.00762255*left[0]+0.014755*above[1]+0.0145382*left[1]+0.0210082*above[2]+0.0203881*left[2]+0.0264968*above[3]+0.0251682*left[3]+0.0313105*above[4]+0.0289294*left[4]+0.0355082*above[5]+0.0317278*left[5]+0.039108*above[6]+0.0336171*left[6]+0.042087*above[7]+0.0346508*left[7]+0.0443857*above[8]+0.0348877*left[8]+0.0459192*above[9]+0.0343977*left[9]+0.0465954*above[10]+0.0332673*left[10]+0.0463435*above[11]+0.0316049*left[11]+0.0451597*above[12]+0.0295518*left[12]+0.0432661*above[13]+0.0273598*left[13]+0.0266326*above[15]+0.0162653*left[15]+0.0410858*above[14]+0.0253614*left[14]);
		temp[15] = (int)(0.00263373*above[-1]+0.00761563*above[0]+0.00756562*left[0]+0.0146623*above[1]+0.0144225*left[1]+0.0208953*above[2]+0.0202092*left[2]+0.0263895*above[3]+0.024919*left[3]+0.0312392*above[4]+0.0286029*left[4]+0.0355074*above[5]+0.0313193*left[5]+0.039215*above[6]+0.0331263*left[6]+0.0423397*above[7]+0.0340833*left[7]+0.0448183*above[8]+0.034255*left[8]+0.0465569*above[9]+0.0337168*left[9]+0.047447*above[10]+0.0325589*left[10]+0.0473951*above[11]+0.0308912*left[11]+0.0463724*above[12]+0.0288534*left[12]+0.0445826*above[13]+0.0266916*left[13]+0.0275425*above[15]+0.015856*left[15]+0.042442*above[14]+0.0247288*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.000376794*above[-1]+0.00105889*above[0]+0.0011146*left[0]+0.00195486*above[1]+0.002223*left[1]+0.00259535*above[2]+0.00336947*left[2]+0.00296962*above[3]+0.00465964*left[3]+0.00310398*above[4]+0.00624092*left[4]+0.00304554*above[5]+0.00833729*left[5]+0.00284977*above[6]+0.0113319*left[6]+0.00257054*above[7]+0.0159608*left[7]+0.00225352*above[8]+0.0238236*left[8]+0.00193321*above[9]+0.0391429*left[9]+0.00163291*above[10]+0.074425*left[10]+0.00136642*above[11]+0.229116*left[11]+0.00114085*above[12]+0.243572*left[12]+0.000961647*above[13]+0.216717*left[13]+0.000519273*above[15]+0.0263407*left[15]+0.000833999*above[14]+0.0611383*left[14]);
		temp[1] = (int)(0.000784865*above[-1]+0.0022084*above[0]+0.00231863*left[0]+0.0040842*above[1]+0.00461449*left[1]+0.00543809*above[2]+0.00696743*left[2]+0.00624657*above[3]+0.0095778*left[3]+0.0065598*above[4]+0.0127161*left[4]+0.00646985*above[5]+0.0167733*left[5]+0.0060871*above[6]+0.022377*left[6]+0.00552072*above[7]+0.030635*left[7]+0.00486532*above[8]+0.0437*left[8]+0.00419418*above[9]+0.0659608*left[9]+0.00355825*above[10]+0.10947*left[10]+0.00298898*above[11]+0.150123*left[11]+0.00250362*above[12]+0.170551*left[12]+0.00211576*above[13]+0.140042*left[13]+0.00114538*above[15]+0.0462868*left[15]+0.00183816*above[14]+0.0987234*left[14]);
		temp[2] = (int)(0.00119654*above[-1]+0.00337294*above[0]+0.00352787*left[0]+0.00625422*above[1]+0.00699903*left[1]+0.00836346*above[2]+0.0105082*left[2]+0.00966285*above[3]+0.0143196*left[3]+0.0102186*above[4]+0.0187719*left[4]+0.0101575*above[5]+0.0243154*left[5]+0.00963563*above[6]+0.0315986*left[6]+0.00881174*above[7]+0.0416055*left[7]+0.00782804*above[8]+0.0558101*left[8]+0.00679884*above[9]+0.0762313*left[9]+0.00580713*above[10]+0.100149*left[10]+0.00490708*above[11]+0.120407*left[11]+0.00413086*above[12]+0.125226*left[12]+0.00350479*above[13]+0.113588*left[13]+0.00190485*above[15]+0.0544807*left[15]+0.0030533*above[14]+0.0931468*left[14]);
		temp[3] = (int)(0.0015853*above[-1]+0.00447905*above[0]+0.00466268*left[0]+0.00833245*above[1]+0.00921464*left[1]+0.0112031*above[2]+0.0137391*left[2]+0.0130389*above[3]+0.0185249*left[3]+0.0139117*above[4]+0.0239181*left[4]+0.0139671*above[5]+0.0303272*left[5]+0.0133905*above[6]+0.0382489*left[6]+0.012378*above[7]+0.0482572*left[7]+0.0111123*above[8]+0.0608476*left[8]+0.00974744*above[9]+0.0755248*left[9]+0.00840135*above[10]+0.0900019*left[10]+0.00715625*above[11]+0.099947*left[11]+0.00606527*above[12]+0.101729*left[12]+0.00517392*above[13]+0.0951685*left[13]+0.00282722*above[15]+0.0527833*left[15]+0.00452433*above[14]+0.0850275*left[14]);
		temp[4] = (int)(0.0019296*above[-1]+0.00546591*above[0]+0.00565986*left[0]+0.0102062*above[1]+0.0111373*left[1]+0.0138073*above[2]+0.0164792*left[2]+0.0162051*above[3]+0.0219644*left[3]+0.0174676*above[4]+0.0279013*left[4]+0.0177418*above[5]+0.0346008*left[5]+0.0172228*above[6]+0.0423519*left[6]+0.0161258*above[7]+0.0513386*left[7]+0.0146615*above[8]+0.0614026*left[8]+0.0130174*above[9]+0.0717181*left[9]+0.0113462*above[10]+0.0805881*left[10]+0.0097619*above[11]+0.0859067*left[11]+0.00834488*above[12]+0.0862092*left[12]+0.00716758*above[13]+0.0820369*left[13]+0.00394365*above[15]+0.0483125*left[15]+0.00629768*above[14]+0.0758971*left[14]);
		temp[5] = (int)(0.00221527*above[-1]+0.0062922*above[0]+0.00647935*left[0]+0.0117954*above[1]+0.0126933*left[1]+0.0160623*above[2]+0.018635*left[2]+0.0190219*above[3]+0.02455*left[3]+0.0207317*above[4]+0.0306886*left[4]+0.0213258*above[5]+0.0372573*left[5]+0.0209895*above[6]+0.0443756*left[6]+0.0199371*above[7]+0.0519931*left[7]+0.0183904*above[8]+0.0597684*left[8]+0.0165586*above[9]+0.0669661*left[9]+0.0146237*above[10]+0.0725172*left[10]+0.0127323*above[11]+0.0753398*left[11]+0.0109967*above[12]+0.0748754*left[12]+0.00952416*above[13]+0.0716798*left[13]+0.00528426*above[15]+0.0432742*left[15]+0.00841695*above[14]+0.0673047*left[14]);
		temp[6] = (int)(0.00243614*above[-1]+0.00693829*above[0]+0.00710539*left[0]+0.0130581*above[1]+0.0138594*left[1]+0.0179001*above[2]+0.0201932*left[2]+0.0213929*above[3]+0.0263086*left[3]+0.0235824*above[4]+0.0323982*left[4]+0.0245806*above[5]+0.0385866*left[5]+0.024548*above[6]+0.044892*left[6]+0.0236796*above[7]+0.0511767*left[7]+0.0221889*above[8]+0.0571012*left[8]+0.0202912*above[9]+0.0621219*left[9]+0.0181872*above[10]+0.0655848*left[10]+0.0160514*above[11]+0.0669316*left[11]+0.0140291*above[12]+0.0659641*left[12]+0.0122687*above[13]+0.0631776*left[13]+0.00687407*above[15]+0.0384927*left[15]+0.0109164*above[14]+0.059636*left[14]);
		temp[7] = (int)(0.00259308*above[-1]+0.00740438*above[0]+0.0075431*left[0]+0.0139884*above[1]+0.0146534*left[1]+0.0192988*above[2]+0.0212004*left[2]+0.0232715*above[3]+0.0273417*left[3]+0.0259431*above[4]+0.0332226*left[4]+0.0274009*above[5]+0.0389166*left[5]+0.0277726*above[6]+0.044398*left[6]+0.02722*above[7]+0.04952*left[7]+0.0259313*above[8]+0.0540091*left[8]+0.0241097*above[9]+0.0574926*left[9]+0.0219593*above[10]+0.0595792*left[10]+0.0196723*above[11]+0.0599829*left[11]+0.0174229*above[12]+0.0586643*left[12]+0.0154035*above[13]+0.056055*left[13]+0.0087268*above[15]+0.0342216*left[15]+0.0138115*above[14]+0.0529592*left[14]);
		temp[8] = (int)(0.00269231*above[-1]+0.00770606*above[0]+0.00781281*left[0]+0.01461*above[1]+0.0151216*left[1]+0.0202776*above[2]+0.0217401*left[2]+0.0246587*above[3]+0.0277867*left[3]+0.0277858*above[4]+0.0333729*left[4]+0.0297248*above[5]+0.0385459*left[5]+0.0305692*above[6]+0.0432713*left[6]+0.0304407*above[7]+0.0474278*left[7]+0.0294897*above[8]+0.0508175*left[8]+0.0278906*above[9]+0.0531985*left[9]+0.0258337*above[10]+0.0543436*left[10]+0.0235137*above[11]+0.0541166*left[11]+0.0211237*above[12]+0.0525534*left[12]+0.0188976*above[13]+0.0500364*left[13]+0.010837*above[15]+0.0305077*left[15]+0.0170873*above[14]+0.047218*left[14]);
		temp[9] = (int)(0.00274343*above[-1]+0.00786916*above[0]+0.00794406*left[0]+0.0149668*above[1]+0.0153258*left[1]+0.0208863*above[2]+0.0219124*left[2]+0.0255957*above[3]+0.0277899*left[3]+0.029129*above[4]+0.0330476*left[4]+0.0315375*above[5]+0.0377227*left[5]+0.0328849*above[6]+0.04179*left[6]+0.0332532*above[7]+0.0451639*left[7]+0.0327486*above[8]+0.0477138*left[8]+0.031505*above[9]+0.0492926*left[9]+0.0296822*above[10]+0.0497777*left[10]+0.0274602*above[11]+0.0491193*left[11]+0.0250361*above[12]+0.0473922*left[12]+0.0226763*above[13]+0.0449416*left[13]+0.0131708*above[15]+0.0273262*left[15]+0.0206856*above[14]+0.0423232*left[14]);
		temp[10] = (int)(0.00275771*above[-1]+0.00792487*above[0]+0.00797055*left[0]+0.0151151*above[1]+0.0153341*left[1]+0.0211949*above[2]+0.0218209*left[2]+0.0261526*above[3]+0.0274916*left[3]+0.0300292*above[4]+0.0324215*left[4]+0.0328683*above[5]+0.0366474*left[5]+0.0347119*above[6]+0.04016*left[6]+0.0356076*above[7]+0.042909*left[7]+0.0356184*above[8]+0.0448181*left[8]+0.0348317*above[9]+0.0458077*left[9]+0.0333648*above[10]+0.0458233*left[10]+0.0313666*above[11]+0.0448678*left[11]+0.029021*above[12]+0.0430371*left[12]+0.0266132*above[13]+0.0406482*left[13]+0.015657*above[15]+0.0246346*left[15]+0.0244932*above[14]+0.0381896*left[14]);
		temp[11] = (int)(0.00274682*above[-1]+0.00790607*above[0]+0.00792658*left[0]+0.0151154*above[1]+0.0152137*left[1]+0.0212829*above[2]+0.021564*left[2]+0.0264172*above[3]+0.0270187*left[3]+0.0305693*above[4]+0.0316447*left[4]+0.0337813*above[5]+0.0354825*left[5]+0.0360823*above[6]+0.0385397*left[6]+0.0374939*above[7]+0.0407966*left[7]+0.0380417*above[8]+0.0422188*left[8]+0.0377674*above[9]+0.0427743*left[9]+0.03674*above[10]+0.0424544*left[10]+0.0350653*above[11]+0.0412948*left[11]+0.032898*above[12]+0.0394038*left[12]+0.0305264*above[13]+0.0370761*left[13]+0.0181806*above[15]+0.0223946*left[15]+0.0283336*above[14]+0.0347509*left[14]);
		temp[12] = (int)(0.00272196*above[-1]+0.00784474*above[0]+0.00784474*left[0]+0.0150278*above[1]+0.0150278*left[1]+0.0212319*above[2]+0.0212319*left[2]+0.0264839*above[3]+0.0264839*left[3]+0.0308457*above[4]+0.0308457*left[4]+0.0343627*above[5]+0.0343627*left[5]+0.0370577*above[6]+0.0370577*left[6]+0.0389354*above[7]+0.0389354*left[7]+0.0399925*above[8]+0.0399925*left[8]+0.0402308*above[9]+0.0402308*left[9]+0.0396721*above[10]+0.0396721*left[10]+0.0383738*above[11]+0.0383738*left[11]+0.0364511*above[12]+0.0364511*left[12]+0.0341815*above[13]+0.0341815*left[13]+0.0205815*above[15]+0.0205815*left[15]+0.031967*above[14]+0.031967*left[14]);
		temp[13] = (int)(0.00269339*above[-1]+0.00777025*above[0]+0.00775465*left[0]+0.014908*above[1]+0.0148332*left[1]+0.0211187*above[2]+0.0209049*left[2]+0.0264443*above[3]+0.0259862*left[3]+0.0309564*above[4]+0.0301359*left[4]+0.0347055*above[5]+0.0334042*left[5]+0.0377131*above[6]+0.0358259*left[6]+0.0399741*above[7]+0.0374238*left[7]+0.0414653*above[8]+0.0382167*left[8]+0.0421575*above[9]+0.0382297*left[9]+0.042032*above[10]+0.037505*left[10]+0.0411007*above[11]+0.036114*left[11]+0.0394375*above[12]+0.0341764*left[12]+0.0373005*above[13]+0.0319563*left[13]+0.0226606*above[15]+0.0191893*left[15]+0.0350995*above[14]+0.0298289*left[14]);
		temp[14] = (int)(0.00267001*above[-1]+0.00770816*above[0]+0.00768204*left[0]+0.0148046*above[1]+0.0146794*left[1]+0.0210108*above[2]+0.0206525*left[2]+0.026379*above[3]+0.0256115*left[3]+0.0309886*above[4]+0.0296131*left[4]+0.0348944*above[5]+0.0327111*left[5]+0.0381187*above[6]+0.0349484*left[6]+0.0406515*above[7]+0.0363596*left[7]+0.0424574*above[8]+0.0369781*left[8]+0.0434865*above[9]+0.0368436*left[9]+0.0436909*above[10]+0.0360112*left[10]+0.0430481*above[11]+0.0345612*left[11]+0.0415978*above[12]+0.0326159*left[12]+0.0395795*above[13]+0.0304307*left[13]+0.0241945*above[15]+0.0182349*left[15]+0.037404*above[14]+0.0283632*left[14]);
		temp[15] = (int)(0.00265896*above[-1]+0.00767888*above[0]+0.00764767*left[0]+0.014756*above[1]+0.0146064*left[1]+0.0209606*above[2]+0.0205325*left[2]+0.02635*above[3]+0.0254328*left[3]+0.0310072*above[4]+0.0293632*left[4]+0.0349895*above[5]+0.0323789*left[5]+0.0383194*above[6]+0.0345266*left[6]+0.0409851*above[7]+0.0358468*left[7]+0.0429454*above[8]+0.0363793*left[8]+0.0441408*above[9]+0.0361712*left[9]+0.0445093*above[10]+0.0352837*left[10]+0.0440111*above[11]+0.0338019*left[11]+0.0426686*above[12]+0.0318498*left[12]+0.0407116*above[13]+0.0296792*left[13]+0.0249582*above[15]+0.0177629*left[15]+0.0385505*above[14]+0.0276393*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.00032318*above[-1]+0.000910597*above[0]+0.000953354*left[0]+0.00168742*above[1]+0.00189306*left[1]+0.00225439*above[2]+0.00284708*left[2]+0.00260175*above[3]+0.00389112*left[3]+0.0027484*above[4]+0.00512565*left[4]+0.00272962*above[5]+0.00669243*left[5]+0.0025882*above[6]+0.00881433*left[6]+0.0023671*above[7]+0.0118793*left[7]+0.00210428*above[8]+0.0166369*left[8]+0.00182991*above[9]+0.0247143*left[9]+0.00156569*above[10]+0.0403893*left[10]+0.00132572*above[11]+0.0762843*left[11]+0.0011184*above[12]+0.232115*left[12]+0.000950785*above[13]+0.249103*left[13]+0.000517907*above[15]+0.0606003*left[15]+0.000829569*above[14]+0.228304*left[14]);
		temp[1] = (int)(0.000677653*above[-1]+0.0019112*above[0]+0.00199697*left[0]+0.00354651*above[1]+0.0039589*left[1]+0.00474894*above[2]+0.00593667*left[2]+0.00549761*above[3]+0.00807758*left[3]+0.00582929*above[4]+0.0105724*left[4]+0.00581399*above[5]+0.0136802*left[5]+0.00553772*above[6]+0.0177891*left[6]+0.00508807*above[7]+0.0235359*left[7]+0.0045437*above[8]+0.0320532*left[8]+0.00396836*above[9]+0.0455451*left[9]+0.00340894*above[10]+0.0684997*left[10]+0.00289677*above[11]+0.113172*left[11]+0.00245126*above[12]+0.155901*left[12]+0.00208903*above[13]+0.180666*left[13]+0.00114076*above[15]+0.0956629*left[15]+0.00182585*above[14]+0.159453*left[14]);
		temp[2] = (int)(0.00104284*above[-1]+0.00294537*above[0]+0.00306844*left[0]+0.00547684*above[1]+0.00606831*left[1]+0.00735886*above[2]+0.00906064*left[2]+0.00855861*above[3]+0.0122473*left[3]+0.00912639*above[4]+0.0158801*left[4]+0.00916083*above[5]+0.0202805*left[5]+0.00878554*above[6]+0.0258935*left[6]+0.00812913*above[7]+0.0333775*left[7]+0.00731005*above[8]+0.0437467*left[8]+0.00642703*above[9]+0.0585362*left[9]+0.00555516*above[10]+0.079879*left[10]+0.00474674*above[11]+0.105271*left[11]+0.00403591*above[12]+0.127993*left[12]+0.00345277*above[13]+0.137449*left[13]+0.00189278*above[15]+0.0861027*left[15]+0.00302591*above[14]+0.134439*left[14]);
		temp[3] = (int)(0.00139768*above[-1]+0.00395474*above[0]+0.00410464*left[0]+0.00737297*above[1]+0.00809305*left[1]+0.00994961*above[2]+0.012019*left[2]+0.0116403*above[3]+0.0161145*left[3]+0.0125025*above[4]+0.0206558*left[4]+0.0126532*above[5]+0.0259662*left[5]+0.0122429*above[6]+0.0324449*left[6]+0.0114325*above[7]+0.0405946*left[7]+0.0103748*above[8]+0.0510162*left[8]+0.00920223*above[9]+0.0642591*left[9]+0.00801958*above[10]+0.0799232*left[10]+0.00690363*above[11]+0.0958854*left[11]+0.00590767*above[12]+0.108099*left[12]+0.00508048*above[13]+0.113638*left[13]+0.00279977*above[15]+0.0738925*left[15]+0.00446865*above[14]+0.11263*left[14]);
		temp[4] = (int)(0.00172356*above[-1]+0.00488695*above[0]+0.00505064*left[0]+0.00913832*above[1]+0.00992422*left[1]+0.0123937*above[2]+0.0146495*left[2]+0.0145993*above[3]+0.0194642*left[3]+0.0158131*above[4]+0.0246361*left[4]+0.0161588*above[5]+0.0304469*left[5]+0.0157997*above[6]+0.0371919*left[6]+0.0149162*above[7]+0.0451589*left[7]+0.0136862*above[8]+0.0545472*left[8]+0.0122702*above[9]+0.0652322*left[9]+0.010802*above[10]+0.0764434*left[10]+0.00938493*above[11]+0.0865647*left[11]+0.00809571*above[12]+0.0936*left[12]+0.00700777*above[13]+0.0963879*left[13]+0.00388748*above[15]+0.0627814*left[15]+0.00619216*above[14]+0.0953884*left[14]);
		temp[5] = (int)(0.00200665*above[-1]+0.00570229*above[0]+0.00586656*left[0]+0.0106975*above[1]+0.0114858*left[1]+0.0145871*above[2]+0.0168472*left[2]+0.0173112*above[3]+0.0221743*left[3]+0.018924*above[4]+0.0277059*left[4]+0.0195446*above[5]+0.0336625*left[5]+0.0193351*above[6]+0.0402284*left[6]+0.0184804*above[7]+0.0475122*left[7]+0.0171705*above[8]+0.0554635*left[8]+0.0155855*above[9]+0.0637494*left[9]+0.0138832*above[10]+0.0716478*left[10]+0.0121933*above[11]+0.0780992*left[11]+0.0106188*above[12]+0.0820187*left[12]+0.00926353*above[13]+0.0830001*left[13]+0.00517961*above[15]+0.053574*left[15]+0.00823033*above[14]+0.0815323*left[14]);
		temp[6] = (int)(0.00223887*above[-1]+0.00637666*above[0]+0.00653014*left[0]+0.0120024*above[1]+0.0127386*left[1]+0.0164579*above[2]+0.0185667*left[2]+0.0196825*above[3]+0.0242115*left[3]+0.021724*above[4]+0.0298743*left[4]+0.0226896*above[5]+0.035716*left[5]+0.0227278*above[6]+0.0418336*left[6]+0.0220137*above[7]+0.0482245*left[7]+0.0207353*above[8]+0.0547363*left[8]+0.0190798*above[9]+0.061016*left[9]+0.0172209*above[10]+0.0665018*left[10]+0.0153107*above[11]+0.0705084*left[11]+0.0134786*above[12]+0.0724348*left[12]+0.0118635*above[13]+0.0722064*left[13]+0.00669384*above[15]+0.0460768*left[15]+0.0106068*above[14]+0.0703365*left[14]);
		temp[7] = (int)(0.0024178*above[-1]+0.0069016*above[0]+0.00703594*left[0]+0.013033*above[1]+0.0136773*left[1]+0.01797*above[2]+0.019814*left[2]+0.0216565*above[3]+0.0256113*left[3]+0.024135*above[4]+0.031234*left[4]+0.0254966*above[5]+0.0367959*left[5]+0.0258682*above[6]+0.0423387*left[6]+0.0254038*above[7]+0.0478083*left[7]+0.0242756*above[8]+0.053032*left[8]+0.0226644*above[9]+0.0577078*left[9]+0.0207482*above[10]+0.0614283*left[10]+0.0186937*above[11]+0.0637568*left[11]+0.0166533*above[12]+0.0643649*left[12]+0.0148025*above[13]+0.0633459*left[13]+0.00843655*above[15]+0.0399752*left[15]+0.0133269*above[14]+0.0612174*left[14]);
		temp[8] = (int)(0.00254578*above[-1]+0.00728224*above[0]+0.00739245*left[0]+0.0137949*above[1]+0.0143232*left[1]+0.0191212*above[2]+0.0206331*left[2]+0.0232153*above[3]+0.0264549*left[3]+0.0261166*above[4]+0.0319225*left[4]+0.0279005*above[5]+0.0371181*left[5]+0.0286692*above[6]+0.0420552*left[6]+0.0285481*above[7]+0.0466663*left[7]+0.0276831*above[8]+0.0507974*left[8]+0.0262354*above[9]+0.0542152*left[9]+0.0243745*above[10]+0.0566379*left[10]+0.0222703*above[11]+0.0577951*left[11]+0.0200907*above[12]+0.0575226*left[12]+0.0180458*above[13]+0.0560016*left[13]+0.0103963*above[15]+0.0349885*left[15]+0.0163684*above[14]+0.0537393*left[14]);
		temp[9] = (int)(0.00262869*above[-1]+0.00753404*above[0]+0.00761821*left[0]+0.0143132*above[1]+0.0147167*left[1]+0.0199377*above[2]+0.0210919*left[2]+0.0243753*above[3]+0.0268474*left[3]+0.0276661*above[4]+0.032093*left[4]+0.0298729*above[5]+0.0368923*left[5]+0.0310741*above[6]+0.0412492*left[6]+0.0313643*above[7]+0.0451037*left[7]+0.0308564*above[8]+0.0483354*left[8]+0.0296822*above[9]+0.050777*left[9]+0.0279897*above[10]+0.0522434*left[10]+0.0259393*above[11]+0.0525773*left[11]+0.0237037*above[12]+0.0517173*left[12]+0.0215215*above[13]+0.0498901*left[13]+0.0125369*above[15]+0.0308992*left[15]+0.0196713*above[14]+0.0475826*left[14]);
		temp[10] = (int)(0.0026747*above[-1]+0.00767929*above[0]+0.00773795*left[0]+0.0146273*above[1]+0.0149085*left[1]+0.0204663*above[2]+0.0212705*left[2]+0.0251802*above[3]+0.0269025*left[3]+0.0288134*above[4]+0.0318969*left[4]+0.0314204*above[5]+0.0363078*left[5]+0.0330596*above[6]+0.0401406*left[6]+0.0337967*above[7]+0.0433525*left[7]+0.03371*above[8]+0.0458591*left[8]+0.0328957*above[9]+0.0475496*left[9]+0.0314705*above[10]+0.0483119*left[10]+0.029573*above[11]+0.0480668*left[11]+0.0273682*above[12]+0.0468183*left[12]+0.0251139*above[13]+0.044815*left[13]+0.0147908*above[15]+0.0275481*left[15]+0.0231295*above[14]+0.0425182*left[14]);
		temp[11] = (int)(0.00269309*above[-1]+0.00774389*above[0]+0.0077793*left[0]+0.0147841*above[1]+0.0149539*left[1]+0.0207669*above[2]+0.0212525*left[2]+0.025693*above[3]+0.0267329*left[3]+0.0296135*above[4]+0.0314754*left[4]+0.0325789*above[5]+0.0355306*left[5]+0.0346333*above[6]+0.0389111*left[6]+0.0358179*above[7]+0.0415935*left[7]+0.0361787*above[8]+0.0435265*left[8]+0.035775*above[9]+0.0446455*left[9]+0.0346866*above[10]+0.0448923*left[10]+0.0330207*above[11]+0.0442417*left[11]+0.0309236*above[12]+0.0427401*left[12]+0.028661*above[13]+0.0406445*left[13]+0.0170546*above[15]+0.0248248*left[15]+0.0265851*above[14]+0.0383894*left[14]);
		temp[12] = (int)(0.00269339*above[-1]+0.00775465*above[0]+0.00777025*left[0]+0.0148332*above[1]+0.014908*left[1]+0.0209049*above[2]+0.0211187*left[2]+0.0259862*above[3]+0.0264443*left[3]+0.0301359*above[4]+0.0309564*left[4]+0.0334042*above[5]+0.0347055*left[5]+0.0358259*above[6]+0.0377131*left[6]+0.0374238*above[7]+0.0399741*left[7]+0.0382167*above[8]+0.0414653*left[8]+0.0382297*above[9]+0.0421575*left[9]+0.037505*above[10]+0.042032*left[10]+0.036114*above[11]+0.0411007*left[11]+0.0341764*above[12]+0.0394375*left[12]+0.0319563*above[13]+0.0373005*left[13]+0.0191893*above[15]+0.0226606*left[15]+0.0298289*above[14]+0.0350995*left[14]);
		temp[13] = (int)(0.00268462*above[-1]+0.00773717*above[0]+0.00773717*left[0]+0.0148223*above[1]+0.0148223*left[1]+0.0209441*above[2]+0.0209441*left[2]+0.0261325*above[3]+0.0261325*left[3]+0.0304536*above[4]+0.0304536*left[4]+0.0339592*above[5]+0.0339592*left[5]+0.0366789*above[6]+0.0366789*left[6]+0.0386229*above[7]+0.0386229*left[7]+0.0397895*above[8]+0.0397895*left[8]+0.0401755*above[9]+0.0401755*left[9]+0.0397898*above[10]+0.0397898*left[10]+0.0386697*above[11]+0.0386697*left[11]+0.0369062*above[12]+0.0369062*left[12]+0.0347556*above[13]+0.0347556*left[13]+0.021024*above[15]+0.021024*left[15]+0.0326072*above[14]+0.0326072*left[14]);
		temp[14] = (int)(0.00267468*above[-1]+0.00771392*above[0]+0.0077032*left[0]+0.0147934*above[1]+0.014742*left[1]+0.0209411*above[2]+0.0207941*left[2]+0.0261965*above[3]+0.0258815*left[3]+0.0306315*above[4]+0.030067*left[4]+0.0342999*above[5]+0.0334039*left[5]+0.0372288*above[6]+0.0359277*left[6]+0.0394208*above[7]+0.0376594*left[7]+0.0408605*above[8]+0.0386117*left[8]+0.041525*above[9]+0.0387986*left[9]+0.0413988*above[10]+0.0382468*left[10]+0.0404926*above[11]+0.0370092*left[11]+0.0388738*above[12]+0.0351872*left[12]+0.0367897*above[13]+0.0330346*left[13]+0.0223678*above[15]+0.0199215*left[15]+0.0346371*above[14]+0.0309262*left[14]);
		temp[15] = (int)(0.00266973*above[-1]+0.0077022*above[0]+0.00768639*left[0]+0.0147783*above[1]+0.0147025*left[1]+0.0209379*above[2]+0.020721*left[2]+0.0262249*above[3]+0.0257601*left[3]+0.0307142*above[4]+0.0298811*left[4]+0.0344605*above[5]+0.0331377*left[5]+0.03749*above[6]+0.0355689*left[6]+0.0398019*above[7]+0.0372003*left[7]+0.0413744*above[8]+0.0380516*left[8]+0.0421754*above[9]+0.0381449*left[9]+0.0421774*above[10]+0.037515*left[10]+0.0413778*above[11]+0.0362222*left[11]+0.0398323*above[12]+0.0343727*left[12]+0.0377832*above[13]+0.0322193*left[13]+0.0230257*above[15]+0.0193993*left[15]+0.0356303*above[14]+0.03013*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.0002867*above[-1]+0.000809266*above[0]+0.000844128*left[0]+0.00150355*above[1]+0.00167114*left[1]+0.00201747*above[2]+0.00249996*left[2]+0.00234219*above[3]+0.00338943*left[3]+0.00249239*above[4]+0.00441481*left[4]+0.00249628*above[5]+0.00567639*left[5]+0.00238877*above[6]+0.00732166*left[6]+0.00220579*above[7]+0.00958849*left[7]+0.00197999*above[8]+0.0128954*left[8]+0.00173823*above[9]+0.0180521*left[9]+0.00150068*above[10]+0.026808*left[10]+0.0012812*above[11]+0.0437131*left[11]+0.00108869*above[12]+0.0820786*left[12]+0.000931043*above[13]+0.244128*left[13]+0.000510242*above[15]+0.231335*left[15]+0.000815766*above[14]+0.27795*left[14]);
		temp[1] = (int)(0.000605415*above[-1]+0.00171022*above[0]+0.00178108*left[0]+0.00318095*above[1]+0.00352152*left[1]+0.00427605*above[2]+0.00525609*left[2]+0.00497678*above[3]+0.00710174*left[3]+0.00531224*above[4]+0.00920543*left[4]+0.00533923*above[5]+0.0117573*left[5]+0.00512877*above[6]+0.0150271*left[6]+0.00475467*above[7]+0.0194313*left[7]+0.00428489*above[8]+0.0256631*left[8]+0.00377619*above[9]+0.034963*left[9]+0.00327195*above[10]+0.049746*left[10]+0.00280262*above[11]+0.0749414*left[11]+0.00238837*above[12]+0.123826*left[12]+0.00204731*above[13]+0.175919*left[13]+0.00112466*above[15]+0.159091*left[15]+0.00179678*above[14]+0.223567*left[14]);
		temp[2] = (int)(0.00094038*above[-1]+0.00265954*above[0]+0.00276312*left[0]+0.00495497*above[1]+0.00545265*left[1]+0.00667947*above[2]+0.00811058*left[2]+0.00780382*above[3]+0.0109022*left[3]+0.00836906*above[4]+0.0140297*left[4]+0.00845697*above[5]+0.0177426*left[5]+0.0081713*above[6]+0.0223748*left[6]+0.00762164*above[7]+0.0284054*left[7]+0.00691091*above[8]+0.0365602*left[8]+0.00612703*above[9]+0.0479762*left[9]+0.005339*above[10]+0.0644111*left[10]+0.00459687*above[11]+0.0884147*left[11]+0.00393515*above[12]+0.118299*left[12]+0.00338568*above[13]+0.149288*left[13]+0.00186688*above[15]+0.120179*left[15]+0.00297912*above[14]+0.170726*left[14]);
		temp[3] = (int)(0.00127291*above[-1]+0.00360535*above[0]+0.00373434*left[0]+0.0067316*above[1]+0.00735115*left[1]+0.00910711*above[2]+0.0108872*left[2]+0.0106927*above[3]+0.0145397*left[3]+0.011537*above[4]+0.0185418*left[4]+0.01174*above[5]+0.0231638*left[5]+0.0114304*above[6]+0.0287377*left[6]+0.0107473*above[7]+0.0356919*left[7]+0.00982443*above[8]+0.0445871*left[8]+0.00877963*above[9]+0.0561229*left[9]+0.00770849*above[10]+0.071006*left[10]+0.00668322*above[11]+0.089075*left[11]+0.00575614*above[12]+0.108522*left[12]+0.00497718*above[13]+0.125632*left[13]+0.00275844*above[15]+0.0936672*left[15]+0.00439496*above[14]+0.136338*left[14]);
		temp[4] = (int)(0.00158612*above[-1]+0.00450028*above[0]+0.00464478*left[0]+0.00842372*above[1]+0.00911751*left[1]+0.0114444*above[2]+0.013436*left[2]+0.0135149*above[3]+0.017811*left[3]+0.014687*above[4]+0.0224822*left[4]+0.0150698*above[5]+0.027706*left[5]+0.0148067*above[6]+0.0337637*left[6]+0.0140564*above[7]+0.0409665*left[7]+0.0129762*above[8]+0.0496359*left[8]+0.0117092*above[9]+0.0600221*left[9]+0.0103765*above[10]+0.0720632*left[10]+0.00907386*above[11]+0.0850351*left[11]+0.00787449*above[12]+0.0972664*left[12]+0.00685133*above[13]+0.106694*left[13]+0.0038213*above[15]+0.075359*left[15]+0.00607652*above[14]+0.111468*left[14]);
		temp[5] = (int)(0.00186655*above[-1]+0.00530597*above[0]+0.00545523*left[0]+0.00995927*above[1]+0.0106756*left[1]+0.0135934*above[2]+0.015648*left[2]+0.0161555*above[3]+0.0205797*left[3]+0.0176967*above[4]+0.0256977*left[4]+0.0183266*above[5]+0.0312229*left[5]+0.0181921*above[6]+0.0373644*left[6]+0.0174596*above[7]+0.0443044*left[7]+0.0162996*above[8]+0.0521559*left[8]+0.0148734*above[9]+0.0608715*left[9]+0.0133234*above[10]+0.070103*left[10]+0.0117681*above[11]+0.0790607*left[11]+0.0103038*above[12]+0.0865271*left[12]+0.00903101*above[13]+0.0913767*left[13]+0.00507495*above[15]+0.062098*left[15]+0.00805142*above[14]+0.092942*left[14]);
		temp[6] = (int)(0.00210516*above[-1]+0.00599605*above[0]+0.00614006*left[0]+0.0112871*above[1]+0.011978*left[1]+0.0154807*above[2]+0.0174608*left[2]+0.0185227*above[3]+0.0227799*left[3]+0.0204613*above[4]+0.028138*left[4]+0.0213998*above[5]+0.0337129*left[5]+0.021478*above[6]+0.0396431*left[6]+0.0208584*above[7]+0.0460064*left[7]+0.0197128*above[8]+0.05278*left[8]+0.0182111*above[9]+0.0597798*left[9]+0.01651*above[10]+0.0665957*left[10]+0.0147474*above[11]+0.0725652*left[11]+0.0130424*above[12]+0.0768642*left[12]+0.0115265*above[13]+0.0789466*left[13]+0.00653262*above[15]+0.0521492*left[15]+0.0103369*above[14]+0.0787322*left[14]);
		temp[7] = (int)(0.00229761*above[-1]+0.00655715*above[0]+0.0066878*left[0]+0.0123793*above[1]+0.0130059*left[1]+0.0170621*above[2]+0.0188569*left[2]+0.0205547*above[3]+0.0244085*left[3]+0.0229021*above[4]+0.0298353*left[4]+0.0241968*above[5]+0.0352743*left[5]+0.0245641*above[6]+0.0408087*left[6]+0.0241522*above[7]+0.0464483*left[7]+0.0231232*above[8]+0.0521025*left[8]+0.0216438*above[9]+0.0575499*left[9]+0.0198761*above[10]+0.0624195*left[10]+0.017971*above[11]+0.0662126*left[11]+0.0160672*above[12]+0.0684033*left[12]+0.0143288*above[13]+0.0687914*left[13]+0.00819657*above[15]+0.0444754*left[15]+0.0129328*above[14]+0.0675897*left[14]);
		temp[8] = (int)(0.0024439*above[-1]+0.00698803*above[0]+0.00709969*left[0]+0.0132302*above[1]+0.0137657*left[1]+0.0183226*above[2]+0.0198556*left[2]+0.0222219*above[3]+0.0255105*left[3]+0.0249711*above[4]+0.0308774*left[4]+0.0266509*above[5]+0.0360609*left[5]+0.0273674*above[6]+0.0411075*left[6]+0.0272476*above[7]+0.0459918*left[7]+0.0264339*above[8]+0.0506026*left[8]+0.0250793*above[9]+0.0547331*left[9]+0.0233407*above[10]+0.0580864*left[10]+0.0213728*above[11]+0.0603126*left[11]+0.0193282*above[12]+0.0610975*left[12]+0.0174019*above[13]+0.0604447*left[13]+0.0100525*above[15]+0.0384334*left[15]+0.0158134*above[14]+0.0587061*left[14]);
		temp[9] = (int)(0.00254754*above[-1]+0.00729756*above[0]+0.00738719*left[0]+0.0138534*above[1]+0.0142831*left[1]+0.0192733*above[2]+0.0205033*left[2]+0.0235252*above[3]+0.0261622*left[3]+0.0266525*above[4]+0.0313831*left[4]+0.0287247*above[5]+0.0362474*left[5]+0.0298284*above[6]+0.0407817*left[6]+0.0300652*above[7]+0.044949*left[7]+0.0295516*above[8]+0.048645*left[8]+0.0284178*above[9]+0.0517001*left[9]+0.0268048*above[10]+0.0538957*left[10]+0.0248609*above[11]+0.0550016*left[11]+0.0227445*above[12]+0.0548497*left[12]+0.0206768*above[13]+0.0535612*left[13]+0.0120638*above[15]+0.0336101*left[15]+0.0189192*above[14]+0.0515462*left[14]);
		temp[10] = (int)(0.00261459*above[-1]+0.00750207*above[0]+0.00756893*left[0]+0.0142769*above[1]+0.0145974*left[1]+0.0199463*above[2]+0.0208636*left[2]+0.0244919*above[3]+0.0264576*left[3]+0.0279597*above[4]+0.0314836*left[4]+0.030411*above[5]+0.0360087*left[5]+0.0319142*above[6]+0.0400512*left[6]+0.0325457*above[7]+0.0435775*left[7]+0.0323929*above[8]+0.0465034*left[8]+0.0315571*above[9]+0.0487028*left[9]+0.0301548*above[10]+0.0500266*left[10]+0.0283182*above[11]+0.0503363*left[11]+0.0262016*above[12]+0.0495624*left[12]+0.0240456*above[13]+0.0478918*left[13]+0.0141663*above[15]+0.0297352*left[15]+0.0221502*above[14]+0.0457515*left[14]);
		temp[11] = (int)(0.00265267*above[-1]+0.00762255*above[0]+0.00766776*left[0]+0.0145382*above[1]+0.014755*left[1]+0.0203881*above[2]+0.0210082*left[2]+0.0251682*above[3]+0.0264968*left[3]+0.0289294*above[4]+0.0313105*left[4]+0.0317278*above[5]+0.0355082*left[5]+0.0336171*above[6]+0.039108*left[6]+0.0346508*above[7]+0.042087*left[7]+0.0348877*above[8]+0.0443857*left[8]+0.0343977*above[9]+0.0459192*left[9]+0.0332673*above[10]+0.0465954*left[10]+0.0316049*above[11]+0.0463435*left[11]+0.0295518*above[12]+0.0451597*left[12]+0.0273598*above[13]+0.0432661*left[13]+0.0162653*above[15]+0.0266326*left[15]+0.0253614*above[14]+0.0410858*left[14]);
		temp[12] = (int)(0.00267001*above[-1]+0.00768204*above[0]+0.00770816*left[0]+0.0146794*above[1]+0.0148046*left[1]+0.0206525*above[2]+0.0210108*left[2]+0.0256115*above[3]+0.026379*left[3]+0.0296131*above[4]+0.0309886*left[4]+0.0327111*above[5]+0.0348944*left[5]+0.0349484*above[6]+0.0381187*left[6]+0.0363596*above[7]+0.0406515*left[7]+0.0369781*above[8]+0.0424574*left[8]+0.0368436*above[9]+0.0434865*left[9]+0.0360112*above[10]+0.0436909*left[10]+0.0345612*above[11]+0.0430481*left[11]+0.0326159*above[12]+0.0415978*left[12]+0.0304307*above[13]+0.0395795*left[13]+0.0182349*above[15]+0.0241945*left[15]+0.0283632*above[14]+0.037404*left[14]);
		temp[13] = (int)(0.00267468*above[-1]+0.0077032*above[0]+0.00771392*left[0]+0.014742*above[1]+0.0147934*left[1]+0.0207941*above[2]+0.0209411*left[2]+0.0258815*above[3]+0.0261965*left[3]+0.030067*above[4]+0.0306315*left[4]+0.0334039*above[5]+0.0342999*left[5]+0.0359277*above[6]+0.0372288*left[6]+0.0376594*above[7]+0.0394208*left[7]+0.0386117*above[8]+0.0408605*left[8]+0.0387986*above[9]+0.041525*left[9]+0.0382468*above[10]+0.0413988*left[10]+0.0370092*above[11]+0.0404926*left[11]+0.0351872*above[12]+0.0388738*left[12]+0.0330346*above[13]+0.0367897*left[13]+0.0199215*above[15]+0.0223678*left[15]+0.0309262*above[14]+0.0346371*left[14]);
		temp[14] = (int)(0.0026738*above[-1]+0.00770603*above[0]+0.00770603*left[0]+0.0147628*above[1]+0.0147628*left[1]+0.0208608*above[2]+0.0208608*left[2]+0.0260307*above[3]+0.0260307*left[3]+0.0303399*above[4]+0.0303399*left[4]+0.033842*above[5]+0.033842*left[5]+0.0365686*above[6]+0.0365686*left[6]+0.0385314*above[7]+0.0385314*left[7]+0.0397292*above[8]+0.0397292*left[8]+0.0401574*above[9]+0.0401574*left[9]+0.0398215*above[10]+0.0398215*left[10]+0.0387528*above[11]+0.0387528*left[11]+0.0370354*above[12]+0.0370354*left[12]+0.0349194*above[13]+0.0349194*left[13]+0.0211507*above[15]+0.0211507*left[15]+0.0327903*above[14]+0.0327903*left[14]);
		temp[15] = (int)(0.00267282*above[-1]+0.0077057*above[0]+0.00770071*left[0]+0.0147694*above[1]+0.0147455*left[1]+0.0208875*above[2]+0.020819*left[2]+0.0260948*above[3]+0.0259482*left[3]+0.0304612*above[4]+0.0301984*left[4]+0.0340405*above[5]+0.0336233*left[5]+0.0368625*above[6]+0.0362566*left[6]+0.038935*above[7]+0.0381146*left[7]+0.0402501*above[8]+0.0392025*left[8]+0.0407949*above[9]+0.0395243*left[9]+0.0405642*above[10]+0.0390948*left[10]+0.0395789*above[11]+0.0379545*left[11]+0.0379144*above[12]+0.0361946*left[12]+0.0358185*above[13]+0.0340662*left[13]+0.0217388*above[15]+0.0205969*left[15]+0.0336812*above[14]+0.0319492*left[14]);
		memcpy(dst, temp, bs);
		dst += stride;
		temp[0] = (int)(0.000261947*above[-1]+0.000740021*above[0]+0.000770563*left[0]+0.00137657*above[1]+0.00152336*left[1]+0.00185085*above[2]+0.00227324*left[2]+0.00215474*above[3]+0.00307059*left[3]+0.0023008*above[4]+0.00397868*left[4]+0.00231349*above[5]+0.00507944*left[5]+0.0022234*above[6]+0.00648914*left[6]+0.00206236*above[7]+0.00838785*left[7]+0.00185968*above[8]+0.0110772*left[8]+0.0016399*above[9]+0.0151041*left[9]+0.00142179*above[10]+0.0215568*left[10]+0.00121855*above[11]+0.0328669*left[11]+0.00103899*above[12]+0.0555194*left[12]+0.000891019*above[13]+0.112236*left[13]+0.000489699*above[15]+0.332491*left[15]+0.000782237*above[14]+0.361886*left[14]);
		temp[1] = (int)(0.000564426*above[-1]+0.00159556*above[0]+0.00165924*left[0]+0.00297077*above[1]+0.00327678*left[1]+0.00400044*above[2]+0.00488069*left[2]+0.0046671*above[3]+0.0065742*left[3]+0.00499643*above[4]+0.00848531*left[4]+0.00503904*above[5]+0.0107758*left[5]+0.0048587*above[6]+0.0136692*left[6]+0.00452228*above[7]+0.0175004*left[7]+0.00409208*above[8]+0.0228085*left[8]+0.00362081*above[9]+0.0305232*left[9]+0.00314944*above[10]+0.0423661*left[10]+0.00270732*above[11]+0.0618234*left[11]+0.00231444*above[12]+0.0965391*left[12]+0.0019891*above[13]+0.16836*left[13]+0.00109561*above[15]+0.203248*left[15]+0.00174892*above[14]+0.25136*left[14]);
		temp[2] = (int)(0.000886281*above[-1]+0.00250796*above[0]+0.00260262*left[0]+0.00467642*above[1]+0.00513122*left[1]+0.00631274*above[2]+0.00762026*left[2]+0.00738963*above[3]+0.0102191*left[3]+0.0079441*above[4]+0.0131088*left[4]+0.00805046*above[5]+0.01651*left[5]+0.00780335*above[6]+0.0207128*left[6]+0.00730343*above[7]+0.0261282*left[7]+0.006646*above[8]+0.0333733*left[8]+0.00591332*above[9]+0.0434175*left[9]+0.00517079*above[10]+0.0578345*left[10]+0.00446665*above[11]+0.079166*left[11]+0.00383491*above[12]+0.111264*left[12]+0.0033075*above[13]+0.15223*left[13]+0.00182837*above[15]+0.139855*left[15]+0.0029154*above[14]+0.190516*left[14]);
		temp[3] = (int)(0.00120963*above[-1]+0.00342751*above[0]+0.00354722*left[0]+0.00640339*above[1]+0.00697833*left[1]+0.00867191*above[2]+0.0103237*left[2]+0.0101964*above[3]+0.0137655*left[3]+0.0110218*above[4]+0.017519*left[4]+0.0112406*above[5]+0.0218323*left[5]+0.0109721*above[6]+0.0270102*left[6]+0.0103453*above[7]+0.0334497*left[7]+0.00948538*above[8]+0.0416854*left[8]+0.0085029*above[9]+0.0524407*left[9]+0.00748859*above[10]+0.0666375*left[10]+0.00651179*above[11]+0.0851984*left[11]+0.0056236*above[12]+0.107702*left[12]+0.00487362*above[13]+0.131107*left[13]+0.00270745*above[15]+0.104576*left[15]+0.00431055*above[14]+0.148792*left[14]);
		temp[4] = (int)(0.00151815*above[-1]+0.00430847*above[0]+0.00444464*left[0]+0.00806765*above[1]+0.00872145*left[1]+0.0109677*above[2]+0.0128446*left[2]+0.012964*above[3]+0.0170133*left[3]+0.0141058*above[4]+0.0214552*left[4]+0.014496*above[5]+0.0264158*left[5]+0.0142694*above[6]+0.0321697*left[6]+0.0135754*above[7]+0.0390338*left[7]+0.012562*above[8]+0.0473691*left[8]+0.0113644*above[9]+0.0575489*left[9]+0.0100974*above[10]+0.0698317*left[10]+0.00885264*above[11]+0.0839893*left[11]+0.00770102*above[12]+0.0987507*left[12]+0.00671423*above[13]+0.111738*left[13]+0.00375301*above[15]+0.0820644*left[15]+0.00596388*above[14]+0.119855*left[14]);
		temp[5] = (int)(0.00179841*above[-1]+0.00511273*above[0]+0.00525562*left[0]+0.00959801*above[1]+0.0102839*left[1]+0.013104*above[2]+0.0150714*left[2]+0.0155808*above[3]+0.0198188*left[3]+0.0170783*above[4]+0.0247478*left[4]+0.0177021*above[5]+0.0300789*left[5]+0.0175928*above[6]+0.0360307*left[6]+0.016909*above[7]+0.0428136*left[7]+0.0158128*above[8]+0.0506044*left[8]+0.0144573*above[9]+0.0594815*left[9]+0.0129779*above[10]+0.0692868*left[10]+0.0114877*above[11]+0.0794182*left[11]+0.010079*above[12]+0.0886659*left[12]+0.00884995*above[13]+0.0955796*left[13]+0.00498285*above[15]+0.0665598*left[15]+0.00790053*above[14]+0.098841*left[14]);
		temp[6] = (int)(0.00204088*above[-1]+0.00581272*above[0]+0.00595289*left[0]+0.0109415*above[1]+0.0116141*left[1]+0.0150061*above[2]+0.016934*left[2]+0.0179549*above[3]+0.0221019*left[3]+0.0198363*above[4]+0.027321*left[4]+0.020752*above[5]+0.032776*left[5]+0.0208383*above[6]+0.0386228*left[6]+0.0202526*above[7]+0.044973*left[7]+0.0191605*above[8]+0.0518613*left[8]+0.0177241*above[9]+0.0591859*left[9]+0.016093*above[10]+0.0666236*left[10]+0.0143988*above[11]+0.0735496*left[11]+0.0127551*above[12]+0.0790544*left[12]+0.0112894*above[13]+0.0823577*left[13]+0.00640866*above[15]+0.0552925*left[15]+0.0101356*above[14]+0.0830485*left[14]);
		temp[7] = (int)(0.00224033*above[-1]+0.00639273*above[0]+0.00652213*left[0]+0.0120664*above[1]+0.0126872*left[1]+0.0166257*above[2]+0.0184041*left[2]+0.0200215*above[3]+0.0238423*left[3]+0.0223*above[4]+0.0291804*left[4]+0.0235542*above[5]+0.0345656*left[5]+0.023909*above[6]+0.0400993*left[6]+0.0235108*above[7]+0.0458202*left[7]+0.0225177*above[8]+0.0516767*left[8]+0.0210908*above[9]+0.0574896*left[9]+0.019386*above[10]+0.0629145*left[10]+0.0175471*above[11]+0.0674337*left[11]+0.0157066*above[12]+0.0704309*left[12]+0.0140226*above[13]+0.0715559*left[13]+0.00803147*above[15]+0.0467887*left[15]+0.0126673*above[14]+0.0708525*left[14]);
		temp[8] = (int)(0.00239568*above[-1]+0.00684856*above[0]+0.0069613*left[0]+0.0129619*above[1]+0.0135026*left[1]+0.0179416*above[2]+0.01949*left[2]+0.0217453*above[3]+0.0250685*left[3]+0.0244172*above[4]+0.0303911*left[4]+0.0260403*above[5]+0.0355725*left[5]+0.0267229*above[6]+0.0406744*left[6]+0.0265929*above[7]+0.0456911*left[7]+0.0257921*above[8]+0.050533*left[8]+0.0244706*above[9]+0.0550088*left[9]+0.0227803*above[10]+0.0588193*left[10]+0.02087*above[11]+0.0615795*left[11]+0.0188855*above[12]+0.0629032*left[12]+0.0170145*above[13]+0.0627016*left[13]+0.00983654*above[15]+0.0401957*left[15]+0.0154697*above[14]+0.0612411*left[14]);
		temp[9] = (int)(0.00250934*above[-1]+0.00718607*above[0]+0.00727859*left[0]+0.0136361*above[1]+0.0140798*left[1]+0.0189582*above[2]+0.0202282*left[2]+0.02312*above[3]+0.0258438*left[3]+0.0261659*above[4]+0.031056*left[4]+0.0281686*above[5]+0.0359547*left[5]+0.0292183*above[6]+0.040577*left[6]+0.0294204*above[7]+0.0448972*left[7]+0.0288934*above[8]+0.0488197*left[8]+0.0277677*above[9]+0.0521757*left[9]+0.0261819*above[10]+0.0547323*left[10]+0.02428*above[11]+0.0562253*left[11]+0.0222144*above[12]+0.0564336*left[12]+0.0201984*above[13]+0.0554233*left[13]+0.0117881*above[15]+0.0349911*left[15]+0.018485*above[14]+0.0535624*left[14]);
		temp[10] = (int)(0.00258641*above[-1]+0.00741884*above[0]+0.00748981*left[0]+0.0141119*above[1]+0.0144522*left[1]+0.0197007*above[2]+0.0206745*left[2]+0.024165*above[3]+0.0262527*left[3]+0.0275517*above[4]+0.0312964*left[4]+0.029925*above[5]+0.035879*left[5]+0.0313578*above[6]+0.0400251*left[6]+0.0319317*above[7]+0.0437057*left[7]+0.0317386*above[8]+0.0468377*left[8]+0.0308829*above[9]+0.0492882*left[9]+0.0294816*above[10]+0.0508912*left[10]+0.0276654*above[11]+0.0514795*left[11]+0.0255842*above[12]+0.0509465*left[12]+0.0234713*above[13]+0.0494471*left[13]+0.0138246*above[15]+0.0308438*left[15]+0.0216173*above[14]+0.047389*left[14]);
		temp[11] = (int)(0.00263373*above[-1]+0.00756562*above[0]+0.00761563*left[0]+0.0144225*above[1]+0.0146623*left[1]+0.0202092*above[2]+0.0208953*left[2]+0.024919*above[3]+0.0263895*left[3]+0.0286029*above[4]+0.0312392*left[4]+0.0313193*above[5]+0.0355074*left[5]+0.0331263*above[6]+0.039215*left[6]+0.0340833*above[7]+0.0423397*left[7]+0.034255*above[8]+0.0448183*left[8]+0.0337168*above[9]+0.0465569*left[9]+0.0325589*above[10]+0.047447*left[10]+0.0308912*above[11]+0.0473951*left[11]+0.0288534*above[12]+0.0463724*left[12]+0.0266916*above[13]+0.0445826*left[13]+0.015856*above[15]+0.0275425*left[15]+0.0247288*above[14]+0.042442*left[14]);
		temp[12] = (int)(0.00265896*above[-1]+0.00764767*above[0]+0.00767888*left[0]+0.0146064*above[1]+0.014756*left[1]+0.0205325*above[2]+0.0209606*left[2]+0.0254328*above[3]+0.02635*left[3]+0.0293632*above[4]+0.0310072*left[4]+0.0323789*above[5]+0.0349895*left[5]+0.0345266*above[6]+0.0383194*left[6]+0.0358468*above[7]+0.0409851*left[7]+0.0363793*above[8]+0.0429454*left[8]+0.0361712*above[9]+0.0441408*left[9]+0.0352837*above[10]+0.0445093*left[10]+0.0338019*above[11]+0.0440111*left[11]+0.0318498*above[12]+0.0426686*left[12]+0.0296792*above[13]+0.0407116*left[13]+0.0177629*above[15]+0.0249582*left[15]+0.0276393*above[14]+0.0385505*left[14]);
		temp[13] = (int)(0.00266973*above[-1]+0.00768639*above[0]+0.0077022*left[0]+0.0147025*above[1]+0.0147783*left[1]+0.020721*above[2]+0.0209379*left[2]+0.0257601*above[3]+0.0262249*left[3]+0.0298811*above[4]+0.0307142*left[4]+0.0331377*above[5]+0.0344605*left[5]+0.0355689*above[6]+0.03749*left[6]+0.0372003*above[7]+0.0398019*left[7]+0.0380516*above[8]+0.0413744*left[8]+0.0381449*above[9]+0.0421754*left[9]+0.037515*above[10]+0.0421774*left[10]+0.0362222*above[11]+0.0413778*left[11]+0.0343727*above[12]+0.0398323*left[12]+0.0322193*above[13]+0.0377832*left[13]+0.0193993*above[15]+0.0230257*left[15]+0.03013*above[14]+0.0356303*left[14]);
		temp[14] = (int)(0.00267282*above[-1]+0.00770071*above[0]+0.0077057*left[0]+0.0147455*above[1]+0.0147694*left[1]+0.020819*above[2]+0.0208875*left[2]+0.0259482*above[3]+0.0260948*left[3]+0.0301984*above[4]+0.0304612*left[4]+0.0336233*above[5]+0.0340405*left[5]+0.0362566*above[6]+0.0368625*left[6]+0.0381146*above[7]+0.038935*left[7]+0.0392025*above[8]+0.0402501*left[8]+0.0395243*above[9]+0.0407949*left[9]+0.0390948*above[10]+0.0405642*left[10]+0.0379545*above[11]+0.0395789*left[11]+0.0361946*above[12]+0.0379144*left[12]+0.0340662*above[13]+0.0358185*left[13]+0.0205969*above[15]+0.0217388*left[15]+0.0319492*above[14]+0.0336812*left[14]);
		temp[15] = (int)(0.00267315*above[-1]+0.00770415*above[0]+0.00770415*left[0]+0.0147592*above[1]+0.0147592*left[1]+0.0208558*above[2]+0.0208558*left[2]+0.0260246*above[3]+0.0260246*left[3]+0.0303332*above[4]+0.0303332*left[4]+0.0338352*above[5]+0.0338352*left[5]+0.0365625*above[6]+0.0365625*left[6]+0.038527*above[7]+0.038527*left[7]+0.0397273*above[8]+0.0397273*left[8]+0.0401589*above[9]+0.0401589*left[9]+0.0398268*above[10]+0.0398268*left[10]+0.0387621*above[11]+0.0387621*left[11]+0.0370482*above[12]+0.0370482*left[12]+0.0349347*above[13]+0.0349347*left[13]+0.0211621*above[15]+0.0211621*left[15]+0.0328069*above[14]+0.0328069*left[14]);	  
		memcpy(dst, temp, bs);
		dst += stride;
		break;	  
	default:
		break;
  }//switch
}

static INLINE void prdct_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  uint8_t *temp;
  int r;
  temp = (uint8_t*) malloc(bs);
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c)
    	temp[c] = (above[c]*((r+1)^3)+left[r]*((c+1)^3))/(((r+1)^3)+((c+1)^3));
    memcpy(dst, temp, bs);
    dst += stride;
  }//for
}

#else

static INLINE void tm_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                const uint8_t *above, const uint8_t *left) {
  int r, c;
  int ytop_left = above[-1];

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++)
      dst[c] = clip_pixel(left[r] + above[c] - ytop_left);
    dst += stride;
  }
}
#endif  // CONFIG_ALT_INTRA

static INLINE void dc_128_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  int r;
  (void)above;
  (void)left;

  for (r = 0; r < bs; r++) {
    memset(dst, 128, bs);
    dst += stride;
  }
}

static INLINE void dc_left_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                     const uint8_t *above,
                                     const uint8_t *left) {
  int i, r, expected_dc, sum = 0;
  (void)above;

  for (i = 0; i < bs; i++) sum += left[i];
  expected_dc = (sum + (bs >> 1)) / bs;

  for (r = 0; r < bs; r++) {
    memset(dst, expected_dc, bs);
    dst += stride;
  }
}

static INLINE void dc_top_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                    const uint8_t *above, const uint8_t *left) {
  int i, r, expected_dc, sum = 0;
  (void)left;

  for (i = 0; i < bs; i++) sum += above[i];
  expected_dc = (sum + (bs >> 1)) / bs;

  for (r = 0; r < bs; r++) {
    memset(dst, expected_dc, bs);
    dst += stride;
  }
}

static INLINE void dc_predictor(uint8_t *dst, ptrdiff_t stride, int bs,
                                const uint8_t *above, const uint8_t *left) {
  int i, r, expected_dc, sum = 0;
  const int count = 2 * bs;

  for (i = 0; i < bs; i++) {
    sum += above[i];
    sum += left[i];
  }

  expected_dc = (sum + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    memset(dst, expected_dc, bs);
    dst += stride;
  }
}

void aom_he_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                            const uint8_t *above, const uint8_t *left) {
  const int H = above[-1];
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];

  memset(dst + stride * 0, AVG3(H, I, J), 2);
  memset(dst + stride * 1, AVG3(I, J, K), 2);
}

void aom_ve_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                            const uint8_t *above, const uint8_t *left) {
  const int H = above[-1];
  const int I = above[0];
  const int J = above[1];
  const int K = above[2];
  (void)left;

  dst[0] = AVG3(H, I, J);
  dst[1] = AVG3(I, J, K);
  memcpy(dst + stride * 1, dst, 2);
}

void aom_d207_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int L = left[3];
  (void)above;
  DST(0, 0) = AVG2(I, J);
  DST(0, 1) = AVG2(J, K);
  DST(1, 0) = AVG3(I, J, K);
  DST(1, 1) = AVG3(J, K, L);
}

void aom_d63_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                             const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)left;
  DST(0, 0) = AVG2(A, B);
  DST(1, 0) = AVG2(B, C);
  DST(0, 1) = AVG3(A, B, C);
  DST(1, 1) = AVG3(B, C, D);
}

void aom_d63f_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)left;
  DST(0, 0) = AVG2(A, B);
  DST(1, 0) = AVG2(B, C);
  DST(0, 1) = AVG3(A, B, C);
  DST(1, 1) = AVG3(B, C, D);
}

void aom_d45_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                             const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)stride;
  (void)left;
  DST(0, 0) = AVG3(A, B, C);
  DST(1, 0) = DST(0, 1) = AVG3(B, C, D);
  DST(1, 1) = AVG3(C, D, D);
}

void aom_d45e_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)stride;
  (void)left;

  DST(0, 0) = AVG3(A, B, C);
  DST(1, 0) = DST(0, 1) = AVG3(B, C, D);
  DST(1, 1) = AVG3(C, D, D);
}

void aom_d117_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  DST(0, 0) = AVG2(X, A);
  DST(1, 0) = AVG2(A, B);
  DST(0, 1) = AVG3(I, X, A);
  DST(1, 1) = AVG3(X, A, B);
}

void aom_d135_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  (void)stride;
  DST(0, 1) = AVG3(X, I, J);
  DST(1, 1) = DST(0, 0) = AVG3(A, X, I);
  DST(1, 0) = AVG3(B, A, X);
}

void aom_d153_predictor_2x2_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int X = above[-1];
  const int A = above[0];

  DST(0, 0) = AVG2(I, X);
  DST(0, 1) = AVG2(J, I);
  DST(1, 0) = AVG3(I, X, A);
  DST(1, 1) = AVG3(J, I, X);
}

void aom_he_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                            const uint8_t *above, const uint8_t *left) {
  const int H = above[-1];
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int L = left[3];

  memset(dst + stride * 0, AVG3(H, I, J), 4);
  memset(dst + stride * 1, AVG3(I, J, K), 4);
  memset(dst + stride * 2, AVG3(J, K, L), 4);
  memset(dst + stride * 3, AVG3(K, L, L), 4);
}

void aom_ve_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                            const uint8_t *above, const uint8_t *left) {
  const int H = above[-1];
  const int I = above[0];
  const int J = above[1];
  const int K = above[2];
  const int L = above[3];
  const int M = above[4];
  (void)left;

  dst[0] = AVG3(H, I, J);
  dst[1] = AVG3(I, J, K);
  dst[2] = AVG3(J, K, L);
  dst[3] = AVG3(K, L, M);
  memcpy(dst + stride * 1, dst, 4);
  memcpy(dst + stride * 2, dst, 4);
  memcpy(dst + stride * 3, dst, 4);
}

void aom_d207_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int L = left[3];
  (void)above;
  DST(0, 0) = AVG2(I, J);
  DST(2, 0) = DST(0, 1) = AVG2(J, K);
  DST(2, 1) = DST(0, 2) = AVG2(K, L);
  DST(1, 0) = AVG3(I, J, K);
  DST(3, 0) = DST(1, 1) = AVG3(J, K, L);
  DST(3, 1) = DST(1, 2) = AVG3(K, L, L);
  DST(3, 2) = DST(2, 2) = DST(0, 3) = DST(1, 3) = DST(2, 3) = DST(3, 3) = L;
}

void aom_d63_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                             const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  const int E = above[4];
  const int F = above[5];
  const int G = above[6];
  (void)left;
  DST(0, 0) = AVG2(A, B);
  DST(1, 0) = DST(0, 2) = AVG2(B, C);
  DST(2, 0) = DST(1, 2) = AVG2(C, D);
  DST(3, 0) = DST(2, 2) = AVG2(D, E);
  DST(3, 2) = AVG2(E, F);  // differs from vp8

  DST(0, 1) = AVG3(A, B, C);
  DST(1, 1) = DST(0, 3) = AVG3(B, C, D);
  DST(2, 1) = DST(1, 3) = AVG3(C, D, E);
  DST(3, 1) = DST(2, 3) = AVG3(D, E, F);
  DST(3, 3) = AVG3(E, F, G);  // differs from vp8
}

void aom_d63f_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  const int E = above[4];
  const int F = above[5];
  const int G = above[6];
  const int H = above[7];
  (void)left;
  DST(0, 0) = AVG2(A, B);
  DST(1, 0) = DST(0, 2) = AVG2(B, C);
  DST(2, 0) = DST(1, 2) = AVG2(C, D);
  DST(3, 0) = DST(2, 2) = AVG2(D, E);
  DST(3, 2) = AVG3(E, F, G);

  DST(0, 1) = AVG3(A, B, C);
  DST(1, 1) = DST(0, 3) = AVG3(B, C, D);
  DST(2, 1) = DST(1, 3) = AVG3(C, D, E);
  DST(3, 1) = DST(2, 3) = AVG3(D, E, F);
  DST(3, 3) = AVG3(F, G, H);
}

void aom_d45_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                             const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  const int E = above[4];
  const int F = above[5];
  const int G = above[6];
  const int H = above[7];
  (void)stride;
  (void)left;
  DST(0, 0) = AVG3(A, B, C);
  DST(1, 0) = DST(0, 1) = AVG3(B, C, D);
  DST(2, 0) = DST(1, 1) = DST(0, 2) = AVG3(C, D, E);
  DST(3, 0) = DST(2, 1) = DST(1, 2) = DST(0, 3) = AVG3(D, E, F);
  DST(3, 1) = DST(2, 2) = DST(1, 3) = AVG3(E, F, G);
  DST(3, 2) = DST(2, 3) = AVG3(F, G, H);
  DST(3, 3) = H;  // differs from vp8
}

void aom_d45e_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  const int E = above[4];
  const int F = above[5];
  const int G = above[6];
  const int H = above[7];
  (void)stride;
  (void)left;
  DST(0, 0) = AVG3(A, B, C);
  DST(1, 0) = DST(0, 1) = AVG3(B, C, D);
  DST(2, 0) = DST(1, 1) = DST(0, 2) = AVG3(C, D, E);
  DST(3, 0) = DST(2, 1) = DST(1, 2) = DST(0, 3) = AVG3(D, E, F);
  DST(3, 1) = DST(2, 2) = DST(1, 3) = AVG3(E, F, G);
  DST(3, 2) = DST(2, 3) = AVG3(F, G, H);
  DST(3, 3) = AVG3(G, H, H);
}

void aom_d117_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  DST(0, 0) = DST(1, 2) = AVG2(X, A);
  DST(1, 0) = DST(2, 2) = AVG2(A, B);
  DST(2, 0) = DST(3, 2) = AVG2(B, C);
  DST(3, 0) = AVG2(C, D);

  DST(0, 3) = AVG3(K, J, I);
  DST(0, 2) = AVG3(J, I, X);
  DST(0, 1) = DST(1, 3) = AVG3(I, X, A);
  DST(1, 1) = DST(2, 3) = AVG3(X, A, B);
  DST(2, 1) = DST(3, 3) = AVG3(A, B, C);
  DST(3, 1) = AVG3(B, C, D);
}

void aom_d135_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int L = left[3];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)stride;
  DST(0, 3) = AVG3(J, K, L);
  DST(1, 3) = DST(0, 2) = AVG3(I, J, K);
  DST(2, 3) = DST(1, 2) = DST(0, 1) = AVG3(X, I, J);
  DST(3, 3) = DST(2, 2) = DST(1, 1) = DST(0, 0) = AVG3(A, X, I);
  DST(3, 2) = DST(2, 1) = DST(1, 0) = AVG3(B, A, X);
  DST(3, 1) = DST(2, 0) = AVG3(C, B, A);
  DST(3, 0) = AVG3(D, C, B);
}

void aom_d153_predictor_4x4_c(uint8_t *dst, ptrdiff_t stride,
                              const uint8_t *above, const uint8_t *left) {
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int L = left[3];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];

  DST(0, 0) = DST(2, 1) = AVG2(I, X);
  DST(0, 1) = DST(2, 2) = AVG2(J, I);
  DST(0, 2) = DST(2, 3) = AVG2(K, J);
  DST(0, 3) = AVG2(L, K);

  DST(3, 0) = AVG3(A, B, C);
  DST(2, 0) = AVG3(X, A, B);
  DST(1, 0) = DST(3, 1) = AVG3(I, X, A);
  DST(1, 1) = DST(3, 2) = AVG3(J, I, X);
  DST(1, 2) = DST(3, 3) = AVG3(K, J, I);
  DST(1, 3) = AVG3(L, K, J);
}

#if CONFIG_AOM_HIGHBITDEPTH
static INLINE void highbd_d207_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  int r, c;
  (void)above;
  (void)bd;

  // First column.
  for (r = 0; r < bs - 1; ++r) {
    dst[r * stride] = AVG2(left[r], left[r + 1]);
  }
  dst[(bs - 1) * stride] = left[bs - 1];
  dst++;

  // Second column.
  for (r = 0; r < bs - 2; ++r) {
    dst[r * stride] = AVG3(left[r], left[r + 1], left[r + 2]);
  }
  dst[(bs - 2) * stride] = AVG3(left[bs - 2], left[bs - 1], left[bs - 1]);
  dst[(bs - 1) * stride] = left[bs - 1];
  dst++;

  // Rest of last row.
  for (c = 0; c < bs - 2; ++c) dst[(bs - 1) * stride + c] = left[bs - 1];

  for (r = bs - 2; r >= 0; --r) {
    for (c = 0; c < bs - 2; ++c)
      dst[r * stride + c] = dst[(r + 1) * stride + c - 2];
  }
}

static INLINE void highbd_d207e_predictor(uint16_t *dst, ptrdiff_t stride,
                                          int bs, const uint16_t *above,
                                          const uint16_t *left, int bd) {
  int r, c;
  (void)above;
  (void)bd;

  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = c & 1 ? AVG3(left[(c >> 1) + r], left[(c >> 1) + r + 1],
                            left[(c >> 1) + r + 2])
                     : AVG2(left[(c >> 1) + r], left[(c >> 1) + r + 1]);
    }
    dst += stride;
  }
}

static INLINE void highbd_d63_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                        const uint16_t *above,
                                        const uint16_t *left, int bd) {
  int r, c;
  (void)left;
  (void)bd;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = r & 1 ? AVG3(above[(r >> 1) + c], above[(r >> 1) + c + 1],
                            above[(r >> 1) + c + 2])
                     : AVG2(above[(r >> 1) + c], above[(r >> 1) + c + 1]);
    }
    dst += stride;
  }
}

#define highbd_d63e_predictor highbd_d63_predictor

static INLINE void highbd_d45_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                        const uint16_t *above,
                                        const uint16_t *left, int bd) {
  int r, c;
  (void)left;
  (void)bd;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = r + c + 2 < bs * 2
                   ? AVG3(above[r + c], above[r + c + 1], above[r + c + 2])
                   : above[bs * 2 - 1];
    }
    dst += stride;
  }
}

static INLINE void highbd_d45e_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  int r, c;
  (void)left;
  (void)bd;
  for (r = 0; r < bs; ++r) {
    for (c = 0; c < bs; ++c) {
      dst[c] = AVG3(above[r + c], above[r + c + 1],
                    above[r + c + 1 + (r + c + 2 < bs * 2)]);
    }
    dst += stride;
  }
}

static INLINE void highbd_d117_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  int r, c;
  (void)bd;

  // first row
  for (c = 0; c < bs; c++) dst[c] = AVG2(above[c - 1], above[c]);
  dst += stride;

  // second row
  dst[0] = AVG3(left[0], above[-1], above[0]);
  for (c = 1; c < bs; c++) dst[c] = AVG3(above[c - 2], above[c - 1], above[c]);
  dst += stride;

  // the rest of first col
  dst[0] = AVG3(above[-1], left[0], left[1]);
  for (r = 3; r < bs; ++r)
    dst[(r - 2) * stride] = AVG3(left[r - 3], left[r - 2], left[r - 1]);

  // the rest of the block
  for (r = 2; r < bs; ++r) {
    for (c = 1; c < bs; c++) dst[c] = dst[-2 * stride + c - 1];
    dst += stride;
  }
}

static INLINE void highbd_d135_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  int r, c;
  (void)bd;
  dst[0] = AVG3(left[0], above[-1], above[0]);
  for (c = 1; c < bs; c++) dst[c] = AVG3(above[c - 2], above[c - 1], above[c]);

  dst[stride] = AVG3(above[-1], left[0], left[1]);
  for (r = 2; r < bs; ++r)
    dst[r * stride] = AVG3(left[r - 2], left[r - 1], left[r]);

  dst += stride;
  for (r = 1; r < bs; ++r) {
    for (c = 1; c < bs; c++) dst[c] = dst[-stride + c - 1];
    dst += stride;
  }
}

static INLINE void highbd_d153_predictor(uint16_t *dst, ptrdiff_t stride,
                                         int bs, const uint16_t *above,
                                         const uint16_t *left, int bd) {
  int r, c;
  (void)bd;
  dst[0] = AVG2(above[-1], left[0]);
  for (r = 1; r < bs; r++) dst[r * stride] = AVG2(left[r - 1], left[r]);
  dst++;

  dst[0] = AVG3(left[0], above[-1], above[0]);
  dst[stride] = AVG3(above[-1], left[0], left[1]);
  for (r = 2; r < bs; r++)
    dst[r * stride] = AVG3(left[r - 2], left[r - 1], left[r]);
  dst++;

  for (c = 0; c < bs - 2; c++)
    dst[c] = AVG3(above[c - 1], above[c], above[c + 1]);
  dst += stride;

  for (r = 1; r < bs; ++r) {
    for (c = 0; c < bs - 2; c++) dst[c] = dst[-stride + c - 2];
    dst += stride;
  }
}

static INLINE void highbd_v_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                      const uint16_t *above,
                                      const uint16_t *left, int bd) {
  int r;
  (void)left;
  (void)bd;
  for (r = 0; r < bs; r++) {
    memcpy(dst, above, bs * sizeof(uint16_t));
    dst += stride;
  }
}

static INLINE void highbd_h_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                      const uint16_t *above,
                                      const uint16_t *left, int bd) {
  int r;
  (void)above;
  (void)bd;
  for (r = 0; r < bs; r++) {
    aom_memset16(dst, left[r], bs);
    dst += stride;
  }
}

void aom_highbd_d207_predictor_2x2_c(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const int I = left[0];
  const int J = left[1];
  const int K = left[2];
  const int L = left[3];
  (void)above;
  (void)bd;
  DST(0, 0) = AVG2(I, J);
  DST(0, 1) = AVG2(J, K);
  DST(1, 0) = AVG3(I, J, K);
  DST(1, 1) = AVG3(J, K, L);
}

void aom_highbd_d63_predictor_2x2_c(uint16_t *dst, ptrdiff_t stride,
                                    const uint16_t *above, const uint16_t *left,
                                    int bd) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)left;
  (void)bd;
  DST(0, 0) = AVG2(A, B);
  DST(1, 0) = AVG2(B, C);
  DST(0, 1) = AVG3(A, B, C);
  DST(1, 1) = AVG3(B, C, D);
}

void aom_highbd_d45e_predictor_2x2_c(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const int A = above[0];
  const int B = above[1];
  const int C = above[2];
  const int D = above[3];
  (void)stride;
  (void)left;
  (void)bd;
  DST(0, 0) = AVG3(A, B, C);
  DST(1, 0) = DST(0, 1) = AVG3(B, C, D);
  DST(1, 1) = AVG3(C, D, D);
}

void aom_highbd_d117_predictor_2x2_c(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const int I = left[0];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  (void)bd;
  DST(0, 0) = AVG2(X, A);
  DST(1, 0) = AVG2(A, B);
  DST(0, 1) = AVG3(I, X, A);
  DST(1, 1) = AVG3(X, A, B);
}

void aom_highbd_d135_predictor_2x2_c(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const int I = left[0];
  const int J = left[1];
  const int X = above[-1];
  const int A = above[0];
  const int B = above[1];
  (void)bd;
  DST(0, 1) = AVG3(X, I, J);
  DST(1, 1) = DST(0, 0) = AVG3(A, X, I);
  DST(1, 0) = AVG3(B, A, X);
}

void aom_highbd_d153_predictor_2x2_c(uint16_t *dst, ptrdiff_t stride,
                                     const uint16_t *above,
                                     const uint16_t *left, int bd) {
  const int I = left[0];
  const int J = left[1];
  const int X = above[-1];
  const int A = above[0];
  (void)bd;
  DST(0, 0) = AVG2(I, X);
  DST(0, 1) = AVG2(J, I);
  DST(1, 0) = AVG3(I, X, A);
  DST(1, 1) = AVG3(J, I, X);
}

#if CONFIG_ALT_INTRA
static INLINE void highbd_paeth_predictor(uint16_t *dst, ptrdiff_t stride,
                                          int bs, const uint16_t *above,
                                          const uint16_t *left, int bd) {
  int r, c;
  const uint16_t ytop_left = above[-1];
  (void)bd;

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++)
      dst[c] = paeth_predictor_single(left[r], above[c], ytop_left);
    dst += stride;
  }
}

static INLINE void highbd_smooth_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16_t below_pred = left[bs - 1];   // estimated by bottom-left pixel
  const uint16_t right_pred = above[bs - 1];  // estimated by top-right pixel
  const int arr_index = (int)lround(log2(bs)) - 1;
  const double *const fwd_weights = sm_weights_fwd[arr_index];
  const double *const rev_weights = sm_weights_rev[arr_index];
  const double scale = 2.0 * bs;
  int r;
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c) {
      const int pixels[] = { above[c], below_pred, left[r], right_pred };
      const double weights[] = { fwd_weights[r], rev_weights[r], fwd_weights[c],
                                 rev_weights[c] };
      double this_pred = 0;
      int i;
      for (i = 0; i < 4; ++i) {
        this_pred += weights[i] * pixels[i];
      }
      dst[c] = clip_pixel_highbd(lround(this_pred / scale), bd);
    }
    dst += stride;
  }
}

static INLINE void highbd_least_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16_t below_pred = left[bs - 1];   // estimated by bottom-left pixel
  const uint16_t right_pred = above[bs - 1];  // estimated by top-right pixel
  const int arr_index = (int)lround(log2(bs)) - 1;
  const double *const fwd_weights = sm_weights_fwd[arr_index];
  const double *const rev_weights = sm_weights_rev[arr_index];
  const double scale = 2.0 * bs;
  int r;
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c) {
      const int pixels[] = { above[c], below_pred, left[r], right_pred };
      const double weights[] = { fwd_weights[r], rev_weights[r], fwd_weights[c],
                                 rev_weights[c] };
      double this_pred = 0;
      int i;
      for (i = 0; i < 4; ++i) {
        this_pred += weights[i] * pixels[i];
      }
      dst[c] = clip_pixel_highbd(lround(this_pred / scale), bd);
    }
    dst += stride;
  }
}

static INLINE void highbd_wcalic_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16_t below_pred = left[bs - 1];   // estimated by bottom-left pixel
  const uint16_t right_pred = above[bs - 1];  // estimated by top-right pixel
  const int arr_index = (int)lround(log2(bs)) - 1;
  const double *const fwd_weights = sm_weights_fwd[arr_index];
  const double *const rev_weights = sm_weights_rev[arr_index];
  const double scale = 2.0 * bs;
  int r;
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c) {
      const int pixels[] = { above[c], below_pred, left[r], right_pred };
      const double weights[] = { fwd_weights[r], rev_weights[r], fwd_weights[c],
                                 rev_weights[c] };
      double this_pred = 0;
      int i;
      for (i = 0; i < 4; ++i) {
        this_pred += weights[i] * pixels[i];
      }
      dst[c] = clip_pixel_highbd(lround(this_pred / scale), bd);
    }
    dst += stride;
  }
}

static INLINE void highbd_isle_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16_t below_pred = left[bs - 1];   // estimated by bottom-left pixel
  const uint16_t right_pred = above[bs - 1];  // estimated by top-right pixel
  const int arr_index = (int)lround(log2(bs)) - 1;
  const double *const fwd_weights = sm_weights_fwd[arr_index];
  const double *const rev_weights = sm_weights_rev[arr_index];
  const double scale = 2.0 * bs;
  int r;
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c) {
      const int pixels[] = { above[c], below_pred, left[r], right_pred };
      const double weights[] = { fwd_weights[r], rev_weights[r], fwd_weights[c],
                                 rev_weights[c] };
      double this_pred = 0;
      int i;
      for (i = 0; i < 4; ++i) {
        this_pred += weights[i] * pixels[i];
      }
      dst[c] = clip_pixel_highbd(lround(this_pred / scale), bd);
    }
    dst += stride;
  }
}

static INLINE void highbd_prdct_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  const uint16_t below_pred = left[bs - 1];   // estimated by bottom-left pixel
  const uint16_t right_pred = above[bs - 1];  // estimated by top-right pixel
  const int arr_index = (int)lround(log2(bs)) - 1;
  const double *const fwd_weights = sm_weights_fwd[arr_index];
  const double *const rev_weights = sm_weights_rev[arr_index];
  const double scale = 2.0 * bs;
  int r;
  for (r = 0; r < bs; ++r) {
    int c;
    for (c = 0; c < bs; ++c) {
      const int pixels[] = { above[c], below_pred, left[r], right_pred };
      const double weights[] = { fwd_weights[r], rev_weights[r], fwd_weights[c],
                                 rev_weights[c] };
      double this_pred = 0;
      int i;
      for (i = 0; i < 4; ++i) {
        this_pred += weights[i] * pixels[i];
      }
      dst[c] = clip_pixel_highbd(lround(this_pred / scale), bd);
    }
    dst += stride;
  }
}

#else
static INLINE void highbd_tm_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                       const uint16_t *above,
                                       const uint16_t *left, int bd) {
  int r, c;
  int ytop_left = above[-1];
  (void)bd;

  for (r = 0; r < bs; r++) {
    for (c = 0; c < bs; c++)
      dst[c] = clip_pixel_highbd(left[r] + above[c] - ytop_left, bd);
    dst += stride;
  }
}
#endif  // CONFIG_ALT_INTRA

static INLINE void highbd_dc_128_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  int r;
  (void)above;
  (void)left;

  for (r = 0; r < bs; r++) {
    aom_memset16(dst, 128 << (bd - 8), bs);
    dst += stride;
  }
}

static INLINE void highbd_dc_left_predictor(uint16_t *dst, ptrdiff_t stride,
                                            int bs, const uint16_t *above,
                                            const uint16_t *left, int bd) {
  int i, r, expected_dc, sum = 0;
  (void)above;
  (void)bd;

  for (i = 0; i < bs; i++) sum += left[i];
  expected_dc = (sum + (bs >> 1)) / bs;

  for (r = 0; r < bs; r++) {
    aom_memset16(dst, expected_dc, bs);
    dst += stride;
  }
}

static INLINE void highbd_dc_top_predictor(uint16_t *dst, ptrdiff_t stride,
                                           int bs, const uint16_t *above,
                                           const uint16_t *left, int bd) {
  int i, r, expected_dc, sum = 0;
  (void)left;
  (void)bd;

  for (i = 0; i < bs; i++) sum += above[i];
  expected_dc = (sum + (bs >> 1)) / bs;

  for (r = 0; r < bs; r++) {
    aom_memset16(dst, expected_dc, bs);
    dst += stride;
  }
}

static INLINE void highbd_dc_predictor(uint16_t *dst, ptrdiff_t stride, int bs,
                                       const uint16_t *above,
                                       const uint16_t *left, int bd) {
  int i, r, expected_dc, sum = 0;
  const int count = 2 * bs;
  (void)bd;

  for (i = 0; i < bs; i++) {
    sum += above[i];
    sum += left[i];
  }

  expected_dc = (sum + (count >> 1)) / count;

  for (r = 0; r < bs; r++) {
    aom_memset16(dst, expected_dc, bs);
    dst += stride;
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

// This serves as a wrapper function, so that all the prediction functions
// can be unified and accessed as a pointer array. Note that the boundary
// above and left are not necessarily used all the time.
#define intra_pred_sized(type, size)                        \
  void aom_##type##_predictor_##size##x##size##_c(          \
      uint8_t *dst, ptrdiff_t stride, const uint8_t *above, \
      const uint8_t *left) {                                \
    type##_predictor(dst, stride, size, above, left);       \
  }

#if CONFIG_AOM_HIGHBITDEPTH
#define intra_pred_highbd_sized(type, size)                        \
  void aom_highbd_##type##_predictor_##size##x##size##_c(          \
      uint16_t *dst, ptrdiff_t stride, const uint16_t *above,      \
      const uint16_t *left, int bd) {                              \
    highbd_##type##_predictor(dst, stride, size, above, left, bd); \
  }

/* clang-format off */
#if CONFIG_TX64X64
#define intra_pred_allsizes(type) \
  intra_pred_sized(type, 2) \
  intra_pred_sized(type, 4) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32) \
  intra_pred_sized(type, 64) \
  intra_pred_highbd_sized(type, 4) \
  intra_pred_highbd_sized(type, 8) \
  intra_pred_highbd_sized(type, 16) \
  intra_pred_highbd_sized(type, 32) \
  intra_pred_highbd_sized(type, 64)

#define intra_pred_above_4x4(type) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32) \
  intra_pred_sized(type, 64) \
  intra_pred_highbd_sized(type, 4) \
  intra_pred_highbd_sized(type, 8) \
  intra_pred_highbd_sized(type, 16) \
  intra_pred_highbd_sized(type, 32) \
  intra_pred_highbd_sized(type, 64)
#else  // CONFIG_TX64X64
#define intra_pred_allsizes(type) \
  intra_pred_sized(type, 2) \
  intra_pred_sized(type, 4) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32) \
  intra_pred_highbd_sized(type, 2) \
  intra_pred_highbd_sized(type, 4) \
  intra_pred_highbd_sized(type, 8) \
  intra_pred_highbd_sized(type, 16) \
  intra_pred_highbd_sized(type, 32)

#define intra_pred_above_4x4(type) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32) \
  intra_pred_highbd_sized(type, 4) \
  intra_pred_highbd_sized(type, 8) \
  intra_pred_highbd_sized(type, 16) \
  intra_pred_highbd_sized(type, 32)
#endif  // CONFIG_TX64X64

#else

#if CONFIG_TX64X64
#define intra_pred_allsizes(type) \
  intra_pred_sized(type, 2) \
  intra_pred_sized(type, 4) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32) \
  intra_pred_sized(type, 64)

#define intra_pred_above_4x4(type) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32) \
  intra_pred_sized(type, 64)
#else  // CONFIG_TX64X64
#define intra_pred_allsizes(type) \
  intra_pred_sized(type, 2) \
  intra_pred_sized(type, 4) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32)

#define intra_pred_above_4x4(type) \
  intra_pred_sized(type, 8) \
  intra_pred_sized(type, 16) \
  intra_pred_sized(type, 32)
#endif  // CONFIG_TX64X64
#endif  // CONFIG_AOM_HIGHBITDEPTH

intra_pred_above_4x4(d207)
intra_pred_above_4x4(d63)
intra_pred_above_4x4(d45)
intra_pred_allsizes(d207e)
intra_pred_allsizes(d63e)
intra_pred_above_4x4(d45e)
intra_pred_above_4x4(d117)
intra_pred_above_4x4(d135)
intra_pred_above_4x4(d153)
intra_pred_allsizes(v)
intra_pred_allsizes(h)
#if CONFIG_ALT_INTRA
intra_pred_allsizes(paeth)
intra_pred_allsizes(smooth)
intra_pred_allsizes(least)
intra_pred_allsizes(wcalic)
intra_pred_allsizes(isle)
intra_pred_allsizes(prdct)
#else
intra_pred_allsizes(tm)
#endif  // CONFIG_ALT_INTRA
intra_pred_allsizes(dc_128)
intra_pred_allsizes(dc_left)
intra_pred_allsizes(dc_top)
intra_pred_allsizes(dc)
/* clang-format on */
#undef intra_pred_allsizes
