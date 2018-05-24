#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "./av1_rtcd.h"
#include "av1/encoder/dwt.h"

// Note: block length must be even for this implementation
static void analysis_53_row(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int n;
  tran_low_t r, *a, *b;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    *a++ = (r = *x++) << 1;
    *b++ = *x - ((r + x[1] + 1) >> 1);
    x++;
  }
  *a = (r = *x++) << 1;
  *b = *x - r;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ += (r + (*b) + 1) >> 1;
    r = *b++;
  }
}

static void analysis_53_col(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int n;
  tran_low_t r, *a, *b;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  while (--n) {
    *a++ = (r = *x++);
    *b++ = (((*x) << 1) - (r + x[1]) + 2) >> 2;
    x++;
  }
  *a = (r = *x++);
  *b = (*x - r + 1) >> 1;

  n = length >> 1;
  b = highpass;
  a = lowpass;
  r = *highpass;
  while (n--) {
    *a++ += (r + (*b) + 1) >> 1;
    r = *b++;
  }
}

static void dyadic_analyze_53(int levels, int width, int height,
                              tran_low_t *x, int pitch_x,
                              tran_low_t *c, int pitch_c,
                              int dwt_scale_bits) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      c[i * pitch_c + j] = x[i * pitch_x + j] << dwt_scale_bits;
    }
  }
  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(tran_low_t));
      analysis_53_row(nw, buffer, &c[i * pitch_c], &c[i * pitch_c] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = c[i * pitch_c + j];
      analysis_53_col(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i];
    }
  }
}

static void dyadic_analyze_53_uint8_input(int levels, int width, int height,
                                          uint8_t *x, int pitch_x,
                                          tran_low_t *c, int pitch_c,
                                          int dwt_scale_bits, int hbd) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];

  if (hbd) {
    uint16_t *x16 = CONVERT_TO_SHORTPTR(x);
    for (i = 0; i < height; i++) {
      for (j = 0; j < width; j++) {
        c[i * pitch_c + j] = x16[i * pitch_x + j] << dwt_scale_bits;
      }
    }
  } else {
    for (i = 0; i < height; i++) {
      for (j = 0; j < width; j++) {
        c[i * pitch_c + j] = x[i * pitch_x + j] << dwt_scale_bits;
      }
    }
  }

  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(tran_low_t));
      analysis_53_row(nw, buffer, &c[i * pitch_c], &c[i * pitch_c] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = c[i * pitch_c + j];
      analysis_53_col(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i];
    }
  }
}

static void analysis_26_row(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int i, n;
  tran_low_t r, s, *a, *b;
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    r = *x++;
    s = *x++;
    *a++ = r + s;
    *b++ = r - s;
  }
  n = length >> 1;
  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ -= (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b -= (r - *a + 4) >> 3;
  }
}

static void analysis_26_col(int length, tran_low_t *x,
                            tran_low_t *lowpass, tran_low_t *highpass) {
  int i, n;
  tran_low_t r, s, *a, *b;
  a = lowpass;
  b = highpass;
  for (i = length >> 1; i; i--) {
    r = *x++;
    s = *x++;
    *a++ = (r + s + 1) >> 1;
    *b++ = (r - s + 1) >> 1;
  }
  n = length >> 1;
  if (n >= 4) {
    a = lowpass;
    b = highpass;
    r = *lowpass;
    while (--n) {
      *b++ -= (r - a[1] + 4) >> 3;
      r = *a++;
    }
    *b -= (r - *a + 4) >> 3;
  }
}

static void dyadic_analyze_26(int levels, int width, int height,
                              int16_t *x, int pitch_x,
                              tran_low_t *c, int pitch_c,
                              int dwt_scale_bits) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  tran_low_t buffer[2 * DWT_MAX_LENGTH];
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      c[i * pitch_c + j] = x[i * pitch_x + j] << dwt_scale_bits;
    }
  }
  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &c[i * pitch_c], nw * sizeof(tran_low_t));
      analysis_26_row(nw, buffer, &c[i * pitch_c], &c[i * pitch_c] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = c[i * pitch_c + j];
      analysis_26_col(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        c[i * pitch_c + j] = buffer[i];
    }
  }
}

static void analysis_97(int length, double *x,
                        double *lowpass, double *highpass) {
  static const double a_predict1 = -1.586134342;
  static const double a_update1 = -0.05298011854;
  static const double a_predict2 = 0.8829110762;
  static const double a_update2 = 0.4435068522;
  static const double s_low = 1.149604398;
  static const double s_high = 1/1.149604398;
  int i;
  double y[DWT_MAX_LENGTH];
  // Predict 1
  for (i = 1; i < length - 2; i += 2) {
    x[i] += a_predict1 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] += 2 * a_predict1 * x[length - 2];
  // Update 1
  for (i = 2; i < length; i += 2) {
    x[i] += a_update1 * (x[i - 1] + x[i + 1]);
  }
  x[0] += 2 * a_update1 * x[1];
  // Predict 2
  for (i = 1; i < length - 2; i += 2) {
    x[i] += a_predict2 * (x[i - 1] + x[i + 1]);
  }
  x[length - 1] += 2 * a_predict2 * x[length - 2];
  // Update 2
  for (i = 2; i < length; i += 2) {
    x[i] += a_update2 * (x[i - 1] + x[i + 1]);
  }
  x[0] += 2 * a_update2 * x[1];
  memcpy(y, x, sizeof(*y) * length);
  // Scale and pack
  for (i = 0; i < length / 2; i++) {
    lowpass[i] = y[2 * i] * s_low;
    highpass[i] = y[2 * i + 1] * s_high;
  }
}

static void dyadic_analyze_97(int levels, int width, int height,
                              int16_t *x, int pitch_x,
                              tran_low_t *c, int pitch_c,
                              int dwt_scale_bits) {
  int lv, i, j, nh, nw, hh = height, hw = width;
  double buffer[2 * DWT_MAX_LENGTH];
  double y[DWT_MAX_LENGTH * DWT_MAX_LENGTH];
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      y[i * DWT_MAX_LENGTH + j] = x[i * pitch_x + j] << dwt_scale_bits;
    }
  }
  for (lv = 0; lv < levels; lv++) {
    nh = hh;
    hh = (hh + 1) >> 1;
    nw = hw;
    hw = (hw + 1) >> 1;
    if ((nh < 2) || (nw < 2)) return;
    for (i = 0; i < nh; i++) {
      memcpy(buffer, &y[i * DWT_MAX_LENGTH], nw * sizeof(*buffer));
      analysis_97(nw, buffer, &y[i * DWT_MAX_LENGTH],
                  &y[i * DWT_MAX_LENGTH] + hw);
    }
    for (j = 0; j < nw; j++) {
      for (i = 0; i < nh; i++)
        buffer[i + nh] = y[i * DWT_MAX_LENGTH + j];
      analysis_97(nh, buffer + nh, buffer, buffer + hh);
      for (i = 0; i < nh; i++)
        y[i * DWT_MAX_LENGTH + j] = buffer[i];
    }
  }
  for (i = 0; i < height; i++) {
    for (j = 0; j < width; j++) {
      c[i * pitch_c + j] = round(y[i * DWT_MAX_LENGTH + j]);
    }
  }
}

void av1_fdwt8x8_c(tran_low_t *input, tran_low_t *output, int stride) {
#if DWT_TYPE == 26
  dyadic_analyze_26(4, 8, 8, input, stride, output, 8, 2);
#elif DWT_TYPE == 97
  dyadic_analyze_97(4, 8, 8, input, stride, output, 8, 2);
#elif DWT_TYPE == 53
  dyadic_analyze_53(4, 8, 8, input, stride, output, 8, 2);
#endif
}

void av1_fdwt8x8_uint8_input_c(uint8_t *input, tran_low_t *output, int stride,
                               int hbd) {
#if DWT_TYPE == 26
  dyadic_analyze_26(4, 8, 8, input, stride, output, 8, 2);
#elif DWT_TYPE == 97
  dyadic_analyze_97(4, 8, 8, input, stride, output, 8, 2);
#elif DWT_TYPE == 53
  dyadic_analyze_53_uint8_input(4, 8, 8, input, stride, output, 8, 2, hbd);
#endif
}

int av1_haar_ac_sad(tran_low_t *output, int bw, int bh, int stride) {
  int acsad = 0;

  for (int r = 0; r < bh; ++r)
    for (int c = 0; c < bw; ++c) {
      if (r >= bh / 2 || c >= bw / 2)
        acsad += abs(output[r * stride + c]);
    }
  return acsad;
}

uint64_t av1_dct_ac_sad(tran_low_t *output, int bw, int bh, int stride) {
  uint64_t acsad = 0;

  for (int r = 0; r < bh; ++r)
    for (int c = 0; c < bw; ++c) {
      if (r > 0 || c > 0)
        acsad += abs(output[r * stride + c]);
    }

  return acsad;
}

uint32_t av1_variance(uint8_t *input, int bw, int bh, int stride) {
  int sum = 0;
  uint32_t sse = 0;

  for (int r = 0; r < bh; ++r)
    for (int c = 0; c < bw; ++c) {
      sum += input[r * stride + c];
      sse += input[r * stride + c] * input[r * stride + c];
    }
  return sse - (uint32_t)(((int64_t)sum * sum)/(bw * bh));
}

int av1_haar_ac_sad_8x8_uint8_input(uint8_t *input, int stride, int hbd) {
  tran_low_t output[64];

  av1_fdwt8x8_uint8_input_c(input, output, stride, hbd);
  return av1_haar_ac_sad(output, 8, 8, 8);
}
