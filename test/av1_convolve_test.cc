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

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./av1_rtcd.h"
#include "./aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "av1/common/filter.h"
#include "av1/common/convolve.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_ports/mem.h"
#include <algorithm>

using libaom_test::ACMRandom;

namespace {

static void filter_block1d_horiz_c(const uint8_t *src_ptr, int src_stride,
                            const int16_t *filter, int tap,
                            uint8_t *dst_ptr, int dst_stride,
                            int w, int h) {
  src_ptr -= tap / 2 - 1;
  for (int r = 0; r < h; ++r) {
    for (int c = 0; c < w; ++c) {
      int sum = 0;
      for (int i = 0; i < tap; ++i) {
        sum += src_ptr[c + i] * filter[i];
      }
      dst_ptr[c] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
    src_ptr += src_stride;
    dst_ptr += dst_stride;
  }
}

static void filter_block1d_vert_c(const uint8_t *src_ptr, int src_stride,
                           const int16_t *filter, int tap,
                           uint8_t *dst_ptr, int dst_stride,
                           int w, int h) {
  src_ptr -= (tap / 2 - 1) * src_stride;
  for (int r = 0; r < h; ++r) {
    for (int c = 0; c < w; ++c) {
      int sum = 0;
      for (int i = 0; i < tap; ++i) {
        sum += src_ptr[c + i * src_stride] * filter[i];
      }
      dst_ptr[c] = clip_pixel(ROUND_POWER_OF_TWO(sum, FILTER_BITS));
    }
    src_ptr += src_stride;
    dst_ptr += dst_stride;
  }
}

static void filter_block2d_c(const uint8_t *src_ptr, const unsigned int src_stride,
                             uint8_t *dst_ptr, unsigned int dst_stride,
                      unsigned int w, unsigned int h,
                      const int16_t *filter_horiz, int tap_horiz,
                      const int16_t *filter_vert, int tap_vert) {
  if (tap_horiz >= tap_vert) {
    uint8_t *temp_ptr = new uint8_t[w * (h + tap_vert - 1)];
    const unsigned int temp_stride = w;
    filter_block1d_horiz_c(src_ptr - (tap_vert/2 - 1) * src_stride, src_stride, filter_horiz, tap_horiz, temp_ptr, temp_stride, w, h + tap_vert - 1);
    filter_block1d_vert_c(temp_ptr + (tap_vert/2 - 1) * temp_stride, temp_stride, filter_vert, tap_vert, dst_ptr, dst_stride, w, h);
    delete[] temp_ptr;
  } else {
    uint8_t *temp_ptr = new uint8_t[(w + tap_horiz - 1) * h];
    const unsigned int temp_stride = w + tap_horiz - 1;
    filter_block1d_vert_c(src_ptr - (tap_horiz/2 - 1), src_stride, filter_vert, tap_vert, temp_ptr, temp_stride, w + tap_horiz - 1, h);
    filter_block1d_horiz_c(temp_ptr + (tap_horiz/2 - 1), temp_stride, filter_horiz, tap_horiz, dst_ptr, dst_stride, w, h);
    delete[] temp_ptr;
  }
}

static int match(const uint8_t* out, int out_stride, const uint8_t* ref_out, int ref_out_stride, int w, int h) {
  for (int r = 0; r < h; ++r) {
    for (int c = 0; c < w; ++c) {
      if (out[r * out_stride + c] != ref_out[r * ref_out_stride + c])
        return 0;
    }
  }
  return 1;
}

void setup_convolve() {
#if HAVE_SSSE3 && CONFIG_RUNTIME_CPU_DETECT
  av1_convolve_horiz = av1_convolve_horiz_c;
  av1_convolve_vert = av1_convolve_vert_c;
#endif
}

typedef void (*Convolve8Func)(const uint8_t *src, ptrdiff_t src_stride,
                             uint8_t *dst, ptrdiff_t dst_stride,
                             const int16_t *filter_x, int filter_x_stride,
                             const int16_t *filter_y, int filter_y_stride,
                             int w, int h);

typedef void (*ConvolveFunc)(const uint8_t *src, int src_stride, uint8_t *dst,
                             int dst_stride, int w, int h,
                             const InterpFilterParams filter_params,
                             const int subpel_y_q4, int y_step_q4, int avg);

struct ConvolveFunctions {
  ConvolveFunctions(Convolve8Func h8, Convolve8Func h8_avg, Convolve8Func v8, Convolve8Func v8_avg, ConvolveFunc hf, ConvolveFunc vf)
      : h8_(h8), h8_avg_(h8_avg), v8_(v8), v8_avg_(v8_avg), hf_(hf), vf_(vf) {}
  Convolve8Func h8_;
  Convolve8Func h8_avg_;
  Convolve8Func v8_;
  Convolve8Func v8_avg_;
  ConvolveFunc hf_;
  ConvolveFunc vf_;
};

class Av1ConvolveTest : public ::testing::TestWithParam<ConvolveFunctions> {
 public:
  virtual void SetUp() {
    // Force input_ to be unaligned, output to be 16 byte aligned.
    input_ = reinterpret_cast<uint8_t *>( aom_memalign(kDataAlignment, kInputBufferSize));
    output_ = reinterpret_cast<uint8_t *>( aom_memalign(kDataAlignment, kOutputBufferSize));
    ref_output_ = reinterpret_cast<uint8_t *>( aom_memalign(kDataAlignment, kOutputBufferSize));
    ConvolveFunctions cfs = GetParam();
    aom_convolve8_horiz = cfs.h8_;
    aom_convolve8_avg_horiz = cfs.h8_avg_;
    aom_convolve8_vert = cfs.v8_;
    aom_convolve8_avg_vert = cfs.v8_avg_;
    av1_convolve_horiz = cfs.hf_;
    av1_convolve_vert = cfs.vf_;
  }
  virtual void TearDown() {
    aom_free(input_);
    aom_free(output_);
    aom_free(ref_output_);
  }
  virtual uint8_t* input(int w, int h, int& stride) {
    stride = w + MAX_FILTER_TAP - 1;
    int offset = MAX_FILTER_TAP/2 - 1;
    for (int r = 0; r < h + MAX_FILTER_TAP - 1; ++r) {
      for (int c = 0; c < w + MAX_FILTER_TAP - 1; ++c) {
        input_[r * stride + c] = r + c;
      }
    }
    return input_ + offset * stride + offset;
  }
  virtual uint8_t* output(int w, int h, int& stride) {
    stride = w;
    return output_;
  }
  virtual uint8_t* ref_output(int w, int h, int& stride) {
    stride = w;
    return ref_output_;
  }

 protected:
  static const int kDataAlignment = 16;
  static const int kOuterBlockSize = MAX_SB_SIZE;
  static const int kInputStride = kOuterBlockSize;
  static const int kOutputStride = kOuterBlockSize;
  static const int kInputBufferSize = kOuterBlockSize * kOuterBlockSize;
  static const int kOutputBufferSize = kOuterBlockSize * kOuterBlockSize;
  uint8_t *input_;
  uint8_t *output_;
  uint8_t *ref_output_;
};

TEST_P(Av1ConvolveTest, av1_convolve) {
  // define parameter
  int w = 16;
  int h = 16;
  int subpel_x_q4 = 8;
  int subpel_y_q4 = 8;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int ref_idx = 0;
  InterpFilter interp_filter[4] = {MULTITAP_SHARP, EIGHTTAP_REGULAR, EIGHTTAP_REGULAR, EIGHTTAP_REGULAR,};

  int in_stride, out_stride, ref_out_stride;
  uint8_t* in = input(w, h, in_stride);
  uint8_t* out = output(w, h, out_stride);
  uint8_t* ref_out = ref_output(w, h, ref_out_stride);

  InterpFilterParams param_vert = av1_get_interp_filter_params(interp_filter[0]);
  InterpFilterParams param_horiz = av1_get_interp_filter_params(interp_filter[1]);
  const int16_t *filter_vert = av1_get_interp_filter_subpel_kernel(param_vert, subpel_y_q4);
  const int16_t *filter_horiz = av1_get_interp_filter_subpel_kernel(param_horiz, subpel_x_q4);

  filter_block2d_c(in, in_stride, ref_out, ref_out_stride, w, h, filter_horiz, param_horiz.taps, filter_vert, param_vert.taps);

  av1_convolve(in, in_stride, out, out_stride, w, h, interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4, y_step_q4, ref_idx);

  EXPECT_EQ(match(out, out_stride, ref_out, ref_out_stride, w, h), 1);
};

ConvolveFunctions convolve_functions_c(aom_convolve8_horiz_c, aom_convolve8_avg_horiz_c, aom_convolve8_vert_c, aom_convolve8_avg_vert_c, av1_convolve_horiz_c, av1_convolve_vert_c);

ConvolveFunctions convolve_functions_ls[] = {convolve_functions_c};

INSTANTIATE_TEST_CASE_P(C, Av1ConvolveTest, ::testing::ValuesIn(convolve_functions_ls));

TEST(AV1ConvolveTest, av1_convolve8) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4] = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR,
                                    EIGHTTAP_REGULAR, EIGHTTAP_REGULAR };
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter[0]);
#else
  InterpFilter interp_filter = EIGHTTAP_REGULAR;
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter);
#endif
  int filter_size = filter_params.taps;
  int filter_center = filter_size / 2 - 1;
  uint8_t src[12 * 12];
  int src_stride = filter_size;
  uint8_t dst[1] = { 0 };
  uint8_t dst1[1] = { 0 };
  int dst_stride = 1;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int subpel_x_q4 = 3;
  int subpel_y_q4 = 2;
  int avg = 0;

  int w = 1;
  int h = 1;

  setup_convolve();

  for (int i = 0; i < filter_size * filter_size; i++) {
    src[i] = rnd.Rand16() % (1 << 8);
  }

  av1_convolve(src + src_stride * filter_center + filter_center, src_stride,
               dst, dst_stride, w, h, interp_filter, subpel_x_q4, x_step_q4,
               subpel_y_q4, y_step_q4, avg);

  const int16_t *x_filter =
      av1_get_interp_filter_subpel_kernel(filter_params, subpel_x_q4);
  const int16_t *y_filter =
      av1_get_interp_filter_subpel_kernel(filter_params, subpel_y_q4);

  aom_convolve8_c(src + src_stride * filter_center + filter_center, src_stride,
                  dst1, dst_stride, x_filter, 16, y_filter, 16, w, h);
  EXPECT_EQ(dst[0], dst1[0]);
}
TEST(AV1ConvolveTest, av1_convolve) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4] = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR,
                                    EIGHTTAP_REGULAR, EIGHTTAP_REGULAR };
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter[0]);
#else
  InterpFilter interp_filter = EIGHTTAP_REGULAR;
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter);
#endif
  int filter_size = filter_params.taps;
  int filter_center = filter_size / 2 - 1;
  uint8_t src[12 * 12];
  int src_stride = filter_size;
  uint8_t dst[1] = { 0 };
  int dst_stride = 1;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int avg = 0;
  int w = 1;
  int h = 1;

  int subpel_x_q4;
  int subpel_y_q4;

  ASSERT_LE(filter_size, 12);
  setup_convolve();

  for (int i = 0; i < static_cast<int>(sizeof(src) / sizeof(src[0])); i++) {
    src[i] = rnd.Rand16() % (1 << 8);
  }

  for (subpel_x_q4 = 0; subpel_x_q4 < SUBPEL_SHIFTS; subpel_x_q4++) {
    for (subpel_y_q4 = 0; subpel_y_q4 < SUBPEL_SHIFTS; subpel_y_q4++) {
      av1_convolve(src + src_stride * filter_center + filter_center, src_stride,
                   dst, dst_stride, w, h, interp_filter, subpel_x_q4, x_step_q4,
                   subpel_y_q4, y_step_q4, avg);

      const int16_t *x_filter =
          av1_get_interp_filter_subpel_kernel(filter_params, subpel_x_q4);
      const int16_t *y_filter =
          av1_get_interp_filter_subpel_kernel(filter_params, subpel_y_q4);

      int temp[12];
      int dst_ref = 0;
      for (int r = 0; r < filter_size; r++) {
        temp[r] = 0;
        for (int c = 0; c < filter_size; c++) {
          temp[r] += x_filter[c] * src[r * filter_size + c];
        }
        temp[r] = clip_pixel(ROUND_POWER_OF_TWO(temp[r], FILTER_BITS));
        dst_ref += temp[r] * y_filter[r];
      }
      dst_ref = clip_pixel(ROUND_POWER_OF_TWO(dst_ref, FILTER_BITS));
      EXPECT_EQ(dst[0], dst_ref);
    }
  }
}

#if CONFIG_EXT_INTERP && CONFIG_DUAL_FILTER
TEST(AV1ConvolveTest, av1_convolve_vert_first) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  InterpFilter interp_filter[4] = { EIGHTTAP_REGULAR, MULTITAP_SHARP,
                                    EIGHTTAP_REGULAR, MULTITAP_SHARP };
  InterpFilterParams filter_params_x =
      av1_get_interp_filter_params(interp_filter[1]);
  InterpFilterParams filter_params_y =
      av1_get_interp_filter_params(interp_filter[0]);
  int filter_size_x = filter_params_x.taps;
  int filter_size_y = filter_params_y.taps;
  int filter_center_x = filter_size_x / 2 - 1;
  int filter_center_y = filter_size_y / 2 - 1;
  uint8_t src[12 * 12];
  int src_stride = filter_size_x;
  uint8_t dst[1] = { 0 };
  int dst_stride = 1;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int avg = 0;
  int w = 1;
  int h = 1;

  int subpel_x_q4;
  int subpel_y_q4;

  ASSERT_LE(filter_size_x, 12);
  ASSERT_LE(filter_size_y, 12);
  setup_convolve();

  for (int i = 0; i < static_cast<int>(sizeof(src) / sizeof(src[0])); i++) {
    src[i] = rnd.Rand16() % (1 << 8);
  }

  for (subpel_x_q4 = 1; subpel_x_q4 < SUBPEL_SHIFTS; subpel_x_q4++) {
    for (subpel_y_q4 = 1; subpel_y_q4 < SUBPEL_SHIFTS; subpel_y_q4++) {
      av1_convolve(src + src_stride * filter_center_y + filter_center_x,
                   src_stride, dst, dst_stride, w, h, interp_filter,
                   subpel_x_q4, x_step_q4, subpel_y_q4, y_step_q4, avg);

      const int16_t *x_filter =
          av1_get_interp_filter_subpel_kernel(filter_params_x, subpel_x_q4);
      const int16_t *y_filter =
          av1_get_interp_filter_subpel_kernel(filter_params_y, subpel_y_q4);

      int temp[12];
      int dst_ref = 0;
      for (int c = 0; c < filter_size_x; c++) {
        temp[c] = 0;
        for (int r = 0; r < filter_size_y; r++) {
          temp[c] += y_filter[r] * src[r * filter_size_x + c];
        }
        temp[c] = clip_pixel(ROUND_POWER_OF_TWO(temp[c], FILTER_BITS));
        dst_ref += temp[c] * x_filter[c];
      }
      dst_ref = clip_pixel(ROUND_POWER_OF_TWO(dst_ref, FILTER_BITS));
      EXPECT_EQ(dst[0], dst_ref);
    }
  }
}
#endif

TEST(AV1ConvolveTest, av1_convolve_avg) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4] = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR,
                                    EIGHTTAP_REGULAR, EIGHTTAP_REGULAR };
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter[0]);
#else
  InterpFilter interp_filter = EIGHTTAP_REGULAR;
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter);
#endif
  int filter_size = filter_params.taps;
  int filter_center = filter_size / 2 - 1;
  uint8_t src0[12 * 12];
  uint8_t src1[12 * 12];
  int src_stride = filter_size;
  uint8_t dst0[1] = { 0 };
  uint8_t dst1[1] = { 0 };
  uint8_t dst[1] = { 0 };
  int dst_stride = 1;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int avg = 0;

  int w = 1;
  int h = 1;

  int subpel_x_q4;
  int subpel_y_q4;

  setup_convolve();

  for (int i = 0; i < filter_size * filter_size; i++) {
    src0[i] = rnd.Rand16() % (1 << 8);
    src1[i] = rnd.Rand16() % (1 << 8);
  }

  int offset = filter_size * filter_center + filter_center;

  for (subpel_x_q4 = 0; subpel_x_q4 < SUBPEL_SHIFTS; subpel_x_q4++) {
    for (subpel_y_q4 = 0; subpel_y_q4 < SUBPEL_SHIFTS; subpel_y_q4++) {
      avg = 0;
      av1_convolve(src0 + offset, src_stride, dst0, dst_stride, w, h,
                   interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                   y_step_q4, avg);
      avg = 0;
      av1_convolve(src1 + offset, src_stride, dst1, dst_stride, w, h,
                   interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                   y_step_q4, avg);

      avg = 0;
      av1_convolve(src0 + offset, src_stride, dst, dst_stride, w, h,
                   interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                   y_step_q4, avg);
      avg = 1;
      av1_convolve(src1 + offset, src_stride, dst, dst_stride, w, h,
                   interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                   y_step_q4, avg);

      EXPECT_EQ(dst[0], ROUND_POWER_OF_TWO(dst0[0] + dst1[0], 1));
    }
  }
}

#if CONFIG_AOM_HIGHBITDEPTH
TEST(AV1ConvolveTest, av1_highbd_convolve) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4] = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR,
                                    EIGHTTAP_REGULAR, EIGHTTAP_REGULAR };
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter[0]);
#else
  InterpFilter interp_filter = EIGHTTAP_REGULAR;
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter);
#endif
  int filter_size = filter_params.taps;
  int filter_center = filter_size / 2 - 1;
  uint16_t src[12 * 12];
  int src_stride = filter_size;
  uint16_t dst[1] = { 0 };
  int dst_stride = 1;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int avg = 0;
  int bd = 10;
  int w = 1;
  int h = 1;

  int subpel_x_q4;
  int subpel_y_q4;

  for (int i = 0; i < filter_size * filter_size; i++) {
    src[i] = rnd.Rand16() % (1 << bd);
  }

  for (subpel_x_q4 = 0; subpel_x_q4 < SUBPEL_SHIFTS; subpel_x_q4++) {
    for (subpel_y_q4 = 0; subpel_y_q4 < SUBPEL_SHIFTS; subpel_y_q4++) {
      av1_highbd_convolve(
          CONVERT_TO_BYTEPTR(src + src_stride * filter_center + filter_center),
          src_stride, CONVERT_TO_BYTEPTR(dst), dst_stride, w, h, interp_filter,
          subpel_x_q4, x_step_q4, subpel_y_q4, y_step_q4, avg, bd);

      const int16_t *x_filter =
          av1_get_interp_filter_subpel_kernel(filter_params, subpel_x_q4);
      const int16_t *y_filter =
          av1_get_interp_filter_subpel_kernel(filter_params, subpel_y_q4);

      int temp[12];
      int dst_ref = 0;
      for (int r = 0; r < filter_size; r++) {
        temp[r] = 0;
        for (int c = 0; c < filter_size; c++) {
          temp[r] += x_filter[c] * src[r * filter_size + c];
        }
        temp[r] =
            clip_pixel_highbd(ROUND_POWER_OF_TWO(temp[r], FILTER_BITS), bd);
        dst_ref += temp[r] * y_filter[r];
      }
      dst_ref = clip_pixel_highbd(ROUND_POWER_OF_TWO(dst_ref, FILTER_BITS), bd);
      EXPECT_EQ(dst[0], dst_ref);
    }
  }
}

TEST(AV1ConvolveTest, av1_highbd_convolve_avg) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
#if CONFIG_DUAL_FILTER
  InterpFilter interp_filter[4] = { EIGHTTAP_REGULAR, EIGHTTAP_REGULAR,
                                    EIGHTTAP_REGULAR, EIGHTTAP_REGULAR };
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter[0]);
#else
  InterpFilter interp_filter = EIGHTTAP_REGULAR;
  InterpFilterParams filter_params =
      av1_get_interp_filter_params(interp_filter);
#endif
  int filter_size = filter_params.taps;
  int filter_center = filter_size / 2 - 1;
  uint16_t src0[12 * 12];
  uint16_t src1[12 * 12];
  int src_stride = filter_size;
  uint16_t dst0[1] = { 0 };
  uint16_t dst1[1] = { 0 };
  uint16_t dst[1] = { 0 };
  int dst_stride = 1;
  int x_step_q4 = 16;
  int y_step_q4 = 16;
  int avg = 0;
  int bd = 10;

  int w = 1;
  int h = 1;

  int subpel_x_q4;
  int subpel_y_q4;

  for (int i = 0; i < filter_size * filter_size; i++) {
    src0[i] = rnd.Rand16() % (1 << bd);
    src1[i] = rnd.Rand16() % (1 << bd);
  }

  for (subpel_x_q4 = 0; subpel_x_q4 < SUBPEL_SHIFTS; subpel_x_q4++) {
    for (subpel_y_q4 = 0; subpel_y_q4 < SUBPEL_SHIFTS; subpel_y_q4++) {
      int offset = filter_size * filter_center + filter_center;

      avg = 0;
      av1_highbd_convolve(CONVERT_TO_BYTEPTR(src0 + offset), src_stride,
                          CONVERT_TO_BYTEPTR(dst0), dst_stride, w, h,
                          interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                          y_step_q4, avg, bd);
      avg = 0;
      av1_highbd_convolve(CONVERT_TO_BYTEPTR(src1 + offset), src_stride,
                          CONVERT_TO_BYTEPTR(dst1), dst_stride, w, h,
                          interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                          y_step_q4, avg, bd);

      avg = 0;
      av1_highbd_convolve(CONVERT_TO_BYTEPTR(src0 + offset), src_stride,
                          CONVERT_TO_BYTEPTR(dst), dst_stride, w, h,
                          interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                          y_step_q4, avg, bd);
      avg = 1;
      av1_highbd_convolve(CONVERT_TO_BYTEPTR(src1 + offset), src_stride,
                          CONVERT_TO_BYTEPTR(dst), dst_stride, w, h,
                          interp_filter, subpel_x_q4, x_step_q4, subpel_y_q4,
                          y_step_q4, avg, bd);

      EXPECT_EQ(dst[0], ROUND_POWER_OF_TWO(dst0[0] + dst1[0], 1));
    }
  }
}
#endif  // CONFIG_AOM_HIGHBITDEPTH

#define CONVOLVE_SPEED_TEST 0
#if CONVOLVE_SPEED_TEST
#define highbd_convolve_speed(func, block_size, frame_size)                  \
  TEST(AV1ConvolveTest, func##_speed_##block_size##_##frame_size) {          \
    ACMRandom rnd(ACMRandom::DeterministicSeed());                           \
    InterpFilter interp_filter = EIGHTTAP;                                   \
    InterpFilterParams filter_params =                                       \
        av1_get_interp_filter_params(interp_filter);                         \
    int filter_size = filter_params.tap;                                     \
    int filter_center = filter_size / 2 - 1;                                 \
    DECLARE_ALIGNED(16, uint16_t,                                            \
                    src[(frame_size + 7) * (frame_size + 7)]) = { 0 };       \
    int src_stride = frame_size + 7;                                         \
    DECLARE_ALIGNED(16, uint16_t, dst[frame_size * frame_size]) = { 0 };     \
    int dst_stride = frame_size;                                             \
    int x_step_q4 = 16;                                                      \
    int y_step_q4 = 16;                                                      \
    int subpel_x_q4 = 8;                                                     \
    int subpel_y_q4 = 6;                                                     \
    int bd = 10;                                                             \
                                                                             \
    int w = block_size;                                                      \
    int h = block_size;                                                      \
                                                                             \
    const int16_t *filter_x =                                                \
        av1_get_interp_filter_kernel(filter_params, subpel_x_q4);            \
    const int16_t *filter_y =                                                \
        av1_get_interp_filter_kernel(filter_params, subpel_y_q4);            \
                                                                             \
    for (int i = 0; i < src_stride * src_stride; i++) {                      \
      src[i] = rnd.Rand16() % (1 << bd);                                     \
    }                                                                        \
                                                                             \
    int offset = filter_center * src_stride + filter_center;                 \
    int row_offset = 0;                                                      \
    int col_offset = 0;                                                      \
    for (int i = 0; i < 100000; i++) {                                       \
      int src_total_offset = offset + col_offset * src_stride + row_offset;  \
      int dst_total_offset = col_offset * dst_stride + row_offset;           \
      func(CONVERT_TO_BYTEPTR(src + src_total_offset), src_stride,           \
           CONVERT_TO_BYTEPTR(dst + dst_total_offset), dst_stride, filter_x, \
           x_step_q4, filter_y, y_step_q4, w, h, bd);                        \
      if (offset + w + w < frame_size) {                                     \
        row_offset += w;                                                     \
      } else {                                                               \
        row_offset = 0;                                                      \
        col_offset += h;                                                     \
      }                                                                      \
      if (col_offset + h >= frame_size) {                                    \
        col_offset = 0;                                                      \
      }                                                                      \
    }                                                                        \
  }

#define lowbd_convolve_speed(func, block_size, frame_size)                  \
  TEST(AV1ConvolveTest, func##_speed_l_##block_size##_##frame_size) {       \
    ACMRandom rnd(ACMRandom::DeterministicSeed());                          \
    InterpFilter interp_filter = EIGHTTAP;                                  \
    InterpFilterParams filter_params =                                      \
        av1_get_interp_filter_params(interp_filter);                        \
    int filter_size = filter_params.tap;                                    \
    int filter_center = filter_size / 2 - 1;                                \
    DECLARE_ALIGNED(16, uint8_t, src[(frame_size + 7) * (frame_size + 7)]); \
    int src_stride = frame_size + 7;                                        \
    DECLARE_ALIGNED(16, uint8_t, dst[frame_size * frame_size]);             \
    int dst_stride = frame_size;                                            \
    int x_step_q4 = 16;                                                     \
    int y_step_q4 = 16;                                                     \
    int subpel_x_q4 = 8;                                                    \
    int subpel_y_q4 = 6;                                                    \
    int bd = 8;                                                             \
                                                                            \
    int w = block_size;                                                     \
    int h = block_size;                                                     \
                                                                            \
    const int16_t *filter_x =                                               \
        av1_get_interp_filter_kernel(filter_params, subpel_x_q4);           \
    const int16_t *filter_y =                                               \
        av1_get_interp_filter_kernel(filter_params, subpel_y_q4);           \
                                                                            \
    for (int i = 0; i < src_stride * src_stride; i++) {                     \
      src[i] = rnd.Rand16() % (1 << bd);                                    \
    }                                                                       \
                                                                            \
    int offset = filter_center * src_stride + filter_center;                \
    int row_offset = 0;                                                     \
    int col_offset = 0;                                                     \
    for (int i = 0; i < 100000; i++) {                                      \
      func(src + offset, src_stride, dst, dst_stride, filter_x, x_step_q4,  \
           filter_y, y_step_q4, w, h);                                      \
      if (offset + w + w < frame_size) {                                    \
        row_offset += w;                                                    \
      } else {                                                              \
        row_offset = 0;                                                     \
        col_offset += h;                                                    \
      }                                                                     \
      if (col_offset + h >= frame_size) {                                   \
        col_offset = 0;                                                     \
      }                                                                     \
    }                                                                       \
  }

// This experiment shows that when frame size is 64x64
// aom_highbd_convolve8_sse2 and aom_convolve8_sse2's speed are similar.
// However when frame size becomes 1024x1024
// aom_highbd_convolve8_sse2 is around 50% slower than aom_convolve8_sse2
// we think the bottleneck is from memory IO
highbd_convolve_speed(aom_highbd_convolve8_sse2, 8, 64);
highbd_convolve_speed(aom_highbd_convolve8_sse2, 16, 64);
highbd_convolve_speed(aom_highbd_convolve8_sse2, 32, 64);
highbd_convolve_speed(aom_highbd_convolve8_sse2, 64, 64);

lowbd_convolve_speed(aom_convolve8_sse2, 8, 64);
lowbd_convolve_speed(aom_convolve8_sse2, 16, 64);
lowbd_convolve_speed(aom_convolve8_sse2, 32, 64);
lowbd_convolve_speed(aom_convolve8_sse2, 64, 64);

highbd_convolve_speed(aom_highbd_convolve8_sse2, 8, 1024);
highbd_convolve_speed(aom_highbd_convolve8_sse2, 16, 1024);
highbd_convolve_speed(aom_highbd_convolve8_sse2, 32, 1024);
highbd_convolve_speed(aom_highbd_convolve8_sse2, 64, 1024);

lowbd_convolve_speed(aom_convolve8_sse2, 8, 1024);
lowbd_convolve_speed(aom_convolve8_sse2, 16, 1024);
lowbd_convolve_speed(aom_convolve8_sse2, 32, 1024);
lowbd_convolve_speed(aom_convolve8_sse2, 64, 1024);
#endif  // CONVOLVE_SPEED_TEST
}  // namespace
