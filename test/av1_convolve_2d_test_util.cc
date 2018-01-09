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

#include "test/av1_convolve_2d_test_util.h"

#include "aom_ports/aom_timer.h"
#include "av1/common/common_data.h"
#include "av1/common/convolve.h"

using std::tr1::tuple;
using std::tr1::make_tuple;

namespace libaom_test {

const int kMaxSize = 128 + 32;  // padding

namespace AV1Convolve2D {

::testing::internal::ParamGenerator<Convolve2DParam> BuildParams(
    convolve_2d_func filter, int has_subx, int has_suby, int is_compound) {
  return ::testing::Values(make_tuple(filter, has_subx, has_suby, is_compound));
}

AV1Convolve2DTest::~AV1Convolve2DTest() {}
void AV1Convolve2DTest::SetUp() { rnd_.Reset(ACMRandom::DeterministicSeed()); }

void AV1Convolve2DTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1Convolve2DTest::RunCheckOutput(convolve_2d_func test_impl) {
  const int w = kMaxSize, h = kMaxSize;
  const int num_iters = 1;
  const int has_subx = GET_PARAM(1);
  const int has_suby = GET_PARAM(2);
  const int is_compound = GET_PARAM(3);
  int i, j, k;
  int hfilter, vfilter, subx, suby;
  uint8_t input[kMaxSize * kMaxSize];
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output2[MAX_SB_SQUARE]);

  (void)is_compound;

  for (i = 0; i < num_iters; ++i) {
    for (j = 0; j < h; ++j)
      for (k = 0; k < w; ++k) input[j * w + k] = rnd_.Rand8();

    for (int block_idx = BLOCK_4X4; block_idx < BLOCK_SIZES_ALL; ++block_idx) {
      const int out_w = block_size_wide[block_idx];
      const int out_h = block_size_high[block_idx];
      for (hfilter = EIGHTTAP_REGULAR; hfilter < INTERP_FILTERS_ALL;
           ++hfilter) {
        for (vfilter = EIGHTTAP_REGULAR; vfilter < INTERP_FILTERS_ALL;
             ++vfilter) {
          InterpFilterParams filter_params_x =
              av1_get_interp_filter_params((InterpFilter)hfilter);
          InterpFilterParams filter_params_y =
              av1_get_interp_filter_params((InterpFilter)vfilter);
          for (int do_average = 0; do_average <= 1; ++do_average) {
            ConvolveParams conv_params1 = get_conv_params_no_round(
                0, do_average, 0, output, MAX_SB_SIZE, 1);
            ConvolveParams conv_params2 = get_conv_params_no_round(
                0, do_average, 0, output2, MAX_SB_SIZE, 1);

            const int subx_range = has_subx ? 16 : 1;
            const int suby_range = has_suby ? 16 : 1;
            for (subx = 0; subx < subx_range; ++subx) {
              for (suby = 0; suby < suby_range; ++suby) {
                memset(output, 0, sizeof(output));
                memset(output2, 0, sizeof(output2));
                // Choose random locations within the source block
                const int offset_r = 3 + rnd_.PseudoUniform(h - out_h - 7);
                const int offset_c = 3 + rnd_.PseudoUniform(w - out_w - 7);
                av1_convolve_2d_c(input + offset_r * w + offset_c, w, NULL, 0,
                                  out_w, out_h, &filter_params_x,
                                  &filter_params_y, subx, suby, &conv_params1);
                test_impl(input + offset_r * w + offset_c, w, NULL, 0, out_w,
                          out_h, &filter_params_x, &filter_params_y, subx, suby,
                          &conv_params2);

                for (j = 0; j < out_h; ++j) {
                  for (k = 0; k < out_w; ++k) {
                    int idx = j * MAX_SB_SIZE + k;
                    ASSERT_EQ(output[idx], output2[idx])
                        << out_w << "x" << out_h << " Pixel mismatch at index "
                        << idx << " = (" << j << ", " << k
                        << "), sub pixel offset = (" << suby << ", " << subx
                        << ")";
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void AV1Convolve2DTest::RunSpeedTest(convolve_2d_func test_impl) {
  const int w = kMaxSize, h = kMaxSize;
  int i, j;
  const int has_subx = GET_PARAM(1);
  const int has_suby = GET_PARAM(2);
  const int is_compound = GET_PARAM(3);
  (void)is_compound;

  uint8_t input[kMaxSize * kMaxSize];
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output[MAX_SB_SQUARE]);

  for (i = 0; i < h; ++i)
    for (j = 0; j < w; ++j) input[i * w + j] = rnd_.Rand8();

  int hfilter = EIGHTTAP_REGULAR, vfilter = EIGHTTAP_REGULAR;
  int subx = 0, suby = 0;

  InterpFilterParams filter_params_x =
      av1_get_interp_filter_params((InterpFilter)hfilter);
  InterpFilterParams filter_params_y =
      av1_get_interp_filter_params((InterpFilter)vfilter);
  const int do_average = 0;
  ConvolveParams conv_params2 =
      get_conv_params_no_round(0, do_average, 0, output, MAX_SB_SIZE, 1);
  int x;

  for (int block_idx = BLOCK_4X4; block_idx < BLOCK_SIZES_ALL; ++block_idx) {
    const int out_w = block_size_wide[block_idx];
    const int out_h = block_size_high[block_idx];
    const int num_loops = 1000000000 / (out_w + out_h);
    aom_usec_timer timer;
    aom_usec_timer_start(&timer);

    for (x = 0; x < num_loops; ++x)
      test_impl(input, w, NULL, 0, out_w, out_h, &filter_params_x,
                &filter_params_y, subx, suby, &conv_params2);

    aom_usec_timer_mark(&timer);
    const int elapsed_time = static_cast<int>(aom_usec_timer_elapsed(&timer));
    printf("%d,%d convolve %3dx%-3d: %7.2f ns\n", has_subx, has_suby, out_w,
           out_h, 1000.0 * elapsed_time / num_loops);
  }
}

#if CONFIG_JNT_COMP
AV1JntConvolve2DTest::~AV1JntConvolve2DTest() {}
void AV1JntConvolve2DTest::SetUp() {
  rnd_.Reset(ACMRandom::DeterministicSeed());
}

void AV1JntConvolve2DTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1JntConvolve2DTest::RunCheckOutput(convolve_2d_func test_impl) {
  const int w = kMaxSize, h = kMaxSize;
  const int num_iters = 1;
  const int has_subx = GET_PARAM(1);
  const int has_suby = GET_PARAM(2);
  const int is_compound = GET_PARAM(3);
  int i, j, k, l, m;
  int hfilter, vfilter, subx, suby;
  uint8_t input[kMaxSize * kMaxSize];
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output2[MAX_SB_SQUARE]);

  (void)is_compound;

  for (i = 0; i < num_iters; ++i) {
    for (j = 0; j < h; ++j)
      for (k = 0; k < w; ++k) input[j * w + k] = rnd_.Rand8();

    for (int block_idx = BLOCK_4X4; block_idx < BLOCK_SIZES_ALL; ++block_idx) {
      const int out_w = block_size_wide[block_idx];
      const int out_h = block_size_high[block_idx];
      for (hfilter = EIGHTTAP_REGULAR; hfilter < INTERP_FILTERS_ALL;
           ++hfilter) {
        for (vfilter = EIGHTTAP_REGULAR; vfilter < INTERP_FILTERS_ALL;
             ++vfilter) {
          InterpFilterParams filter_params_x =
              av1_get_interp_filter_params((InterpFilter)hfilter);
          InterpFilterParams filter_params_y =
              av1_get_interp_filter_params((InterpFilter)vfilter);
          for (int do_average = 0; do_average <= 1; ++do_average) {
            ConvolveParams conv_params1 = get_conv_params_no_round(
                0, do_average, 0, output, MAX_SB_SIZE, 1);
            ConvolveParams conv_params2 = get_conv_params_no_round(
                0, do_average, 0, output2, MAX_SB_SIZE, 1);

            // Test special case where jnt_comp_avg is not used
            conv_params1.use_jnt_comp_avg = 0;
            conv_params2.use_jnt_comp_avg = 0;

            const int subx_range = has_subx ? 16 : 1;
            const int suby_range = has_suby ? 16 : 1;
            for (subx = 0; subx < subx_range; ++subx) {
              for (suby = 0; suby < suby_range; ++suby) {
                memset(output, 0, sizeof(output));
                memset(output2, 0, sizeof(output2));
                // Choose random locations within the source block
                const int offset_r = 3 + rnd_.PseudoUniform(h - out_h - 7);
                const int offset_c = 3 + rnd_.PseudoUniform(w - out_w - 7);
                av1_jnt_convolve_2d_c(input + offset_r * w + offset_c, w, NULL,
                                      0, out_w, out_h, &filter_params_x,
                                      &filter_params_y, subx, suby,
                                      &conv_params1);
                test_impl(input + offset_r * w + offset_c, w, NULL, 0, out_w,
                          out_h, &filter_params_x, &filter_params_y, subx, suby,
                          &conv_params2);

                for (j = 0; j < out_h; ++j) {
                  for (k = 0; k < out_w; ++k) {
                    int idx = j * MAX_SB_SIZE + k;
                    ASSERT_EQ(output[idx], output2[idx])
                        << "Mismatch at unit tests for av1_jnt_convolve_2d\n"
                        << out_w << "x" << out_h << " Pixel mismatch at index "
                        << idx << " = (" << j << ", " << k
                        << "), sub pixel offset = (" << suby << ", " << subx
                        << ")";
                  }
                }
              }
            }

            // Test different combination of fwd and bck offset weights
            for (l = 0; l < 2; ++l) {
              for (m = 0; m < 4; ++m) {
                conv_params1.use_jnt_comp_avg = 1;
                conv_params2.use_jnt_comp_avg = 1;
                conv_params1.fwd_offset = quant_dist_lookup_table[l][m][0];
                conv_params1.bck_offset = quant_dist_lookup_table[l][m][1];
                conv_params2.fwd_offset = quant_dist_lookup_table[l][m][0];
                conv_params2.bck_offset = quant_dist_lookup_table[l][m][1];

                for (subx = 0; subx < subx_range; ++subx) {
                  for (suby = 0; suby < suby_range; ++suby) {
                    memset(output, 0, sizeof(output));
                    memset(output2, 0, sizeof(output2));
                    // Choose random locations within the source block
                    const int offset_r = 3 + rnd_.PseudoUniform(h - out_h - 7);
                    const int offset_c = 3 + rnd_.PseudoUniform(w - out_w - 7);
                    av1_jnt_convolve_2d_c(input + offset_r * w + offset_c, w,
                                          NULL, 0, out_w, out_h,
                                          &filter_params_x, &filter_params_y,
                                          subx, suby, &conv_params1);
                    test_impl(input + offset_r * w + offset_c, w, NULL, 0,
                              out_w, out_h, &filter_params_x, &filter_params_y,
                              subx, suby, &conv_params2);

                    for (j = 0; j < out_h; ++j) {
                      for (k = 0; k < out_w; ++k) {
                        int idx = j * MAX_SB_SIZE + k;
                        ASSERT_EQ(output[idx], output2[idx])
                            << "Mismatch at unit tests for "
                               "av1_jnt_convolve_2d\n"
                            << out_w << "x" << out_h
                            << " Pixel mismatch at index " << idx << " = (" << j
                            << ", " << k << "), sub pixel offset = (" << suby
                            << ", " << subx << ")";
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
#endif  // CONFIG_JNT_COMP
}  // namespace AV1Convolve2D

namespace AV1HighbdConvolve2D {

::testing::internal::ParamGenerator<HighbdConvolve2DParam> BuildParams(
    highbd_convolve_2d_func filter) {
  return ::testing::Values(make_tuple(8, filter), make_tuple(10, filter),
                           make_tuple(12, filter));
}

AV1HighbdConvolve2DTest::~AV1HighbdConvolve2DTest() {}
void AV1HighbdConvolve2DTest::SetUp() {
  rnd_.Reset(ACMRandom::DeterministicSeed());
}

void AV1HighbdConvolve2DTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1HighbdConvolve2DTest::RunCheckOutput(
    highbd_convolve_2d_func test_impl) {
  const int w = kMaxSize, h = kMaxSize;
  const int num_iters = 1;
  const int bd = GET_PARAM(0);
  int i, j, k;
  int hfilter, vfilter, subx, suby;
  uint16_t input[kMaxSize * kMaxSize];
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output2[MAX_SB_SQUARE]);

  for (i = 0; i < num_iters; ++i) {
    for (j = 0; j < h; ++j)
      for (k = 0; k < w; ++k)
        input[j * w + k] = rnd_.Rand16() & ((1 << bd) - 1);

    for (int block_idx = BLOCK_4X4; block_idx < BLOCK_SIZES_ALL; ++block_idx) {
      const int out_w = block_size_wide[block_idx];
      const int out_h = block_size_high[block_idx];
      for (hfilter = EIGHTTAP_REGULAR; hfilter < INTERP_FILTERS_ALL;
           ++hfilter) {
        for (vfilter = EIGHTTAP_REGULAR; vfilter < INTERP_FILTERS_ALL;
             ++vfilter) {
          InterpFilterParams filter_params_x =
              av1_get_interp_filter_params((InterpFilter)hfilter);
          InterpFilterParams filter_params_y =
              av1_get_interp_filter_params((InterpFilter)vfilter);
          ConvolveParams conv_params1 =
              get_conv_params_no_round(0, 0, 0, output, MAX_SB_SIZE, 1);
          ConvolveParams conv_params2 =
              get_conv_params_no_round(0, 0, 0, output2, MAX_SB_SIZE, 1);

          for (subx = 0; subx < 16; ++subx) {
            for (suby = 0; suby < 16; ++suby) {
              memset(output, 0, sizeof(output));
              memset(output2, 0, sizeof(output2));
              // Choose random locations within the source block
              const int offset_r = 3 + rnd_.PseudoUniform(h - out_h - 7);
              const int offset_c = 3 + rnd_.PseudoUniform(w - out_w - 7);
              av1_highbd_convolve_2d_c(input + offset_r * w + offset_c, w,
                                       output, MAX_SB_SIZE, out_w, out_h,
                                       &filter_params_x, &filter_params_y, subx,
                                       suby, &conv_params1, bd);
              test_impl(input + offset_r * w + offset_c, w, output2,
                        MAX_SB_SIZE, out_w, out_h, &filter_params_x,
                        &filter_params_y, subx, suby, &conv_params2, bd);

              for (j = 0; j < out_h; ++j) {
                for (k = 0; k < out_w; ++k) {
                  int idx = j * MAX_SB_SIZE + k;
                  ASSERT_EQ(output[idx], output2[idx])
                      << out_w << "x" << out_h << " Pixel mismatch at index "
                      << idx << " = (" << j << ", " << k
                      << "), sub pixel offset = (" << suby << ", " << subx
                      << ")";
                }
              }
            }
          }
        }
      }
    }
  }
}

#if CONFIG_JNT_COMP
AV1HighbdJntConvolve2DTest::~AV1HighbdJntConvolve2DTest() {}
void AV1HighbdJntConvolve2DTest::SetUp() {
  rnd_.Reset(ACMRandom::DeterministicSeed());
}

void AV1HighbdJntConvolve2DTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1HighbdJntConvolve2DTest::RunCheckOutput(
    highbd_convolve_2d_func test_impl) {
  const int w = kMaxSize, h = kMaxSize;
  const int num_iters = 1;
  const int bd = GET_PARAM(0);
  int i, j, k, l, m;
  int hfilter, vfilter, subx, suby;
  uint16_t input[kMaxSize * kMaxSize];
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, CONV_BUF_TYPE, output2[MAX_SB_SQUARE]);

  for (i = 0; i < num_iters; ++i) {
    for (j = 0; j < h; ++j)
      for (k = 0; k < w; ++k)
        input[j * w + k] = rnd_.Rand16() & ((1 << bd) - 1);

    for (int block_idx = BLOCK_4X4; block_idx < BLOCK_SIZES_ALL; ++block_idx) {
      const int out_w = block_size_wide[block_idx];
      const int out_h = block_size_high[block_idx];
      for (hfilter = EIGHTTAP_REGULAR; hfilter < INTERP_FILTERS_ALL;
           ++hfilter) {
        for (vfilter = EIGHTTAP_REGULAR; vfilter < INTERP_FILTERS_ALL;
             ++vfilter) {
          InterpFilterParams filter_params_x =
              av1_get_interp_filter_params((InterpFilter)hfilter);
          InterpFilterParams filter_params_y =
              av1_get_interp_filter_params((InterpFilter)vfilter);
          ConvolveParams conv_params1 =
              get_conv_params_no_round(0, 0, 0, output, MAX_SB_SIZE, 1);
          ConvolveParams conv_params2 =
              get_conv_params_no_round(0, 0, 0, output2, MAX_SB_SIZE, 1);

          // Test special case where jnt_comp_avg is not used
          conv_params1.use_jnt_comp_avg = 0;
          conv_params2.use_jnt_comp_avg = 0;

          for (subx = 0; subx < 16; ++subx) {
            for (suby = 0; suby < 16; ++suby) {
              memset(output, 0, sizeof(output));
              memset(output2, 0, sizeof(output2));
              // Choose random locations within the source block
              const int offset_r = 3 + rnd_.PseudoUniform(h - out_h - 7);
              const int offset_c = 3 + rnd_.PseudoUniform(w - out_w - 7);
              av1_highbd_jnt_convolve_2d_c(input + offset_r * w + offset_c, w,
                                           output, MAX_SB_SIZE, out_w, out_h,
                                           &filter_params_x, &filter_params_y,
                                           subx, suby, &conv_params1, bd);
              test_impl(input + offset_r * w + offset_c, w, output2,
                        MAX_SB_SIZE, out_w, out_h, &filter_params_x,
                        &filter_params_y, subx, suby, &conv_params2, bd);

              for (j = 0; j < out_h; ++j) {
                for (k = 0; k < out_w; ++k) {
                  int idx = j * MAX_SB_SIZE + k;
                  ASSERT_EQ(output[idx], output2[idx])
                      << out_w << "x" << out_h << " Pixel mismatch at index "
                      << idx << " = (" << j << ", " << k
                      << "), sub pixel offset = (" << suby << ", " << subx
                      << ")";
                }
              }
            }
          }

          // Test different combination of fwd and bck offset weights
          for (l = 0; l < 2; ++l) {
            for (m = 0; m < 4; ++m) {
              conv_params1.use_jnt_comp_avg = 1;
              conv_params2.use_jnt_comp_avg = 1;
              conv_params1.fwd_offset = quant_dist_lookup_table[l][m][0];
              conv_params1.bck_offset = quant_dist_lookup_table[l][m][1];
              conv_params2.fwd_offset = quant_dist_lookup_table[l][m][0];
              conv_params2.bck_offset = quant_dist_lookup_table[l][m][1];

              for (subx = 0; subx < 16; ++subx) {
                for (suby = 0; suby < 16; ++suby) {
                  memset(output, 0, sizeof(output));
                  memset(output2, 0, sizeof(output2));
                  // Choose random locations within the source block
                  const int offset_r = 3 + rnd_.PseudoUniform(h - out_h - 7);
                  const int offset_c = 3 + rnd_.PseudoUniform(w - out_w - 7);
                  av1_highbd_jnt_convolve_2d_c(
                      input + offset_r * w + offset_c, w, output, MAX_SB_SIZE,
                      out_w, out_h, &filter_params_x, &filter_params_y, subx,
                      suby, &conv_params1, bd);
                  test_impl(input + offset_r * w + offset_c, w, output2,
                            MAX_SB_SIZE, out_w, out_h, &filter_params_x,
                            &filter_params_y, subx, suby, &conv_params2, bd);

                  for (j = 0; j < out_h; ++j) {
                    for (k = 0; k < out_w; ++k) {
                      int idx = j * MAX_SB_SIZE + k;
                      ASSERT_EQ(output[idx], output2[idx])
                          << out_w << "x" << out_h
                          << " Pixel mismatch at index " << idx << " = (" << j
                          << ", " << k << "), sub pixel offset = (" << suby
                          << ", " << subx << ")";
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
#endif  // CONFIG_JNT_COMP
}  // namespace AV1HighbdConvolve2D
}  // namespace libaom_test
