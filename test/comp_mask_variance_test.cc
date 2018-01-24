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

#include <cstdlib>
#include <new>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/util.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "aom/aom_codec.h"
#include "aom/aom_integer.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/mem.h"
#include "aom_ports/aom_timer.h"
#include "av1/common/reconinter.h"

namespace AV1CompMaskVariance {
typedef void (*comp_mask_pred_func)(uint8_t *comp_pred, const uint8_t *pred,
                                    int width, int height, const uint8_t *ref,
                                    int ref_stride, const uint8_t *mask,
                                    int mask_stride, int invert_mask);
typedef void (*comp_mask_up_pred_func)(uint8_t *comp_pred, const uint8_t *pred,
                                       int width, int height, int subpel_x_q3,
                                       int subpel_y_q3, const uint8_t *ref,
                                       int ref_stride, const uint8_t *mask,
                                       int mask_stride, int invert_mask);

typedef std::tr1::tuple<comp_mask_pred_func, BLOCK_SIZE> CompMaskPredParam;
typedef std::tr1::tuple<comp_mask_up_pred_func, BLOCK_SIZE> CompMaskUpPredParam;

class AV1CompMaskVarianceTest
    : public ::testing::TestWithParam<CompMaskPredParam> {
 public:
  ~AV1CompMaskVarianceTest();
  void SetUp();

  void TearDown();

 protected:
  void RunCheckOutput(comp_mask_pred_func test_impl, BLOCK_SIZE bsize);
  void RunSpeedTest(comp_mask_pred_func test_impl, BLOCK_SIZE bsize);

  libaom_test::ACMRandom rnd_;
};

AV1CompMaskVarianceTest::~AV1CompMaskVarianceTest() { ; }

void AV1CompMaskVarianceTest::SetUp() {
  rnd_.Reset(libaom_test::ACMRandom::DeterministicSeed());
}

void AV1CompMaskVarianceTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1CompMaskVarianceTest::RunCheckOutput(comp_mask_pred_func test_impl,
                                             BLOCK_SIZE bsize) {
  DECLARE_ALIGNED(16, uint8_t, comp_pred1[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, comp_pred2[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, pred[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, ref[MAX_SB_SQUARE]);

  // init
  av1_init_wedge_masks();
  const int w = block_size_wide[bsize];
  const int h = block_size_high[bsize];
  memset(comp_pred1, 0, MAX_SB_SQUARE);
  memset(comp_pred2, 0, MAX_SB_SQUARE);
  for (int i = 0; i < MAX_SB_SQUARE; ++i) {
    pred[i] = rnd_.Rand8();
    ref[i] = rnd_.Rand8();
  }
  int wedge_types = (1 << get_wedge_bits_lookup(bsize));
  for (int wedge_index = 0; wedge_index < wedge_types; ++wedge_index) {
    const uint8_t *mask = av1_get_contiguous_soft_mask(wedge_index, 1, bsize);

    for (int inv = 0; inv < 2; ++inv) {
      aom_comp_mask_pred_c(comp_pred1, pred, w, h, ref, MAX_SB_SIZE, mask, w,
                           inv);
      test_impl(comp_pred2, pred, w, h, ref, MAX_SB_SIZE, mask, w, inv);
      // check result
      for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
          int idx = i * w + j;
          ASSERT_EQ(comp_pred1[idx], comp_pred2[idx])
              << w << "x" << h << " Pixel mismatch at index " << idx << " = ("
              << i << ", " << j << "), wedge " << wedge_index << " inv " << inv;
        }
      }
    }
  }
}

void AV1CompMaskVarianceTest::RunSpeedTest(comp_mask_pred_func test_impl,
                                           BLOCK_SIZE bsize) {
  DECLARE_ALIGNED(16, uint8_t, comp_pred[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, pred[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, ref[MAX_SB_SQUARE]);

  // init
  av1_init_wedge_masks();
  const int w = block_size_wide[bsize];
  const int h = block_size_high[bsize];
  memset(comp_pred, 0, MAX_SB_SQUARE);
  for (int i = 0; i < MAX_SB_SQUARE; ++i) {
    pred[i] = rnd_.Rand8();
    ref[i] = rnd_.Rand8();
  }

  int wedge_types = (1 << get_wedge_bits_lookup(bsize));
  int wedge_index = wedge_types / 2;
  const uint8_t *mask = av1_get_contiguous_soft_mask(wedge_index, 1, bsize);
  const int num_loops = 1000000000 / (w + h);

  comp_mask_pred_func funcs[2] = { aom_comp_mask_pred_c, test_impl };
  double elapsed_time[2] = { 0 };
  for (int i = 0; i < 2; ++i) {
    aom_usec_timer timer;
    aom_usec_timer_start(&timer);
    comp_mask_pred_func func = funcs[i];
    for (int j = 0; j < num_loops; ++j) {
      func(comp_pred, pred, w, h, ref, MAX_SB_SIZE, mask, w, 0);
    }
    aom_usec_timer_mark(&timer);
    double time = static_cast<double>(aom_usec_timer_elapsed(&timer));
    elapsed_time[i] = 1000.0 * time / num_loops;
  }
  printf("comp_mask_pred %3dx%-3d: %7.2f/%7.2f ns", w, h, elapsed_time[0],
         elapsed_time[1]);
  printf(" (%3.2f)\n", elapsed_time[0] / elapsed_time[1]);
}

TEST_P(AV1CompMaskVarianceTest, CheckOutput) {
  RunCheckOutput(GET_PARAM(0), GET_PARAM(1));
}

TEST_P(AV1CompMaskVarianceTest, DISABLED_Speed) {
  RunSpeedTest(GET_PARAM(0), GET_PARAM(1));
}

const CompMaskPredParam kArrayCompMaskPred_ssse3[] = {
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_8X8),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_8X16),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_16X8),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_16X16),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_16X32),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_32X16),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_32X32),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_8X32),
  testing::make_tuple(&aom_comp_mask_pred_ssse3, BLOCK_32X8),
};

INSTANTIATE_TEST_CASE_P(SSSE3, AV1CompMaskVarianceTest,
                        ::testing::ValuesIn(kArrayCompMaskPred_ssse3));

class AV1CompMaskUpVarianceTest
    : public ::testing::TestWithParam<CompMaskUpPredParam> {
 public:
  ~AV1CompMaskUpVarianceTest();
  void SetUp();

  void TearDown();

 protected:
  void RunCheckOutput(comp_mask_up_pred_func test_impl, BLOCK_SIZE bsize);
  void RunSpeedTest(comp_mask_up_pred_func test_impl, BLOCK_SIZE bsize);
  void RunSpeedTestSub(comp_mask_up_pred_func test_impl, BLOCK_SIZE bsize,
                       int havSub);

  libaom_test::ACMRandom rnd_;
};

AV1CompMaskUpVarianceTest::~AV1CompMaskUpVarianceTest() { ; }

void AV1CompMaskUpVarianceTest::SetUp() {
  rnd_.Reset(libaom_test::ACMRandom::DeterministicSeed());
}

void AV1CompMaskUpVarianceTest::TearDown() { libaom_test::ClearSystemState(); }

void AV1CompMaskUpVarianceTest::RunCheckOutput(comp_mask_up_pred_func test_impl,
                                               BLOCK_SIZE bsize) {
  DECLARE_ALIGNED(16, uint8_t, comp_pred1[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, comp_pred2[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, pred[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, ref[MAX_SB_SQUARE]);

  // init
  av1_init_wedge_masks();
  memset(comp_pred1, 0, MAX_SB_SQUARE);
  memset(comp_pred2, 0, MAX_SB_SQUARE);
  for (int i = 0; i < MAX_SB_SQUARE; ++i) {
    pred[i] = rnd_.Rand8();
    ref[i] = rnd_.Rand8();
  }
  const int w = block_size_wide[bsize];
  const int h = block_size_high[bsize];
  int wedge_types = (1 << get_wedge_bits_lookup(bsize));

  // loop through subx and suby
  for (int sub = 0; sub < 16 * 16; ++sub) {
    int subx = sub & 0xf;
    int suby = (sub >> 4);

    for (int wedge_index = 0; wedge_index < wedge_types; ++wedge_index) {
      const uint8_t *mask = av1_get_contiguous_soft_mask(wedge_index, 1, bsize);

      for (int inv = 0; inv < 2; ++inv) {
        aom_comp_mask_upsampled_pred_c(comp_pred1, pred, w, h, subx, suby, ref,
                                       MAX_SB_SIZE, mask, w, inv);
        test_impl(comp_pred2, pred, w, h, subx, suby, ref, MAX_SB_SIZE, mask, w,
                  inv);
        // check result
        for (int i = 0; i < h; ++i) {
          for (int j = 0; j < w; ++j) {
            int idx = i * w + j;
            ASSERT_EQ(comp_pred1[idx], comp_pred2[idx])
                << w << "x" << h << " Pixel mismatch at index " << idx << " = ("
                << i << ", " << j << "), wedge " << wedge_index << " inv "
                << inv << "sub (" << subx << "," << suby << ")";
          }
        }
      }
    }
  }
}

void AV1CompMaskUpVarianceTest::RunSpeedTestSub(
    comp_mask_up_pred_func test_impl, BLOCK_SIZE bsize, int havSub) {
  DECLARE_ALIGNED(16, uint8_t, comp_pred[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, pred[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(16, uint8_t, ref[MAX_SB_SQUARE]);

  // init data
  av1_init_wedge_masks();
  const int w = block_size_wide[bsize];
  const int h = block_size_high[bsize];
  memset(comp_pred, 0, MAX_SB_SQUARE);
  for (int i = 0; i < MAX_SB_SQUARE; ++i) {
    pred[i] = rnd_.Rand8();
    ref[i] = rnd_.Rand8();
  }

  const int subx = havSub ? 3 : 0;
  const int suby = havSub ? 4 : 0;

  int wedge_types = (1 << get_wedge_bits_lookup(bsize));
  int wedge_index = wedge_types / 2;
  const uint8_t *mask = av1_get_contiguous_soft_mask(wedge_index, 1, bsize);

  const int num_loops = 1000000000 / (w + h);
  comp_mask_up_pred_func funcs[2] = { &aom_comp_mask_upsampled_pred_c,
                                      test_impl };
  double elapsed_time[2] = { 0 };
  for (int i = 0; i < 2; ++i) {
    aom_usec_timer timer;
    aom_usec_timer_start(&timer);

    comp_mask_up_pred_func func = funcs[i];
    for (int j = 0; j < num_loops; ++j) {
      func(comp_pred, pred, w, h, subx, suby, ref, MAX_SB_SIZE, mask, w, 0);
    }

    aom_usec_timer_mark(&timer);
    double time = static_cast<double>(aom_usec_timer_elapsed(&timer));
    elapsed_time[i] = 1000.0 * time / num_loops;
  }
  printf("CompMask[%d] %3dx%-3d:%7.2f/%7.2fns", havSub, w, h, elapsed_time[0],
         elapsed_time[1]);
  printf("(%3.2f)\n", elapsed_time[0] / elapsed_time[1]);
}

void AV1CompMaskUpVarianceTest::RunSpeedTest(comp_mask_up_pred_func test_impl,
                                             BLOCK_SIZE bsize) {
  RunSpeedTestSub(test_impl, bsize, 0);  // could skip upsample
  RunSpeedTestSub(test_impl, bsize, 1);
}

TEST_P(AV1CompMaskUpVarianceTest, CheckOutput) {
  RunCheckOutput(GET_PARAM(0), GET_PARAM(1));
}

TEST_P(AV1CompMaskUpVarianceTest, DISABLED_Speed) {
  RunSpeedTest(GET_PARAM(0), GET_PARAM(1));
}

const CompMaskUpPredParam kArrayCompMaskUpPred_ssse3[] = {
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_8X8),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_8X16),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_16X8),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_16X16),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_16X32),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_32X16),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_32X32),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_8X32),
  testing::make_tuple(&aom_comp_mask_upsampled_pred_ssse3, BLOCK_32X8),
};

INSTANTIATE_TEST_CASE_P(SSSE3, AV1CompMaskUpVarianceTest,
                        ::testing::ValuesIn(kArrayCompMaskUpPred_ssse3));
}  // namespace AV1CompMaskVariance
