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

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "./aom_config.h"
#include "./av1_rtcd.h"
#include "aom_ports/aom_timer.h"
#include "aom_ports/mem.h"
#include "av1/common/idct.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/scan.h"
#include "av1/common/txb_common.h"
#include "av1/decoder/decodeframe.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

namespace {
using libaom_test::ACMRandom;

typedef void (*DequantTxbFunc)(const uint8_t *const levels,
                               const int8_t *const signs,
                               const int16_t *const dequant,
                               const int16_t *const scan, const int bwl,
                               const int height, const int eob, const int shift,
                               tran_low_t *const tcoeffs);

class DecodeTxbTest : public ::testing::TestWithParam<DequantTxbFunc> {
 public:
  DecodeTxbTest() : dequant_txb_func_(GetParam()) {}

  virtual ~DecodeTxbTest() {}

  virtual void SetUp() {
    signs_ = reinterpret_cast<int8_t *>(
        aom_memalign(16, sizeof(*signs_) * MAX_TX_SQUARE));
    ASSERT_TRUE(signs_ != NULL);
    tcoeffs_ref_ = reinterpret_cast<tran_low_t *>(
        aom_memalign(16, sizeof(*tcoeffs_ref_) * MAX_TX_SQUARE));
    ASSERT_TRUE(tcoeffs_ref_ != NULL);
    tcoeffs_ = reinterpret_cast<tran_low_t *>(
        aom_memalign(16, sizeof(*tcoeffs_) * MAX_TX_SQUARE));
    ASSERT_TRUE(tcoeffs_ != NULL);
  }

  virtual void TearDown() {
    aom_free(signs_);
    aom_free(tcoeffs_ref_);
    aom_free(tcoeffs_);
    libaom_test::ClearSystemState();
  }

  void DequantTxbRun() {
    const int kNumTests = 100;
    int result = 0;
    int16_t dequant[2];
    const int qindex = rnd_.Rand8();
    const int delta = 0;
    aom_bit_depth_t bit_depth = AOM_BITS_8;

    dequant[0] = av1_dc_quant_QTX(qindex, delta, bit_depth);
    dequant[1] = av1_ac_quant_QTX(qindex, delta, bit_depth);

    for (int tx_size = TX_4X4; tx_size < TX_SIZES_ALL; ++tx_size) {
      const int width = tx_size_wide[tx_size];
      const int height = tx_size_high[tx_size];
      const int bwl = tx_size_wide_log2[tx_size];
      const int shift = av1_get_tx_scale((TX_SIZE)tx_size);
      const int16_t *const scan = av1_inter_scan_orders[tx_size][DCT_DCT].scan;

      levels_ = set_levels(levels_buf_, width);
      for (int i = 0; i < kNumTests && !result; ++i) {
        for (int eob = 0; eob < width * height && !result; ++eob) {
          InitDataWithEob(scan, bwl, eob);

          av1_dequant_txb_c(levels_, signs_, dequant, scan, bwl, height, eob,
                            shift, tcoeffs_ref_);
          dequant_txb_func_(levels_, signs_, dequant, scan, bwl, height, eob,
                            shift, tcoeffs_);

          PrintDiff(width, height);

          result = memcmp(tcoeffs_, tcoeffs_ref_,
                          sizeof(*tcoeffs_ref_) * MAX_TX_SQUARE);

          EXPECT_EQ(result, 0)
              << " width " << width << " height " << height << " eob " << eob
              << " shift " << shift << " dequant[0] " << dequant[0]
              << " dequant[1] " << dequant[1];
        }
      }
    }
  }

  void SpeedTestDequantTxbRun() {
    const int kNumTests = 10000000;
    aom_usec_timer timer;
    int16_t dequant[2];
    const int qindex = rnd_.Rand8();
    const int delta = 0;
    aom_bit_depth_t bit_depth = AOM_BITS_8;

    dequant[0] = av1_dc_quant_QTX(qindex, delta, bit_depth);
    dequant[1] = av1_ac_quant_QTX(qindex, delta, bit_depth);

    printf("Note: Only test the largest possible eob case!\n");
    for (int tx_size = TX_4X4; tx_size < TX_SIZES_ALL; ++tx_size) {
      const int width = tx_size_wide[tx_size];
      const int height = tx_size_high[tx_size];
      const int bwl = tx_size_wide_log2[tx_size];
      const int shift = av1_get_tx_scale((TX_SIZE)tx_size);
      const int16_t *const scan = av1_inter_scan_orders[tx_size][DCT_DCT].scan;
      const int eob = width * height - 1;

      levels_ = set_levels(levels_buf_, width);
      InitDataWithEob(scan, bwl, eob);

      aom_usec_timer_start(&timer);
      for (int i = 0; i < kNumTests; ++i) {
        av1_dequant_txb_c(levels_, signs_, dequant, scan, bwl, height, eob,
                          shift, tcoeffs_);
      }
      aom_usec_timer_mark(&timer);

      const int elapsed_time = static_cast<int>(aom_usec_timer_elapsed(&timer));
      printf("dequant_txb_%2dx%2d: %7.1f ms\n", width, height,
             elapsed_time / 1000.0);
    }
  }

 private:
  void InitDataWithEob(const int16_t *const scan, const int bwl,
                       const int eob) {
    memset(levels_buf_, 0, sizeof(levels_buf_));
    memset(tcoeffs_, 0, sizeof(*tcoeffs_) * MAX_TX_SQUARE);

    // Make sure signs_[] has random numbers in zero coefficients' positions.
    for (unsigned int i = 0; i < sizeof(signs_) / sizeof(*signs_); ++i) {
      signs_[i] = static_cast<int8_t>(rnd_.Rand16());
    }

    for (int c = 0; c < eob; ++c) {
      levels_[get_paded_idx(scan[c], bwl)] =
          static_cast<uint8_t>(clamp(rnd_.Rand8(), 0, INT8_MAX));
      signs_[scan[c]] = static_cast<int8_t>(rnd_.Rand16() & 1);
#if !CONFIG_HIGHBITDEPTH
      // TODO(linfengz): Shift right to avoid addition overflow. Use realistic
      // input coefficients instead later.
      tcoeffs_[scan[c]] = rnd_.Rand16() >> 1;
#else
      tcoeffs_[scan[c]] = rnd_.Rand31();
#endif
    }

    memcpy(tcoeffs_ref_, tcoeffs_, sizeof(*tcoeffs_) * MAX_TX_SQUARE);
  }

  void PrintDiff(const int width, const int height) const {
    if (memcmp(tcoeffs_, tcoeffs_ref_, sizeof(*tcoeffs_ref_) * MAX_TX_SQUARE)) {
      for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
          if (tcoeffs_ref_[y * width + x] != tcoeffs_[y * width + x]) {
            printf("tcoeffs_[%d][%d] diff:%6d (ref),%6d (opt)\n", y, x,
                   tcoeffs_ref_[y * width + x], tcoeffs_[y * width + x]);
            break;
          }
        }
      }
    }
  }

  DequantTxbFunc dequant_txb_func_;
  ACMRandom rnd_;
  uint8_t levels_buf_[TX_PAD_2D];
  uint8_t *levels_;
  int8_t *signs_;
  tran_low_t *tcoeffs_ref_;
  tran_low_t *tcoeffs_;
};

TEST_P(DecodeTxbTest, BitExact) { DequantTxbRun(); }

TEST_P(DecodeTxbTest, DISABLED_Speed) { SpeedTestDequantTxbRun(); }

#if HAVE_SSE2 && !CONFIG_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(SSE2, DecodeTxbTest,
                        ::testing::Values(av1_dequant_txb_sse2));
#endif

#if HAVE_SSE4_1 && CONFIG_HIGHBITDEPTH
INSTANTIATE_TEST_CASE_P(SSE4_1, DecodeTxbTest,
                        ::testing::Values(av1_dequant_txb_sse4_1));
#endif
}  // namespace
