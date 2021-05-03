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

#include <memory>
#include <ostream>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "aom_ports/mem.h"
#include "aom_dsp/ssim.h"
#include "av1/common/blockd.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/end_to_end_test.h"
#include "test/util.h"
#include "test/y4m_video_source.h"
#include "test/yuv_video_source.h"

namespace endtoend_test {
namespace {

class EndToEndPsnrTest
    : public ::libaom_test::CodecTestWith4Params<
          libaom_test::TestMode, TestVideoParam, int, aom_tune_metric>,
      public EndToEndTest {
 public:
  EndToEndPsnrTest() : EndToEndTest(GET_PARAM(0)) {
    encoding_mode_ = GET_PARAM(1);
    test_video_param_ = GET_PARAM(2);
    cpu_used_ = GET_PARAM(3);
    metric_ = GET_PARAM(4);
  }

 protected:
  ~EndToEndPsnrTest() override {}

  virtual void SetUp() {
    InitializeConfig(encoding_mode_);
    if (encoding_mode_ == ::libaom_test::kOnePassGood ||
        encoding_mode_ == ::libaom_test::kTwoPassGood) {
      cfg_.g_lag_in_frames = 5;
    } else if (encoding_mode_ == ::libaom_test::kRealTime) {
      cfg_.rc_buf_sz = 1000;
      cfg_.rc_buf_initial_sz = 500;
      cfg_.rc_buf_optimal_sz = 600;
    }
    if (metric_ == AOM_TUNE_PSNR) init_flags_ = AOM_CODEC_USE_PSNR;
  }
};

class EndToEndPsnrTestLarge : public EndToEndPsnrTest {};

class EndToEndPsnrAllIntraTestLarge : public EndToEndPsnrTest {};

class EndToEndPsnrAllIntraTest : public EndToEndPsnrTest {};

TEST_P(EndToEndPsnrTestLarge, PsnrMetricTest) { DoTest(); }

TEST_P(EndToEndPsnrTest, PsnrMetricTest) { DoTest(); }

TEST_P(EndToEndPsnrAllIntraTestLarge, PsnrMetricTest) { DoTest(); }

TEST_P(EndToEndPsnrAllIntraTest, PsnrMetricTest) { DoTest(); }

AV1_INSTANTIATE_TEST_SUITE(EndToEndPsnrTestLarge,
                           ::testing::ValuesIn(kEncodingModeVectors),
                           ::testing::ValuesIn(kTestVectors),
                           ::testing::ValuesIn(kCpuUsedVectors),
                           ::testing::Values(AOM_TUNE_PSNR));  // metric

AV1_INSTANTIATE_TEST_SUITE(EndToEndPsnrTest,
                           ::testing::Values(::libaom_test::kTwoPassGood),
                           ::testing::Values(kTestVectors[2]),  // 444
                           ::testing::Values(3),                // cpu_used
                           ::testing::Values(AOM_TUNE_PSNR));   // metric

AV1_INSTANTIATE_TEST_SUITE(EndToEndPsnrAllIntraTestLarge,
                           ::testing::Values(::libaom_test::kAllIntra),
                           ::testing::ValuesIn(kTestVectors),
                           ::testing::Values(2, 4, 6, 8),      // cpu_used
                           ::testing::Values(AOM_TUNE_PSNR));  // metric

AV1_INSTANTIATE_TEST_SUITE(EndToEndPsnrAllIntraTest,
                           ::testing::Values(::libaom_test::kAllIntra),
                           ::testing::Values(kTestVectors[0]),  // 420
                           ::testing::Values(6),                // cpu_used
                           ::testing::Values(AOM_TUNE_PSNR));   // metric
}  // namespace
}  // namespace endtoend_test
