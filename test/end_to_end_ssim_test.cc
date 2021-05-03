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

class EndToEndSsimTest
    : public ::libaom_test::CodecTestWith4Params<
          libaom_test::TestMode, TestVideoParam, int, aom_tune_metric>,
      public EndToEndTest {
 public:
  EndToEndSsimTest() : EndToEndTest(GET_PARAM(0)) {
    encoding_mode_ = GET_PARAM(1);
    test_video_param_ = GET_PARAM(2);
    cpu_used_ = GET_PARAM(3);
    metric_ = GET_PARAM(4);
  }

 protected:
  ~EndToEndSsimTest() override {}

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
  }

  void CalculateFrameLevelSSIM(const aom_image_t *img_src,
                               const aom_image_t *img_enc,
                               aom_bit_depth_t bit_depth,
                               unsigned int input_bit_depth) override {
    // When quality metric is psnr, do not evaluate ssim.
    if (GET_PARAM(4) == AOM_TUNE_PSNR) return;
    double frame_ssim;
    double plane_ssim[MAX_MB_PLANE] = { 0.0, 0.0, 0.0 };
    int crop_widths[PLANE_TYPES];
    int crop_heights[PLANE_TYPES];
    crop_widths[PLANE_TYPE_Y] = img_src->d_w;
    crop_heights[PLANE_TYPE_Y] = img_src->d_h;
    // Width of UV planes calculated based on chroma_shift values.
    crop_widths[PLANE_TYPE_UV] =
        img_src->x_chroma_shift == 1 ? (img_src->w + 1) >> 1 : img_src->w;
    crop_heights[PLANE_TYPE_UV] =
        img_src->y_chroma_shift == 1 ? (img_src->h + 1) >> 1 : img_src->h;
    nframes_++;

#if CONFIG_AV1_HIGHBITDEPTH
    uint8_t is_hbd = bit_depth > AOM_BITS_8;
    if (is_hbd) {
      // HBD ssim calculation.
      uint8_t shift = bit_depth - input_bit_depth;
      for (int i = AOM_PLANE_Y; i < MAX_MB_PLANE; ++i) {
        const int is_uv = i > AOM_PLANE_Y;
        plane_ssim[i] = aom_highbd_ssim2(
            CONVERT_TO_BYTEPTR(img_src->planes[i]),
            CONVERT_TO_BYTEPTR(img_enc->planes[i]),
            img_src->stride[is_uv] >> is_hbd, img_enc->stride[is_uv] >> is_hbd,
            crop_widths[is_uv], crop_heights[is_uv], input_bit_depth, shift);
      }
      frame_ssim = plane_ssim[AOM_PLANE_Y] * .8 +
                   .1 * (plane_ssim[AOM_PLANE_U] + plane_ssim[AOM_PLANE_V]);
      // Accumulate to find sequence level ssim value.
      ssim_ += frame_ssim;
      return;
    }
#else
    (void)bit_depth;
    (void)input_bit_depth;
#endif  // CONFIG_AV1_HIGHBITDEPTH

    // LBD ssim calculation.
    for (int i = AOM_PLANE_Y; i < MAX_MB_PLANE; ++i) {
      const int is_uv = i > AOM_PLANE_Y;
      plane_ssim[i] = aom_ssim2(img_src->planes[i], img_enc->planes[i],
                                img_src->stride[is_uv], img_enc->stride[is_uv],
                                crop_widths[is_uv], crop_heights[is_uv]);
    }
    frame_ssim = plane_ssim[AOM_PLANE_Y] * .8 +
                 .1 * (plane_ssim[AOM_PLANE_U] + plane_ssim[AOM_PLANE_V]);
    // Accumulate to find sequence level ssim value.
    ssim_ += frame_ssim;
  }
};

class EndToEndSsimAllIntraTestLarge : public EndToEndSsimTest {};

class EndToEndSsimAllIntraTest : public EndToEndSsimTest {};

TEST_P(EndToEndSsimAllIntraTestLarge, SsimMetricTest) { DoTest(); }

TEST_P(EndToEndSsimAllIntraTest, SsimMetricTest) { DoTest(); }

AV1_INSTANTIATE_TEST_SUITE(EndToEndSsimAllIntraTestLarge,
                           ::testing::Values(::libaom_test::kAllIntra),
                           ::testing::ValuesIn(kTestVectors),
                           ::testing::Values(2, 4, 6, 8),      // cpu_used
                           ::testing::Values(AOM_TUNE_SSIM));  // metric

AV1_INSTANTIATE_TEST_SUITE(EndToEndSsimAllIntraTest,
                           ::testing::Values(::libaom_test::kAllIntra),
                           ::testing::Values(kTestVectors[0]),  // 420
                           ::testing::Values(6),                // cpu_used
                           ::testing::Values(AOM_TUNE_SSIM));   // metric
}  // namespace
}  // namespace endtoend_test
