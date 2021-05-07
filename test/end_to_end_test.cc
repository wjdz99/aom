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
#include "test/util.h"
#include "test/y4m_video_source.h"
#include "test/yuv_video_source.h"

namespace {

const unsigned int kWidth = 160;
const unsigned int kHeight = 90;
const unsigned int kFramerate = 50;
const unsigned int kFrames = 10;
const int kBitrate = 500;
const unsigned int kCqLevel = 18;
// List of psnr thresholds for speed settings 0-8 and 4 encoding modes
const double kPsnrThreshold[][4] = {
  { 35.7, 44.4, 39.5, 41.9 }, { 35.7, 44.4, 39.5, 41.9 },
  { 35.7, 44.4, 39.4, 41.9 }, { 35.7, 44.4, 39.1, 41.8 },
  { 35.6, 44.4, 39.1, 41.8 }, { 35.0, 44.3, 38.7, 41.8 },
  { 35.0, 44.3, 38.7, 41.3 }, { 35.0, 44.3, 38.7, 40.8 },
  { 35.0, 44.3, 38.7, 40.8 }
};

// List of ssim thresholds for speed settings 0-8 with allintra encoding mode.
const double kSsimThreshold[] = { 83.4, 83.4, 83.4, 83.3, 83.3,
                                  83.0, 82.3, 81.1, 81.1 };

typedef struct {
  const char *filename;
  unsigned int input_bit_depth;
  aom_img_fmt fmt;
  aom_bit_depth_t bit_depth;
  unsigned int profile;
} TestVideoParam;

std::ostream &operator<<(std::ostream &os, const TestVideoParam &test_arg) {
  return os << "TestVideoParam { filename:" << test_arg.filename
            << " input_bit_depth:" << test_arg.input_bit_depth
            << " fmt:" << test_arg.fmt << " bit_depth:" << test_arg.bit_depth
            << " profile:" << test_arg.profile << " }";
}

const TestVideoParam kTestVectors[] = {
  { "park_joy_90p_8_420.y4m", 8, AOM_IMG_FMT_I420, AOM_BITS_8, 0 },
  { "park_joy_90p_8_422.y4m", 8, AOM_IMG_FMT_I422, AOM_BITS_8, 2 },
  { "park_joy_90p_8_444.y4m", 8, AOM_IMG_FMT_I444, AOM_BITS_8, 1 },
#if CONFIG_AV1_HIGHBITDEPTH
  { "park_joy_90p_10_420.y4m", 10, AOM_IMG_FMT_I42016, AOM_BITS_10, 0 },
  { "park_joy_90p_10_422.y4m", 10, AOM_IMG_FMT_I42216, AOM_BITS_10, 2 },
  { "park_joy_90p_10_444.y4m", 10, AOM_IMG_FMT_I44416, AOM_BITS_10, 1 },
  { "park_joy_90p_12_420.y4m", 12, AOM_IMG_FMT_I42016, AOM_BITS_12, 2 },
  { "park_joy_90p_12_422.y4m", 12, AOM_IMG_FMT_I42216, AOM_BITS_12, 2 },
  { "park_joy_90p_12_444.y4m", 12, AOM_IMG_FMT_I44416, AOM_BITS_12, 2 },
#endif
};

// Encoding modes tested
const libaom_test::TestMode kEncodingModeVectors[] = {
  ::libaom_test::kTwoPassGood,
  ::libaom_test::kOnePassGood,
  ::libaom_test::kRealTime,
};

// Speed settings tested
const int kCpuUsedVectors[] = { 1, 2, 3, 5, 6 };

int is_extension_y4m(const char *filename) {
  const char *dot = strrchr(filename, '.');
  if (!dot || dot == filename)
    return 0;
  else
    return !strcmp(dot, ".y4m");
}

class EndToEndTest
    : public ::libaom_test::CodecTestWith4Params<libaom_test::TestMode,
                                                 TestVideoParam, int, bool>,
      public ::libaom_test::EncoderTest {
 protected:
  EndToEndTest()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        test_video_param_(GET_PARAM(2)), cpu_used_(GET_PARAM(3)),
        enable_ssim_(GET_PARAM(4)), nframes_(0), psnr_(0.0), ssim_(0.0) {}

  virtual ~EndToEndTest() {}

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
    if (enable_ssim_ == 0) init_flags_ = AOM_CODEC_USE_PSNR;
  }

  virtual void BeginPassHook(unsigned int) {
    nframes_ = 0;
    psnr_ = 0.0;
    ssim_ = 0.0;
  }

  virtual void PSNRPktHook(const aom_codec_cx_pkt_t *pkt) {
    psnr_ += pkt->data.psnr.psnr[0];
    if (enable_ssim_ == 0) nframes_++;
  }

  virtual void calc_ssim_frame_level(const aom_image_t *img_src,
                                     const aom_image_t *img_enc,
                                     aom_bit_depth_t bd, unsigned int in_bd) {
    // When quality metric is psnr, do not evaluate ssim.
    if (enable_ssim_ == 0) return;
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
    uint8_t is_hbd = bd > AOM_BITS_8;
    if (is_hbd) {
      // HBD ssim calculation.
      uint8_t shift = bd - in_bd;
      for (int i = AOM_PLANE_Y; i < MAX_MB_PLANE; ++i) {
        const int is_uv = i > AOM_PLANE_Y;
        plane_ssim[i] = aom_highbd_ssim2(
            CONVERT_TO_BYTEPTR(img_src->planes[i]),
            CONVERT_TO_BYTEPTR(img_enc->planes[i]),
            img_src->stride[is_uv] >> is_hbd, img_enc->stride[is_uv] >> is_hbd,
            crop_widths[is_uv], crop_heights[is_uv], in_bd, shift);
      }
      frame_ssim = plane_ssim[AOM_PLANE_Y] * .8 +
                   .1 * (plane_ssim[AOM_PLANE_U] + plane_ssim[AOM_PLANE_V]);
      // Accumulate to find sequence level ssim value.
      ssim_ += frame_ssim;
      return;
    }
#else
    (void)bd;
    (void)in_bd;
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

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AV1E_SET_FRAME_PARALLEL_DECODING, 1);
      encoder->Control(AV1E_SET_TILE_COLUMNS, 4);
      encoder->Control(AOME_SET_CPUUSED, cpu_used_);
      if (enable_ssim_ == 1) encoder->Control(AOME_SET_TUNING, AOM_TUNE_SSIM);
      // Test screen coding tools at cpu_used = 1 && encoding mode is two-pass.
      if (cpu_used_ == 1 && encoding_mode_ == ::libaom_test::kTwoPassGood)
        encoder->Control(AV1E_SET_TUNE_CONTENT, AOM_CONTENT_SCREEN);
      else
        encoder->Control(AV1E_SET_TUNE_CONTENT, AOM_CONTENT_DEFAULT);
      if (encoding_mode_ == ::libaom_test::kOnePassGood ||
          encoding_mode_ == ::libaom_test::kTwoPassGood) {
        encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
        encoder->Control(AOME_SET_ARNR_MAXFRAMES, 7);
        encoder->Control(AOME_SET_ARNR_STRENGTH, 5);
      } else if (encoding_mode_ == ::libaom_test::kAllIntra) {
        encoder->Control(AOME_SET_CQ_LEVEL, kCqLevel);
      }
    }
  }

  double GetAveragePsnr() const {
    if (nframes_) return psnr_ / nframes_;
    return 0.0;
  }

  double GetAverageSsim() const {
    if (nframes_) return 100 * pow(ssim_ / nframes_, 8.0);
    return 0.0;
  }

  double GetPsnrThreshold() {
    return kPsnrThreshold[cpu_used_][encoding_mode_];
  }

  double GetSsimThreshold() { return kSsimThreshold[cpu_used_]; }

  void DoTest() {
    cfg_.rc_target_bitrate = kBitrate;
    cfg_.g_error_resilient = 0;
    cfg_.g_profile = test_video_param_.profile;
    cfg_.g_input_bit_depth = test_video_param_.input_bit_depth;
    cfg_.g_bit_depth = test_video_param_.bit_depth;
    if (cfg_.g_bit_depth > 8) init_flags_ |= AOM_CODEC_USE_HIGHBITDEPTH;

    std::unique_ptr<libaom_test::VideoSource> video;
    if (is_extension_y4m(test_video_param_.filename)) {
      video.reset(new libaom_test::Y4mVideoSource(test_video_param_.filename, 0,
                                                  kFrames));
    } else {
      video.reset(new libaom_test::YUVVideoSource(
          test_video_param_.filename, test_video_param_.fmt, kWidth, kHeight,
          kFramerate, 1, 0, kFrames));
    }
    ASSERT_TRUE(video.get() != NULL);
    ASSERT_NO_FATAL_FAILURE(RunLoop(video.get()));
    if (enable_ssim_ == 0) {
      const double psnr = GetAveragePsnr();
      EXPECT_GT(psnr, GetPsnrThreshold())
          << "encoding mode = " << encoding_mode_
          << ", cpu used = " << cpu_used_ << ", enable ssim = " << enable_ssim_;
    } else {
      const double ssim = GetAverageSsim();
      EXPECT_GT(ssim, GetSsimThreshold())
          << "encoding mode = " << encoding_mode_
          << ", cpu used = " << cpu_used_ << ", enable ssim = " << enable_ssim_;
    }
  }

 private:
  const libaom_test::TestMode encoding_mode_;
  const TestVideoParam test_video_param_;
  const int cpu_used_;
  const bool enable_ssim_;
  unsigned int nframes_;
  double psnr_;
  double ssim_;
};

class EndToEndTestLarge : public EndToEndTest {};

class EndToEndAllIntraTestLarge : public EndToEndTest {};

class EndToEndAllIntraTest : public EndToEndTest {};

TEST_P(EndToEndTestLarge, QualityMetricTest) { DoTest(); }

TEST_P(EndToEndTest, QualityMetricTest) { DoTest(); }

TEST_P(EndToEndAllIntraTestLarge, QualityMetricTest) { DoTest(); }

TEST_P(EndToEndAllIntraTest, QualityMetricTest) { DoTest(); }

AV1_INSTANTIATE_TEST_SUITE(EndToEndTestLarge,
                           ::testing::ValuesIn(kEncodingModeVectors),
                           ::testing::ValuesIn(kTestVectors),
                           ::testing::ValuesIn(kCpuUsedVectors),
                           ::testing::Values(0));  // enable_ssim

AV1_INSTANTIATE_TEST_SUITE(EndToEndTest,
                           ::testing::Values(::libaom_test::kTwoPassGood),
                           ::testing::Values(kTestVectors[2]),  // 444
                           ::testing::Values(3),                // cpu_used
                           ::testing::Values(0));               // enable_ssim

AV1_INSTANTIATE_TEST_SUITE(EndToEndAllIntraTestLarge,
                           ::testing::Values(::libaom_test::kAllIntra),
                           ::testing::ValuesIn(kTestVectors),
                           ::testing::Values(2, 4, 6, 8),  // cpu_used
                           ::testing::Values(0, 1));       // enable_ssim

AV1_INSTANTIATE_TEST_SUITE(EndToEndAllIntraTest,
                           ::testing::Values(::libaom_test::kAllIntra),
                           ::testing::Values(kTestVectors[0]),  // 420
                           ::testing::Values(6),                // cpu_used
                           ::testing::Values(0, 1));            // enable_ssim
}  // namespace
