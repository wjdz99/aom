/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
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
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/util.h"
#include "test/y4m_video_source.h"
#include "test/yuv_video_source.h"

namespace endtoend_test {
namespace {

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

class EndToEndTest : public ::libaom_test::EncoderTest {
 public:
  explicit EndToEndTest(const ::libaom_test::CodecFactory *codec)
      : EncoderTest(codec), nframes_(0), psnr_(0.0), ssim_(0.0) {}

 protected:
  ~EndToEndTest() override {}

  void BeginPassHook(unsigned int) override {
    nframes_ = 0;
    psnr_ = 0.0;
    ssim_ = 0.0;
  }

  void PSNRPktHook(const aom_codec_cx_pkt_t *pkt) override {
    psnr_ += pkt->data.psnr.psnr[0];
    if (metric_ == AOM_TUNE_PSNR) nframes_++;
  }

  void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                          ::libaom_test::Encoder *encoder) override {
    if (video->frame() == 0) {
      encoder->Control(AV1E_SET_FRAME_PARALLEL_DECODING, 1);
      encoder->Control(AV1E_SET_TILE_COLUMNS, 4);
      encoder->Control(AOME_SET_CPUUSED, cpu_used_);
      if (metric_ == AOM_TUNE_SSIM)
        encoder->Control(AOME_SET_TUNING, AOM_TUNE_SSIM);
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
    // InitializeConfig(encoding_mode_);
    cfg_.rc_target_bitrate = kBitrate;
    cfg_.g_error_resilient = 0;
    cfg_.g_profile = test_video_param_.profile;
    cfg_.g_input_bit_depth = test_video_param_.input_bit_depth;
    cfg_.g_bit_depth = test_video_param_.bit_depth;
    if (cfg_.g_bit_depth > 8) init_flags_ |= AOM_CODEC_USE_HIGHBITDEPTH;
    std::unique_ptr<libaom_test::VideoSource> video(
        new libaom_test::Y4mVideoSource(test_video_param_.filename, 0,
                                        kFrames));
    ASSERT_TRUE(video.get() != NULL);
    ASSERT_NO_FATAL_FAILURE(RunLoop(video.get()));
    if (metric_ == AOM_TUNE_PSNR) {
      const double psnr = GetAveragePsnr();
      EXPECT_GT(psnr, GetPsnrThreshold())
          << "encoding mode = " << encoding_mode_
          << ", cpu used = " << cpu_used_ << ", metric = " << metric_;
    } else {
      const double ssim = GetAverageSsim();
      EXPECT_GT(ssim, GetSsimThreshold())
          << "encoding mode = " << encoding_mode_
          << ", cpu used = " << cpu_used_ << ", metric = " << metric_;
    }
  }

  libaom_test::TestMode encoding_mode_;
  TestVideoParam test_video_param_;
  int cpu_used_;
  aom_tune_metric metric_;
  unsigned int nframes_;
  double psnr_;
  double ssim_;
};

}  // namespace
}  // namespace endtoend_test
