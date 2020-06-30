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

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/util.h"
namespace {

class AltRefForcedKeyTestLarge
    : public ::libaom_test::CodecTestWith2Params<libaom_test::TestMode, int>,
      public ::libaom_test::EncoderTest {
 protected:
  AltRefForcedKeyTestLarge()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        cpu_used_(GET_PARAM(2)), forced_kf_frame_num_(1), frame_num_(0) {}
  virtual ~AltRefForcedKeyTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);
    cfg_.rc_end_usage = AOM_VBR;
    cfg_.g_threads = 0;
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, cpu_used_);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
#if CONFIG_AV1_ENCODER
      // override test default for tile columns if necessary.
      if (GET_PARAM(0) == &libaom_test::kAV1) {
        encoder->Control(AV1E_SET_TILE_COLUMNS, 6);
      }
#endif
    }
    frame_flags_ =
        (video->frame() == forced_kf_frame_num_) ? AOM_EFLAG_FORCE_KF : 0;
  }

  virtual void FramePktHook(const aom_codec_cx_pkt_t *pkt) {
    if (frame_num_ == forced_kf_frame_num_) {
      ASSERT_EQ(pkt->data.frame.flags & AOM_FRAME_IS_KEY,
                static_cast<aom_codec_frame_flags_t>(AOM_FRAME_IS_KEY))
          << "Frame #" << frame_num_ << " isn't a keyframe!";
    }
    ++frame_num_;
  }

  ::libaom_test::TestMode encoding_mode_;
  int cpu_used_;
  unsigned int forced_kf_frame_num_;
  unsigned int frame_num_;
};

TEST_P(AltRefForcedKeyTestLarge, Frame1IsKey) {
  const aom_rational timebase = { 1, 30 };
  const int lag_values[] = { 3, 15, 25, -1 };

  forced_kf_frame_num_ = 1;
  for (int i = 0; lag_values[i] != -1; ++i) {
    frame_num_ = 0;
    cfg_.g_lag_in_frames = lag_values[i];
    libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                       timebase.den, timebase.num, 0, 30);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }
}

TEST_P(AltRefForcedKeyTestLarge, ForcedFrameIsKey) {
  const aom_rational timebase = { 1, 30 };
  const int lag_values[] = { 3, 15, 25, -1 };

  for (int i = 0; lag_values[i] != -1; ++i) {
    frame_num_ = 0;
    forced_kf_frame_num_ = lag_values[i] - 1;
    cfg_.g_lag_in_frames = lag_values[i];
    libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                       timebase.den, timebase.num, 0, 30);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }
}

AV1_INSTANTIATE_TEST_CASE(AltRefForcedKeyTestLarge,
                          ::testing::Values(::libaom_test::kOnePassGood),
                          ::testing::Values(2, 5));

typedef struct {
  const unsigned int min_kf_dist;
  const unsigned int max_kf_dist;
} kfIntervalParam;

const kfIntervalParam kfTestParams[] = {
  { 0, 10 }, { 10, 10 }, { 0, 30 }, { 30, 30 }
};

std::ostream &operator<<(std::ostream &os, const kfIntervalParam &test_arg) {
  return os << "kfIntervalParam { min_kf_dist:" << test_arg.min_kf_dist
            << " max_kf_dist:" << test_arg.max_kf_dist << " }";
}

// This class is used to check the presence of altref frame.
class AltrefFramePresenceTestLarge
    : public ::libaom_test::CodecTestWith3Params<libaom_test::TestMode,
                                                 kfIntervalParam, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  AltrefFramePresenceTestLarge()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        kf_dist_param_(GET_PARAM(2)), end_usage_check_(GET_PARAM(3)) {}
  virtual ~AltrefFramePresenceTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = end_usage_check_;
    cfg_.g_threads = 1;
    cfg_.kf_min_dist = kf_dist_param_.min_kf_dist;
    cfg_.kf_max_dist = kf_dist_param_.max_kf_dist;
    cfg_.g_lag_in_frames = 19;
  }

  virtual bool DoDecode() const { return 1; }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, 5);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    if (is_arf_frame_present_ != 1 && AOM_CODEC_OK == res_dec) {
      aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_ALTREF_FRAME_PRESENT,
                                    &is_arf_frame_present_);
    }
    return AOM_CODEC_OK == res_dec;
  }

  ::libaom_test::TestMode encoding_mode_;
  const kfIntervalParam kf_dist_param_;
  int is_arf_frame_present_;
  aom_rc_mode end_usage_check_;
};

TEST_P(AltrefFramePresenceTestLarge, AltrefFrameEncodePresenceTest) {
  is_arf_frame_present_ = 0;
  libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                     cfg_.g_timebase.den, cfg_.g_timebase.num,
                                     0, 150);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  ASSERT_EQ(is_arf_frame_present_, 1);
}

AV1_INSTANTIATE_TEST_CASE(AltrefFramePresenceTestLarge,
                          ::testing::Values(::libaom_test::kOnePassGood,
                                            ::libaom_test::kTwoPassGood),
                          ::testing::ValuesIn(kfTestParams),
                          ::testing::Values(AOM_Q, AOM_VBR, AOM_CBR, AOM_CQ));

}  // namespace
