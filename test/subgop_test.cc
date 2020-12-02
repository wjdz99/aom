/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
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
#include "test/y4m_video_source.h"
#include "test/i420_video_source.h"
#include "test/util.h"
#include "aom/aom_codec.h"

#include "av1/encoder/encoder.h"
#include "av1/encoder/subgop.h"

// Silence compiler warning for unused static functions
static void yuvconfig2image(aom_image_t *img, const YV12_BUFFER_CONFIG *yv12,
                            void *user_priv) AOM_UNUSED;
static aom_codec_err_t image2yuvconfig(const aom_image_t *img,
                                       YV12_BUFFER_CONFIG *yv12) AOM_UNUSED;
#include "av1/av1_iface_common.h"

#define MAX_SUBGOP_CODES 3

namespace {
// Default config
extern "C" const char subgop_config_str_def[];
// An enhanced config where the last subgop uses a shorter dist to arf
extern "C" const char subgop_config_str_enh[];
// A config that honors temporally scalable prediction structure, i.e.
// no frame is coded with references at higher pyramid depths.
extern "C" const char subgop_config_str_ts[];
// An asymmetrical config where the hierarchical frames are not exactly
// dyadic, but slightly skewed.
extern "C" const char subgop_config_str_asym[];
// low delay config without references
extern "C" const char subgop_config_str_ld[];

typedef enum {
  DEFAULT,
  ENHANCE,
  ASYMMETRIC,
  TEMPORAL_SCALABLE,
  LOW_DELAY,
} subgop_config_tag;

typedef struct {
  const char *preset_tag;
  const char *preset_str;
} subgop_config_str_preset_map_type;

const subgop_config_str_preset_map_type subgop_config_str_preset_map[] = {
  { "def", subgop_config_str_def },   { "enh", subgop_config_str_enh },
  { "asym", subgop_config_str_asym }, { "ts", subgop_config_str_ts },
  { "ld", subgop_config_str_ld },
};

typedef struct {
  const char *subgop_str;
  const char *input_file;
  int min_gf_interval;
  int max_gf_interval;
  int frame_w;
  int frame_h;
  int cpu_used;
} SubgopTestParams;

int is_extension_y4m(const char *filename) {
  const char *dot = strrchr(filename, '.');
  if (!dot || dot == filename)
    return 0;
  else
    return !strcmp(dot, ".y4m");
}

static const SubgopTestParams SubGopTestVectors[] = {
  // Default sub-gop config
  { subgop_config_str_preset_map[DEFAULT].preset_tag,
    "hantro_collage_w352h288.yuv", 0, 16, 352, 288, 3 },
  { subgop_config_str_preset_map[DEFAULT].preset_tag,
    "pixel_capture_w320h240.yuv", 0, 16, 320, 240, 3 },

  { subgop_config_str_preset_map[ENHANCE].preset_tag, "niklas_640_480_30.yuv",
    0, 15, 640, 480, 5 },
  { subgop_config_str_preset_map[ENHANCE].preset_tag, "paris_352_288_30.y4m", 0,
    6, 352, 288, 3 },

  { subgop_config_str_preset_map[ASYMMETRIC].preset_tag,
    "pixel_capture_w320h240.yuv", 0, 16, 320, 240, 5 },

  { subgop_config_str_preset_map[TEMPORAL_SCALABLE].preset_tag,
    "hantro_collage_w352h288.yuv", 0, 16, 352, 288, 5 },

  // TODO(vishnu) : Enable ld config
  // { subgop_config_str_preset_map[LOW_DELAY].preset_tag,
  // "paris_352_288_30.y4m",
  //   0, 16, 352, 288, 5 },
  // { subgop_config_str_preset_map[LOW_DELAY].preset_tag,
  // "desktop1.320_180.yuv",
  //   0, 16, 320, 180, 3 },

  // TODO(vishnu) : Add non-default subgop config
};

std::ostream &operator<<(std::ostream &os, const SubgopTestParams &test_arg) {
  return os << "SubgopTestParams { sub_gop_config:" << test_arg.subgop_str
            << " source_file:" << test_arg.input_file
            << " min_gf_interval:" << test_arg.min_gf_interval
            << " max_gf_interval:" << test_arg.max_gf_interval
            << " frame_width:" << test_arg.frame_w
            << " frame_height:" << test_arg.frame_h
            << " cpu_used:" << test_arg.cpu_used << " }";
}
// This class is used to validate the subgop config in a gop.
class SubGopTest
    : public ::libaom_test::CodecTestWith2Params<SubgopTestParams, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  SubGopTest()
      : EncoderTest(GET_PARAM(0)), subgop_test_params_(GET_PARAM(1)),
        rc_end_usage_(GET_PARAM(2)) {}
  virtual ~SubGopTest() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(::libaom_test::kOnePassGood);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.g_threads = 1;
    cfg_.rc_end_usage = rc_end_usage_;
    // Note: kf_min_dist, kf_max_dist, g_lag_in_frames are configurable
    // parameters
    cfg_.kf_min_dist = 65;
    cfg_.kf_max_dist = 65;
    cfg_.g_lag_in_frames = 35;
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, subgop_test_params_.cpu_used);
      encoder->Control(AV1E_SET_SUBGOP_CONFIG_STR,
                       subgop_test_params_.subgop_str);
      encoder->Control(AV1E_SET_MIN_GF_INTERVAL,
                       subgop_test_params_.min_gf_interval);
      encoder->Control(AV1E_SET_MAX_GF_INTERVAL,
                       subgop_test_params_.max_gf_interval);
    }
  }

  SubgopTestParams subgop_test_params_;
  aom_rc_mode rc_end_usage_;
};

TEST_P(SubGopTest, SubGopTest) {
  if (!is_extension_y4m(subgop_test_params_.input_file)) {
    libaom_test::I420VideoSource video(
        subgop_test_params_.input_file, subgop_test_params_.frame_w,
        subgop_test_params_.frame_h, cfg_.g_timebase.den, cfg_.g_timebase.num,
        0, 200);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  } else {
    ::libaom_test::Y4mVideoSource video(subgop_test_params_.input_file, 0, 200);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  }
}

AV1_INSTANTIATE_TEST_SUITE(SubGopTest, ::testing::ValuesIn(SubGopTestVectors),
                           ::testing::Values(AOM_Q, AOM_VBR, AOM_CQ, AOM_CBR));

}  // namespace
