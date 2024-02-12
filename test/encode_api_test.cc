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

#include <cassert>
#include <cstdlib>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "config/aom_config.h"

#include "aom/aomcx.h"
#include "aom/aom_encoder.h"
#include "aom/aom_image.h"

namespace {

#if CONFIG_REALTIME_ONLY
const int kUsage = AOM_USAGE_REALTIME;
#else
const int kUsage = AOM_USAGE_GOOD_QUALITY;
#endif

TEST(EncodeAPI, InvalidParams) {
  uint8_t buf[1] = { 0 };
  aom_image_t img;
  aom_codec_ctx_t enc;
  aom_codec_enc_cfg_t cfg;

  EXPECT_EQ(&img, aom_img_wrap(&img, AOM_IMG_FMT_I420, 1, 1, 1, buf));

  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_init(nullptr, nullptr, nullptr, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_init(&enc, nullptr, nullptr, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_encode(nullptr, nullptr, 0, 0, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM, aom_codec_encode(nullptr, &img, 0, 0, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM, aom_codec_destroy(nullptr));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_config_default(nullptr, nullptr, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_config_default(nullptr, &cfg, 0));
  EXPECT_NE(aom_codec_error(nullptr), nullptr);

  aom_codec_iface_t *iface = aom_codec_av1_cx();
  SCOPED_TRACE(aom_codec_iface_name(iface));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_init(nullptr, iface, nullptr, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_init(&enc, iface, nullptr, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_enc_config_default(iface, &cfg, 3));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_config_default(iface, &cfg, kUsage));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface, &cfg, 0));
  EXPECT_EQ(nullptr, aom_codec_get_global_headers(nullptr));

  aom_fixed_buf_t *glob_headers = aom_codec_get_global_headers(&enc);
  EXPECT_NE(glob_headers->buf, nullptr);
  if (glob_headers) {
    free(glob_headers->buf);
    free(glob_headers);
  }
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_encode(&enc, nullptr, 0, 0, 0));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));
}

TEST(EncodeAPI, InvalidControlId) {
  aom_codec_iface_t *iface = aom_codec_av1_cx();
  aom_codec_ctx_t enc;
  aom_codec_enc_cfg_t cfg;
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_config_default(iface, &cfg, kUsage));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface, &cfg, 0));
  EXPECT_EQ(AOM_CODEC_ERROR, aom_codec_control(&enc, -1, 0));
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM, aom_codec_control(&enc, 0, 0));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));
}

TEST(EncodeAPI, SetSFrameOnFirstFrame) {
  constexpr int kWidth = 2;
  constexpr int kHeight = 128;
  unsigned char kBuffer[kWidth * kHeight * 3] = { 0 };
  aom_image_t img;
  ASSERT_EQ(aom_img_wrap(&img, AOM_IMG_FMT_I420, kWidth, kHeight, 1, kBuffer),
            &img);

  aom_codec_iface_t *iface = aom_codec_av1_cx();
  aom_codec_enc_cfg_t cfg;
  ASSERT_EQ(aom_codec_enc_config_default(iface, &cfg, kUsage), AOM_CODEC_OK);
  cfg.g_w = kWidth;
  cfg.g_h = kHeight;

  aom_codec_ctx_t enc;
  ASSERT_EQ(aom_codec_enc_init(&enc, iface, &cfg, 0), AOM_CODEC_OK);
  // One of these aom_codec_encode() calls should fail.
  if (aom_codec_encode(&enc, &img, 0, 1, AOM_EFLAG_SET_S_FRAME) ==
      AOM_CODEC_OK) {
    EXPECT_NE(aom_codec_encode(&enc, NULL, 0, 0, 0), AOM_CODEC_OK);
  }
  EXPECT_EQ(aom_codec_destroy(&enc), AOM_CODEC_OK);
}

TEST(EncodeAPI, MonochromeInProfiles) {
  aom_codec_iface_t *iface = aom_codec_av1_cx();
  aom_codec_enc_cfg_t cfg;
  ASSERT_EQ(AOM_CODEC_OK, aom_codec_enc_config_default(iface, &cfg, kUsage));
  cfg.g_w = 128;
  cfg.g_h = 128;
  cfg.monochrome = 1;
  aom_codec_ctx_t enc;

  // Test Profile 0
  cfg.g_profile = 0;
  ASSERT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface, &cfg, 0));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));

  // Test Profile 1
  cfg.g_profile = 1;
  ASSERT_EQ(AOM_CODEC_INVALID_PARAM, aom_codec_enc_init(&enc, iface, &cfg, 0));

  // Test Profile 3
  cfg.g_profile = 2;
  ASSERT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface, &cfg, 0));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));
}

class AV1Encoder {
 public:
  explicit AV1Encoder(int speed) : speed_(speed) {}
  ~AV1Encoder();

  void Configure(unsigned int threads, unsigned int width, unsigned int height,
                 aom_rc_mode end_usage, unsigned int usage);
  void Encode(bool key_frame);

 private:
  // Flushes the encoder. Should be called after all the Encode() calls.
  void Flush();

  const int speed_;
  bool initialized_ = false;
  aom_codec_enc_cfg_t cfg_;
  aom_codec_ctx_t enc_;
  int frame_index_ = 0;
};

AV1Encoder::~AV1Encoder() {
  if (initialized_) {
    Flush();
    EXPECT_EQ(aom_codec_destroy(&enc_), AOM_CODEC_OK);
  }
}

void AV1Encoder::Configure(unsigned int threads, unsigned int width,
                           unsigned int height, aom_rc_mode end_usage,
                           unsigned int usage) {
  if (!initialized_) {
    aom_codec_iface_t *const iface = aom_codec_av1_cx();
    ASSERT_EQ(aom_codec_enc_config_default(iface, &cfg_, usage), AOM_CODEC_OK);
    cfg_.g_threads = threads;
    cfg_.g_w = width;
    cfg_.g_h = height;
    cfg_.g_forced_max_frame_width = cfg_.g_w;
    cfg_.g_forced_max_frame_height = cfg_.g_h;
    cfg_.g_timebase.num = 1;
    cfg_.g_timebase.den = 1000 * 1000;  // microseconds
    cfg_.g_pass = AOM_RC_ONE_PASS;
    cfg_.g_lag_in_frames = 0;
    cfg_.rc_end_usage = end_usage;
    cfg_.rc_min_quantizer = 2;
    cfg_.rc_max_quantizer = 58;
    ASSERT_EQ(aom_codec_enc_init(&enc_, iface, &cfg_, 0), AOM_CODEC_OK);
    ASSERT_EQ(aom_codec_control(&enc_, AOME_SET_CPUUSED, speed_), AOM_CODEC_OK);
    initialized_ = true;
    return;
  }

  ASSERT_EQ(usage, cfg_.g_usage);
  cfg_.g_threads = threads;
  cfg_.g_w = width;
  cfg_.g_h = height;
  cfg_.rc_end_usage = end_usage;
  ASSERT_EQ(aom_codec_enc_config_set(&enc_, &cfg_), AOM_CODEC_OK)
      << aom_codec_error_detail(&enc_);
}

void AV1Encoder::Encode(bool key_frame) {
  assert(initialized_);
  // TODO(wtc): Support high bit depths and other YUV formats.
  aom_image_t *const image =
      CreateGrayImage(AOM_IMG_FMT_I420, cfg_.g_w, cfg_.g_h);
  ASSERT_NE(image, nullptr);
  const aom_enc_frame_flags_t flags = key_frame ? AOM_EFLAG_FORCE_KF : 0;
  ASSERT_EQ(aom_codec_encode(&enc_, image, frame_index_, 1, flags),
            AOM_CODEC_OK);
  frame_index_++;
  const aom_codec_cx_pkt_t *pkt;
  aom_codec_iter_t iter = nullptr;
  while ((pkt = aom_codec_get_cx_data(&enc_, &iter)) != nullptr) {
    ASSERT_EQ(pkt->kind, AOM_CODEC_CX_FRAME_PKT);
    if (key_frame) {
      ASSERT_EQ(pkt->data.frame.flags & AOM_FRAME_IS_KEY, AOM_FRAME_IS_KEY);
    }
  }
  aom_img_free(image);
}

void AV1Encoder::Flush() {
  bool got_data;
  do {
    ASSERT_EQ(aom_codec_encode(&enc_, nullptr, 0, 0, 0), AOM_CODEC_OK);
    got_data = false;
    const aom_codec_cx_pkt_t *pkt;
    aom_codec_iter_t iter = nullptr;
    while ((pkt = aom_codec_get_cx_data(&enc_, &iter)) != nullptr) {
      ASSERT_EQ(pkt->kind, AOM_CODEC_CX_FRAME_PKT);
      got_data = true;
    }
  } while (got_data);
}

TEST(EncodeAPI, Buganizer314858909) {
  AV1Encoder encoder(7);

  encoder.Configure(6, 1582, 750, AOM_CBR, AOM_USAGE_REALTIME);

  // Encode a frame.
  encoder.Encode(false);

  encoder.Configure(0, 1582, 23, AOM_CBR, AOM_USAGE_REALTIME);

  // Encode a frame..
  encoder.Encode(false);

  encoder.Configure(16, 1542, 363, AOM_CBR, AOM_USAGE_REALTIME);

  // Encode a frame..
  encoder.Encode(false);
}

// This test covers a possible use case where the change of frame sizes and
// thread numbers happens before and after the first frame coding.
TEST(EncodeAPI, Buganizer310455204) {
  AV1Encoder encoder(7);

  encoder.Configure(0, 1915, 503, AOM_VBR, AOM_USAGE_REALTIME);

  encoder.Configure(4, 1, 1, AOM_VBR, AOM_USAGE_REALTIME);

  encoder.Configure(6, 559, 503, AOM_CBR, AOM_USAGE_REALTIME);

  // Encode a frame.
  encoder.Encode(false);

  // Increase the number of threads.
  encoder.Configure(16, 1915, 503, AOM_CBR, AOM_USAGE_REALTIME);

  // Encode a frame.
  encoder.Encode(false);
}

#if !CONFIG_REALTIME_ONLY
TEST(EncodeAPI, AllIntraMode) {
  aom_codec_iface_t *iface = aom_codec_av1_cx();
  aom_codec_ctx_t enc;
  aom_codec_enc_cfg_t cfg;
  EXPECT_EQ(AOM_CODEC_OK,
            aom_codec_enc_config_default(iface, &cfg, AOM_USAGE_ALL_INTRA));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface, &cfg, 0));
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));

  // Set g_lag_in_frames to a nonzero value. This should cause
  // aom_codec_enc_init() to fail.
  EXPECT_EQ(AOM_CODEC_OK,
            aom_codec_enc_config_default(iface, &cfg, AOM_USAGE_ALL_INTRA));
  cfg.g_lag_in_frames = 1;
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM, aom_codec_enc_init(&enc, iface, &cfg, 0));

  // Set kf_max_dist to a nonzero value. This should cause aom_codec_enc_init()
  // to fail.
  EXPECT_EQ(AOM_CODEC_OK,
            aom_codec_enc_config_default(iface, &cfg, AOM_USAGE_ALL_INTRA));
  cfg.kf_max_dist = 1;
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM, aom_codec_enc_init(&enc, iface, &cfg, 0));
}
#endif

}  // namespace
