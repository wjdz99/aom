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

/*
   Reproduces https://crbug.com/aomedia/3376. Emulates the command line:

   ./aomenc --cpu-used=6 --threads=10 --cq-level=14 --passes=1 --limit=1 \
     --lag-in-frames=0 --end-usage=q --deltaq-mode=3 --min-q=0 --max-q=63 \
     -o output.av1 niklas_1280_720_30.y4m
*/
TEST(EncodeAPI, DeltaQMode3MultiThread) {
  constexpr int kWidth = 1280;
  constexpr int kHeight = 720;
  // Dummy buffer of neutral gray samples.
  constexpr size_t kBufferSize = kWidth * kHeight + kWidth * kHeight / 2;
  std::vector<unsigned char> buffer(kBufferSize,
                                    static_cast<unsigned char>(128));

  aom_image_t img;
  EXPECT_EQ(&img, aom_img_wrap(&img, AOM_IMG_FMT_I420, kWidth, kHeight, 1,
                               buffer.data()));

  aom_codec_iface_t *iface = aom_codec_av1_cx();
  aom_codec_enc_cfg_t cfg;
  const unsigned int usage = AOM_USAGE_GOOD_QUALITY;
  EXPECT_EQ(aom_codec_enc_config_default(iface, &cfg, usage), AOM_CODEC_OK);
  cfg.g_w = kWidth;
  cfg.g_h = kHeight;
  cfg.g_threads = 10;
  cfg.rc_end_usage = AOM_Q;
  cfg.g_profile = 0;
  cfg.g_bit_depth = AOM_BITS_8;
  cfg.g_input_bit_depth = 8;
  cfg.g_lag_in_frames = 0;
  cfg.rc_min_quantizer = 0;
  cfg.rc_max_quantizer = 63;
  cfg.g_pass = AOM_RC_ONE_PASS;
  cfg.g_limit = 1;
  aom_codec_ctx_t enc;
  EXPECT_EQ(aom_codec_enc_init(&enc, iface, &cfg, 0), AOM_CODEC_OK);
  EXPECT_EQ(aom_codec_control(&enc, AOME_SET_CPUUSED, 6), AOM_CODEC_OK);
  EXPECT_EQ(aom_codec_control(&enc, AOME_SET_CQ_LEVEL, 14), AOM_CODEC_OK);
  EXPECT_EQ(aom_codec_control(&enc, AV1E_SET_DELTAQ_MODE, 3), AOM_CODEC_OK);
  EXPECT_EQ(aom_codec_set_option(&enc, "passes", "1"), AOM_CODEC_OK);
  EXPECT_EQ(aom_codec_control(&enc, AV1E_SET_COLOR_RANGE, AOM_CR_STUDIO_RANGE),
            AOM_CODEC_OK);
  EXPECT_EQ(aom_codec_encode(&enc, &img, 0, 1, 0), AOM_CODEC_OK);
  aom_codec_iter_t iter = nullptr;
  const aom_codec_cx_pkt_t *pkt = aom_codec_get_cx_data(&enc, &iter);
  EXPECT_NE(pkt, nullptr);
  EXPECT_EQ(pkt->kind, AOM_CODEC_CX_FRAME_PKT);
  // pkt->data.frame.flags is 0x1f0011.
  EXPECT_EQ(pkt->data.frame.flags, 0x1f0011u);
  EXPECT_EQ(pkt->data.frame.flags & AOM_FRAME_IS_KEY, AOM_FRAME_IS_KEY);
  pkt = aom_codec_get_cx_data(&enc, &iter);
  EXPECT_EQ(pkt, nullptr);

  // Flush encoder
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_encode(&enc, nullptr, 0, 1, 0));
  iter = nullptr;
  pkt = aom_codec_get_cx_data(&enc, &iter);
  EXPECT_EQ(pkt, nullptr);

  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));
}
#endif

}  // namespace
