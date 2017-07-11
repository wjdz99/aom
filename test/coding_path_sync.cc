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
#include "test/acm_random.h"

#include "./aom_config.h"
#include "aom_ports/mem.h"  // ROUND_POWER_OF_TWO
#include "aom/aomcx.h"
#include "aom/aomdx.h"
#include "aom/aom_encoder.h"
#include "aom/aom_decoder.h"

#include <vector>

using libaom_test::ACMRandom;
namespace {

struct CompressedSource {
  CompressedSource(int seed) : rnd_(seed) {
    frame_count_ = 0;
    aom_codec_iface_t *algo = &aom_codec_av1_cx_algo;

    aom_codec_enc_cfg_t cfg;
    aom_codec_enc_config_default(algo, &cfg, 0);

    cfg.g_w = WIDTH;
    cfg.g_h = HEIGHT;
    cfg.g_lag_in_frames = 0;

    aom_codec_enc_init(&enc_, algo, &cfg, 0);
  }

  ~CompressedSource() { aom_codec_destroy(&enc_); }

  const aom_codec_cx_pkt_t *readFrame() {
    uint8_t buf[WIDTH * HEIGHT * 3 / 2] = { 0 };
    for (size_t i = 0; i < sizeof buf; ++i) buf[i] = rnd_.Rand8();

    aom_image_t img;
    aom_img_wrap(&img, AOM_IMG_FMT_I420, WIDTH, HEIGHT, 0, buf);

    aom_codec_encode(&enc_, &img, frame_count_++, 1, 0, 0);

    aom_codec_iter_t iter = NULL;
    return aom_codec_get_cx_data(&enc_, &iter);
  }

 private:
  ACMRandom rnd_;
  aom_codec_ctx_t enc_;
  int frame_count_;
  static const int WIDTH = 4;
  static const int HEIGHT = 4;
};

// lowers an aom_image_t to a easily comparable/printable form
std::vector<short> serialize(const aom_image_t *img) {
  const int w_uv = ROUND_POWER_OF_TWO(img->d_w, img->x_chroma_shift);
  const int h_uv = ROUND_POWER_OF_TWO(img->d_h, img->y_chroma_shift);
  const int w[] = { (int)img->d_w, w_uv, w_uv };
  const int h[] = { (int)img->d_h, h_uv, h_uv };

  std::vector<short> bytes;
  bytes.reserve(img->d_w * img->d_h * 3);
  for (int plane = 0; plane < 3; ++plane)
    for (int r = 0; r < h[plane]; ++r)
      for (int c = 0; c < w[plane]; ++c)
        bytes.push_back(img->planes[plane][r * img->stride[plane] + c]);

  return bytes;
}

struct Decoder {
  Decoder(int allowLowbitdepth) {
    aom_codec_iface_t *algo = &aom_codec_av1_dx_algo;

    aom_codec_dec_cfg cfg = { 0 };
    cfg.allow_lowbitdepth = allowLowbitdepth;

    aom_codec_dec_init(&dec_, algo, &cfg, 0);
  }

  ~Decoder() { aom_codec_destroy(&dec_); }

  std::vector<short> decode(const aom_codec_cx_pkt_t *pkt) {
    aom_codec_decode(&dec_, (uint8_t *)pkt->data.frame.buf, pkt->data.frame.sz,
                     NULL, 0);

    aom_codec_iter_t iter = NULL;
    return serialize(aom_codec_get_frame(&dec_, &iter));
  }

 private:
  aom_codec_ctx_t dec_;
};

// Try to reveal a mismatch between LBD and HBD coding paths.
TEST(CodingPathSync, SearchForHbdLbdMismatch) {
  // TODO: re-enable this test once the LBD transforms match the HBD transforms
  return;

  for (int seed = 0; seed < 10; ++seed) {
    CompressedSource enc(seed);
    Decoder dec_HBD(0);
    Decoder dec_LBD(1);

    const aom_codec_cx_pkt_t *frame = enc.readFrame();
    ASSERT_EQ(dec_LBD.decode(frame), dec_HBD.decode(frame));
  }
}

}  // namespace
