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

#include <cstdio>
#include <cstdlib>
#include <string>
#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/util.h"
#include "test/md5_helper.h"
#include "aom_mem/aom_mem.h"

namespace {
<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
class TileIndependenceTest
    : public ::libaom_test::EncoderTest,
      public ::libaom_test::CodecTestWith2Params<int, int> {
=======
class TileIndependenceTest : public ::libaom_test::EncoderTest,
                             public ::libaom_test::CodecTestWithParam<int> {
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
 protected:
  TileIndependenceTest()
      : EncoderTest(GET_PARAM(0)), md5_fw_order_(), md5_inv_order_(),
<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
        n_tile_cols_(GET_PARAM(1)), n_tile_rows_(GET_PARAM(2)) {
=======
        n_tiles_(GET_PARAM(1)) {
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
    init_flags_ = AOM_CODEC_USE_PSNR;
    aom_codec_dec_cfg_t cfg = aom_codec_dec_cfg_t();
    cfg.w = 704;
    cfg.h = 144;
    cfg.threads = 1;
    fw_dec_ = codec_->CreateDecoder(cfg, 0);
    inv_dec_ = codec_->CreateDecoder(cfg, 0);
    inv_dec_->Control(AV1_INVERT_TILE_DECODE_ORDER, 1);
<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)

#if CONFIG_AV1 && CONFIG_EXT_TILE
    if (fw_dec_->IsAV1() && inv_dec_->IsAV1()) {
      fw_dec_->Control(AV1_SET_DECODE_TILE_ROW, -1);
      fw_dec_->Control(AV1_SET_DECODE_TILE_COL, -1);
      inv_dec_->Control(AV1_SET_DECODE_TILE_ROW, -1);
      inv_dec_->Control(AV1_SET_DECODE_TILE_COL, -1);
    }
#endif
=======
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
  }

  virtual ~TileIndependenceTest() {
    delete fw_dec_;
    delete inv_dec_;
  }

  virtual void SetUp() {
    InitializeConfig();
    SetMode(libaom_test::kTwoPassGood);
  }

  virtual void PreEncodeFrameHook(libaom_test::VideoSource *video,
                                  libaom_test::Encoder *encoder) {
    if (video->frame() == 1) {
<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
      encoder->Control(AV1E_SET_TILE_COLUMNS, n_tile_cols_);
      encoder->Control(AV1E_SET_TILE_ROWS, n_tile_rows_);
      SetCpuUsed(encoder);
=======
      encoder->Control(AV1E_SET_TILE_COLUMNS, n_tiles_);
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
    }
  }

<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
  virtual void SetCpuUsed(libaom_test::Encoder *encoder) {
    static const int kCpuUsed = 3;
    encoder->Control(AOME_SET_CPUUSED, kCpuUsed);
  }

=======
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
  void UpdateMD5(::libaom_test::Decoder *dec, const aom_codec_cx_pkt_t *pkt,
                 ::libaom_test::MD5 *md5) {
    const aom_codec_err_t res = dec->DecodeFrame(
        reinterpret_cast<uint8_t *>(pkt->data.frame.buf), pkt->data.frame.sz);
    if (res != AOM_CODEC_OK) {
      abort_ = true;
      ASSERT_EQ(AOM_CODEC_OK, res);
    }
    const aom_image_t *img = dec->GetDxData().Next();
    md5->Add(img);
  }

  virtual void FramePktHook(const aom_codec_cx_pkt_t *pkt) {
    UpdateMD5(fw_dec_, pkt, &md5_fw_order_);
    UpdateMD5(inv_dec_, pkt, &md5_inv_order_);
  }

<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
  void DoTest() {
    const aom_rational timebase = { 33333333, 1000000000 };
    cfg_.g_timebase = timebase;
    cfg_.rc_target_bitrate = 500;
    cfg_.g_lag_in_frames = 12;
    cfg_.rc_end_usage = AOM_VBR;

    libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 704, 576,
                                       timebase.den, timebase.num, 0, 5);
    ASSERT_NO_FATAL_FAILURE(RunLoop(&video));

    const char *md5_fw_str = md5_fw_order_.Get();
    const char *md5_inv_str = md5_inv_order_.Get();
    ASSERT_STREQ(md5_fw_str, md5_inv_str);
  }

=======
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
  ::libaom_test::MD5 md5_fw_order_, md5_inv_order_;
  ::libaom_test::Decoder *fw_dec_, *inv_dec_;

 private:
  int n_tile_cols_;
  int n_tile_rows_;
};

// run an encode with 2 or 4 tiles, and do the decode both in normal and
// inverted tile ordering. Ensure that the MD5 of the output in both cases
// is identical. If so, tiles are considered independent and the test passes.
<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
TEST_P(TileIndependenceTest, MD5Match) { DoTest(); }
=======
TEST_P(TileIndependenceTest, MD5Match) {
  const aom_rational timebase = { 33333333, 1000000000 };
  cfg_.g_timebase = timebase;
  cfg_.rc_target_bitrate = 500;
  cfg_.g_lag_in_frames = 25;
  cfg_.rc_end_usage = AOM_VBR;
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)

<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
class TileIndependenceTestLarge : public TileIndependenceTest {
  virtual void SetCpuUsed(libaom_test::Encoder *encoder) {
    static const int kCpuUsed = 0;
    encoder->Control(AOME_SET_CPUUSED, kCpuUsed);
  }
};
=======
  libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 704, 144,
                                     timebase.den, timebase.num, 0, 30);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)

TEST_P(TileIndependenceTestLarge, MD5Match) { DoTest(); }

<<<<<<< HEAD   (005ff8 Merge "warped_motion: Fix ubsan warning for signed integer o)
#if CONFIG_EC_ADAPT
// TODO(thdavies): EC_ADAPT does not support tiles
#else
#if CONFIG_EXT_TILE
AV1_INSTANTIATE_TEST_CASE(TileIndependenceTest, ::testing::Values(1, 2, 32),
                          ::testing::Values(1, 2, 32));
AV1_INSTANTIATE_TEST_CASE(TileIndependenceTestLarge,
                          ::testing::Values(1, 2, 32),
                          ::testing::Values(1, 2, 32));
#else
AV1_INSTANTIATE_TEST_CASE(TileIndependenceTest, ::testing::Values(0, 1),
                          ::testing::Values(0, 1));
AV1_INSTANTIATE_TEST_CASE(TileIndependenceTestLarge, ::testing::Values(0, 1),
                          ::testing::Values(0, 1));
#endif  // CONFIG_EXT_TILE
#endif  // CONFIG_EC_ADAPT
=======
  // could use ASSERT_EQ(!memcmp(.., .., 16) here, but this gives nicer
  // output if it fails. Not sure if it's helpful since it's really just
  // a MD5...
  ASSERT_STREQ(md5_fw_str, md5_inv_str);
}
#if CONFIG_EC_ADAPT
// TODO(thdavies): EC_ADAPT does not support tiles
INSTANTIATE_TEST_CASE_P(
    DISABLED_AV1, TileIndependenceTest,
    ::testing::Combine(
        ::testing::Values(
            static_cast<const libaom_test::CodecFactory *>(&libaom_test::kAV1)),
        ::testing::Range(0, 2, 1)));
#else
AV1_INSTANTIATE_TEST_CASE(TileIndependenceTest, ::testing::Range(0, 2, 1));
#endif

>>>>>>> BRANCH (5bf37c Use --enable-daala_ec by default.)
}  // namespace
