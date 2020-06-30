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

#include <ostream>

#include "aom/aom_codec.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/y4m_video_source.h"
#include "test/util.h"

namespace {
typedef struct {
  // log2(number of tile rows)
  const unsigned int tile_rows;
  // log2(number of tile columns)
  const unsigned int tile_cols;
} uniformTileConfigParam;

const uniformTileConfigParam uniformTileConfigParams[] = {
  { 0, 0 }, { 0, 2 }, { 2, 0 }, { 1, 2 }, { 2, 2 }
};

typedef struct {
  // number of tile widths
  const unsigned int tile_width_count;
  // list of tile widths
  int tile_widths[64];
  // number of tile heights
  unsigned int tile_height_count;
  // list of tile heights
  int tile_heights[64];
} nonUniformTileConfigParam;

const nonUniformTileConfigParam nonUniformTileConfigParams[] = {
  { 1, { 3 }, 1, { 3 } },
  { 2, { 1, 2 }, 2, { 1, 2 } },
  { 3, { 2, 3, 4 }, 2, { 2, 3 } },
};

// This class is used to validate tile configuration for uniform spacing.
class UniformTileConfigTestLarge
    : public ::libaom_test::CodecTestWith3Params<
          libaom_test::TestMode, uniformTileConfigParam, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  UniformTileConfigTestLarge()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        tile_config_param_(GET_PARAM(2)), end_usage_check_(GET_PARAM(3)) {
    tile_config_violated_ = false;
  }
  virtual ~UniformTileConfigTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = end_usage_check_;
    cfg_.g_threads = 1;
    cfg_.g_lag_in_frames = 19;
  }

  virtual bool DoDecode() const { return 1; }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AV1E_SET_TILE_COLUMNS, tile_config_param_.tile_cols);
      encoder->Control(AV1E_SET_TILE_ROWS, tile_config_param_.tile_rows);
      encoder->Control(AOME_SET_CPUUSED, 5);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    if (AOM_CODEC_OK == res_dec) {
      aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
      aom_tile_info tile_info;
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_TILE_INFO, &tile_info);
      if (tile_info.tile_columns != (int)tile_config_param_.tile_cols ||
          tile_info.tile_rows != (int)tile_config_param_.tile_rows) {
        tile_config_violated_ = true;
      }
    }
    return AOM_CODEC_OK == res_dec;
  }

  ::libaom_test::TestMode encoding_mode_;
  const uniformTileConfigParam tile_config_param_;
  bool tile_config_violated_;
  aom_rc_mode end_usage_check_;
};

// This class is used to validate tile configuration for uniform spacing.
class NonUniformTileConfigTestLarge
    : public ::libaom_test::CodecTestWith3Params<
          libaom_test::TestMode, nonUniformTileConfigParam, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  NonUniformTileConfigTestLarge()
      : EncoderTest(GET_PARAM(0)), encoding_mode_(GET_PARAM(1)),
        tile_config_param_(GET_PARAM(2)), end_usage_check_(GET_PARAM(3)) {
    tile_config_violated_ = false;
  }
  virtual ~NonUniformTileConfigTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(encoding_mode_);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = end_usage_check_;
    cfg_.g_threads = 1;
    cfg_.g_lag_in_frames = 19;
    cfg_.tile_width_count = tile_config_param_.tile_width_count;
    memcpy(cfg_.tile_widths, tile_config_param_.tile_widths,
           sizeof(int) * tile_config_param_.tile_width_count);
    cfg_.tile_height_count = tile_config_param_.tile_height_count;
    memcpy(cfg_.tile_heights, tile_config_param_.tile_heights,
           sizeof(int) * tile_config_param_.tile_height_count);
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
    if (AOM_CODEC_OK == res_dec) {
      aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
      aom_tile_info tile_info;
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_TILE_INFO, &tile_info);

      // check validity of tile cols
      int i, j = 0;
      for (i = 0; i < tile_info.cols - 1;
           i++, j = (j + 1) % (int)tile_config_param_.tile_width_count) {
        if (tile_config_param_.tile_widths[j] != tile_info.tile_widths[i]) {
          tile_config_violated_ = true;
        }
      }
      if (tile_config_param_.tile_widths[j] < tile_info.tile_widths[i]) {
        tile_config_violated_ = true;
      }

      // check validity of tile rows
      j = 0;
      for (i = 0; i < tile_info.rows - 1;
           i++, j = (j + 1) % (int)tile_config_param_.tile_height_count) {
        if (tile_config_param_.tile_heights[j] != tile_info.tile_heights[i]) {
          tile_config_violated_ = true;
        }
      }
      if (tile_config_param_.tile_heights[j] < tile_info.tile_heights[i]) {
        tile_config_violated_ = true;
      }
    }
    return AOM_CODEC_OK == res_dec;
  }

  ::libaom_test::TestMode encoding_mode_;
  const nonUniformTileConfigParam tile_config_param_;
  bool tile_config_violated_;
  aom_rc_mode end_usage_check_;
};

TEST_P(UniformTileConfigTestLarge, UniformTileConfigTest) {
  ::libaom_test::Y4mVideoSource video("niklas_1280_720_30.y4m", 0, 6);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  ASSERT_EQ(tile_config_violated_, false);
}

TEST_P(NonUniformTileConfigTestLarge, NonUniformTileConfigTest) {
  ::libaom_test::Y4mVideoSource video("niklas_1280_720_30.y4m", 0, 1);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
  ASSERT_EQ(tile_config_violated_, false);
}

AV1_INSTANTIATE_TEST_CASE(UniformTileConfigTestLarge,
                          ::testing::Values(::libaom_test::kOnePassGood,
                                            ::libaom_test::kTwoPassGood),
                          ::testing::ValuesIn(uniformTileConfigParams),
                          ::testing::Values(AOM_Q, AOM_VBR, AOM_CBR, AOM_CQ));

AV1_INSTANTIATE_TEST_CASE(NonUniformTileConfigTestLarge,
                          ::testing::Values(::libaom_test::kOnePassGood,
                                            ::libaom_test::kTwoPassGood),
                          ::testing::ValuesIn(nonUniformTileConfigParams),
                          ::testing::Values(AOM_Q, AOM_VBR, AOM_CBR, AOM_CQ));
}  // namespace
