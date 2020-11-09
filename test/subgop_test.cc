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

#include "av1/common/av1_common_int.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/subgop.h"

#define DEBUG 1
namespace {
typedef struct {
  const char *sub_gop_config;
  libaom_test::TestMode encoding_mode;
} SubgopTestParams;

static const SubgopTestParams InputParams[] = {
  { "16:2:16F1/8F2/4U3/2U4/1V5/2S/3V5/4S/6U4/5V5/6S/7V5/8R2/"
    "12U3/10U4/9V5/10S/11V5/12S/14U4/13V5/14S/15V5/16R1,"
    "16:0:16F1/8F2/4U3/2U4/1V5/2S/3V5/4S/6U4/5V5/6S/7V5/8R2/"
    "12U3/10U4/9V5/10S/11V5/12S/14U4/13V5/14S/15V5/16R1,"
    "16:1:13F1/7F2/4U3/2U4/1V5/2S/3V5/4S/5V4/6V5/7R2/"
    "10U3/8V4/9V5/10S/11V4/12V5/13R1/16U4/14V4/15V5/16S",
    ::libaom_test::kOnePassGood },

  { "6:0:6U1/4U2/2U3/1V4/2S/3V4/4S/5V4/6S,"
    "6:1:5U1/3U2/1V3/2V3/3S/4V3/5S/6V3,"
    "8:0:8F1/4U2/2U3/1V4/2S/3V4/4S/6U3/5V4/6S/7V4/8R1,"
    "8:1:7F1/3U2/1V4/2V4/3S/5U3/4V4/5S/6V4/7R1/8V4,"
    "10:0:10F1/7U2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7S/9U4/8V5/9S/10R1,"
    "10:1:8F1/4U2/2U3/1V4/2S/3V4/4S/6U3/5V4/6S/7V4/8R1/10U4/9V5/10S,"
    "11:0:11F1/7U2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7S/9U4/8V5/9S/10V5/"
    "11R1,"
    "11:1:9F1/4U2/2U3/1V4/2S/3V4/4S/7U3/5V4/6V5/7S/8V5/9R1/11U4/10V5/"
    "11S,"
    "12:0:12F1/7U2/3U3/1V5/2V5/3S/5U4/4V5/"
    "5S/6V5/7S/9U3/8V5/9S/11U4/10V5/11S/12R1,"
    "12:1:10F1/7U2/3U3/1V5/2V5/3S/5U4/4V5/"
    "5S/6V5/7S/9U4/8V5/9S/10R1/12U4/11V5/12S,"
    "13:0:13F1/7F2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7R2/"
    "10U3/8V5/9V5/10S/12U4/11V5/12S/13R1,"
    "13:1:10F1/6F2/3U3/1V5/2V5/3S/5U4/4V5/5S/6R2/"
    "8U4/7V5/8S/9V5/10R1/12U4/11V5/12S/13V5,"
    "14:0:14F1/7F2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7R2/"
    "10U3/8V5/9V5/10S/12U4/11V5/12S/13V5/14R1,"
    "14:1:11F1/7F2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7R2/"
    "9U4/8V5/9S/10V5/11R1/13U4/12V5/13S/14V5,"
    "15:0:15F1/7F2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7R2/"
    "11U3/9U4/8V5/9S/10V5/11S/13U4/12V5/13S/14V5/15R1,"
    "15:1:12F1/7F2/3U3/1V5/2V5/3S/5U4/4V5/5S/6V5/7R2/"
    "10U3/8V5/9V5/10S/11V5/12R1/15U4/13V5/14V5/15S,"
    "16:2:16F1/8F2/4U3/2U4/1V5/2S/3V5/4S/6U4/5V5/6S/7V5/8R2/"
    "12U3/10U4/9V5/10S/11V5/12S/14U4/13V5/14S/15V5/16R1,"
    "16:0:16F1/8F2/4U3/2U4/1V5/2S/3V5/4S/6U4/5V5/6S/7V5/8R2/"
    "12U3/10U4/9V5/10S/11V5/12S/14U4/13V5/14S/15V5/16R1,"
    "16:1:13F1/7F2/4U3/2U4/1V5/2S/3V5/4S/5V4/6V5/7R2/"
    "10U3/8V4/9V5/10S/11V4/12V5/13R1/16U4/14V4/15V5/16S",
    ::libaom_test::kOnePassGood },

  { "16:0:16F1/10F2/5U3/3U4/1V5/2V5/3S/"
    "4V5/5S/8U4/6V5/7V5/8S/9V5/10R2/"
    "13U3/11V5/12V5/13S/14V5/15V5/16R1,"
    "16:1:13F1/7F2/4U3/2U4/1V5/2S/3V5/4S/"
    "5V4/6V5/7R2/10U3/8V4/9V5/10S/11V4/12V5/"
    "13R1/16U4/14V4/15V5/16S",
    ::libaom_test::kOnePassGood },

  { "16:0:1V5/2V4/3V5/4V3/5V5/6V4/7V5/8V2/"
    "9V5/10V4/11V5/12V3/13V5/14V4/15V5/16V1,"
    "16:1:1V5/2V4/3V5/4V3/5V5/6V4/7V5/8V2/"
    "9V5/10V4/11V5/12V3/13V5/14V4/15V5/16V5",
    ::libaom_test::kOnePassGood },

};

std::ostream &operator<<(std::ostream &os, const SubgopTestParams &test_arg) {
  return os << "SubgopTestParams { sub_gop_config:" << test_arg.sub_gop_config
            << " encoding_mode:" << test_arg.encoding_mode << " }";
}

// This class is used to check the frame type in a gop.
class SubGopTestLarge
    : public ::libaom_test::CodecTestWith2Params<SubgopTestParams, aom_rc_mode>,
      public ::libaom_test::EncoderTest {
 protected:
  SubGopTestLarge()
      : EncoderTest(GET_PARAM(0)), subgop_test_params_(GET_PARAM(1)),
        rc_end_usage_(GET_PARAM(2)) {
    InitSubgop();
  }
  virtual ~SubGopTestLarge() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(subgop_test_params_.encoding_mode);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.g_threads = 1;
    // Note: kf_min_dist, kf_max_dist, g_lag_in_frames are configurable
    // parameters
    cfg_.kf_min_dist = 65;
    cfg_.kf_max_dist = 65;
    cfg_.g_lag_in_frames = 35;
  }

  virtual bool DoDecode() const { return 1; }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, 5);
      encoder->Control(AV1E_ENABLE_SUBGOP_STATS, enable_subgop_stats_);
      encoder->Control(AV1E_SET_SUBGOP_CONFIG_STR,
                       subgop_test_params_.sub_gop_config);
      av1_process_subgop_config_set(subgop_test_params_.sub_gop_config,
                                    &user_cfg_set_);
    }
  }

  virtual void PreDecodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Decoder *decoder) {
    aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
    if (video->frame() == 0)
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AV1D_ENABLE_SUBGOP_STATS,
                                    enable_subgop_stats_);
  }

  void InitSubgop() {
    user_cfg_set_ = { 0 };
    subgop_data_.num_steps = MAX_SUBGOP_STATS_SIZE;
    ResetSubgop();
    is_prev_frame_key_ = 0;
    frames_from_key_ = 0;
    frame_num_ = 0;
    // Note: This test case is not tested for 'CONFIG_REALTIME_ONLY'
    enable_subgop_stats_ = 1;
  }

  void ResetSubgop() {
    subgop_info_.is_user_specified = 0;
    subgop_info_.frames_to_key = 0;
    subgop_info_.gf_interval = 0;
    subgop_info_.size = 0;
    subgop_info_.code = SUBGOP_IN_GOP_GENERIC;

    for (int idx = 0; idx < MAX_SUBGOP_STATS_SIZE; idx++) {
      subgop_data_.step[idx].disp_frame_idx = -1;
      subgop_data_.step[idx].show_existing_frame = -1;
      subgop_data_.step[idx].show_frame = -1;
      subgop_data_.step[idx].is_filtered = -1;
    }
    subgop_data_.num_steps = 0;
    subgop_data_.step_idx_enc = 0;
    subgop_data_.step_idx_dec = 0;

    subgop_code_test_ = SUBGOP_IN_GOP_GENERIC;
    subgop_size_ = 0;
    frame_num_in_subgop_ = 0;
  }

  void DetermineSubgopCode(libaom_test::Encoder *encoder) {
    int gf_interval, frames_to_key;
    encoder->Control(AV1E_GET_FRAME_TYPE, &frame_type_test_);
    if (frame_type_test_ == KEY_FRAME) {
      is_prev_frame_key_ = 1;
      return;
    }
    frames_to_key = subgop_info_.frames_to_key;
    gf_interval = subgop_info_.gf_interval;
    if ((frames_to_key - 1) <= gf_interval + 2)
      subgop_code_test_ = SUBGOP_IN_GOP_LAST;
    else if (is_prev_frame_key_)
      subgop_code_test_ = SUBGOP_IN_GOP_FIRST;
    else
      subgop_code_test_ = SUBGOP_IN_GOP_GENERIC;
    is_prev_frame_key_ = 0;
    subgop_size_ = gf_interval;
  }

  virtual bool HandleEncodeResult(libaom_test::Encoder *encoder) {
    if (!frame_num_in_subgop_) {
      encoder->Control(AV1E_GET_SUB_GOP_CONFIG, &subgop_info_);
      DetermineSubgopCode(encoder);
      subgop_data_.num_steps = subgop_info_.subgop_cfg.num_steps;
    }
    frame_num_in_subgop_++;
    if (subgop_info_.is_user_specified)
      encoder->Control(AV1E_GET_FRAME_INFO, &subgop_data_);
    return 1;
  }

  void FillTestSubgopConfig() {
    int filtered_frames[REF_FRAMES] = { 0 }, buf_idx = 0;
    if (frame_type_test_ == KEY_FRAME) return;
    test_subgop_.num_frames = subgop_info_.size;
    test_subgop_.num_steps = subgop_data_.num_steps;
    test_subgop_.subgop_in_gop_code = subgop_info_.code;
    // Populating the filter-type of out-of-order frames appropriately for all
    // steps in sub-gop
    for (int idx = 0; idx < subgop_data_.num_steps; idx++) {
      test_subgop_.step[idx].disp_frame_idx =
          subgop_data_.step[idx].disp_frame_idx - frames_from_key_;
      if (subgop_data_.step[idx].is_filtered)
        filtered_frames[buf_idx++] =
            subgop_data_.step[idx].disp_frame_idx - frames_from_key_;
      else {
        for (char ref_frame = 0; ref_frame < buf_idx; ref_frame++) {
          if (test_subgop_.step[idx].disp_frame_idx ==
              filtered_frames[ref_frame])
            subgop_data_.step[idx].is_filtered = 1;
        }
      }
    }
    for (int idx = 0; idx < subgop_data_.num_steps; idx++) {
      FRAME_TYPE_CODE frame_type_code;
      int show_existing_frame = subgop_data_.step[idx].show_existing_frame;
      int show_frame = subgop_data_.step[idx].show_frame;
      int is_filtered = subgop_data_.step[idx].is_filtered;
      if (show_existing_frame == 0) {
        if (show_frame == 0)
          frame_type_code = (is_filtered == 1) ? FRAME_TYPE_OOO_FILTERED
                                               : FRAME_TYPE_OOO_UNFILTERED;
        else if (show_frame == 1)
          frame_type_code = FRAME_TYPE_INO_VISIBLE;
      } else if (show_existing_frame == 1) {
        if (show_frame == 1)
          frame_type_code = (is_filtered == 1) ? FRAME_TYPE_INO_REPEAT
                                               : FRAME_TYPE_INO_SHOWEXISTING;
      }
      test_subgop_.step[idx].type_code = frame_type_code;
    }
  }

  SubGOPCfg *DetermineSubgopConfig() {
    SubGOPCfg *subgop_cfg = user_cfg_set_.config;
    for (int idx = 0; idx < user_cfg_set_.num_configs; idx++) {
      if (subgop_cfg[idx].num_frames == subgop_size_) {
        if (ValidateSubgopCode()) return &subgop_cfg[idx];
      }
    }
    return NULL;
  }

  // Validates frametype(along with temporal filtering), frame coding order
  bool ValidateSubgopFrametype() {
    SubGOPCfg *subgop_cfg = &subgop_info_.subgop_cfg;
    for (int idx = 0; idx < subgop_cfg->num_steps; idx++) {
      EXPECT_EQ(subgop_cfg->step[idx].disp_frame_idx,
                test_subgop_.step[idx].disp_frame_idx)
          << "Error:display_index doesn't match";
      EXPECT_EQ(subgop_cfg->step[idx].type_code,
                test_subgop_.step[idx].type_code)
          << "Error:frame type doesn't match";
#if DEBUG
      printf("DisplayIdx=%d\tFrameType=%c\n",
             subgop_data_.step[idx].disp_frame_idx,
             test_subgop_.step[idx].type_code);
#endif
    }
    return 1;
  }

  bool ValidateSubgopCode() {
    EXPECT_EQ(subgop_code_test_, subgop_info_.code)
        << "Error:subgop code doesn't match";
#if DEBUG
    printf("Correct subgop code:%d selected\n", subgop_code_test_);
#endif
    return 1;
  }

  void ValidateSubgopConfig() {
    SubGOPCfg *subgop_cfg_test;
    subgop_cfg_test = DetermineSubgopConfig();
    if (subgop_cfg_test) {
      EXPECT_EQ(subgop_size_, subgop_cfg_test->num_frames)
          << "Error:subgop config selection wrong";
#if DEBUG
      printf("###User defined subgop config selected\n");
      printf(
          "###Following structure: Subgop_Code:%d, Size:%d, Steps:%d matched "
          "with Encoder\n",
          test_subgop_.subgop_in_gop_code, test_subgop_.num_frames,
          test_subgop_.num_steps);
#endif
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    if (AOM_CODEC_OK == res_dec) {
      aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
      if (subgop_info_.is_user_specified)
        AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_FRAME_INFO,
                                      &subgop_data_);
      if (frame_num_in_subgop_ == subgop_info_.size) {
        // Validation of user-specified sub-gop structure adoption in encoder
        // path. Validation of sub-gop structure propagation to decoder.
        if (subgop_info_.is_user_specified) {
          FillTestSubgopConfig();
          ValidateSubgopCode();
          ValidateSubgopConfig();
          ValidateSubgopFrametype();
        }
#if DEBUG
        else {
          printf("\n###Encoder defined subgop config selected\n");
          printf("###Following is the structure : Size:%d\n",
                 subgop_info_.size);
          ValidateSubgopConfig();
          if (frame_type_test_ == KEY_FRAME)
            printf("FrameNumber=%d\tFrameType=KEY_FRAME\n", frame_num_);
        }
#endif
        frames_from_key_ += subgop_info_.size;
        if (frame_type_test_ == KEY_FRAME) frames_from_key_ = 0;
        ResetSubgop();
      }
      frame_num_++;
    }
    return AOM_CODEC_OK == res_dec;
  }

  const SubgopTestParams subgop_test_params_;
  SubGOPSetCfg user_cfg_set_;
  SubGOPCfg test_subgop_;
  SUBGOP_INFO subgop_info_;
  SUBGOP_DATA subgop_data_;
  SUBGOP_IN_GOP_CODE subgop_code_test_;
  FRAME_TYPE frame_type_test_;
  aom_rc_mode rc_end_usage_;
  int subgop_size_;
  bool is_prev_frame_key_;
  int frames_from_key_;
  unsigned int frame_num_in_subgop_;
  unsigned int frame_num_;
  unsigned int enable_subgop_stats_;
};

TEST_P(SubGopTestLarge, SubGopTest) {
  libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                     cfg_.g_timebase.den, cfg_.g_timebase.num,
                                     0, 200);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
}

AV1_INSTANTIATE_TEST_SUITE(SubGopTestLarge, ::testing::ValuesIn(InputParams),
                           ::testing::Values(AOM_Q, AOM_VBR, AOM_CQ, AOM_CBR));

}  // namespace