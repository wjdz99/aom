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
        rc_end_usage_(GET_PARAM(2)) {
    InitSubgop();
  }
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

  // check if subgop_config_str is a preset tag
  void GetSubGOPConfigStr() {
    int num_preset_configs = sizeof(subgop_config_str_preset_map) /
                             sizeof(*subgop_config_str_preset_map);
    for (int p = 0; p < num_preset_configs; ++p) {
      if (!strcmp(subgop_test_params_.subgop_str,
                  subgop_config_str_preset_map[p].preset_tag)) {
        subgop_test_params_.subgop_str =
            subgop_config_str_preset_map[p].preset_str;
        break;
      }
    }
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      encoder->Control(AOME_SET_CPUUSED, subgop_test_params_.cpu_used);
      encoder->Control(AV1E_ENABLE_SUBGOP_STATS, enable_subgop_stats_);
      GetSubGOPConfigStr();
      encoder->Control(AV1E_SET_SUBGOP_CONFIG_STR,
                       subgop_test_params_.subgop_str);
      av1_process_subgop_config_set(subgop_test_params_.subgop_str,
                                    &user_cfg_set_);
      encoder->Control(AV1E_SET_MIN_GF_INTERVAL,
                       subgop_test_params_.min_gf_interval);
      encoder->Control(AV1E_SET_MAX_GF_INTERVAL,
                       subgop_test_params_.max_gf_interval);
    }
  }

  virtual bool DoDecode() const { return 1; }

  virtual void PreDecodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Decoder *decoder) {
    aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
    if (video->frame() == 0)
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AV1D_ENABLE_SUBGOP_STATS,
                                    enable_subgop_stats_);
  }

  void InitSubgop() {
    memset(&user_cfg_set_, 0, sizeof(user_cfg_set_));
    subgop_data_.num_steps = MAX_SUBGOP_STATS_SIZE;
    ResetSubgop();
    is_first_frame_in_subgop_key_ = 0;
    frames_from_key_ = 0;
    frame_num_ = 0;
    // TODO(any): Extend this unit test for 'CONFIG_REALTIME_ONLY'
    enable_subgop_stats_ = 1;
    memset(&subgop_last_step_, 0, sizeof(subgop_last_step_));
  }

  void ResetSubgop() {
    subgop_info_.is_user_specified = 0;
    subgop_info_.frames_to_key = 0;
    subgop_info_.gf_interval = 0;
    subgop_info_.size = 0;
    subgop_info_.pos_code = SUBGOP_IN_GOP_GENERIC;

    for (int idx = 0; idx < MAX_SUBGOP_STATS_SIZE; idx++) {
      subgop_data_.step[idx].disp_frame_idx = -1;
      subgop_data_.step[idx].show_existing_frame = -1;
      subgop_data_.step[idx].show_frame = -1;
      subgop_data_.step[idx].is_filtered = -1;
      subgop_data_.step[idx].pyramid_level = 0;
      subgop_data_.step[idx].qindex = 0;
      subgop_data_.step[idx].refresh_frame_flags = 0;
      subgop_data_.step[idx].num_references = -1;
      memset(subgop_data_.step[idx].ref_frame_pyr_level, 0,
             sizeof(subgop_data_.step[idx].ref_frame_pyr_level));
      memset(subgop_data_.step[idx].ref_frame_disp_order, -1,
             sizeof(subgop_data_.step[idx].ref_frame_disp_order));
      memset(subgop_data_.step[idx].is_valid_ref_frame, 0,
             sizeof(subgop_data_.step[idx].is_valid_ref_frame));
      memset(subgop_data_.step[idx].ref_frame_map, 0,
             sizeof(subgop_data_.step[idx].ref_frame_map));
    }
    memset(display_order_test_, -1, sizeof(display_order_test_));
    subgop_data_.num_steps = 0;
    subgop_data_.step_idx_enc = 0;
    subgop_data_.step_idx_dec = 0;

    subgop_code_test_ = SUBGOP_IN_GOP_GENERIC;
    subgop_size_ = 0;
    frame_num_in_subgop_ = 0;
  }

  void DetermineSubgopCode(libaom_test::Encoder *encoder) {
    encoder->Control(AV1E_GET_FRAME_TYPE, &frame_type_test_);
    if (frame_type_test_ == KEY_FRAME) {
      is_first_frame_in_subgop_key_ = 1;
      return;
    }
    const int is_last_subgop =
        subgop_info_.frames_to_key <= (subgop_info_.gf_interval + 1);
    const int is_first_subgop = is_first_frame_in_subgop_key_;
    if (is_last_subgop)
      subgop_code_test_ = SUBGOP_IN_GOP_LAST;
    else if (is_first_subgop)
      subgop_code_test_ = SUBGOP_IN_GOP_FIRST;
    else
      subgop_code_test_ = SUBGOP_IN_GOP_GENERIC;
    subgop_size_ = subgop_info_.gf_interval;
  }

  virtual bool HandleEncodeResult(::libaom_test::VideoSource *video,
                                  libaom_test::Encoder *encoder) {
    (void)video;
    // Capturing the subgop info at the start of subgop.
    if (!frame_num_in_subgop_) {
      encoder->Control(AV1E_GET_SUB_GOP_CONFIG, &subgop_info_);
      DetermineSubgopCode(encoder);
      // Validation of user specified subgop structure adoption in encoder path.
      ValidateSubgopConfig();
      subgop_data_.num_steps = subgop_info_.num_steps;
    }
    if (subgop_info_.is_user_specified)
      encoder->Control(AV1E_GET_FRAME_INFO, &subgop_data_);
    return 1;
  }

  void FillTestSubgopConfig() {
    int filtered_frames[REF_FRAMES] = { 0 }, buf_idx = 0;
    if (frame_type_test_ == KEY_FRAME) return;

    subgop_cfg_test_.num_frames = subgop_info_.size;
    subgop_cfg_test_.num_steps = subgop_data_.num_steps;
    subgop_cfg_test_.subgop_in_gop_code = subgop_info_.pos_code;
    // Populating the filter-type of out-of-order frames appropriately for all
    // steps in sub-gop.
    for (int idx = 0; idx < subgop_data_.num_steps; idx++) {
      subgop_cfg_test_.step[idx].disp_frame_idx =
          subgop_data_.step[idx].disp_frame_idx - frames_from_key_;
      subgop_cfg_test_.step[idx].pyr_level =
          subgop_data_.step[idx].pyramid_level;
      subgop_cfg_test_.step[idx].num_references =
          subgop_data_.step[idx].num_references;
      for (int ref = 0; ref < INTER_REFS_PER_FRAME; ref++) {
        subgop_cfg_test_.step[idx].references[ref] =
            subgop_data_.step[idx].ref_frame_pyr_level[ref];
        display_order_test_[idx][ref] =
            subgop_data_.step[idx].ref_frame_disp_order[ref];
      }
      if (subgop_data_.step[idx].is_filtered) {
        filtered_frames[buf_idx++] =
            subgop_data_.step[idx].disp_frame_idx - frames_from_key_;
      } else {
        for (int ref_frame = 0; ref_frame < buf_idx; ref_frame++) {
          if (subgop_cfg_test_.step[idx].disp_frame_idx ==
              filtered_frames[ref_frame])
            subgop_data_.step[idx].is_filtered = 1;
        }
      }
    }
    // Calculating frame type code for all the steps in subgop.
    for (int idx = 0; idx < subgop_data_.num_steps; idx++) {
      FRAME_TYPE_CODE frame_type_code = FRAME_TYPE_INO_VISIBLE;
      int show_existing_frame = subgop_data_.step[idx].show_existing_frame;
      int show_frame = subgop_data_.step[idx].show_frame;
      int is_filtered = subgop_data_.step[idx].is_filtered;

      assert(show_existing_frame >= 0);
      assert(show_frame >= 0);
      assert(frame_type_code != 0);
      if (show_existing_frame == 0) {
        if (show_frame == 0)
          frame_type_code = (is_filtered == 1) ? FRAME_TYPE_OOO_FILTERED
                                               : FRAME_TYPE_OOO_UNFILTERED;
        else if (show_frame == 1)
          frame_type_code = (is_filtered == 1) ? FRAME_TYPE_INO_REPEAT
                                               : FRAME_TYPE_INO_VISIBLE;
      } else if (show_existing_frame == 1) {
        frame_type_code = (is_filtered == 1) ? FRAME_TYPE_INO_REPEAT
                                             : FRAME_TYPE_INO_SHOWEXISTING;
      }
      subgop_cfg_test_.step[idx].type_code = frame_type_code;
    }
  }

  SubGOPCfg *DetermineSubgopConfig() {
    SubGOPCfg *subgop_cfg = user_cfg_set_.config;
    SubGOPCfg *viable_subgop_cfg[MAX_SUBGOP_CODES] = { NULL };
    unsigned int cfg_count = 0;

    for (int idx = 0; idx < user_cfg_set_.num_configs; idx++) {
      if (subgop_cfg[idx].num_frames == subgop_size_) {
        if (subgop_cfg[idx].subgop_in_gop_code == subgop_code_test_)
          return &subgop_cfg[idx];
        else
          viable_subgop_cfg[cfg_count++] = &subgop_cfg[idx];
        assert(cfg_count < MAX_SUBGOP_CODES);
      }
    }
    subgop_code_test_ = SUBGOP_IN_GOP_GENERIC;
    for (unsigned int cfg = 0; cfg < cfg_count; cfg++) {
      if (viable_subgop_cfg[cfg]->subgop_in_gop_code == subgop_code_test_)
        return viable_subgop_cfg[cfg];
    }
    return NULL;
  }

  // Validates frametype(along with temporal filtering), frame coding order
  bool ValidateSubgopFrametype() {
    for (int idx = 0; idx < subgop_cfg_ref_->num_steps; idx++) {
      EXPECT_EQ(subgop_cfg_ref_->step[idx].disp_frame_idx,
                subgop_cfg_test_.step[idx].disp_frame_idx)
          << "Error:display_index doesn't match";
      EXPECT_EQ(subgop_cfg_ref_->step[idx].type_code,
                subgop_cfg_test_.step[idx].type_code)
          << "Error:frame type doesn't match";
    }
    return 1;
  }

  // Validates Pyramid level with user config
  void ValidatePyramidLevel() {
    int max_pyramid_level = 0;
    for (int idx = 0; idx < subgop_cfg_ref_->num_steps; idx++) {
      if (max_pyramid_level < subgop_cfg_ref_->step[idx].pyr_level)
        max_pyramid_level = subgop_cfg_ref_->step[idx].pyr_level;
    }
    for (int idx = 0; idx < subgop_cfg_ref_->num_steps; idx++) {
      int ref_pyramid_level =
          (subgop_cfg_ref_->step[idx].pyr_level == max_pyramid_level)
              ? 6
              : subgop_cfg_ref_->step[idx].pyr_level;
      EXPECT_EQ(subgop_cfg_test_.step[idx].pyr_level, ref_pyramid_level)
          << "Error:pyramid level doesn't match";
    }
  }

  // Validates Pyramid level along with qindex assignment
  void ValidatePyramidLevelQIndex() {
    int level_qindex[REF_FRAMES] = { 0 };
    int pyramid_level;
    for (int idx = 0; idx < subgop_cfg_ref_->num_steps; idx++) {
      pyramid_level = subgop_cfg_test_.step[idx].pyr_level;
      if (!level_qindex[pyramid_level])
        level_qindex[pyramid_level] = subgop_data_.step[idx].qindex;
      else if (!subgop_data_.step[idx].show_existing_frame &&
               !subgop_data_.step[idx].is_filtered)
        EXPECT_EQ(level_qindex[pyramid_level], subgop_data_.step[idx].qindex)
            << "Error:qindex in a pyramid level doesn't match";
    }
    pyramid_level = 1;
    for (int idx = 1; idx < REF_FRAMES; idx++) {
      if (level_qindex[pyramid_level]) {
        EXPECT_LT(level_qindex[pyramid_level - 1], level_qindex[pyramid_level])
            << "Error:qindex should be higher in hierarchical pyramid level";
        pyramid_level++;
      }
    }
  }

  // Validates reference buffer refresh
  void ValidateRefBufRefresh() {
    int start_idx = 0;
    SubGOPStepData *prev_step_data = &subgop_last_step_;
    if (is_first_frame_in_subgop_key_) {
      start_idx = 1;
      prev_step_data = &subgop_data_.step[0];
    }

    for (int idx = start_idx; idx < subgop_cfg_ref_->num_steps; idx++) {
      SubGOPStepData *curr_step_data = &subgop_data_.step[idx];
      int ref_count = 0;
      int refresh_frame_flags = curr_step_data->refresh_frame_flags;
      // Validates user-defined refresh_flag with decoder
      if (subgop_cfg_ref_->step[idx].refresh != -1 &&
          !curr_step_data->show_existing_frame)
        EXPECT_EQ(subgop_cfg_ref_->step[idx].refresh, refresh_frame_flags)
            << "Error: refresh flag mismatch";
      // Validates reference picture management w.r.t refresh_flags
      if (refresh_frame_flags && !curr_step_data->show_existing_frame) {
        for (int mask = refresh_frame_flags; mask; mask >>= 1) {
          if (mask & 1)
            EXPECT_EQ(curr_step_data->disp_frame_idx,
                      (int)curr_step_data->ref_frame_map[ref_count])
                << "Error: reference buffer refresh failed";
          else
            EXPECT_EQ(prev_step_data->ref_frame_map[ref_count],
                      curr_step_data->ref_frame_map[ref_count])
                << "Error: reference buffer refresh failed";
          assert(ref_count < REF_FRAMES);
          ref_count++;
        }
      }

      for (int ref_idx = ref_count; ref_idx < REF_FRAMES; ref_idx++)
        EXPECT_EQ(prev_step_data->ref_frame_map[ref_idx],
                  curr_step_data->ref_frame_map[ref_idx])
            << "Error: reference buffer refresh failed";
      prev_step_data = curr_step_data;
    }
  }

  void ValidateRefFrames() {
    int start_idx = is_first_frame_in_subgop_key_;
    for (int idx = start_idx; idx < subgop_cfg_ref_->num_steps; idx++) {
      unsigned int *ref_frame_map =
          (idx > 0) ? subgop_data_.step[idx - 1].ref_frame_map
                    : subgop_last_step_.ref_frame_map;
      if (!subgop_data_.step[idx].show_existing_frame)
        EXPECT_EQ(subgop_cfg_ref_->step[idx].num_references,
                  subgop_cfg_test_.step[idx].num_references)
            << "Error:Reference frames count doesn't match";
      // To validate the count of ref_frames and their management with user
      // config.
      for (int ref = 0; ref < subgop_cfg_test_.step[idx].num_references;
           ref++) {
        if (subgop_data_.step[idx].is_valid_ref_frame[ref]) {
          EXPECT_EQ(subgop_cfg_ref_->step[idx].references[ref],
                    subgop_cfg_test_.step[idx].references[ref])
              << "Error:Reference frame level doesn't match";
          for (int buf = 0; buf < REF_FRAMES; buf++) {
            if (display_order_test_[idx][ref] == (int8_t)ref_frame_map[buf])
              break;
            EXPECT_NE(buf, REF_FRAMES - 1)
                << "Error:Ref frame isn't part of ref_picture_buf";
          }
        }
      }
    }
  }

  bool IsInputSubgopCfgUsed() {
    int num_ooo_frames_ref = 0;
    int num_ooo_frames_test = 0;
    // Encoder may choose to opt out input sub-gop config when use_altref is 0
    // for the given sub-gop
    for (int idx = 0; idx < subgop_cfg_ref_->num_steps; idx++) {
      // Count number out-of-order frames in reference config
      num_ooo_frames_ref +=
          subgop_cfg_ref_->step[idx].type_code == FRAME_TYPE_OOO_FILTERED;
      num_ooo_frames_ref +=
          subgop_cfg_ref_->step[idx].type_code == FRAME_TYPE_OOO_UNFILTERED;
      // Count number out-of-order frames in test config
      num_ooo_frames_test +=
          subgop_cfg_test_.step[idx].type_code == FRAME_TYPE_OOO_FILTERED;
      num_ooo_frames_test +=
          subgop_cfg_test_.step[idx].type_code == FRAME_TYPE_OOO_UNFILTERED;
    }
    return num_ooo_frames_ref == num_ooo_frames_test;
  }

  void ValidateSubgopConfig() {
    if (frame_type_test_ == KEY_FRAME) return;
    subgop_cfg_ref_ = DetermineSubgopConfig();
    if (subgop_cfg_ref_) {
      EXPECT_EQ(subgop_size_, subgop_cfg_ref_->num_frames)
          << "Error:subgop config selection wrong";
      subgop_info_.is_user_specified = 1;
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    if (AOM_CODEC_OK != res_dec) return 0;
    aom_codec_ctx_t *ctx_dec = decoder->GetDecoder();
    frame_num_in_subgop_++;
    int is_last_frame_in_subgop = (frame_num_in_subgop_ == subgop_info_.size);

    if (subgop_info_.is_user_specified ||
        is_last_frame_in_subgop)  // To collect last step info of subgop in
                                  // encoder defined config
      AOM_CODEC_CONTROL_TYPECHECKED(ctx_dec, AOMD_GET_FRAME_INFO,
                                    &subgop_data_);
    if (is_last_frame_in_subgop) {
      // Validation of sub-gop structure propagation to decoder.
      if (subgop_info_.is_user_specified) {
        FillTestSubgopConfig();
        if ((subgop_info_.is_user_specified = IsInputSubgopCfgUsed())) {
          EXPECT_EQ(subgop_code_test_, subgop_info_.pos_code)
              << "Error:subgop code doesn't match";
          ValidateSubgopFrametype();
          ValidatePyramidLevel();
          if (rc_end_usage_ == AOM_Q) ValidatePyramidLevelQIndex();
          ValidateRefBufRefresh();
          ValidateRefFrames();
        }
      }
      frames_from_key_ += subgop_info_.size;
      if (frame_type_test_ == KEY_FRAME) {
        frames_from_key_ = 0;
      } else {
        is_first_frame_in_subgop_key_ = 0;
        // To collect last step info of subgop.
        assert(subgop_data_.step_idx_dec >= 0);
        memcpy(&subgop_last_step_,
               &subgop_data_.step[subgop_data_.step_idx_dec - 1],
               sizeof(subgop_last_step_));
      }
      ResetSubgop();
    }
    frame_num_++;
    return AOM_CODEC_OK == res_dec;
  }

  SubgopTestParams subgop_test_params_;
  SubGOPSetCfg user_cfg_set_;
  SubGOPCfg subgop_cfg_test_;
  SubGOPCfg *subgop_cfg_ref_;
  SubGOPInfo subgop_info_;
  SubGOPData subgop_data_;
  SubGOPStepData subgop_last_step_;
  SUBGOP_IN_GOP_CODE subgop_code_test_;
  FRAME_TYPE frame_type_test_;
  aom_rc_mode rc_end_usage_;
  int8_t display_order_test_[MAX_SUBGOP_SIZE][REF_FRAMES];
  int subgop_size_;
  bool is_first_frame_in_subgop_key_;
  int frames_from_key_;
  unsigned int frame_num_in_subgop_;
  unsigned int frame_num_;
  unsigned int enable_subgop_stats_;
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
