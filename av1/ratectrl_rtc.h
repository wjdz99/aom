/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_RATECTRL_RTC_H_
#define AOM_AV1_RATECTRL_RTC_H_

#ifdef __cplusplus
#include <cstddef>
#include <cstdint>
#include <memory>
#else
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#endif  // __cplusplus

struct AV1_COMP;

struct AomAV1LoopfilterLevel {
  int filter_level[2];
  int filter_level_u;
  int filter_level_v;
};

struct AomAV1CdefInfo {
  int cdef_strength_y;
  int cdef_strength_uv;
  int damping;
};

struct AomAV1SegmentationData {
  const uint8_t *segmentation_map;
  size_t segmentation_map_size;
  const int *delta_q;
  size_t delta_q_size;
};

typedef enum AomFrameType { kAomKeyFrame, kAomInterFrame } AomFrameType;

struct AomAV1FrameParamsRTC {
  AomFrameType frame_type;
  int spatial_layer_id;
  int temporal_layer_id;
};

typedef enum {
  kAomFrameDropDecisionOk,    // Frame is encoded.
  kAomFrameDropDecisionDrop,  // Frame is dropped.
} AomFrameDropDecision;

// These constants come from AV1 spec.
#define kAomMaxLayers 32
#define kAomMaxTemporalLayers 8
#define kAomMaxSpatialLayers 4

struct AomAV1RateControlRtcConfig {
#ifdef __cplusplus
  AomAV1RateControlRtcConfig();
#endif

  int width;
  int height;
  // Flag indicating if the content is screen or not.
  bool is_screen;
  // 0-63
  int max_quantizer;
  int min_quantizer;
  int64_t target_bandwidth;
  int64_t buf_initial_sz;
  int64_t buf_optimal_sz;
  int64_t buf_sz;
  int undershoot_pct;
  int overshoot_pct;
  int max_intra_bitrate_pct;
  int max_inter_bitrate_pct;
  int frame_drop_thresh;
  int max_consec_drop_ms;
  double framerate;
  int layer_target_bitrate[kAomMaxLayers];
  int ts_rate_decimator[kAomMaxTemporalLayers];
  int aq_mode;
  // Number of spatial layers
  int ss_number_layers;
  // Number of temporal layers
  int ts_number_layers;
  int max_quantizers[kAomMaxLayers];
  int min_quantizers[kAomMaxLayers];
  int scaling_factor_num[kAomMaxSpatialLayers];
  int scaling_factor_den[kAomMaxSpatialLayers];
};

#ifdef __cplusplus
namespace aom {

using AV1LoopfilterLevel = AomAV1LoopfilterLevel;
using AV1CdefInfo = AomAV1CdefInfo;
using AV1SegmentationData = AomAV1SegmentationData;
using AV1FrameParamsRTC = AomAV1FrameParamsRTC;
using AV1RateControlRtcConfig = AomAV1RateControlRtcConfig;

using FrameType = AomFrameType;
constexpr FrameType kKeyFrame = kAomKeyFrame;
constexpr FrameType kInterFrame = kAomInterFrame;

using FrameDropDecision = AomFrameDropDecision;
constexpr FrameDropDecision kFrameDropDecisionOk = kAomFrameDropDecisionOk;
constexpr FrameDropDecision kFrameDropDecisionDrop = kAomFrameDropDecisionDrop;

class AV1RateControlRTC {
 public:
  static std::unique_ptr<AV1RateControlRTC> Create(
      const AV1RateControlRtcConfig &cfg);
  ~AV1RateControlRTC();

  bool UpdateRateControl(const AV1RateControlRtcConfig &rc_cfg);
  // GetQP() needs to be called after ComputeQP() to get the latest QP
  int GetQP() const;
  // GetLoopfilterLevel() needs to be called after ComputeQP()
  AV1LoopfilterLevel GetLoopfilterLevel() const;
  // GetCdefInfo() needs to be called after ComputeQP()
  AV1CdefInfo GetCdefInfo() const;
  // Returns the segmentation map used for cyclic refresh, based on 4x4 blocks.
  bool GetSegmentationData(AV1SegmentationData *segmentation_data) const;
  // ComputeQP returns the QP if the frame is not dropped (kOk return),
  // otherwise it returns kDrop and subsequent GetQP and PostEncodeUpdate
  // are not to be called (av1_rc_postencode_update_drop_frame is already
  // called via ComputeQP if drop is decided).
  FrameDropDecision ComputeQP(const AV1FrameParamsRTC &frame_params);
  // Feedback to rate control with the size of current encoded frame
  void PostEncodeUpdate(uint64_t encoded_frame_size);

 private:
  AV1RateControlRTC() = default;
  bool InitRateControl(const AV1RateControlRtcConfig &cfg);
  AV1_COMP *cpi_;
  int initial_width_;
  int initial_height_;
};
}  // namespace aom
#endif  // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif
void *av1_ratecontrol_rtc_create(
    const struct AomAV1RateControlRtcConfig *rc_cfg);
void av1_ratecontrol_rtc_destroy(void *controller);
bool av1_ratecontrol_rtc_update(
    void *controller, const struct AomAV1RateControlRtcConfig *rc_cfg);
int av1_ratecontrol_rtc_get_qp(void *controller);

struct AomAV1LoopfilterLevel av1_ratecontrol_rtc_get_loop_filter_level(
    void *controller);
AomFrameDropDecision av1_ratecontrol_rtc_compute_qp(
    void *controller, const struct AomAV1FrameParamsRTC *frame_params);

void av1_ratecontrol_rtc_post_encode_update(void *controller,
                                            uint64_t encoded_frame_size);

bool av1_ratecontrol_rtc_get_segmentation(
    void *controller, struct AomAV1SegmentationData *segmentation_data);

struct AomAV1CdefInfo av1_ratecontrol_rtc_get_cdef_info(void *controller);

void av1_ratecontrol_rtc_init_ratecontrol_config(
    struct AomAV1RateControlRtcConfig *config);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_RATECTRL_RTC_H_
