/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
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

#include <cstdint>
#include <memory>

#include "aom/aomcx_defines.h"
#include "aom_ports/mem_bytes.h"
#include "av1/common/blockd_defines.h"

struct AV1_COMP;

namespace aom {

struct AV1RateControlRtcConfig {
 public:
  AV1RateControlRtcConfig();

  int width;
  int height;
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
  double framerate;
  int layer_target_bitrate[AOM_MAX_LAYERS];
  int ts_rate_decimator[AOM_MAX_TS_LAYERS];
  int aq_mode;
  // Number of spatial layers
  int ss_number_layers;
  // Number of temporal layers
  int ts_number_layers;
  int max_quantizers[AOM_MAX_LAYERS];
  int min_quantizers[AOM_MAX_LAYERS];
  int scaling_factor_num[AOM_MAX_SS_LAYERS];
  int scaling_factor_den[AOM_MAX_SS_LAYERS];
};

struct AV1FrameParamsRTC {
  FRAME_TYPE frame_type;
  int spatial_layer_id;
  int temporal_layer_id;
};

class AV1RateControlRTC {
 public:
  static std::unique_ptr<AV1RateControlRTC> Create(
      const AV1RateControlRtcConfig &cfg);
  ~AV1RateControlRTC();

  void UpdateRateControl(const AV1RateControlRtcConfig &rc_cfg);
  // GetQP() needs to be called after ComputeQP() to get the latest QP
  int GetQP() const;
  signed char *GetCyclicRefreshMap() const;
  int *GetDeltaQ() const;
  void ComputeQP(const AV1FrameParamsRTC &frame_params);
  // Feedback to rate control with the size of current encoded frame
  void PostEncodeUpdate(uint64_t encoded_frame_size);

 private:
  AV1RateControlRTC() = default;
  void InitRateControl(const AV1RateControlRtcConfig &cfg);
  AV1_COMP *cpi_;
  int initial_width_;
  int initial_height_;
};

}  // namespace aom

#endif  // AOM_AV1_RATECTRL_RTC_H_
