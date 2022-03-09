/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AV1_RATECTRL_QMODE_H_
#define AOM_AV1_RATECTRL_QMODE_H_

#include "av1/ratectrl_qmode_interface.h"

namespace aom {

#define MAX_FIRSTPASS_ANALYSIS_FRAMES 150
typedef enum region_types {
  STABLE_REGION = 0,
  HIGH_VAR_REGION = 1,
  SCENECUT_REGION = 2,
  BLENDING_REGION = 3,
} REGION_TYPES;

typedef struct regions {
  int start;
  int last;
  double avg_noise_var;
  double avg_cor_coeff;
  double avg_sr_fr_ratio;
  double avg_intra_err;
  double avg_coded_err;
  REGION_TYPES type;
} REGIONS;

class AV1RateControlQMode : public AV1RateControlQModeInterface {
 public:
  AV1RateControlQMode() = default;
  ~AV1RateControlQMode() = default;
  void SetRcParam(const RateControlParam &rc_param) override;
  GopChunkList DetermineGopInfo(const FirstpassInfo &firstpass_info) override;
  virtual std::vector<FrameEncodeParameters> GetGopEncodeInfo(
      const TplGopStats &tpl_stats_list) override;

 private:
  RateControlParam rc_param_;
};
}  // namespace aom

#endif
