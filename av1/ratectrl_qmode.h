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

#include <vector>
#include "av1/encoder/firstpass.h"

struct FirstPassInfo {
  FIRSTPASS_STATS *stats_list;
  int frame_count;
};

struct TplDepStats {
  int64_t intra_cost;
  int64_t inter_cost;
  int64_t srcrf_dist;
  int64_t recrf_dist;
  int64_t cmp_recrf_dist[2];
  int64_t srcrf_rate;
  int64_t recrf_rate;
  int64_t srcrf_sse;
  int64_t cmp_recrf_rate[2];
  int64_t mc_dep_rate;
  int64_t mc_dep_dist;
  int_mv mv[INTER_REFS_PER_FRAME];
  int ref_frame_index[2];
  int64_t pred_error[INTER_REFS_PER_FRAME];
};

struct GopInfo {
  GF_GROUP *gop_list;
  int gop_count;
};

struct TplInfo {
  TplDepStats *tpl_stats_list;
  int tpl_stats_count;
};

struct GopEncodeParameters {
  std::vector<int> q_index_list;
  std::vector<int> rdmult_list;
  int gop_frame_count;
};

struct GopEncodeInfo {
  std::vector<GopEncodeParameters> gop_encode_parameters_list;
  int gop_count;
};

namespace aom {
class AV1RateControlQMode {
 public:
  AV1RateControlQMode() = default;
  ~AV1RateControlQMode() {}

  void SetBaseQ(int base_q) { base_q = base_q; }

  virtual void SetFirstPassInfo(const FirstPassInfo &first_pass_info) = 0;
  virtual void SetTplInfo(const GopInfo &gop_info, const TplInfo &tpl_info) = 0;
  virtual void DetermineGopInfo() = 0;
  virtual GopEncodeInfo *GetGopEncodeInfo() = 0;

 private:
  int base_q_;
  FirstPassInfo first_pass_info_;
  TplInfo tpl_info_;
  GopEncodeInfo gop_encode_info_;
};  // class AV1RateCtrlQMode
}  // namespace aom

#endif  // AOM_AV1_RATECTRL_QMODE_H_
