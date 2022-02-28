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

namespace aom {

struct FirstPassInfo {
  std::vector<FIRSTPASS_STATS> stats_list;
  int frame_count;
};

struct TplDepStats {
  BLOCK_SIZE block_size;
  int row;
  int col;
  int64_t intra_cost;
  int64_t inter_cost;
  int_mv mv[2];
  int ref_frame_index[2];
};

struct GopInfo {
  std::vector<GF_GROUP> gop_list;
  int gop_count;
};

struct TplInfo {
  std::vector<TplDepStats> tpl_stats_list;
  int tpl_stats_count;
};

struct FrameEncodeParameters {
  int q_index;
  int rdmult;
};

class AV1RateControlQMode {
 public:
  AV1RateControlQMode() = default;
  ~AV1RateControlQMode() {}

  void SetBaseQ(int base_q) { base_q = base_q; }

  virtual void SetFirstPassInfo(const FirstPassInfo &first_pass_info) = 0;
  virtual void SetTplInfo(const TplInfo &tpl_info) = 0;
  virtual void DetermineGopInfo() = 0;
  virtual std::vector<FrameEncodeParameters> GetGopEncodeInfo() = 0;

 private:
  int base_q_;
  FirstPassInfo first_pass_info_;
  TplInfo tpl_info_;
};  // class AV1RateCtrlQMode
}  // namespace aom

#endif  // AOM_AV1_RATECTRL_QMODE_H_
