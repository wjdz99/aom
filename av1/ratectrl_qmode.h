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

struct RateControlParam {
  int max_gop_length;
  int min_gop_length;
  int max_ref_frames;
  int base_q_index;
};

struct TplBlockStats {
  BLOCK_SIZE block_size;
  int row;
  int col;
  int64_t intra_cost;
  int64_t inter_cost;
  int_mv mv[2];
  int ref_frame_index[2];
};

typedef std::vector<GF_GROUP> GopChunkList;

struct FrameEncodeParameters {
  int q_index;
  int rdmult;
};

typedef std::vector<FIRSTPASS_STATS> FirstpassInfo;
typedef std::vector<TplBlockStats> TplFrameStats;
typedef std::vector<TplFrameStats> TplGopStats;

class AV1RateControlQMode {
 public:
  AV1RateControlQMode() = default;
  ~AV1RateControlQMode() {}

  virtual void SetBaseQ(RateControlParam &rc_param) = 0;
  // Accept firstpass and tpl info from the encoder, call the following private
  // functions and return q index and rdmult.
  virtual std::vector<FrameEncodeParameters> GetGopEncodeInfo(
      const FirstpassInfo &firstpass_stats_list,
      const TplGopStats &tpl_stats_list) = 0;

 private:
  virtual void SetFirstPassInfo() = 0;
  virtual void SetTplInfo() = 0;
  virtual void DetermineGopInfo() = 0;
};  // class AV1RateCtrlQMode
}  // namespace aom

#endif  // AOM_AV1_RATECTRL_QMODE_H_
