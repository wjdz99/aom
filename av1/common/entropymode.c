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

#include "aom_mem/aom_mem.h"

#include "av1/common/reconinter.h"
#include "av1/common/scan.h"
#include "av1/common/onyxc_int.h"
#include "av1/common/seg_common.h"
#include "av1/common/txb_common.h"

static const aom_cdf_prob default_newmv_cdf[NEWMV_MODE_CONTEXTS][CDF_SIZE(2)] =
    {{AOM_CDF2(22867)}, {AOM_CDF2(16903)}, {AOM_CDF2(16348)}, {AOM_CDF2(9064)},
     {AOM_CDF2(13406)}, {AOM_CDF2(16384)}, {AOM_CDF2(5797)}};

static const aom_cdf_prob default_zeromv_cdf[GLOBALMV_MODE_CONTEXTS][CDF_SIZE(
    2)] = {{AOM_CDF2(5102)}, {AOM_CDF2(2722)}};

static const aom_cdf_prob default_refmv_cdf[REFMV_MODE_CONTEXTS][CDF_SIZE(2)] =
    {{AOM_CDF2(23625)}, {AOM_CDF2(23874)}, {AOM_CDF2(19037)},
     {AOM_CDF2(28714)}, {AOM_CDF2(25015)}, {AOM_CDF2(20098)},
     {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}};

static const aom_cdf_prob default_drl_cdf[DRL_MODE_CONTEXTS][CDF_SIZE(2)] = {
    {AOM_CDF2(12855)}, {AOM_CDF2(24603)}, {AOM_CDF2(17531)}};

static const aom_cdf_prob
    default_inter_compound_mode_cdf[INTER_MODE_CONTEXTS][CDF_SIZE(
        INTER_COMPOUND_MODES)] = {
        {AOM_CDF8(6157, 10180, 11778, 13207, 14565, 15826, 26575)},
        {AOM_CDF8(7899, 14050, 15695, 17210, 18638, 19788, 28015)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(8431, 15199, 16702, 18090, 19700, 20974, 27587)},
        {AOM_CDF8(11229, 14980, 18483, 21906, 22942, 23754, 25633)},
        {AOM_CDF8(15078, 18573, 20518, 22597, 23083, 23446, 30444)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)},
        {AOM_CDF8(8770, 14754, 17595, 20389, 22041, 23250, 25253)},
        {AOM_CDF8(14772, 19854, 21533, 23300, 24067, 24621, 30182)},
        {AOM_CDF8(11750, 20202, 21647, 23235, 24808, 25969, 28957)}};

static const aom_cdf_prob default_compound_type_cdf[BLOCK_SIZES_ALL][CDF_SIZE(
    COMPOUND_TYPES - 1)] = {
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(21407)},
    {AOM_CDF2(12047)}, {AOM_CDF2(9485)},  {AOM_CDF2(8726)},  {AOM_CDF2(8780)},
    {AOM_CDF2(8195)},  {AOM_CDF2(6716)},  {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(11122)}, {AOM_CDF2(6407)},
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}};

#if WEDGE_IDX_ENTROPY_CODING
static const aom_cdf_prob default_wedge_idx_cdf[BLOCK_SIZES_ALL][CDF_SIZE(16)] =
    {{AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(1736, 3962, 6502, 8300, 10398, 11946, 15439, 18225, 20017,
                21534, 23013, 24359, 26693, 28657, 30864)},
     {AOM_CDF16(650, 3050, 5810, 6451, 6851, 7002, 7377, 14449, 15933, 16873,
                17786, 18683, 22189, 25514, 29450)},
     {AOM_CDF16(2449, 3617, 4773, 7171, 7981, 8267, 9026, 14519, 17756, 20717,
                23732, 26545, 28665, 29982, 31429)},
     {AOM_CDF16(1353, 3526, 6005, 7431, 9479, 11128, 14422, 17233, 18942, 20476,
                22097, 23607, 25903, 27999, 30581)},
     {AOM_CDF16(990, 3387, 6406, 7385, 8070, 8266, 8879, 13990, 15606, 16805,
                18118, 19281, 22200, 25546, 29279)},
     {AOM_CDF16(2115, 3798, 5559, 7816, 9079, 9559, 10644, 13985, 16743, 19455,
                22347, 24793, 26834, 28711, 31078)},
     {AOM_CDF16(1407, 3490, 5974, 7662, 9589, 11972, 15676, 19022, 20367, 21817,
                23436, 24861, 26512, 28421, 30805)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(158, 1033, 2223, 2337, 2373, 2396, 2438, 21914, 22629, 23230,
                23931, 24724, 26382, 28141, 30501)},
     {AOM_CDF16(871, 1059, 1247, 2190, 2259, 2287, 2341, 22605, 24133, 25869,
                27642, 29373, 30348, 31117, 32029)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)},
     {AOM_CDF16(2048, 4096, 6144, 8192, 10240, 12288, 14336, 16384, 18432,
                20480, 22528, 24576, 26624, 28672, 30720)}};
#endif

static const aom_cdf_prob default_interintra_cdf[BLOCK_SIZE_GROUPS][CDF_SIZE(
    2)] = {
    {AOM_CDF2(16384)}, {AOM_CDF2(25985)}, {AOM_CDF2(28575)}, {AOM_CDF2(30709)}};

static const aom_cdf_prob
    default_interintra_mode_cdf[BLOCK_SIZE_GROUPS][CDF_SIZE(INTERINTRA_MODES)] =
        {{AOM_CDF4(8192, 16384, 24576)},
         {AOM_CDF4(2645, 13446, 26605)},
         {AOM_CDF4(4800, 14102, 25659)},
         {AOM_CDF4(6909, 12661, 24547)}};

static const aom_cdf_prob
    default_wedge_interintra_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)] = {
        {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
        {AOM_CDF2(19447)}, {AOM_CDF2(24708)}, {AOM_CDF2(25314)},
        {AOM_CDF2(25855)}, {AOM_CDF2(27604)}, {AOM_CDF2(28241)},
        {AOM_CDF2(27242)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
        {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
        {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
        {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
        {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}};

static const aom_cdf_prob default_motion_mode_cdf[BLOCK_SIZES_ALL][CDF_SIZE(
    MOTION_MODES)] = {{AOM_CDF3(10923, 21845)}, {AOM_CDF3(10923, 21845)},
                      {AOM_CDF3(10923, 21845)}, {AOM_CDF3(3415, 21646)},
                      {AOM_CDF3(3366, 22137)},  {AOM_CDF3(3774, 23118)},
                      {AOM_CDF3(15626, 24001)}, {AOM_CDF3(3885, 21057)},
                      {AOM_CDF3(11291, 22736)}, {AOM_CDF3(24541, 27336)},
                      {AOM_CDF3(22883, 26905)}, {AOM_CDF3(23841, 27627)},
                      {AOM_CDF3(27434, 29899)}, {AOM_CDF3(30337, 31293)},
                      {AOM_CDF3(31525, 32009)}, {AOM_CDF3(32395, 32523)},
                      {AOM_CDF3(10923, 21845)}, {AOM_CDF3(10923, 21845)},
                      {AOM_CDF3(27289, 30883)}, {AOM_CDF3(27770, 31192)},
                      {AOM_CDF3(29523, 31392)}, {AOM_CDF3(30903, 31907)},
                      {AOM_CDF3(10923, 21845)}, {AOM_CDF3(10923, 21845)}};

static const aom_cdf_prob default_obmc_cdf[BLOCK_SIZES_ALL][CDF_SIZE(2)] = {
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(5685)},
    {AOM_CDF2(8289)},  {AOM_CDF2(7864)},  {AOM_CDF2(12554)}, {AOM_CDF2(15690)},
    {AOM_CDF2(16679)}, {AOM_CDF2(25400)}, {AOM_CDF2(24702)}, {AOM_CDF2(22966)},
    {AOM_CDF2(24310)}, {AOM_CDF2(31079)}, {AOM_CDF2(31090)}, {AOM_CDF2(31528)},
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(20613)}, {AOM_CDF2(19892)},
    {AOM_CDF2(24920)}, {AOM_CDF2(25293)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}};

static const aom_cdf_prob default_delta_q_cdf[CDF_SIZE(DELTA_Q_PROBS + 1)] = {
  AOM_CDF4(28160, 32120, 32677)
};
#if CONFIG_EXT_DELTA_Q
static const aom_cdf_prob default_delta_lf_multi_cdf[FRAME_LF_COUNT][CDF_SIZE(
    DELTA_LF_PROBS + 1)] = { { AOM_CDF4(28160, 32120, 32677) },
                             { AOM_CDF4(28160, 32120, 32677) },
                             { AOM_CDF4(28160, 32120, 32677) },
                             { AOM_CDF4(28160, 32120, 32677) } };
static const aom_cdf_prob default_delta_lf_cdf[CDF_SIZE(DELTA_LF_PROBS + 1)] = {
  AOM_CDF4(28160, 32120, 32677)
};
#endif

static const aom_cdf_prob default_intra_inter_cdf[INTRA_INTER_CONTEXTS]
                                                 [CDF_SIZE(2)] = {
                                                     {AOM_CDF2(898)},
                                                     {AOM_CDF2(16561)},
                                                     {AOM_CDF2(17968)},
                                                     {AOM_CDF2(25223)}};

static const aom_cdf_prob default_comp_inter_cdf[COMP_INTER_CONTEXTS][CDF_SIZE(
    2)] = {{AOM_CDF2(25955)},
           {AOM_CDF2(23417)},
           {AOM_CDF2(11827)},
           {AOM_CDF2(10281)},
           {AOM_CDF2(2791)}};

static const aom_cdf_prob default_comp_ref_type_cdf[COMP_REF_TYPE_CONTEXTS]
                                                   [CDF_SIZE(2)] = {
                                                       {AOM_CDF2(375)},
                                                       {AOM_CDF2(766)},
                                                       {AOM_CDF2(7274)},
                                                       {AOM_CDF2(4489)},
                                                       {AOM_CDF2(24144)}};
static const aom_cdf_prob
    default_uni_comp_ref_cdf[UNI_COMP_REF_CONTEXTS][UNIDIR_COMP_REFS -
                                                    1][CDF_SIZE(2)] = {
        {{AOM_CDF2(1273)}, {AOM_CDF2(4)}, {AOM_CDF2(4)}},
        {{AOM_CDF2(8734)}, {AOM_CDF2(4)}, {AOM_CDF2(4)}},
        {{AOM_CDF2(30435)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}}};

static const aom_cdf_prob
    default_comp_ref_cdf[REF_CONTEXTS][FWD_REFS - 1][CDF_SIZE(2)] = {
        {{AOM_CDF2(4522)}, {AOM_CDF2(11366)}, {AOM_CDF2(551)}},
        {{AOM_CDF2(19041)}, {AOM_CDF2(24905)}, {AOM_CDF2(10977)}},
        {{AOM_CDF2(30769)}, {AOM_CDF2(31308)}, {AOM_CDF2(28293)}}};

static const aom_cdf_prob
    default_comp_bwdref_cdf[REF_CONTEXTS][BWD_REFS - 1][CDF_SIZE(2)] = {
        {{AOM_CDF2(1277)}, {AOM_CDF2(674)}},
        {{AOM_CDF2(16888)}, {AOM_CDF2(15334)}},
        {{AOM_CDF2(30932)}, {AOM_CDF2(31259)}}};

static const aom_cdf_prob default_single_ref_cdf[REF_CONTEXTS][SINGLE_REFS - 1]
                                                [CDF_SIZE(2)] = {
                                                    {{AOM_CDF2(5987)},
                                                     {AOM_CDF2(2650)},
                                                     {AOM_CDF2(4432)},
                                                     {AOM_CDF2(10344)},
                                                     {AOM_CDF2(167)},
                                                     {AOM_CDF2(1927)}},
                                                    {{AOM_CDF2(17478)},
                                                     {AOM_CDF2(19202)},
                                                     {AOM_CDF2(18636)},
                                                     {AOM_CDF2(27724)},
                                                     {AOM_CDF2(5534)},
                                                     {AOM_CDF2(15369)}},
                                                    {{AOM_CDF2(29199)},
                                                     {AOM_CDF2(30660)},
                                                     {AOM_CDF2(30719)},
                                                     {AOM_CDF2(31843)},
                                                     {AOM_CDF2(29092)},
                                                     {AOM_CDF2(29893)}}};

// TODO(huisu): tune these cdfs
const aom_cdf_prob
    default_palette_y_size_cdf[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)] = {
      { AOM_CDF7(12288, 19408, 24627, 26662, 28499, 30667) },
      { AOM_CDF7(12288, 19408, 24627, 26662, 28499, 30667) },
      { AOM_CDF7(12288, 19408, 24627, 26662, 28499, 30667) },
      { AOM_CDF7(2815, 4570, 9416, 10875, 13782, 19863) },
      { AOM_CDF7(12032, 14948, 22187, 23138, 24756, 27635) },
      { AOM_CDF7(14847, 20167, 25433, 26751, 28278, 30119) },
      { AOM_CDF7(18816, 25574, 29030, 29877, 30656, 31506) },
      { AOM_CDF7(23039, 27333, 30220, 30708, 31070, 31826) },
      { AOM_CDF7(12543, 20838, 27455, 28762, 29763, 31546) },
    };

const aom_cdf_prob
    default_palette_uv_size_cdf[PALATTE_BSIZE_CTXS][CDF_SIZE(PALETTE_SIZES)] = {
      { AOM_CDF7(20480, 29888, 32453, 32715, 32751, 32766) },
      { AOM_CDF7(20480, 29888, 32453, 32715, 32751, 32766) },
      { AOM_CDF7(20480, 29888, 32453, 32715, 32751, 32766) },
      { AOM_CDF7(11135, 23641, 31056, 31998, 32496, 32668) },
      { AOM_CDF7(9984, 21999, 29192, 30645, 31640, 32402) },
      { AOM_CDF7(7552, 16614, 24880, 27283, 29254, 31203) },
      { AOM_CDF7(11391, 18656, 23727, 26058, 27788, 30278) },
      { AOM_CDF7(8576, 13585, 17632, 20884, 23948, 27152) },
      { AOM_CDF7(9216, 14276, 19043, 22689, 25799, 28712) },
    };

// When palette mode is enabled, following probability tables indicate the
// probabilities to code the "is_palette" bit (i.e. the bit that indicates
// if this block uses palette mode or DC_PRED mode).
const aom_cdf_prob default_palette_y_mode_cdf[PALATTE_BSIZE_CTXS]
                                             [PALETTE_Y_MODE_CONTEXTS]
                                             [CDF_SIZE(2)] = {
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                               { { AOM_CDF2(128 * 240) },
                                                 { AOM_CDF2(128 * 180) },
                                                 { AOM_CDF2(128 * 100) } },
                                             };

const aom_cdf_prob
    default_palette_uv_mode_cdf[PALETTE_UV_MODE_CONTEXTS][CDF_SIZE(2)] = {
      { AOM_CDF2(128 * 253) }, { AOM_CDF2(128 * 229) }
    };

const aom_cdf_prob default_palette_y_color_index_cdf
    [PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][CDF_SIZE(PALETTE_COLORS)] = {
      {
          { AOM_CDF2(29568), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(16384), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(8832), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(28672), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(31872), 0, 0, 0, 0, 0, 0 },
      },
      {
          { AOM_CDF3(28032, 30326), 0, 0, 0, 0, 0 },
          { AOM_CDF3(11647, 27405), 0, 0, 0, 0, 0 },
          { AOM_CDF3(4352, 30659), 0, 0, 0, 0, 0 },
          { AOM_CDF3(23552, 27800), 0, 0, 0, 0, 0 },
          { AOM_CDF3(32256, 32504), 0, 0, 0, 0, 0 },
      },
      {
          { AOM_CDF4(26112, 28374, 30039), 0, 0, 0, 0 },
          { AOM_CDF4(9472, 22576, 27712), 0, 0, 0, 0 },
          { AOM_CDF4(6656, 26138, 29608), 0, 0, 0, 0 },
          { AOM_CDF4(19328, 23791, 28946), 0, 0, 0, 0 },
          { AOM_CDF4(31744, 31984, 32336), 0, 0, 0, 0 },
      },
      {
          { AOM_CDF5(27904, 29215, 30075, 31190), 0, 0, 0 },
          { AOM_CDF5(9728, 22598, 26134, 29425), 0, 0, 0 },
          { AOM_CDF5(2688, 30066, 31058, 31933), 0, 0, 0 },
          { AOM_CDF5(22015, 25039, 27726, 29932), 0, 0, 0 },
          { AOM_CDF5(32383, 32482, 32554, 32660), 0, 0, 0 },
      },
      {
          { AOM_CDF6(24319, 26299, 27486, 28600, 29804), 0, 0 },
          { AOM_CDF6(7935, 18217, 21116, 25440, 28589), 0, 0 },
          { AOM_CDF6(6656, 25016, 27105, 28698, 30399), 0, 0 },
          { AOM_CDF6(19967, 24117, 26550, 28566, 30224), 0, 0 },
          { AOM_CDF6(31359, 31607, 31775, 31977, 32258), 0, 0 },
      },
      {
          { AOM_CDF7(26368, 27768, 28588, 29274, 29997, 30917), 0 },
          { AOM_CDF7(8960, 18260, 20810, 23986, 26627, 28882), 0 },
          { AOM_CDF7(7295, 24111, 25836, 27515, 29033, 30769), 0 },
          { AOM_CDF7(22016, 25208, 27305, 28159, 29221, 30274), 0 },
          { AOM_CDF7(31744, 31932, 32050, 32199, 32335, 32521), 0 },
      },
      {
          { AOM_CDF8(26624, 27872, 28599, 29153, 29633, 30172, 30841) },
          { AOM_CDF8(6655, 17569, 19587, 23345, 25884, 28088, 29678) },
          { AOM_CDF8(3584, 27296, 28429, 29158, 30032, 30780, 31572) },
          { AOM_CDF8(23551, 25855, 27070, 27893, 28597, 29721, 30970) },
          { AOM_CDF8(32128, 32173, 32245, 32337, 32416, 32500, 32609) },
      },
    };

const aom_cdf_prob default_palette_uv_color_index_cdf
    [PALETTE_SIZES][PALETTE_COLOR_INDEX_CONTEXTS][CDF_SIZE(PALETTE_COLORS)] = {
      {
          { AOM_CDF2(29824), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(16384), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(8832), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(30720), 0, 0, 0, 0, 0, 0 },
          { AOM_CDF2(31744), 0, 0, 0, 0, 0, 0 },
      },
      {
          { AOM_CDF3(27648, 30208), 0, 0, 0, 0, 0 },
          { AOM_CDF3(14080, 26563), 0, 0, 0, 0, 0 },
          { AOM_CDF3(5120, 30932), 0, 0, 0, 0, 0 },
          { AOM_CDF3(24448, 27828), 0, 0, 0, 0, 0 },
          { AOM_CDF3(31616, 32219), 0, 0, 0, 0, 0 },
      },
      {
          { AOM_CDF4(25856, 28259, 30584), 0, 0, 0, 0 },
          { AOM_CDF4(11520, 22476, 27944), 0, 0, 0, 0 },
          { AOM_CDF4(8064, 26882, 30308), 0, 0, 0, 0 },
          { AOM_CDF4(19455, 23823, 29134), 0, 0, 0, 0 },
          { AOM_CDF4(30848, 31501, 32174), 0, 0, 0, 0 },
      },
      {
          { AOM_CDF5(26751, 28020, 29541, 31230), 0, 0, 0 },
          { AOM_CDF5(12032, 26045, 30772, 31497), 0, 0, 0 },
          { AOM_CDF5(1280, 32153, 32458, 32560), 0, 0, 0 },
          { AOM_CDF5(23424, 24154, 29201, 29856), 0, 0, 0 },
          { AOM_CDF5(32256, 32402, 32561, 32682), 0, 0, 0 },
      },
      {
          { AOM_CDF6(24576, 26720, 28114, 28950, 31694), 0, 0 },
          { AOM_CDF6(7551, 16613, 20462, 25269, 29077), 0, 0 },
          { AOM_CDF6(6272, 23039, 25623, 28163, 30861), 0, 0 },
          { AOM_CDF6(17024, 18808, 20771, 27941, 29845), 0, 0 },
          { AOM_CDF6(31616, 31936, 32079, 32321, 32546), 0, 0 },
      },
      {
          { AOM_CDF7(23296, 25590, 27833, 29337, 29954, 31229), 0 },
          { AOM_CDF7(7552, 13659, 16570, 21695, 24506, 27701), 0 },
          { AOM_CDF7(6911, 24788, 26284, 27753, 29575, 30872), 0 },
          { AOM_CDF7(17535, 22236, 24457, 26242, 27363, 30191), 0 },
          { AOM_CDF7(30592, 31289, 31745, 31921, 32149, 32321), 0 },
      },
      {
          { AOM_CDF8(22016, 24242, 25141, 27137, 27797, 29331, 30848) },
          { AOM_CDF8(8063, 13564, 16940, 21948, 24568, 25689, 26989) },
          { AOM_CDF8(6528, 27028, 27835, 28741, 30031, 31795, 32285) },
          { AOM_CDF8(18047, 23797, 25444, 26274, 27111, 27929, 30367) },
          { AOM_CDF8(30208, 30628, 31046, 31658, 31762, 32367, 32469) },
      }
    };

static const aom_cdf_prob default_intrabc_cdf[CDF_SIZE(2)] = {AOM_CDF2(30751)};

#define MAX_COLOR_CONTEXT_HASH 8
// Negative values are invalid
static const int palette_color_index_context_lookup[MAX_COLOR_CONTEXT_HASH +
                                                    1] = { -1, -1, 0, -1, -1,
                                                           4,  3,  2, 1 };

static const aom_cdf_prob default_switchable_restore_cdf[CDF_SIZE(
    RESTORE_SWITCHABLE_TYPES)] = { AOM_CDF3(32 * 128, 144 * 128) };

static const aom_cdf_prob default_wiener_restore_cdf[CDF_SIZE(2)] = { AOM_CDF2(
    64 * 128) };

static const aom_cdf_prob default_sgrproj_restore_cdf[CDF_SIZE(2)] = { AOM_CDF2(
    64 * 128) };

#define NUM_PALETTE_NEIGHBORS 3  // left, top-left and top.
int av1_get_palette_color_index_context(const uint8_t *color_map, int stride,
                                        int r, int c, int palette_size,
                                        uint8_t *color_order, int *color_idx) {
  int i;
  // The +10 below should not be needed. But we get a warning "array subscript
  // is above array bounds [-Werror=array-bounds]" without it, possibly due to
  // this (or similar) bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=59124
  int scores[PALETTE_MAX_SIZE + 10];
  const int weights[NUM_PALETTE_NEIGHBORS] = { 2, 1, 2 };
  const int hash_multipliers[NUM_PALETTE_NEIGHBORS] = { 1, 2, 2 };
  int color_index_ctx_hash;
  int color_index_ctx;
  int color_neighbors[NUM_PALETTE_NEIGHBORS];
  int inverse_color_order[PALETTE_MAX_SIZE];
  assert(palette_size <= PALETTE_MAX_SIZE);
  assert(r > 0 || c > 0);

  // Get color indices of neighbors.
  color_neighbors[0] = (c - 1 >= 0) ? color_map[r * stride + c - 1] : -1;
  color_neighbors[1] =
      (c - 1 >= 0 && r - 1 >= 0) ? color_map[(r - 1) * stride + c - 1] : -1;
  color_neighbors[2] = (r - 1 >= 0) ? color_map[(r - 1) * stride + c] : -1;

  for (i = 0; i < PALETTE_MAX_SIZE; ++i) {
    color_order[i] = i;
    inverse_color_order[i] = i;
  }
  memset(scores, 0, PALETTE_MAX_SIZE * sizeof(scores[0]));
  for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
    if (color_neighbors[i] >= 0) {
      scores[color_neighbors[i]] += weights[i];
    }
  }

  // Get the top NUM_PALETTE_NEIGHBORS scores (sorted from large to small).
  for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
    int max = scores[i];
    int max_idx = i;
    int j;
    for (j = i + 1; j < palette_size; ++j) {
      if (scores[j] > max) {
        max = scores[j];
        max_idx = j;
      }
    }
    if (max_idx != i) {
      // Move the score at index 'max_idx' to index 'i', and shift the scores
      // from 'i' to 'max_idx - 1' by 1.
      const int max_score = scores[max_idx];
      const uint8_t max_color_order = color_order[max_idx];
      int k;
      for (k = max_idx; k > i; --k) {
        scores[k] = scores[k - 1];
        color_order[k] = color_order[k - 1];
        inverse_color_order[color_order[k]] = k;
      }
      scores[i] = max_score;
      color_order[i] = max_color_order;
      inverse_color_order[color_order[i]] = i;
    }
  }

  // Get hash value of context.
  color_index_ctx_hash = 0;
  for (i = 0; i < NUM_PALETTE_NEIGHBORS; ++i) {
    color_index_ctx_hash += scores[i] * hash_multipliers[i];
  }
  assert(color_index_ctx_hash > 0);
  assert(color_index_ctx_hash <= MAX_COLOR_CONTEXT_HASH);

  // Lookup context from hash.
  color_index_ctx = palette_color_index_context_lookup[color_index_ctx_hash];
  assert(color_index_ctx >= 0);
  assert(color_index_ctx < PALETTE_COLOR_INDEX_CONTEXTS);

  if (color_idx != NULL) {
    *color_idx = inverse_color_order[color_map[r * stride + c]];
  }
  return color_index_ctx;
}
#undef NUM_PALETTE_NEIGHBORS
#undef MAX_COLOR_CONTEXT_HASH

static const aom_cdf_prob
    default_txfm_partition_cdf[TXFM_PARTITION_CONTEXTS][CDF_SIZE(2)] = {
        {AOM_CDF2(29409)}, {AOM_CDF2(26504)}, {AOM_CDF2(16941)},
        {AOM_CDF2(32764)}, {AOM_CDF2(32764)}, {AOM_CDF2(32760)},
        {AOM_CDF2(21067)}, {AOM_CDF2(14252)}, {AOM_CDF2(16377)},
        {AOM_CDF2(32764)}, {AOM_CDF2(32764)}, {AOM_CDF2(32764)},
        {AOM_CDF2(27078)}, {AOM_CDF2(20933)}, {AOM_CDF2(17656)},
        {AOM_CDF2(32764)}, {AOM_CDF2(32764)}, {AOM_CDF2(32764)},
        {AOM_CDF2(27632)}, {AOM_CDF2(21712)}, {AOM_CDF2(14535)},
        {AOM_CDF2(16384)}};

static const aom_cdf_prob default_skip_mode_cdfs[SKIP_MODE_CONTEXTS][CDF_SIZE(
    2)] = {{AOM_CDF2(32641)}, {AOM_CDF2(20994)}, {AOM_CDF2(10392)}};
static const aom_cdf_prob default_skip_cdfs[SKIP_CONTEXTS][CDF_SIZE(2)] = {
    {AOM_CDF2(31866)}, {AOM_CDF2(16646)}, {AOM_CDF2(5176)}};

static const aom_cdf_prob
    default_compound_idx_cdfs[COMP_INDEX_CONTEXTS][CDF_SIZE(2)] = {
      { AOM_ICDF(24576), AOM_ICDF(32768), 0 },
      { AOM_ICDF(16384), AOM_ICDF(32768), 0 },
      { AOM_ICDF(8192), AOM_ICDF(32768), 0 },
      { AOM_ICDF(24576), AOM_ICDF(32768), 0 },
      { AOM_ICDF(16384), AOM_ICDF(32768), 0 },
      { AOM_ICDF(8192), AOM_ICDF(32768), 0 },
    };

static const aom_cdf_prob
    default_comp_group_idx_cdfs[COMP_GROUP_IDX_CONTEXTS][CDF_SIZE(2)] = {
      { AOM_ICDF(29491), AOM_ICDF(32768), 0 },
      { AOM_ICDF(24576), AOM_ICDF(32768), 0 },
      { AOM_ICDF(16384), AOM_ICDF(32768), 0 },
      { AOM_ICDF(24576), AOM_ICDF(32768), 0 },
      { AOM_ICDF(16384), AOM_ICDF(32768), 0 },
      { AOM_ICDF(13107), AOM_ICDF(32768), 0 },
      { AOM_ICDF(13107), AOM_ICDF(32768), 0 },
    };

static const aom_cdf_prob default_filter_intra_mode_cdf[CDF_SIZE(
    FILTER_INTRA_MODES)] = {AOM_CDF5(8419, 12492, 16432, 29799)};

static const aom_cdf_prob default_filter_intra_cdfs[BLOCK_SIZES_ALL][CDF_SIZE(
    2)] = {
    {AOM_CDF2(4912)},  {AOM_CDF2(6741)},  {AOM_CDF2(6376)},  {AOM_CDF2(8509)},
    {AOM_CDF2(13558)}, {AOM_CDF2(12794)}, {AOM_CDF2(15293)}, {AOM_CDF2(16189)},
    {AOM_CDF2(16293)}, {AOM_CDF2(23403)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)},
    {AOM_CDF2(13070)}, {AOM_CDF2(10338)}, {AOM_CDF2(20525)}, {AOM_CDF2(17844)},
    {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}, {AOM_CDF2(16384)}};

// FIXME(someone) need real defaults here
static const aom_cdf_prob default_seg_tree_cdf[CDF_SIZE(MAX_SEGMENTS)] = {
  AOM_CDF8(4096, 8192, 12288, 16384, 20480, 24576, 28672)
};

static const aom_cdf_prob
    default_segment_pred_cdf[SEG_TEMPORAL_PRED_CTXS][CDF_SIZE(2)] = {
      { AOM_CDF2(128 * 128) }, { AOM_CDF2(128 * 128) }, { AOM_CDF2(128 * 128) }
    };

static const aom_cdf_prob
    default_switchable_interp_cdf[SWITCHABLE_FILTER_CONTEXTS][CDF_SIZE(
        SWITCHABLE_FILTERS)] = {
        {AOM_CDF3(32760, 32764)}, {AOM_CDF3(10923, 21845)},
        {AOM_CDF3(10923, 21845)}, {AOM_CDF3(32760, 32764)},
        {AOM_CDF3(32760, 32764)}, {AOM_CDF3(10923, 21845)},
        {AOM_CDF3(10923, 21845)}, {AOM_CDF3(32760, 32764)},
        {AOM_CDF3(31505, 32696)}, {AOM_CDF3(5926, 32695)},
        {AOM_CDF3(274, 2338)},    {AOM_CDF3(27887, 32488)},
        {AOM_CDF3(30979, 31753)}, {AOM_CDF3(5974, 31893)},
        {AOM_CDF3(389, 640)},     {AOM_CDF3(14879, 18448)}};

#if CONFIG_SPATIAL_SEGMENTATION
static const aom_cdf_prob
    default_spatial_pred_seg_tree_cdf[SPATIAL_PREDICTION_PROBS][CDF_SIZE(
        MAX_SEGMENTS)] = {
      {
          AOM_CDF8(5622, 7893, 16093, 18233, 27809, 28373, 32533),
      },
      {
          AOM_CDF8(14274, 18230, 22557, 24935, 29980, 30851, 32344),
      },
      {
          AOM_CDF8(27527, 28487, 28723, 28890, 32397, 32647, 32679),
      },
    };
#endif

static const aom_cdf_prob default_tx_size_cdf[MAX_TX_CATS][TX_SIZE_CONTEXTS]
                                             [CDF_SIZE(MAX_TX_DEPTH + 1)] = {
#if MAX_TX_DEPTH == 2
                                               { { AOM_CDF2(19968) },
                                                 { AOM_CDF2(19968) },
                                                 { AOM_CDF2(24320) } },
                                               { { AOM_CDF3(12272, 30172) },
                                                 { AOM_CDF3(12272, 30172) },
                                                 { AOM_CDF3(18677, 30848) } },
                                               { { AOM_CDF3(12986, 15180) },
                                                 { AOM_CDF3(12986, 15180) },
                                                 { AOM_CDF3(24302, 25602) } },
                                               { { AOM_CDF3(5782, 11475) },
                                                 { AOM_CDF3(5782, 11475) },
                                                 { AOM_CDF3(16803, 22759) } },
#elif MAX_TX_DEPTH == 3
                                               { { AOM_CDF2(19968) },
                                                 { AOM_CDF2(24320) } },
                                               { { AOM_CDF3(12272, 30172) },
                                                 { AOM_CDF3(18677, 30848) } },
                                               { { AOM_CDF4(12986, 15180,
                                                            32384) },
                                                 { AOM_CDF4(24302, 25602,
                                                            32128) } },
                                               { { AOM_CDF4(5782, 11475,
                                                            24480) },
                                                 { AOM_CDF4(16803, 22759,
                                                            28560) } },
#else
                                               { { AOM_CDF2(19968) },
                                                 { AOM_CDF2(24320) } },
                                               { { AOM_CDF3(12272, 30172) },
                                                 { AOM_CDF3(18677, 30848) } },
                                               { { AOM_CDF4(12986, 15180,
                                                            32384) },
                                                 { AOM_CDF4(24302, 25602,
                                                            32128) } },
                                               { { AOM_CDF5(5782, 11475, 24480,
                                                            32640) },
                                                 { AOM_CDF5(16803, 22759, 28560,
                                                            32640) } },
#endif  // MAX_TX_DEPTH == 2
                                             };

static const aom_cdf_prob default_if_y_mode_cdf[BLOCK_SIZE_GROUPS][CDF_SIZE(
    INTRA_MODES)] = {{AOM_CDF13(21328, 22235, 23002, 23491, 25028, 26012, 27096,
                                27910, 28722, 30009, 30553, 31222)},
                     {AOM_CDF13(17745, 19623, 21566, 22258, 23709, 25228, 26999,
                                28040, 29237, 30497, 30859, 31409)},
                     {AOM_CDF13(20615, 21979, 23831, 24406, 25281, 26030, 27012,
                                27905, 28735, 30190, 30718, 31471)},
                     {AOM_CDF13(24283, 24866, 26273, 26508, 26811, 26996, 27491,
                                28025, 28278, 29334, 29822, 30565)}};

static const aom_cdf_prob
    default_uv_mode_cdf[CFL_ALLOWED_TYPES][INTRA_MODES][CDF_SIZE(
        UV_INTRA_MODES)] = {
        {{AOM_CDF13(23059, 24483, 25740, 25999, 26285, 26756, 27360, 28208,
                    28471, 30163, 30976, 31898)},
         {AOM_CDF13(9887, 27453, 27530, 27618, 27709, 28140, 28216, 28273,
                    28758, 29697, 29958, 31912)},
         {AOM_CDF13(9530, 9599, 29199, 29235, 29258, 29287, 29608, 30104, 30138,
                    30840, 31819, 32076)},
         {AOM_CDF13(14736, 14975, 15377, 25331, 25671, 25784, 25935, 27130,
                    28515, 30654, 31371, 32504)},
         {AOM_CDF13(8487, 8639, 8868, 8906, 25613, 27592, 28315, 28544, 28734,
                    31132, 31626, 32235)},
         {AOM_CDF13(13648, 14854, 14957, 15016, 15810, 29282, 29371, 29532,
                    29959, 31327, 31591, 32577)},
         {AOM_CDF13(14089, 14152, 16764, 16970, 17730, 17840, 28431, 30330,
                    30362, 31707, 32214, 32657)},
         {AOM_CDF13(13165, 13274, 15199, 15829, 15868, 15947, 16323, 29125,
                    29252, 30625, 31941, 32383)},
         {AOM_CDF13(12619, 14046, 14178, 15707, 15884, 16178, 16399, 16744,
                    28304, 30518, 30996, 32444)},
         {AOM_CDF13(19387, 20131, 20951, 21240, 21425, 21720, 22029, 22773,
                    23193, 28830, 30437, 32130)},
         {AOM_CDF13(18570, 19604, 20967, 21251, 21463, 21778, 22143, 23106,
                    23570, 27011, 31454, 32209)},
         {AOM_CDF13(17527, 18942, 20065, 20401, 20682, 21103, 21372, 21992,
                    22462, 25836, 26479, 31964)},
         {AOM_CDF13(13320, 14738, 15810, 15846, 15873, 15919, 15971, 16044,
                    16130, 16522, 16842, 17315)}},
        {{AOM_CDF14(10475, 11362, 12942, 13236, 13894, 14291, 15050, 15889,
                    16251, 20104, 21111, 22460, 24161)},
         {AOM_CDF14(4176, 20208, 20427, 20578, 20781, 21427, 21542, 21776,
                    22419, 24559, 25030, 27209, 29559)},
         {AOM_CDF14(4977, 5077, 20257, 20347, 20461, 20514, 21061, 21858, 21942,
                    23203, 24314, 24709, 26697)},
         {AOM_CDF14(6598, 7046, 7530, 14027, 14407, 14677, 14937, 16704, 18396,
                    21454, 22192, 23284, 24072)},
         {AOM_CDF14(5865, 6226, 6745, 6892, 19132, 20320, 21879, 22284, 22437,
                    25368, 25849, 26440, 26740)},
         {AOM_CDF14(5356, 6894, 7205, 7352, 9374, 20971, 21267, 21612, 21855,
                    24718, 25376, 26410, 27047)},
         {AOM_CDF14(5547, 5695, 7760, 7881, 9274, 9448, 20286, 20820, 20941,
                    23644, 24260, 24866, 25358)},
         {AOM_CDF14(6030, 6170, 7816, 8371, 8607, 8725, 9156, 18617, 18877,
                    21705, 22671, 23354, 24210)},
         {AOM_CDF14(5978, 7430, 7761, 8880, 9198, 9507, 9706, 10374, 19000,
                    22506, 23200, 24756, 25923)},
         {AOM_CDF14(10758, 11407, 12616, 13021, 13382, 13711, 14137, 15122,
                    15674, 20695, 21829, 23387, 24777)},
         {AOM_CDF14(10504, 11531, 12479, 12894, 13324, 13753, 14203, 15174,
                    15869, 20018, 22552, 23601, 25706)},
         {AOM_CDF14(9690, 10274, 12169, 12591, 12949, 13248, 13755, 14878,
                    15389, 19344, 20208, 23077, 24751)},
         {AOM_CDF14(2910, 4964, 7192, 7334, 7441, 7557, 7685, 7987, 8189, 9325,
                    10075, 10759, 29276)}}};

static const aom_cdf_prob default_partition_cdf[PARTITION_CONTEXTS][CDF_SIZE(
    EXT_PARTITION_TYPES)] = {
  // 8x8 -> 4x4 only supports the four legacy partition types
  { AOM_CDF4(25472, 28949, 31052), 0, 0, 0, 0, 0, 0 },
  { AOM_CDF4(18816, 22250, 28783), 0, 0, 0, 0, 0, 0 },
  { AOM_CDF4(18944, 26126, 29188), 0, 0, 0, 0, 0, 0 },
  { AOM_CDF4(15488, 22508, 27077), 0, 0, 0, 0, 0, 0 },
  // 16x16 -> 8x8
  {AOM_CDF10(15785, 20863, 25027, 27388, 28125, 29181, 29856, 30771, 32056)},
  {AOM_CDF10(8061, 11270, 17074, 23589, 24885, 25856, 27350, 29450, 30450)},
  {AOM_CDF10(6026, 12591, 14989, 21211, 22699, 24684, 25632, 26386, 32237)},
  {AOM_CDF10(2767, 5557, 7763, 22513, 24280, 25712, 27219, 28524, 31514)},
  // 32x32 -> 16x16
  {AOM_CDF10(18545, 20901, 23140, 28349, 28748, 29417, 29779, 30329, 31637)},
  {AOM_CDF10(7246, 8523, 11055, 26485, 26993, 27396, 28020, 29038, 29914)},
  {AOM_CDF10(6403, 8783, 9890, 25945, 26592, 27566, 28014, 28388, 32154)},
  {AOM_CDF10(1311, 1798, 2249, 29850, 30154, 30438, 30722, 30996, 32154)},
  // 64x64 -> 32x32
  {AOM_CDF10(17119, 18723, 20188, 29732, 29936, 30291, 30485, 30773, 31675)},
  {AOM_CDF10(4985, 5649, 7120, 29441, 29669, 29850, 30153, 30605, 31037)},
  {AOM_CDF10(5166, 6618, 7221, 29800, 30096, 30547, 30741, 30913, 32343)},
  {AOM_CDF10(894, 1148, 1350, 31790, 31896, 31997, 32092, 32172, 32489)},
  // 128x128 -> 64x64
  { AOM_CDF8(28416, 28705, 28926, 32258, 32402, 32547, 32548) },
  { AOM_CDF8(9216, 9952, 11849, 30134, 30502, 30870, 30871) },
  { AOM_CDF8(7424, 9008, 9528, 30664, 31456, 32248, 32249) },
  { AOM_CDF8(1280, 1710, 2069, 31978, 32193, 32409, 32410) },
};

static const aom_cdf_prob default_intra_ext_tx_cdf
    [EXT_TX_SETS_INTRA][EXT_TX_SIZES][INTRA_MODES][CDF_SIZE(TX_TYPES)] = {
      {
          // FIXME: unused zero positions, from uncoded trivial transform set
          {
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
          },
          {
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
          },
          {
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
          },
          {
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
              { 0 },
          },
      },
      {
          {
              { AOM_CDF7(1024, 28800, 29048, 29296, 30164, 31466) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1152, 25852, 26284, 26717, 28230, 30499) },
              { AOM_CDF7(1024, 2016, 3938, 5860, 29404, 31086) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
          },
          {
              { AOM_CDF7(1024, 28800, 29048, 29296, 30164, 31466) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1152, 25852, 26284, 26717, 28230, 30499) },
              { AOM_CDF7(1024, 2016, 3938, 5860, 29404, 31086) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
          },
          {
              { AOM_CDF7(1024, 28800, 29048, 29296, 30164, 31466) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1152, 25852, 26284, 26717, 28230, 30499) },
              { AOM_CDF7(1024, 2016, 3938, 5860, 29404, 31086) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
          },
          {
              { AOM_CDF7(1024, 28800, 29048, 29296, 30164, 31466) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1152, 25852, 26284, 26717, 28230, 30499) },
              { AOM_CDF7(1024, 2016, 3938, 5860, 29404, 31086) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 27118) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 4109, 5900, 7691, 15528, 27380) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
              { AOM_CDF7(1280, 5216, 6938, 8660, 10167, 15817) },
          },
      },
      {
          {
              { AOM_CDF5(1024, 28800, 29792, 31280) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1152, 25852, 27581, 30174) },
              { AOM_CDF5(1024, 2016, 28924, 30846) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
          },
          {
              { AOM_CDF5(1024, 28800, 29792, 31280) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1152, 25852, 27581, 30174) },
              { AOM_CDF5(1024, 2016, 28924, 30846) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
          },
          {
              { AOM_CDF5(1024, 28800, 29792, 31280) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1152, 25852, 27581, 30174) },
              { AOM_CDF5(1024, 2016, 28924, 30846) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
          },
          {
              { AOM_CDF5(1024, 28800, 29792, 31280) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1152, 25852, 27581, 30174) },
              { AOM_CDF5(1024, 2016, 28924, 30846) },
              { AOM_CDF5(1280, 5216, 6938, 26310) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 4109, 13065, 26611) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
              { AOM_CDF5(1280, 5216, 6938, 13396) },
          },
      },
    };
static const aom_cdf_prob
    default_inter_ext_tx_cdf[EXT_TX_SETS_INTER][EXT_TX_SIZES][CDF_SIZE(
        TX_TYPES)] = {
      { { 0 }, { 0 }, { 0 }, { 0 } },
      { { AOM_CDF16(1280, 1453, 1626, 2277, 2929, 3580, 4232, 16717, 19225,
                    21733, 24241, 26749, 28253, 29758, 31263) },
        { AOM_CDF16(1280, 1453, 1626, 2277, 2929, 3580, 4232, 16717, 19225,
                    21733, 24241, 26749, 28253, 29758, 31263) },
        { AOM_CDF16(1280, 1453, 1626, 2277, 2929, 3580, 4232, 16717, 19225,
                    21733, 24241, 26749, 28253, 29758, 31263) },
        { AOM_CDF16(1280, 1453, 1626, 2277, 2929, 3580, 4232, 16717, 19225,
                    21733, 24241, 26749, 28253, 29758, 31263) } },
      { { AOM_CDF12(1280, 3125, 4970, 17132, 19575, 22018, 24461, 26904, 28370,
                    29836, 31302) },
        { AOM_CDF12(1280, 3125, 4970, 17132, 19575, 22018, 24461, 26904, 28370,
                    29836, 31302) },
        { AOM_CDF12(1280, 3125, 4970, 17132, 19575, 22018, 24461, 26904, 28370,
                    29836, 31302) },
        { AOM_CDF12(1280, 3125, 4970, 17132, 19575, 22018, 24461, 26904, 28370,
                    29836, 31302) } },
      { { AOM_CDF2(1536) },
        { AOM_CDF2(1536) },
        { AOM_CDF2(1536) },
        { AOM_CDF2(1536) } },
    };

static const aom_cdf_prob default_cfl_sign_cdf[CDF_SIZE(CFL_JOINT_SIGNS)] = {
    AOM_CDF8(1353, 2009, 13371, 18261, 26874, 28140, 32297)};

static const aom_cdf_prob
    default_cfl_alpha_cdf[CFL_ALPHA_CONTEXTS][CDF_SIZE(CFL_ALPHABET_SIZE)] = {
        {AOM_CDF16(9237, 21623, 31523, 32518, 32692, 32696, 32700, 32704, 32708,
                   32712, 32716, 32720, 32724, 32728, 32732)},
        {AOM_CDF16(14277, 24037, 28697, 31427, 32270, 32491, 32566, 32634,
                   32662, 32672, 32676, 32680, 32684, 32688, 32692)},
        {AOM_CDF16(11761, 22512, 28501, 31409, 32372, 32570, 32634, 32680,
                   32684, 32688, 32692, 32696, 32700, 32704, 32708)},
        {AOM_CDF16(26819, 31413, 32284, 32552, 32684, 32688, 32692, 32696,
                   32700, 32704, 32708, 32712, 32716, 32720, 32724)},
        {AOM_CDF16(18000, 26298, 29131, 30728, 31484, 32036, 32280, 32451,
                   32513, 32555, 32590, 32615, 32624, 32633, 32651)},
        {AOM_CDF16(15014, 21756, 25903, 27829, 28872, 30279, 31006, 31842,
                   32115, 32430, 32516, 32600, 32621, 32656, 32660)}};

const aom_cdf_prob
    default_kf_y_mode_cdf[KF_MODE_CONTEXTS][KF_MODE_CONTEXTS][CDF_SIZE(
        INTRA_MODES)] = {
        {{AOM_CDF13(15578, 17033, 19368, 20258, 20749, 21211, 21960, 23415,
                    24372, 28339, 29303, 30717)},
         {AOM_CDF13(12157, 18222, 19644, 20424, 20857, 21593, 22057, 23217,
                    24622, 28811, 30339, 31603)},
         {AOM_CDF13(10125, 10843, 22277, 22767, 23047, 23247, 24172, 25692,
                    26239, 29379, 29994, 31684)},
         {AOM_CDF13(14025, 15316, 16339, 18801, 19147, 19583, 20044, 22150,
                    24774, 29562, 30953, 32474)},
         {AOM_CDF13(12137, 13276, 15550, 16437, 18611, 20070, 22429, 25473,
                    26329, 30159, 31019, 32478)}},
        {{AOM_CDF13(9973, 19727, 20967, 21559, 21968, 22893, 23241, 24180,
                    25544, 29085, 30560, 31527)},
         {AOM_CDF13(6026, 24171, 24634, 24973, 25157, 25882, 26004, 26512,
                    27674, 29940, 31307, 31837)},
         {AOM_CDF13(7487, 12853, 20121, 20679, 21029, 21569, 22137, 23423,
                    24461, 27909, 29111, 30559)},
         {AOM_CDF13(8581, 14706, 15416, 17111, 17460, 18217, 18484, 19912,
                    24611, 29076, 31081, 32211)},
         {AOM_CDF13(7666, 14270, 15506, 16232, 17993, 20752, 22021, 24445,
                    25919, 29576, 30947, 31994)}},
        {{AOM_CDF13(12753, 13729, 21279, 21917, 22248, 22527, 23394, 25072,
                    25756, 29525, 30324, 32130)},
         {AOM_CDF13(9901, 13628, 18517, 19238, 19616, 20175, 20753, 22119,
                    23280, 27831, 29290, 30997)},
         {AOM_CDF13(6394, 6712, 25944, 26180, 26298, 26371, 27041, 28335, 28544,
                    30453, 30790, 32130)},
         {AOM_CDF13(10799, 11777, 14879, 17178, 17541, 17927, 18609, 21516,
                    23861, 28894, 30246, 32241)},
         {AOM_CDF13(9378, 10048, 16425, 17186, 18359, 19092, 21499, 25827,
                    26439, 29976, 30661, 32455)}},
        {{AOM_CDF13(12664, 14429, 15519, 18390, 18792, 19359, 19782, 21476,
                    25057, 29640, 30900, 32414)},
         {AOM_CDF13(8370, 13816, 14487, 17008, 17372, 18094, 18384, 19485,
                    25197, 29431, 31098, 32310)},
         {AOM_CDF13(8620, 9881, 15045, 17669, 18047, 18490, 19170, 21608, 24566,
                    29060, 30050, 32033)},
         {AOM_CDF13(7894, 9163, 9607, 16775, 16929, 17104, 17280, 18875, 26825,
                    30404, 31549, 32539)},
         {AOM_CDF13(9055, 10477, 11672, 15314, 16645, 17956, 19127, 22525,
                    25525, 29813, 30900, 32469)}},
        {{AOM_CDF13(12626, 13681, 15942, 16771, 19085, 20984, 23045, 25445,
                    26216, 30130, 31098, 32475)},
         {AOM_CDF13(9670, 13543, 14976, 15788, 17713, 20841, 22172, 24183,
                    25381, 29579, 31127, 32290)},
         {AOM_CDF13(8495, 9012, 17553, 18205, 19328, 20102, 22565, 25873, 26408,
                    29835, 30486, 32267)},
         {AOM_CDF13(10963, 11937, 13218, 16127, 17089, 18174, 19061, 22716,
                    24594, 29455, 30875, 32497)},
         {AOM_CDF13(7848, 8491, 9991, 10615, 15952, 19309, 23424, 28737, 29147,
                    31313, 31778, 32596)}}};

static const aom_cdf_prob default_angle_delta_cdf[DIRECTIONAL_MODES][CDF_SIZE(
    2 * MAX_ANGLE_DELTA + 1)] = {
    {AOM_CDF7(2141, 4968, 7395, 22873, 27109, 30277)},
    {AOM_CDF7(2363, 5692, 8939, 23918, 27284, 30554)},
    {AOM_CDF7(3695, 11020, 13726, 19311, 22986, 31291)},
    {AOM_CDF7(4653, 11210, 15062, 16892, 21592, 28415)},
    {AOM_CDF7(1731, 10875, 14501, 19443, 22626, 28736)},
    {AOM_CDF7(2863, 10428, 12754, 17754, 21673, 30533)},
    {AOM_CDF7(1940, 11043, 15391, 20148, 22301, 28830)},
    {AOM_CDF7(3595, 10437, 12499, 17489, 21114, 30673)}};

static void init_mode_probs(FRAME_CONTEXT *fc) {
  av1_copy(fc->palette_y_size_cdf, default_palette_y_size_cdf);
  av1_copy(fc->palette_uv_size_cdf, default_palette_uv_size_cdf);
  av1_copy(fc->palette_y_color_index_cdf, default_palette_y_color_index_cdf);
  av1_copy(fc->palette_uv_color_index_cdf, default_palette_uv_color_index_cdf);
  av1_copy(fc->kf_y_cdf, default_kf_y_mode_cdf);
  av1_copy(fc->angle_delta_cdf, default_angle_delta_cdf);
  av1_copy(fc->comp_inter_cdf, default_comp_inter_cdf);
  av1_copy(fc->comp_ref_type_cdf, default_comp_ref_type_cdf);
  av1_copy(fc->uni_comp_ref_cdf, default_uni_comp_ref_cdf);
  av1_copy(fc->palette_y_mode_cdf, default_palette_y_mode_cdf);
  av1_copy(fc->palette_uv_mode_cdf, default_palette_uv_mode_cdf);
  av1_copy(fc->comp_ref_cdf, default_comp_ref_cdf);
  av1_copy(fc->comp_bwdref_cdf, default_comp_bwdref_cdf);
  av1_copy(fc->single_ref_cdf, default_single_ref_cdf);
  av1_copy(fc->txfm_partition_cdf, default_txfm_partition_cdf);
  av1_copy(fc->compound_index_cdf, default_compound_idx_cdfs);
  av1_copy(fc->comp_group_idx_cdf, default_comp_group_idx_cdfs);
  av1_copy(fc->newmv_cdf, default_newmv_cdf);
  av1_copy(fc->zeromv_cdf, default_zeromv_cdf);
  av1_copy(fc->refmv_cdf, default_refmv_cdf);
  av1_copy(fc->drl_cdf, default_drl_cdf);
  av1_copy(fc->motion_mode_cdf, default_motion_mode_cdf);
  av1_copy(fc->obmc_cdf, default_obmc_cdf);
  av1_copy(fc->inter_compound_mode_cdf, default_inter_compound_mode_cdf);
  av1_copy(fc->compound_type_cdf, default_compound_type_cdf);
#if WEDGE_IDX_ENTROPY_CODING
  av1_copy(fc->wedge_idx_cdf, default_wedge_idx_cdf);
#endif
  av1_copy(fc->interintra_cdf, default_interintra_cdf);
  av1_copy(fc->wedge_interintra_cdf, default_wedge_interintra_cdf);
  av1_copy(fc->interintra_mode_cdf, default_interintra_mode_cdf);
  av1_copy(fc->seg.pred_cdf, default_segment_pred_cdf);
  av1_copy(fc->seg.tree_cdf, default_seg_tree_cdf);
  av1_copy(fc->filter_intra_cdfs, default_filter_intra_cdfs);
  av1_copy(fc->filter_intra_mode_cdf, default_filter_intra_mode_cdf);
  av1_copy(fc->switchable_restore_cdf, default_switchable_restore_cdf);
  av1_copy(fc->wiener_restore_cdf, default_wiener_restore_cdf);
  av1_copy(fc->sgrproj_restore_cdf, default_sgrproj_restore_cdf);
  av1_copy(fc->y_mode_cdf, default_if_y_mode_cdf);
  av1_copy(fc->uv_mode_cdf, default_uv_mode_cdf);
  av1_copy(fc->switchable_interp_cdf, default_switchable_interp_cdf);
  av1_copy(fc->partition_cdf, default_partition_cdf);
  av1_copy(fc->intra_ext_tx_cdf, default_intra_ext_tx_cdf);
  av1_copy(fc->inter_ext_tx_cdf, default_inter_ext_tx_cdf);
  av1_copy(fc->skip_mode_cdfs, default_skip_mode_cdfs);
  av1_copy(fc->skip_cdfs, default_skip_cdfs);
  av1_copy(fc->intra_inter_cdf, default_intra_inter_cdf);
#if CONFIG_SPATIAL_SEGMENTATION
  for (int i = 0; i < SPATIAL_PREDICTION_PROBS; i++)
    av1_copy(fc->seg.spatial_pred_seg_cdf[i],
             default_spatial_pred_seg_tree_cdf[i]);
#endif
  av1_copy(fc->tx_size_cdf, default_tx_size_cdf);
  av1_copy(fc->delta_q_cdf, default_delta_q_cdf);
#if CONFIG_EXT_DELTA_Q
  av1_copy(fc->delta_lf_cdf, default_delta_lf_cdf);
  av1_copy(fc->delta_lf_multi_cdf, default_delta_lf_multi_cdf);
#endif
  av1_copy(fc->cfl_sign_cdf, default_cfl_sign_cdf);
  av1_copy(fc->cfl_alpha_cdf, default_cfl_alpha_cdf);
  av1_copy(fc->intrabc_cdf, default_intrabc_cdf);
}

void av1_set_default_ref_deltas(int8_t *ref_deltas) {
  assert(ref_deltas != NULL);

  ref_deltas[INTRA_FRAME] = 1;
  ref_deltas[LAST_FRAME] = 0;
  ref_deltas[LAST2_FRAME] = ref_deltas[LAST_FRAME];
  ref_deltas[LAST3_FRAME] = ref_deltas[LAST_FRAME];
  ref_deltas[BWDREF_FRAME] = ref_deltas[LAST_FRAME];
  ref_deltas[GOLDEN_FRAME] = -1;
  ref_deltas[ALTREF2_FRAME] = -1;
  ref_deltas[ALTREF_FRAME] = -1;
}

void av1_set_default_mode_deltas(int8_t *mode_deltas) {
  assert(mode_deltas != NULL);

  mode_deltas[0] = 0;
  mode_deltas[0] = 0;
}

static void set_default_lf_deltas(struct loopfilter *lf) {
  lf->mode_ref_delta_enabled = 1;
  lf->mode_ref_delta_update = 1;

  av1_set_default_ref_deltas(lf->ref_deltas);
  av1_set_default_mode_deltas(lf->mode_deltas);
}

void av1_setup_frame_contexts(AV1_COMMON *cm) {
  // Store the frame context into a special slot (not associated with any
  // reference buffer), so that we can set up cm->pre_fc correctly later
  // This function must ONLY be called when cm->fc has been initialized with
  // default probs, either by av1_setup_past_independence or after manually
  // initializing them
  cm->frame_contexts[FRAME_CONTEXT_DEFAULTS] = *cm->fc;
}

void av1_setup_past_independence(AV1_COMMON *cm) {
  // Reset the segment feature data to the default stats:
  // Features disabled, 0, with delta coding (Default state).
  av1_clearall_segfeatures(&cm->seg);

  cm->current_frame_seg_map = cm->cur_frame->seg_map;

  if (cm->current_frame_seg_map)
    memset(cm->current_frame_seg_map, 0, (cm->mi_rows * cm->mi_cols));

  // reset mode ref deltas
  av1_set_default_ref_deltas(cm->cur_frame->ref_deltas);
  av1_set_default_mode_deltas(cm->cur_frame->mode_deltas);
  set_default_lf_deltas(&cm->lf);

  av1_default_coef_probs(cm);
  init_mode_probs(cm->fc);
  av1_init_mv_probs(cm);
  av1_init_lv_map(cm);
  cm->fc->initialized = 1;
  av1_setup_frame_contexts(cm);

  // prev_mip will only be allocated in encoder.
  if (frame_is_intra_only(cm) && cm->prev_mip)
    memset(cm->prev_mip, 0,
           cm->mi_stride * cm->mi_rows * sizeof(*cm->prev_mip));
}
