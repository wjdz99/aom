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

#ifndef AOM_AV1_ENCODER_PICKLPF_H_
#define AOM_AV1_ENCODER_PICKLPF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "av1/encoder/encoder.h"

struct yv12_buffer_config;
struct AV1_COMP;
int av1_get_max_filter_level(const AV1_COMP *cpi);

/*!\brief Algorithm for AV1 loop filter level selection.
 *
 * \ingroup in_loop_filter
 * This function determines proper filter levels used for in-loop filter
 * (deblock filter).
 *
 * \param[in]    sd             The pointer of frame buffer
 * \param[in]    cpi            Top-level encoder structure
 * \param[in]    method         The method used to select filter levels
 *
 * \par
 * method includes:
 * \arg \c LPF_PICK_FROM_FULL_IMAGE:  Try the full image with different values.
 * \arg \c LPF_PICK_FROM_FULL_IMAGE_NON_DUAL: Try the full image filter search
 * with non-dual filter only.
 * \arg \c LPF_PICK_FROM_SUBIMAGE: Try a small portion of the image with
 * different values.
 * \arg \c LPF_PICK_FROM_Q: Estimate the level based on quantizer and frame type
 * \arg \c LPF_PICK_MINIMAL_LPF: Pick 0 to disable LPF if LPF was enabled last
 * frame
 *
 * \remark Nothing is returned. Instead, filter levels below are stored in the
 * "loopfilter" structure inside "cpi":
 * \arg \c filter_level[0]: the vertical filter level for Y plane
 * \arg \c filter_level[1]: the horizontal filter level for Y plane
 * \arg \c filter_level_u: the filter level for U plane
 * \arg \c filter_level_v: the filter level for V plane
 *
 * \n
 * \b Overview
 * \par
 * The workflow of deblock filter is shown in Fig.1. \n
 * Boundary pixels pass through a non-flatness check, followed by a step that
 * determines smoothness and selects proper types of filters
 * (4-, 6-, 8-, 14-tap filter). \n
 * If non-flatness criteria is not satisfied, the encoder will not apply
 * deblock filtering on these boundary pixels.
 * \image html filter_flow.png "Fig.1. The workflow of deblock filter" width=70%
 *
 * \par
 * The non-flatness is determined by the boundary pixels and thresholds as shown
 * in Fig.2. \n
 * Filtering is applied when \n
 * \f$|p_0-p_1|<thr_1\f$   and   \f$|q_0-q_1|<thr_1\f$   and
 * \f$2*|p_0-q_0|+|p_1-q_1|/2<thr_2\f$ \n
 * \image html filter_thr.png "Fig.2. Non-flatness of pixel boundary" height=40%
 *
 * \par
 * Thresholds ("thr_1" and "thr_2") are determined by the filter level. \n
 * In AV1, for each frame, we employ the four filter levels, based on these
 * observations: \n
 * Luma and chroma planes have different characteristics, including subsampling
 * (different plane size), coding quality (chroma planes are better coded). \n
 * Therefore chroma planes need less deblocking filtering than luma plane. \n
 * In addition, content texture has different spatial characteristics: vertical
 * and horizontal direction may need different level of filtering. \n
 * The selection of these filter levels is described in the following section.
 *
 * \par
 * \b Algorithm
 * \par
 * The encoder selects filter levels given the current frame buffer, and the
 * method. \n
 * By default, "LPF_PICK_FROM_FULL_IMAGE" is used, which should provide
 * the most appropriate filter levels. \n
 * For video on demand (VOD) mode, if speed setting is larger than 5,
 * "LPF_PICK_FROM_FULL_IMAGE_NON_DUAL" is used. \n
 * For real-time mode, if speed setting is larger than 5, "LPF_PICK_FROM_Q" is
 * used.
 *
 * \par
 * "LPF_PICK_FROM_FULL_IMAGE" method: determine filter levels sequentially
 * by a filter level search procedure (function "search_filter_level"). \n
 * The order is: \n
 * First search and determine the filter level for Y plane.
 * Let vertical filter level (filter_level[0]) and the horizontal filter level
 * (filter_level[1]) be equal to it. \n
 * Keep the horizontal filter level the same and search and determine the
 * vertical filter level. \n
 * Search and determine the horizontal filter level. \n
 * Search and determine filter level for U plane. \n
 * Search and determine filter level for V plane.
 *
 * \par
 * Search and determine filter level is fulfilled by function
 * "search_filter_level". \n
 * It starts with a base filter level ("filt_mid") initialized by the
 * corresponding last frame's filter level. \n
 * A filter step ("filter_step") is determined as:
 * filter_step = filt_mid < 16 ? 4 : filt_mid / 4. \n
 * Then a modified binary search strategy is employed to find a proper
 * filter level. \n
 * In each iteration, set filt_low = filt_mid - filter_step,
 * filt_high = filt_mid + filter_step. \n
 * We now have three candidate levels, "filt_mid", "filt_low" and "filt_high".
 * \n
 * Deblock filtering is applied on the current frame with candidate filter
 * levels and the sum of squared error (SSE) between source and filtered frame
 * is computed. \n
 * Set "filt_best" to the filter level of the smallest SSE. If "filter_best"
 * equals to "filt_mid", halve the filter_step. Otherwise, set filt_mid =
 * filt_best. \n
 * Go to the next iteration until "filter_step" is 0. \n
 * Note that in the comparison of SSEs between SSE[filt_low] and SSE[filt_mid],
 * a "bias" is introduced to slightly raise the filter level. \n
 * It is based on the observation that low filter levels tend to yield a smaller
 * SSE and produce a higher PSNR for the current frame, \n
 * while oversmoothing it and degradating the quality for prediction for future
 * frames and leanding to a suboptimal performance overall. \n
 * Function "try_filter_frame" is the referrence for applying deblock filtering
 * with a given filter level and computatition of SSE.
 *
 * \par
 * "LPF_PICK_FROM_FULL_IMAGE_NON_DUAL" method: almost the same as
 * "LPF_PICK_FROM_FULL_IMAGE", \n
 * just without separately searching for appropriate filter levels for vertical
 * and horizontal filters.
 *
 * \par
 * "LPF_PICK_FROM_Q" method: filter levels are determined by the
 * quantization factor (q). \n
 * For 8 bit: \n
 *   Keyframes: filt_guess = q * 0.06699 - 1.60817 \n
 *   Other frames: filt_guess = q * inter_frame_multiplier + 2.48225 \n
 *   inter_frame_multiplier = q > 700 ? 0.04590 : 0.02295 \n
 * For 10 bit and 12 bit: \n
 * filt_guess = q * 0.316206 + 3.87252 \n
 * Then filter_level[0] = filter_level[1] = filter_level_u = filter_level_v =
 * clamp(filt_guess, min_filter_level, max_filter_level) \n
 * Where min_filter_level = 0, max_filter_level = 64 \n
 * The equations were determined by linear fitting using filter levels
 * generated by "LPF_PICK_FROM_FULL_IMAGE" method.
 *
 */
void av1_pick_filter_level(const struct yv12_buffer_config *sd,
                           struct AV1_COMP *cpi, LPF_PICK_METHOD method);

static AOM_INLINE int disable_lf(AV1_COMP *cpi) {
  const AV1EncoderConfig *const oxcf = &cpi->oxcf;
  int disable_filter_rt_screen = 0;
  if (oxcf->tune_cfg.content == AOM_CONTENT_SCREEN &&
      oxcf->q_cfg.aq_mode == CYCLIC_REFRESH_AQ && cpi->sf.rt_sf.skip_lf_screen)
    disable_filter_rt_screen = av1_cyclic_refresh_disable_lf_cdef(cpi);

  if (disable_filter_rt_screen ||
      oxcf->algo_cfg.loopfilter_control == LOOPFILTER_NONE ||
      (oxcf->algo_cfg.loopfilter_control == LOOPFILTER_REFERENCE &&
       cpi->rtc_ref.non_reference_frame)) {
    return 1;
  }
  return 0;
}

static AOM_INLINE void pick_lf_level_from_q(AV1_COMP *cpi) {
  AV1_COMMON *const cm = &cpi->common;
  const SequenceHeader *const seq_params = cm->seq_params;
  struct loopfilter *const lf = &cm->lf;
  const int min_filter_level = 0;
  const int max_filter_level = av1_get_max_filter_level(cpi);
  const int q =
      av1_ac_quant_QTX(cm->quant_params.base_qindex, 0, seq_params->bit_depth);
  // based on tests result for rtc test set
  // 0.04590 boosted or 0.02295 non-booseted in 18-bit fixed point
  const int strength_boost_q_treshold = 0;
  int inter_frame_multiplier =
      (q > strength_boost_q_treshold ||
       (cpi->sf.rt_sf.use_nonrd_pick_mode &&
        cpi->common.width * cpi->common.height > 352 * 288))
          ? 12034
          : 6017;
  // These values were determined by linear fitting the result of the
  // searched level for 8 bit depth:
  // Keyframes: filt_guess = q * 0.06699 - 1.60817
  // Other frames: filt_guess = q * inter_frame_multiplier + 2.48225
  //
  // And high bit depth separately:
  // filt_guess = q * 0.316206 + 3.87252
  int filt_guess;
  switch (seq_params->bit_depth) {
    case AOM_BITS_8:
      filt_guess =
          (cm->current_frame.frame_type == KEY_FRAME)
              ? ROUND_POWER_OF_TWO(q * 17563 - 421574, 18)
              : ROUND_POWER_OF_TWO(q * inter_frame_multiplier + 650707, 18);
      break;
    case AOM_BITS_10:
      filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 4060632, 20);
      break;
    case AOM_BITS_12:
      filt_guess = ROUND_POWER_OF_TWO(q * 20723 + 16242526, 22);
      break;
    default:
      assert(0 &&
             "bit_depth should be AOM_BITS_8, AOM_BITS_10 "
             "or AOM_BITS_12");
      return;
  }
  if (seq_params->bit_depth != AOM_BITS_8 &&
      cm->current_frame.frame_type == KEY_FRAME)
    filt_guess -= 4;
  // TODO(chengchen): retrain the model for Y, U, V filter levels
  lf->filter_level[0] = clamp(filt_guess, min_filter_level, max_filter_level);
  lf->filter_level[1] = clamp(filt_guess, min_filter_level, max_filter_level);
  lf->filter_level_u = clamp(filt_guess, min_filter_level, max_filter_level);
  lf->filter_level_v = clamp(filt_guess, min_filter_level, max_filter_level);
  if (cpi->oxcf.algo_cfg.loopfilter_control == LOOPFILTER_SELECTIVELY &&
      !frame_is_intra_only(cm) && !cpi->rc.high_source_sad) {
    if (cpi->oxcf.tune_cfg.content == AOM_CONTENT_SCREEN) {
      lf->filter_level[0] = 0;
      lf->filter_level[1] = 0;
    } else {
      const int num4x4 = (cm->width >> 2) * (cm->height >> 2);
      const int newmv_thresh = 7;
      const int distance_since_key_thresh = 5;
      if ((cpi->td.rd_counts.newmv_or_intra_blocks * 100 / num4x4) <
              newmv_thresh &&
          cpi->rc.frames_since_key > distance_since_key_thresh) {
        lf->filter_level[0] = 0;
        lf->filter_level[1] = 0;
      }
    }
  }
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_PICKLPF_H_
