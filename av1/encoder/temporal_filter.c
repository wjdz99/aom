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

#include <math.h>
#include <limits.h>

#include "config/aom_config.h"

#include "av1/common/alloccommon.h"
#include "av1/common/av1_common_int.h"
#include "av1/common/odintrin.h"
#include "av1/common/quant_common.h"
#include "av1/common/reconinter.h"
#include "av1/encoder/av1_quantize.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/extend.h"
#include "av1/encoder/firstpass.h"
#include "av1/encoder/mathutils.h"
#include "av1/encoder/mcomp.h"
#include "av1/encoder/ratectrl.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/segmentation.h"
#include "av1/encoder/temporal_filter.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"
#include "aom_ports/aom_timer.h"
#include "aom_ports/mem.h"
#include "aom_ports/system_state.h"
#include "aom_scale/aom_scale.h"

int COUNT = 0;
// NOTE: All `tf` in this file means `temporal filtering`.
// Forward Declaration.
static void tf_determine_block_partition(const MV block_mv, const int block_mse,
                                         MV *subblock_mvs, int *subblock_mses);

// Does motion search for blocks in temporal filtering. This is the first step
// for temporal filtering. More specifically, given a frame to be filtered and
// another frame as reference, this function searches the reference frame to
// find out the most similar block as that from the frame to be filtered. This
// found block will be further used for weighted averaging.
// NOTE: Besides doing motion search for the entire block, this function will
//       also do motion search for each 1/4 sub-block to get more precise
//       predictions. Then, this function will determines whether to use 4
//       sub-blocks to replace the entire block. If we do need to split the
//       entire block, 4 elements in `subblock_mvs` and `subblock_mses` refer to
//       the searched motion vector and search error (MSE) w.r.t. each sub-block
//       respectively. Otherwise, the 4 elements will be the same, all of which
//       are assigned as the searched motion vector and search error (MSE) for
//       the entire block.
// Inputs:
//   cpi: Pointer to the composed information of input video.
//   frame_to_filter: Pointer to the frame to be filtered.
//   ref_frame: Pointer to the reference frame.
//   block_size: Block size used for motion search.
//   mb_row: Row index of the block in the entire frame.
//   mb_col: Column index of the block in the entire frame.
//   ref_mv: Reference motion vector, which is commonly inherited from the
//           motion search result of previous frame.
//   subblock_mvs: Pointer to the motion vectors for 4 sub-blocks.
//   subblock_mses: Pointer to the search errors (MSE) for 4 sub-blocks.
// Returns:
//   Nothing will be returned. Results are saved in `subblock_mvs` and
//   `subblock_mses`.
static void tf_motion_search(AV1_COMP *cpi,
                             const YV12_BUFFER_CONFIG *frame_to_filter,
                             const YV12_BUFFER_CONFIG *ref_frame,
                             const BLOCK_SIZE block_size, const int mb_row,
                             const int mb_col, MV *ref_mv, MV *subblock_mvs,
                             int *subblock_mses) {
  // Frame information
  const int min_frame_size = AOMMIN(cpi->common.width, cpi->common.height);

  // Block information (ONLY Y-plane is used for motion search).
  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  const int mb_pels = mb_height * mb_width;
  const int y_stride = frame_to_filter->y_stride;
  assert(y_stride == ref_frame->y_stride);
  const int y_offset = mb_row * mb_height * y_stride + mb_col * mb_width;

  // Save input state.
  MACROBLOCK *const mb = &cpi->td.mb;
  MACROBLOCKD *const mbd = &mb->e_mbd;
  const struct buf_2d ori_src_buf = mb->plane[0].src;
  const struct buf_2d ori_pre_buf = mbd->plane[0].pre[0];

  // Parameters used for motion search.
  FULLPEL_MOTION_SEARCH_PARAMS full_ms_params;
  SUBPEL_MOTION_SEARCH_PARAMS ms_params;
  const search_site_config search_site_cfg =
      cpi->mv_search_params.search_site_cfg[SS_CFG_LOOKAHEAD];
  const SEARCH_METHODS full_search_method = NSTEP;
  const int step_param = av1_init_search_range(
      AOMMAX(frame_to_filter->y_crop_width, frame_to_filter->y_crop_height));
  const SUBPEL_SEARCH_TYPE subpel_search_type = USE_8_TAPS;
  const int force_integer_mv = cpi->common.features.cur_frame_force_integer_mv;
  const MV_COST_TYPE mv_cost_type =
      min_frame_size >= 720
          ? MV_COST_L1_HDRES
          : (min_frame_size >= 480 ? MV_COST_L1_MIDRES : MV_COST_L1_LOWRES);
  // Starting position for motion search.
  FULLPEL_MV start_mv = get_fullmv_from_mv(ref_mv);
  // Baseline position for motion search (used for rate distortion comparison).
  const MV baseline_mv = kZeroMv;

  // Setup.
  mb->plane[0].src.buf = frame_to_filter->y_buffer + y_offset;
  mb->plane[0].src.stride = y_stride;
  mbd->plane[0].pre[0].buf = ref_frame->y_buffer + y_offset;
  mbd->plane[0].pre[0].stride = y_stride;
  // Unused intermediate results for motion search.
  unsigned int sse, error;
  int distortion;
  int cost_list[5];

  // Do motion search.
  int_mv best_mv;  // Searched motion vector.
  int block_mse = INT_MAX;
  MV block_mv = kZeroMv;

  av1_make_default_fullpel_ms_params(&full_ms_params, cpi, mb, block_size,
                                     &baseline_mv, &search_site_cfg,
                                     /*fine_search_interval=*/0);
  full_ms_params.run_mesh_search = 1;
  full_ms_params.search_method = full_search_method;
  full_ms_params.mv_cost_params.mv_cost_type = mv_cost_type;

  av1_full_pixel_search(start_mv, &full_ms_params, step_param,
                        cond_cost_list(cpi, cost_list), &best_mv.as_fullmv,
                        NULL);

  if (force_integer_mv == 1) {  // Only do full search on the entire block.
    const int mv_row = best_mv.as_mv.row;
    const int mv_col = best_mv.as_mv.col;
    best_mv.as_mv.row = GET_MV_SUBPEL(mv_row);
    best_mv.as_mv.col = GET_MV_SUBPEL(mv_col);
    const int mv_offset = mv_row * y_stride + mv_col;
    error = cpi->fn_ptr[block_size].vf(
        ref_frame->y_buffer + y_offset + mv_offset, y_stride,
        frame_to_filter->y_buffer + y_offset, y_stride, &sse);
    block_mse = DIVIDE_AND_ROUND(error, mb_pels);
    block_mv = best_mv.as_mv;
  } else {  // Do fractional search on the entire block and all sub-blocks.
    av1_make_default_subpel_ms_params(&ms_params, cpi, mb, block_size,
                                      &baseline_mv, cost_list);
    ms_params.forced_stop = EIGHTH_PEL;
    ms_params.var_params.subpel_search_type = subpel_search_type;
    // Since we are merely refining the result from full pixel search, we don't
    // need regularization for subpel search
    ms_params.mv_cost_params.mv_cost_type = MV_COST_NONE;

    MV subpel_start_mv = get_mv_from_fullmv(&best_mv.as_fullmv);
    error = cpi->mv_search_params.find_fractional_mv_step(
        &mb->e_mbd, &cpi->common, &ms_params, subpel_start_mv, &best_mv.as_mv,
        &distortion, &sse, NULL);
    block_mse = DIVIDE_AND_ROUND(error, mb_pels);
    block_mv = best_mv.as_mv;
    *ref_mv = best_mv.as_mv;
    // On 4 sub-blocks.
    const BLOCK_SIZE subblock_size = ss_size_lookup[block_size][1][1];
    const int subblock_height = block_size_high[subblock_size];
    const int subblock_width = block_size_wide[subblock_size];
    const int subblock_pels = subblock_height * subblock_width;
    start_mv = get_fullmv_from_mv(ref_mv);

    int subblock_idx = 0;
    for (int i = 0; i < mb_height; i += subblock_height) {
      for (int j = 0; j < mb_width; j += subblock_width) {
        const int offset = i * y_stride + j;
        mb->plane[0].src.buf = frame_to_filter->y_buffer + y_offset + offset;
        mbd->plane[0].pre[0].buf = ref_frame->y_buffer + y_offset + offset;

        av1_make_default_fullpel_ms_params(&full_ms_params, cpi, mb,
                                           subblock_size, &baseline_mv,
                                           &search_site_cfg,
                                           /*fine_search_interval=*/0);
        full_ms_params.run_mesh_search = 1;
        full_ms_params.search_method = full_search_method;
        full_ms_params.mv_cost_params.mv_cost_type = mv_cost_type;

        av1_full_pixel_search(start_mv, &full_ms_params, step_param,
                              cond_cost_list(cpi, cost_list),
                              &best_mv.as_fullmv, NULL);

        av1_make_default_subpel_ms_params(&ms_params, cpi, mb, subblock_size,
                                          &baseline_mv, cost_list);
        ms_params.forced_stop = EIGHTH_PEL;
        ms_params.var_params.subpel_search_type = subpel_search_type;
        // Since we are merely refining the result from full pixel search, we
        // don't need regularization for subpel search
        ms_params.mv_cost_params.mv_cost_type = MV_COST_NONE;

        subpel_start_mv = get_mv_from_fullmv(&best_mv.as_fullmv);
        error = cpi->mv_search_params.find_fractional_mv_step(
            &mb->e_mbd, &cpi->common, &ms_params, subpel_start_mv,
            &best_mv.as_mv, &distortion, &sse, NULL);
        subblock_mses[subblock_idx] = DIVIDE_AND_ROUND(error, subblock_pels);
        subblock_mvs[subblock_idx] = best_mv.as_mv;
        ++subblock_idx;
      }
    }
  }

  // Restore input state.
  mb->plane[0].src = ori_src_buf;
  mbd->plane[0].pre[0] = ori_pre_buf;

  // Make partition decision.
  tf_determine_block_partition(block_mv, block_mse, subblock_mvs,
                               subblock_mses);

  // Do not pass down the reference motion vector if error is too large.
  const int thresh = (min_frame_size >= 720) ? 12 : 3;
  if (block_mse > (thresh << (mbd->bd - 8))) {
    *ref_mv = kZeroMv;
  }
}

// Determines whether to split the entire block to 4 sub-blocks for filtering.
// In particular, this decision is made based on the comparison between the
// motion search error of the entire block and the errors of all sub-blocks.
// Inputs:
//   block_mv: Motion vector for the entire block (ONLY as reference).
//   block_mse: Motion search error (MSE) for the entire block (ONLY as
//              reference).
//   subblock_mvs: Pointer to the motion vectors for 4 sub-blocks (will be
//                 modified based on the partition decision).
//   subblock_mses: Pointer to the search errors (MSE) for 4 sub-blocks (will
//                  be modified based on the partition decision).
// Returns:
//   Nothing will be returned. Results are saved in `subblock_mvs` and
//   `subblock_mses`.
static void tf_determine_block_partition(const MV block_mv, const int block_mse,
                                         MV *subblock_mvs, int *subblock_mses) {
  int min_subblock_mse = INT_MAX;
  int max_subblock_mse = INT_MIN;
  int sum_subblock_mse = 0;
  for (int i = 0; i < 4; ++i) {
    sum_subblock_mse += subblock_mses[i];
    min_subblock_mse = AOMMIN(min_subblock_mse, subblock_mses[i]);
    max_subblock_mse = AOMMAX(max_subblock_mse, subblock_mses[i]);
  }

  // TODO(any): The following magic numbers may be tuned to improve the
  // performance OR find a way to get rid of these magic numbers.
  if (((block_mse * 15 < sum_subblock_mse * 4) &&
       max_subblock_mse - min_subblock_mse < 48) ||
      ((block_mse * 14 < sum_subblock_mse * 4) &&
       max_subblock_mse - min_subblock_mse < 24)) {  // No split.
    for (int i = 0; i < 4; ++i) {
      subblock_mvs[i] = block_mv;
      subblock_mses[i] = block_mse;
    }
  }
}

// Helper function to determine whether a frame is encoded with high bit-depth.
static INLINE int is_frame_high_bitdepth(const YV12_BUFFER_CONFIG *frame) {
  return (frame->flags & YV12_FLAG_HIGHBITDEPTH) ? 1 : 0;
}

//debug modes: 4 - pixel mses and mvs used in weights,
//3 - hard decision for pixel vs. block mvs
//2 - pixel mvs, use block mses
//1 - block mvs

#define TF_DEBUG_MODE 4 
#define BP_INTERP 1
#define NEIGHBORHOOD_MSE 1 
#define SAVE_IMGS 0 
#define PRINT_STATEMENTS 0 
// Builds predictor for blocks in temporal filtering. This is the second step
// for temporal filtering, which is to construct predictions from all reference
// frames INCLUDING the frame to be filtered itself. These predictors are built
// based on the motion search results (motion vector is set as 0 for the frame
// to be filtered), and will be futher used for weighted averaging.
// Inputs:
//   ref_frame: Pointer to the reference frame (or the frame to be filtered).
//   mbd: Pointer to the block for filtering. Besides containing the subsampling
//        information of all planes, this field also gives the searched motion
//        vector for the entire block, i.e., `mbd->mi[0]->mv[0]`. This vector
//        should be 0 if the `ref_frame` itself is the frame to be filtered.
//   block_size: Size of the block.
//   mb_row: Row index of the block in the entire frame.
//   mb_col: Column index of the block in the entire frame.
//   num_planes: Number of planes in the frame.
//   scale: Scaling factor.
//   subblock_mvs: The motion vectors for each sub-block (row-major order).
//   pred: Pointer to the predictor to build.
// Returns:
//   Nothing will be returned. But the content to which `pred` points will be
//   modified.
static void tf_build_predictor(const YV12_BUFFER_CONFIG *ref_frame,
                               const MACROBLOCKD *mbd,
                               const BLOCK_SIZE block_size, const int mb_row,
                               const int mb_col, const int num_planes,
                               const struct scale_factors *scale,
                               const MV *subblock_mvs, uint8_t *pred,
                               const YV12_BUFFER_CONFIG *frame, int *subblock_mses) {
  // Information of the entire block.
  const int mb_height = block_size_high[block_size];  // Height.
  const int mb_width = block_size_wide[block_size];   // Width.
  const int mb_pels = mb_height * mb_width;           // Number of pixels.
  const int mb_y = mb_height * mb_row;                // Y-coord (Top-left).
  const int mb_x = mb_width * mb_col;                 // X-coord (Top-left).
  const int bit_depth = mbd->bd;                      // Bit depth.
  const int is_intrabc = 0;                           // Is intra-copied?
  const int is_high_bitdepth = is_frame_high_bitdepth(ref_frame);

  // Default interpolation filters.
  const int_interpfilters interp_filters =
      av1_broadcast_interp_filter(MULTITAP_SHARP);
  // Handle Y-plane, U-plane and V-plane (if needed) in sequence.
  int plane_offset = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    const int subsampling_y = mbd->plane[plane].subsampling_y;
    const int subsampling_x = mbd->plane[plane].subsampling_x;
    // Information of each sub-block in current plane.
    const int plane_h = mb_height >> subsampling_y;  // Plane height.
    const int plane_w = mb_width >> subsampling_x;   // Plane width.
    const int plane_y = mb_y >> subsampling_y;       // Y-coord (Top-left).
    const int plane_x = mb_x >> subsampling_x;       // X-coord (Top-left).
    const int h = plane_h >> 1;                      // Sub-block height.
    const int w = plane_w >> 1;                      // Sub-block width.
    const int is_y_plane = (plane == 0);             // Is Y-plane?

    const struct buf_2d ref_buf = { NULL, ref_frame->buffers[plane],
                                    ref_frame->widths[is_y_plane ? 0 : 1],
                                    ref_frame->heights[is_y_plane ? 0 : 1],
                                    ref_frame->strides[is_y_plane ? 0 : 1] };

#if PRINT_STATEMENTS == 1
    if(is_y_plane  && plane_y < 33 && plane_x < 10)
      printf("Plane info: h, w %d %d, y, x %d %d, sbblck h,w %d %d subsamp y,x %d %d mby %d mbx %d\n",
             plane_h, plane_w, plane_y, plane_x, h, w, subsampling_y, subsampling_x, mb_y, mb_x);
#endif
#if SAVE_IMGS == 1
    int *predimg = malloc(mb_height*mb_width*sizeof(int));
    int *origimg = malloc(mb_height*mb_width*sizeof(int));
#endif
    int newerror = 0;
    int to_match = 5;
    // Handle each subblock.
    int subblock_idx = 0;
#if TF_DEBUG_MODE == 4
//use pixel mvs everywhere. Use pixel mses instead of subblock
// for weights.
    int newblock_size = 1;
    for (int i = 0; i < plane_h; i += newblock_size) {
      for (int j = 0; j < plane_w; j += newblock_size) {
        // Choose proper motion vector.
        const MV mv = subblock_mvs[subblock_idx];
        //TODO: build pred interp doesn't behave as expected for 1x1
        // although it's still better than my interp implementation
        uint8_t temppred[4]; 
        assert(mv.row >= INT16_MIN && mv.row <= INT16_MAX &&
               mv.col >= INT16_MIN && mv.col <= INT16_MAX);
        const int y = plane_y + i;
        const int x = plane_x + j;
        // Build predictor for each sub-block on current plane.
#if BP_INTERP == 1
        InterPredParams inter_pred_params;
        av1_init_inter_params(&inter_pred_params, newblock_size, newblock_size, y, x, subsampling_x,
                              subsampling_y, bit_depth, is_high_bitdepth,
                              is_intrabc, scale, &ref_buf, interp_filters);
        inter_pred_params.conv_params = get_conv_params(0, plane, bit_depth);
        av1_enc_build_one_inter_predictor(&temppred,
                                          newblock_size, &mv, &inter_pred_params);
        
#else
        //I think dividing by 8 is same as >> 3 ? but probably won't keep this
        get_subpixels(ref_buf, &temppred, 1, 1, mv.row/8.0, mv.col/8.0, 
                   x, y);
#endif
        pred[plane_offset+i*plane_w+j] = temppred[0];
        int pixelerror = pow((temppred[0] - frame->y_buffer[y*frame->y_stride + x]), 2);
        subblock_mses[subblock_idx] = pixelerror;
        subblock_idx += newblock_size; //probably incorrect for blocks not 1x1
      }
    }
#if NEIGHBORHOOD_MSE == 1
    subblock_idx = 0;
    //TODO: this can be modified to a window almost always
    //centered at the pixel, if we remove blocks everywhere
    //and calculate this after building entire image predictor. 
    // -> corner of a block can be a window if entire image is used
    if (is_y_plane){
    for (int i = 0; i < plane_h; i += 1) {
      for (int j = 0; j < plane_w; j += 1) {
        int pixelerror = 0;
        int countpixels = 0;
        int winsize = 5;//3;
        int starty = AOMMAX(0, i-winsize/2);
        int startx = AOMMAX(0, j-winsize/2);
        int endy = AOMMIN(mb_height-1, starty+winsize);
        int endx = AOMMIN(mb_width-1, startx+winsize);
        //without these, smaller neighborhood at corners of mb:
        //starty = AOMMIN(starty, endy-winsize);
        //startx = AOMMAX(startx,endx-winsize);
        for(int ii = starty; ii < endy; ii++){
          for(int jj = startx; jj < endx; jj++){
            const int y = plane_y + ii;
            const int x = plane_x + jj;
            pixelerror += pow((pred[plane_offset + ii*plane_w + jj] - frame->y_buffer[y*frame->y_stride + x]), 2);
            countpixels++;
          }
        }
        subblock_mses[subblock_idx] = round(1.0*pixelerror/countpixels);
        subblock_idx += 1;

      }
    }}
#endif
#elif TF_DEBUG_MODE == 2
    int newblock_size = 1;
    for (int i = 0; i < plane_h; i += newblock_size) {
      for (int j = 0; j < plane_w; j += newblock_size) {
        // Choose proper motion vector.
        const MV mv = subblock_mvs[subblock_idx];
        subblock_idx += newblock_size;
        uint8_t temppred[4];
        assert(mv.row >= INT16_MIN && mv.row <= INT16_MAX &&
               mv.col >= INT16_MIN && mv.col <= INT16_MAX);
        const int y = plane_y + i;
        const int x = plane_x + j;
        // Build predictor for each sub-block on current plane.
        InterPredParams inter_pred_params;
        av1_init_inter_params(&inter_pred_params, newblock_size, newblock_size, y, x, subsampling_x,
                              subsampling_y, bit_depth, is_high_bitdepth,
                              is_intrabc, scale, &ref_buf, interp_filters);
        inter_pred_params.conv_params = get_conv_params(0, plane, bit_depth);
        av1_enc_build_one_inter_predictor(&temppred,
                                          1, &mv, &inter_pred_params);
        pred[plane_offset+i*plane_w+j] = temppred[0];
      }
    }
#else
    for (int i = 0; i < plane_h; i += h) {
      for (int j = 0; j < plane_w; j += w) {
        // Choose proper motion vector.
        const MV mv = subblock_mvs[subblock_idx++];
        assert(mv.row >= INT16_MIN && mv.row <= INT16_MAX &&
               mv.col >= INT16_MIN && mv.col <= INT16_MAX);

        const int y = plane_y + i;
        const int x = plane_x + j;

        // Build predictior for each sub-block on current plane.
        InterPredParams inter_pred_params;
        av1_init_inter_params(&inter_pred_params, w, h, y, x, subsampling_x,
                              subsampling_y, bit_depth, is_high_bitdepth,
                              is_intrabc, scale, &ref_buf, interp_filters);
        inter_pred_params.conv_params = get_conv_params(0, plane, bit_depth);
        av1_enc_build_one_inter_predictor(&pred[plane_offset + i * plane_w + j],
                                          plane_w, &mv, &inter_pred_params);
      }
    }
#endif
#if TF_DEBUG_MODE > 0
    if(is_y_plane){
#if TF_DEBUG_MODE == 2
    for(int i = 0; i < 4; i++){
      subblock_mses[i] = 0;
    }
#endif
    for (int i = 0; i < plane_h; i += 1) {
      for (int j = 0; j < plane_w; j += 1) {
        const int y = plane_y + i;
        const int x = plane_x + j;
        const int subblock_idx =
          (i >= mb_height / 2) * 2 + (j >= mb_width / 2);
#if TF_DEBUG_MODE == 2
        subblock_mses[subblock_idx] += pow((pred[plane_offset + i*plane_w + j] - frame->y_buffer[y*frame->y_stride + x]), 2);
#else
        newerror += pow((pred[plane_offset + i*plane_w + j] - frame->y_buffer[y*frame->y_stride + x]), 2);
#endif
        //predimg[i*mb_width+j] = pred[plane_offset + i * plane_w + j];
        //origimg[i*mb_width +j] =  frame->y_buffer[y*frame->y_stride + x];
      }
    }
#if TF_DEBUG_MODE == 2
    for(int i = 0; i < 4; i++){
      subblock_mses[i] = round(1.0*subblock_mses[i]/(h*w)); // % subblockheight*subblockwidth
    }
#endif
    }
    if ( plane == 0 && (mb_row == to_match || mb_col == to_match)){
#if SAVE_IMGS == 1
        save_img(predimg, mb_height, mb_width, 0, "predblock");
        save_img(origimg, mb_height, mb_width, 0, "origblock");
#endif
      }
#endif
    plane_offset += mb_pels;
  }
}

// Computes temporal filter weights and accumulators for the frame to be
// filtered. More concretely, the filter weights for all pixels are the same.
// Inputs:
//   mbd: Pointer to the block for filtering, which is ONLY used to get
//        subsampling information of all planes as well as the bit-depth.
//   block_size: Size of the block.
//   num_planes: Number of planes in the frame.
//   pred: Pointer to the well-built predictors.
//   accum: Pointer to the pixel-wise accumulator for filtering.
//   count: Pointer to the pixel-wise counter fot filtering.
// Returns:
//   Nothing will be returned. But the content to which `accum` and `pred`
//   point will be modified.
void tf_apply_temporal_filter_self(const MACROBLOCKD *mbd,
                                   const BLOCK_SIZE block_size,
                                   const int num_planes, const uint8_t *pred,
                                   uint32_t *accum, uint16_t *count) {
  // Block information.
  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  const int mb_pels = mb_height * mb_width;
  const int is_high_bitdepth = is_cur_buf_hbd(mbd);
  const uint16_t *pred16 = CONVERT_TO_SHORTPTR(pred);

  int plane_offset = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    const int subsampling_y = mbd->plane[plane].subsampling_y;
    const int subsampling_x = mbd->plane[plane].subsampling_x;
    const int h = mb_height >> subsampling_y;  // Plane height.
    const int w = mb_width >> subsampling_x;   // Plane width.

    int pred_idx = 0;
    for (int i = 0; i < h; ++i) {
      for (int j = 0; j < w; ++j) {
        const int idx = plane_offset + pred_idx;  // Index with plane shift.
        const int pred_value = is_high_bitdepth ? pred16[idx] : pred[idx];
        accum[idx] += TF_WEIGHT_SCALE * pred_value;
        count[idx] += TF_WEIGHT_SCALE;
        ++pred_idx;
      }
    }
    plane_offset += mb_pels;
  }
}

// Function to compute pixel-wise squared difference between two buffers.
// Inputs:
//   ref: Pointer to reference buffer.
//   ref_offset: Start position of reference buffer for computation.
//   ref_stride: Stride for reference buffer.
//   tgt: Pointer to target buffer.
//   tgt_offset: Start position of target buffer for computation.
//   tgt_stride: Stride for target buffer.
//   height: Height of block for computation.
//   width: Width of block for computation.
//   is_high_bitdepth: Whether the two buffers point to high bit-depth frames.
//   square_diff: Pointer to save the squared differces.
// Returns:
//   Nothing will be returned. But the content to which `square_diff` points
//   will be modified.
static INLINE void compute_square_diff(const uint8_t *ref, const int ref_offset,
                                       const int ref_stride, const uint8_t *tgt,
                                       const int tgt_offset,
                                       const int tgt_stride, const int height,
                                       const int width,
                                       const int is_high_bitdepth,
                                       uint32_t *square_diff) {
  const uint16_t *ref16 = CONVERT_TO_SHORTPTR(ref);
  const uint16_t *tgt16 = CONVERT_TO_SHORTPTR(tgt);

  int ref_idx = 0;
  int tgt_idx = 0;
  int idx = 0;
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const uint16_t ref_value = is_high_bitdepth ? ref16[ref_offset + ref_idx]
                                                  : ref[ref_offset + ref_idx];
      const uint16_t tgt_value = is_high_bitdepth ? tgt16[tgt_offset + tgt_idx]
                                                  : tgt[tgt_offset + tgt_idx];
      const uint32_t diff = (ref_value > tgt_value) ? (ref_value - tgt_value)
                                                    : (tgt_value - ref_value);
      square_diff[idx] = diff * diff;

      ++ref_idx;
      ++tgt_idx;
      ++idx;
    }
    ref_idx += (ref_stride - width);
    tgt_idx += (tgt_stride - width);
  }
}

// Applies temporal filtering.
// Inputs:
//   frame_to_filter: Pointer to the frame to be filtered, which is used as
//                    reference to compute squared differece from the predictor.
//   mbd: Pointer to the block for filtering, which is ONLY used to get
//        subsampling information of all planes.
//   block_size: Size of the block.
//   mb_row: Row index of the block in the entire frame.
//   mb_col: Column index of the block in the entire frame.
//   num_planes: Number of planes in the frame.
//   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
//                 with each plane (in Y, U, V order).
//   subblock_mvs:  Pointer to the motion vectors for 4 sub-blocks.
//   subblock_mses: Pointer to the search errors (MSE) for 4 sub-blocks.
//   q_factor: Quantization factor. This is actually the `q` defined in libaom,
//             which is converted from `qindex`.
//   filter_strength: Filtering strength. This value lies in range [0, 6], where
//                    6 is the maximum filtering strength.
//   pred: Pointer to the well-built predictors.
//   accum: Pointer to the pixel-wise accumulator for filtering.
//   count: Pointer to the pixel-wise counter fot filtering.
// Returns:
//   Nothing will be returned. But the content to which `accum` and `pred`
//   point will be modified.
void av1_apply_temporal_filter_c(
    const YV12_BUFFER_CONFIG *frame_to_filter, const MACROBLOCKD *mbd,
    const BLOCK_SIZE block_size, const int mb_row, const int mb_col,
    const int num_planes, const double *noise_levels, const MV *subblock_mvs,
    const int *subblock_mses, const int q_factor, const int filter_strength,
    const uint8_t *pred, uint32_t *accum, uint16_t *count) {
  // Block information.
  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  const int mb_pels = mb_height * mb_width;
  const int is_high_bitdepth = is_frame_high_bitdepth(frame_to_filter);
  const uint16_t *pred16 = CONVERT_TO_SHORTPTR(pred);
  // Frame information.
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
  const int min_frame_size = AOMMIN(frame_height, frame_width);

  // Allocate memory for pixel-wise squared differences for all planes. They,
  // regardless of the subsampling, are assigned with memory of size `mb_pels`.
  uint32_t *square_diff =
      aom_memalign(16, num_planes * mb_pels * sizeof(uint32_t));
  memset(square_diff, 0, num_planes * mb_pels * sizeof(square_diff[0]));

  int plane_offset = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    // Locate pixel on reference frame.
    const int plane_h = mb_height >> mbd->plane[plane].subsampling_y;
    const int plane_w = mb_width >> mbd->plane[plane].subsampling_x;
    const int frame_stride = frame_to_filter->strides[plane == 0 ? 0 : 1];
    const int frame_offset = mb_row * plane_h * frame_stride + mb_col * plane_w;
    const uint8_t *ref = frame_to_filter->buffers[plane];
    compute_square_diff(ref, frame_offset, frame_stride, pred, plane_offset,
                        plane_w, plane_h, plane_w, is_high_bitdepth,
                        square_diff + plane_offset);
    plane_offset += mb_pels;
  }

  // Get window size for pixel-wise filtering.
  assert(TF_WINDOW_LENGTH % 2 == 1);
  const int half_window = TF_WINDOW_LENGTH >> 1;

  // Handle planes in sequence.
  plane_offset = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    const int subsampling_y = mbd->plane[plane].subsampling_y;
    const int subsampling_x = mbd->plane[plane].subsampling_x;
    const int h = mb_height >> subsampling_y;  // Plane height.
    const int w = mb_width >> subsampling_x;   // Plane width.

    // Perform filtering.
    int pred_idx = 0;
    for (int i = 0; i < h; ++i) {
      for (int j = 0; j < w; ++j) {
        // non-local mean approach
        uint64_t sum_square_diff = 0;
        int num_ref_pixels = 0;

        for (int wi = -half_window; wi <= half_window; ++wi) {
          for (int wj = -half_window; wj <= half_window; ++wj) {
            const int y = CLIP(i + wi, 0, h - 1);  // Y-coord on current plane.
            const int x = CLIP(j + wj, 0, w - 1);  // X-coord on current plane.
            sum_square_diff += square_diff[plane_offset + y * w + x];
            ++num_ref_pixels;
          }
        }

        // Filter U-plane and V-plane using Y-plane. This is because motion
        // search is only done on Y-plane, so the information from Y-plane will
        // be more accurate.
        if (plane != 0) {
          const int ss_y_shift = subsampling_y - mbd->plane[0].subsampling_y;
          const int ss_x_shift = subsampling_x - mbd->plane[0].subsampling_x;
          for (int ii = 0; ii < (1 << ss_y_shift); ++ii) {
            for (int jj = 0; jj < (1 << ss_x_shift); ++jj) {
              const int yy = (i << ss_y_shift) + ii;  // Y-coord on Y-plane.
              const int xx = (j << ss_x_shift) + jj;  // X-coord on Y-plane.
              const int ww = w << ss_x_shift;         // Width of Y-plane.
              sum_square_diff += square_diff[yy * ww + xx];
              ++num_ref_pixels;
            }
          }
        }

        // Scale down the difference for high bit depth input.
        if (mbd->bd > 8) sum_square_diff >>= (mbd->bd - 8) * (mbd->bd - 8);

        // Combine window error and block error, and normalize it.
        const double window_error = (double)sum_square_diff / num_ref_pixels;
        const int subblock_idx = (i >= h / 2) * 2 + (j >= w / 2);
        const double block_error = (double)subblock_mses[subblock_idx];
        const double combined_error =
            (TF_WINDOW_BLOCK_BALANCE_WEIGHT * window_error + block_error) /
            (TF_WINDOW_BLOCK_BALANCE_WEIGHT + 1) / TF_SEARCH_ERROR_NORM_WEIGHT;

        // Decay factors for non-local mean approach.
        // Larger noise -> larger filtering weight.
        const double n_decay = 0.5 + log(2 * noise_levels[plane] + 5.0);
        // Smaller q -> smaller filtering weight.
        const double q_decay =
            CLIP(pow((double)q_factor / TF_Q_DECAY_THRESHOLD, 2), 1e-5, 1);
        // Smaller strength -> smaller filtering weight.
        const double s_decay = CLIP(
            pow((double)filter_strength / TF_STRENGTH_THRESHOLD, 2), 1e-5, 1);
        // Larger motion vector -> smaller filtering weight.
        const MV mv = subblock_mvs[subblock_idx];
        const double distance = sqrt(pow(mv.row, 2) + pow(mv.col, 2));
        const double distance_threshold =
            (double)AOMMAX(min_frame_size * TF_SEARCH_DISTANCE_THRESHOLD, 1);
        const double d_factor = AOMMAX(distance / distance_threshold, 1);

        // Compute filter weight.
        const double scaled_error =
            AOMMIN(combined_error * d_factor / n_decay / q_decay / s_decay, 7);
        const int weight = (int)(exp(-scaled_error) * TF_WEIGHT_SCALE);

        const int idx = plane_offset + pred_idx;  // Index with plane shift.
        const int pred_value = is_high_bitdepth ? pred16[idx] : pred[idx];
        accum[idx] += weight * pred_value;
        count[idx] += weight;

        ++pred_idx;
      }
    }
    plane_offset += mb_pels;
  }

  aom_free(square_diff);
}

// Normalizes the accumulated filtering result to produce the filtered frame.
// Inputs:
//   mbd: Pointer to the block for filtering, which is ONLY used to get
//        subsampling information of all planes.
//   block_size: Size of the block.
//   mb_row: Row index of the block in the entire frame.
//   mb_col: Column index of the block in the entire frame.
//   num_planes: Number of planes in the frame.
//   accum: Pointer to the pre-computed accumulator.
//   count: Pointer to the pre-computed count.
//   result_buffer: Pointer to result buffer.
// Returns:
//   Nothing will be returned. But the content to which `result_buffer` point
//   will be modified.
static void tf_normalize_filtered_frame(
    const MACROBLOCKD *mbd, const BLOCK_SIZE block_size, const int mb_row,
    const int mb_col, const int num_planes, const uint32_t *accum,
    const uint16_t *count, YV12_BUFFER_CONFIG *result_buffer) {
  // Block information.
  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  const int mb_pels = mb_height * mb_width;
  const int is_high_bitdepth = is_frame_high_bitdepth(result_buffer);

  int plane_offset = 0;
  for (int plane = 0; plane < num_planes; ++plane) {
    const int plane_h = mb_height >> mbd->plane[plane].subsampling_y;
    const int plane_w = mb_width >> mbd->plane[plane].subsampling_x;
    const int frame_stride = result_buffer->strides[plane == 0 ? 0 : 1];
    const int frame_offset = mb_row * plane_h * frame_stride + mb_col * plane_w;
    uint8_t *const buf = result_buffer->buffers[plane];
    uint16_t *const buf16 = CONVERT_TO_SHORTPTR(buf);

    int plane_idx = 0;             // Pixel index on current plane (block-base).
    int frame_idx = frame_offset;  // Pixel index on the entire frame.
    for (int i = 0; i < plane_h; ++i) {
      for (int j = 0; j < plane_w; ++j) {
        const int idx = plane_idx + plane_offset;
        const uint16_t rounding = count[idx] >> 1;
        if (is_high_bitdepth) {
          buf16[frame_idx] =
              (uint16_t)OD_DIVU(accum[idx] + rounding, count[idx]);
        } else {
          buf[frame_idx] = (uint8_t)OD_DIVU(accum[idx] + rounding, count[idx]);
        }
        ++plane_idx;
        ++frame_idx;
      }
      frame_idx += (frame_stride - plane_w);
    }
    plane_offset += mb_pels;
  }
}

// Helper function to compute number of blocks on either side of the frame.
static INLINE int get_num_blocks(const int frame_length, const int mb_length) {
  return (frame_length + mb_length - 1) / mb_length;
}

// Helper function to get `q` used for encoding.
static INLINE int get_q(const AV1_COMP *cpi) {
  const FRAME_TYPE frame_type =
      (cpi->common.current_frame.frame_number > 1) ? INTER_FRAME : KEY_FRAME;
  const int q = (int)av1_convert_qindex_to_q(
      cpi->rc.avg_frame_qindex[frame_type], cpi->common.seq_params.bit_depth);
  return q;
}

typedef struct {
  int64_t sum;
  int64_t sse;
} FRAME_DIFF;

#if TF_DEBUG_MODE > 0
void dump_frame(AV1_COMP *cpi, const YV12_BUFFER_CONFIG *frame_buf,
                const int frame_idx, char *fname_addition) {
  AV1_COMMON *const cm = &cpi->common;
  if (frame_buf == NULL) {
    printf("Frame is not ready.\n");
    return;
  }

  int h;
  char file_name[256];
  snprintf(file_name, 256, "/tmp/frame_%d%s.yuv", frame_idx, fname_addition);
  FILE *f_frame = NULL;
  if ((f_frame = fopen(file_name, "w")) == NULL) {
    printf("Unable to open file %s to write.\n", file_name);
    return;
  }
  printf(
      "\nDumping Frame=%d, Filename=%s, encode_update_type[%5d]=%1d,"
      "show_frame=%d, show_existing_frame=%d, source_alt_ref_active=%d, "
      "refresh_alt_ref_frame=%d, "
      "y_stride=%4d, uv_stride=%4d, cm->width=%4d, cm->height=%4d\n\n",
      frame_idx, file_name, cpi->gf_group.index,
      cpi->gf_group.update_type[cpi->gf_group.index], cm->show_frame,
      cm->show_existing_frame, cpi->rc.source_alt_ref_active,
      cpi->refresh_frame.alt_ref_frame, frame_buf->y_stride,
      frame_buf->uv_stride, cm->width, cm->height);

  // --- Y ---
  for (h = 0; h < cm->height; ++h) {
    fwrite(&frame_buf->y_buffer[h * frame_buf->y_stride], sizeof(uint8_t),
           cm->width, f_frame);
  }
  // --- U ---
  for (h = 0; h < (cm->height >> 1); ++h) {
    fwrite(&frame_buf->u_buffer[h * frame_buf->uv_stride], sizeof(uint8_t),
           (cm->width >> 1), f_frame);
  }
  // --- V ---
  for (h = 0; h < (cm->height >> 1); ++h) {
    fwrite(&frame_buf->v_buffer[h * frame_buf->uv_stride], sizeof(uint8_t),
           (cm->width >> 1), f_frame);
  }

  fclose(f_frame);
}
#endif

#if TF_DEBUG_MODE > 1
//TODO: try out FULLPEL_MV ? already defined?
typedef struct LOCALMV{
  double row;
  double col;
} LOCALMV;
// Sobel filter to compute spatial gradient
void spatial_gradient(const YV12_BUFFER_CONFIG *frame, const int x_coord,
                      const int y_coord, const int direction,
                      double *derivative) {
  double values[9];
  double d = 0;
  double filter[9];
  int idx = 0;
  for (int yy = -1; yy <= 1; yy++) {
    for (int xx = -1; xx <= 1; xx++) {
      values[idx++] =
          frame->y_buffer[(y_coord + yy) * frame->y_stride + (x_coord + xx)] /
          255.0;
    }
  }
  // Sobel filters
  double gx[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
  double gy[9] = { -1, -2, -1, 0, 0, 0, 1, 2, 1 };
  //Scharr filters
  //double gx[9] = { -3, 0, 3, -10, 0, 10, -3, 0, 3 };
  //double gy[9] = { -3, -10, -3, 0, 0, 0, 3, 10, 3 };
  if (direction == 0) {  // x direction
    memcpy(filter, gx, sizeof(filter));
  } else {  // y direction
    memcpy(filter, gy, sizeof(filter));
  }
  for (int i = 0; i < 9; i++) {
    d += filter[i] * values[i];
  }
  *derivative = d;
}
void temporal_gradient(const YV12_BUFFER_CONFIG *frame,
                       const YV12_BUFFER_CONFIG *frame2, const int x_coord,
                       const int y_coord, double *derivative, LOCALMV *mv, 
                       //const int y_coord, double *derivative, MV mv, 
                       AV1_COMP *cpi) {
  MACROBLOCK *const mb = &cpi->td.mb;
  MACROBLOCKD *const mbd = &mb->e_mbd;
  const int bit_depth = mbd->bd;                      // Bit depth.
  const int is_intrabc = 0;                           // Is intra-copied?
  const int is_high_bitdepth = is_frame_high_bitdepth(frame2);

  const int w = 2;
  const int h = 2;
  const int y = y_coord;
  const int x = x_coord;
  const int subsampling_x = 0, subsampling_y=0;
  const int_interpfilters interp_filters =
      av1_broadcast_interp_filter(MULTITAP_SHARP);
  uint8_t pred[4];
  //uint8_t *pred = malloc(w*h*sizeof(uint8_t));
  const int plane = 0;
  const int is_y_plane = (plane == 0);             // Is Y-plane?
  const struct buf_2d ref_buf = { NULL, frame2->buffers[plane],
                                  frame2->widths[is_y_plane ? 0 : 1],
                                  frame2->heights[is_y_plane ? 0 : 1],
                                  frame2->strides[is_y_plane ? 0 : 1] };
  struct scale_factors scale;
  av1_setup_scale_factors_for_frame(
      &scale, frame->y_crop_width, frame->y_crop_height,
      frame->y_crop_width, frame->y_crop_height);
  InterPredParams inter_pred_params;
  av1_init_inter_params(&inter_pred_params, w, h, y, x, subsampling_x,
                        subsampling_y, bit_depth, is_high_bitdepth,
                        is_intrabc, &scale, &ref_buf, interp_filters);
  inter_pred_params.conv_params = get_conv_params(0, plane, bit_depth);
  MV newmv = {.row = (int)(mv->row*8), .col = (int)(mv->col*8)};
  av1_enc_build_one_inter_predictor(&pred,w, &newmv, &inter_pred_params);
  
  double values1[4], values2[4];
  int idx = 0;
  /*int status = get_subpixels(frame2,&pred,w,h,*mv, x_coord, y_coord);
  if(status == -1){
    *derivative = 0;
    mv->row = 0;
    mv->col = 0;
    return;
  }*/
  for (int yy = -1; yy < 1; yy++) {
    for (int xx = -1; xx < 1; xx++) {
      values1[idx] =
          frame->y_buffer[(y_coord + yy) * frame->y_stride + (x_coord + xx)] /
          255.0;
      //repeat line above for values2 with frame2 for previous version of grad.
      values2[idx] = pred[idx]/255.0;
      idx++;
    }
  }
  double d1 = 0, d2 = 0;
  double filter[4] = { 0.25,0.25,0.25,0.25 };
  //double filter[4] = { 1, 1, 1, 1 };
  for (int i = 0; i < 4; i++) {
    d1 += filter[i] * values1[i];
    d2 += filter[i] * values2[i];
  }
  *derivative = d1 - d2;
}
// Numerical differentiate over window_size x window_size surrounding (x,y)
// location. Alters ix, iy, it to contain numerical partial derivatives
// TODO(lpartin) if LK is computed over all pixels instead of on corners,
// it would be more efficient to compute gradients on entire image (instead of
// window)
void gradients_over_window(const YV12_BUFFER_CONFIG *frame,
                           const YV12_BUFFER_CONFIG *ref_frame,
                           const int x_coord, const int y_coord,
                           const int window_size, double *ix, double *iy,
                           double *it, LOCALMV *mv, AV1_COMP *cpi) {
                           //double *it, MV mv, AV1_COMP *cpi) {
  // gradient operators need pixel before and after (start at 1)
  const int x_start =
      1 > x_coord - window_size / 2 ? 1 : x_coord - window_size / 2;
  const int y_start =
      1 > y_coord - window_size / 2 ? 1 : y_coord - window_size / 2;
  const int frame_height = frame->y_crop_height;
  const int frame_width = frame->y_crop_width;
  double deriv_x;
  double deriv_y;
  double deriv_t;

  int endx = x_start + window_size < frame_width - 1 ? x_start + window_size
                                                     : frame_width - 1;
  int endy = y_start + window_size < frame_height - 1 ? y_start + window_size
                                                      : frame_height - 1;
  // compute numerical differentiation for every pixel in window
  for (int j = y_start; j < endy; j++) {
    for (int i = x_start; i < endx; i++) {
      temporal_gradient(frame, ref_frame, i, j, &deriv_t, mv, cpi);
      spatial_gradient(frame, i, j, 0, &deriv_x);
      spatial_gradient(frame, i, j, 1, &deriv_y);
      ix[(j - y_start) * window_size + (i - x_start)] = deriv_x;
      iy[(j - y_start) * window_size + (i - x_start)] = deriv_y;
      it[(j - y_start) * window_size + (i - x_start)] = deriv_t;
    }
  }
}
// To compute eigenvalues of 2x2 matrix: Solve for lambda where
// Determinant of matrix - lambda*identity == 0
void eigenvalues_2x2(const double *matrix, double *eig) {
  double a = 1;
  double b = -1 * matrix[0] - matrix[3];
  double c = -matrix[1] * matrix[2] + matrix[1] * matrix[3];
  // quadratic formula
  double discriminant = b * b - 4 * a * c;

  *(eig++) = (-b - sqrt(discriminant)) / (2 * a);
  *(eig) = (-b + sqrt(discriminant)) / (2 * a);
}
void lucas_kanade(const YV12_BUFFER_CONFIG *frame_to_filter,
             const YV12_BUFFER_CONFIG *ref_frame, const int frame_idx,
             const int ref_frame_idx, LOCALMV *mvs, AV1_COMP *cpi,
             const int level, const int max_level,
             const int num_ref_corners, int *ref_corners, 
             const int highres_frame_width, const int highres_frame_height) {
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
#if PRINT_STATEMENTS == 1
  //differential_test(frame_to_filter, ref_frame);
  printf("frame height and frame width %d %d\n", frame_height, frame_width);
  printf("corner count %d\n", num_ref_corners);
#endif
  const int expand_multiplier = pow(2, level);
#if SAVE_IMGS == 1
  int fullimg_size = frame_height*expand_multiplier*frame_width*expand_multiplier;
  int fullimg_size = frame_height*frame_width;
  int *cornersimg = (int*)calloc(fullimg_size,sizeof(int));
  int *img = (int*)calloc(fullimg_size,sizeof(int));
  int *imgx = (int*)calloc(fullimg_size,sizeof(int));
  int *imgy = (int*)calloc(fullimg_size,sizeof(int));
#endif
  //algorithm is sensitive to window size:
  const int n = 15;
  double *i_x = (double*)malloc(n*n*sizeof(double));
  double *i_y= (double*)malloc(n*n*sizeof(double));
  double *i_t = (double*)malloc(n*n*sizeof(double));
  int repeatedcorners=0;
  for (int i = 0; i < num_ref_corners; i++) {
    const int x_coord = ref_corners[i*2]/expand_multiplier;
    const int y_coord = ref_corners[i*2+1]/expand_multiplier;
    int highres_y = ref_corners[i*2+1];
    int highres_x = ref_corners[i*2];
    int mv_idx = highres_y*(highres_frame_width) + highres_x;
    LOCALMV mv_old = mvs[mv_idx];
    mv_old.row = mv_old.row/expand_multiplier;
    mv_old.col = mv_old.col/expand_multiplier;
    gradients_over_window2(frame_to_filter, ref_frame, x_coord,
                               y_coord, n, i_x, i_y, i_t, &mv_old, cpi);
    double Mres1[1] = {0}, Mres2[1] ={0}, Mres3[1] = {0};
    double bres1[1] = {0}, bres2[1] = {0};
    multiply_mat(i_x, i_x, Mres1, 1, n*n, 1);
    multiply_mat(i_x, i_y, Mres2, 1, n*n, 1);
    multiply_mat(i_y, i_y, Mres3, 1, n*n, 1);

    multiply_mat(i_x, i_t, bres1, 1, n*n, 1);
    multiply_mat(i_y, i_t, bres2, 1, n*n, 1);

    double M[4] = { Mres1[0], Mres2[0], Mres2[0], Mres3[0] };
    double b[2] = { -1 * bres1[0], -1 * bres2[0] };
    double u[2];
    linsolve(2, M, 2, b, u);
#if PRINT_STATEMENTS == 1
    int to_match = 10;//ref_corners[i*2];
    if(x_coord == to_match && y_coord == to_match){
      for (int j = 0; j < 3/*n*n*/; j++){
        printf("\ndx=%lf ", i_x[j]);
        printf("dy=%lf ", i_y[j]);
        printf("dt=%lf", i_t[j]);
      }
      printf("\nM matrix %lf, %lf, %lf, %lf", M[0], M[1], M[2], M[3]);
      printf("\nb matrix %lf, %lf\n", b[0], b[1]);
      printf(" MV: ");
      for (int j = 0; j < 2; j++) {
        printf("%lf,", (int)u[j]);
      }
      printf("\n");
    }
#endif
    int magnitude = (int)sqrt(u[0]*u[0] + u[1]*u[1]);
    int img_idx = y_coord*frame_width+ x_coord;
    double eig[2] = {1,1};
    eigenvalues_2x22(M, eig);
    double threshold = 0.01;
    if (eig[0] < threshold) {
      int mult = 1;
      if(level !=0)
        mult = expand_multiplier; //mv doubles when resolution doubles
      LOCALMV mv = { .row = (mult*(u[0] + mv_old.row)), .col = (mult*(u[1]+mv_old.col)) };
#if PRINT_STATEMENTS == 1 
      if(abs(u[0]) >= 1 || abs(u[1]) >= 1){
        printf("x,y %d %d ", x_coord, y_coord);
        printf("old mv %lf %lf ", mv_old.row, mv_old.col);
        printf("new mv %lf %lf & MV %lf %lf\n", mv.row, mv.col, u[0], u[1]);
      }
#endif
      mvs[mv_idx] = mv;
      img_idx = mv_idx;
#if SAVE_IMGS == 1
      if(cornersimg[img_idx] == 1)
        repeatedcorners++;
      cornersimg[img_idx] = 1;
#endif
    }
  }
#if PRINT_STATEMENTS == 1
  printf("repeated corners %d\n", repeatedcorners);
#endif
#if SAVE_IMGS == 1
  if(frame_idx == 0){
    for(int i = 0; i < frame_height*frame_width*pow(2,level)*pow(2,level); i++){
      int dist = sqrt(pow(mvs[i].row*8,2) + pow(mvs[i].col*8,2));
        img[i] = dist;
        imgx[i] = abs(mvs[i].row*8);
        imgy[i] = abs(mvs[i].col*8);
    }
    char l[10];
    sprintf(l, "%d", level);
    char filename[100] = "corners";
    strcat(filename, l);
    int fh = frame_height * expand_multiplier;
    int fw = frame_width * expand_multiplier;
    save_img(cornersimg, fh, fw, frame_idx, filename);
    strcpy(filename, "motionfield");
    strcat(filename, l);
    save_img(img, fh, fw, frame_idx, filename);
    strcpy(filename, "motionfieldx");
    strcat(filename, l);
    save_img(imgx, fh, fw, frame_idx, filename);
    strcpy(filename, "motionfieldy");
    strcat(filename, l);
    save_img(imgy, fh, fw, frame_idx, filename);
  }
  free(img); free(imgx); free(imgy);
#endif
  free(i_t);
  free(i_x); 
  free(i_y); 
}

double convolve(double *filter, int *img, int size){
  double result = 0;
  for(int i = 0; i < size; i++){
    result += filter[i]*img[i];
  }
  return result;
}

void reduce(uint8_t *img, int height, int width, int stride, uint8_t *reduced_img){
  int new_width = width/2;
  int new_height = height/2;
  int window_size = 3;
  double gaussian_filter[9] = {1./16, 1.0/8, 1./16, 1./8, 1./4, 1./8, 1./16, 1./8, 1./16};
  //filter is 3x3 so need prev and forward 1 pixel
  int img_section[9];
  for(int y = 0; y < height-1; y+=2){  //without dummy pixels (y=2; y < height -2; y+=2)
    for(int x = 0; x < width-1; x+=2){
      int i = 0;
      for(int yy = y-window_size/2; yy <= y+window_size/2; yy++){
        for(int xx = x-window_size/2; xx <= x+window_size/2; xx++){
          int yvalue = yy;
          int xvalue = xx;
          //dummy pixels outside the boundary
          //only valid for 3x3 window 
          if(yvalue < 0)
            yvalue = 0;
          if(xvalue < 0)
            xvalue = 0;
          if(yvalue >= height)
            yvalue -= 1;
          if(xvalue >= width)
            xvalue -= 1;
          img_section[i++] = img[yvalue*stride + xvalue];
        }
      }
      reduced_img[(y/2)*new_width + (x/2)] = (uint8_t)convolve(gaussian_filter,img_section, pow(window_size,2));
    }
  }
}

void optflow(AV1_COMP *cpi, const YV12_BUFFER_CONFIG *frame_to_filter,
             const YV12_BUFFER_CONFIG *ref_frame, const int frame_idx, const int ref_frame_idx, MV *mvs) {
  int levels = 3; //use levels = 1 for non-iterative version
  uint8_t *images1[3];
  uint8_t *images2[3];
  images1[0] = frame_to_filter->y_buffer;
  images2[0] = ref_frame->y_buffer;
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
  LOCALMV *localmvs = malloc(frame_height*frame_width*sizeof(LOCALMV));
  //for initializing mvs to zero instead of block MVs
  /*for(int i = 0; i < frame_width*frame_height; i++){
      MV mv = {.row = 0, .col = 0};
      mvs[i] = mv;
  }*/
  YV12_BUFFER_CONFIG *buffers1 = malloc(levels*sizeof(YV12_BUFFER_CONFIG)); 
  YV12_BUFFER_CONFIG *buffers2 = malloc(levels*sizeof(YV12_BUFFER_CONFIG)); 
  buffers1[0] = *frame_to_filter;
  buffers2[0] = *ref_frame;
  int fw = frame_width;
  int fh = frame_height;
  for(int i = 1; i < levels; i++){
    images1[i] = (uint8_t *)calloc(fh/2*fw/2,sizeof(uint8_t));
    images2[i] = (uint8_t *)calloc(fh/2*fw/2,sizeof(uint8_t));
    int stride;// = buffers[i-1]->y_stride;
    if(i == 1)
      stride = frame_to_filter->y_stride;
    else
      stride = fw;
    reduce(images1[i-1], fh, fw, stride, images1[i]);
    reduce(images2[i-1], fh, fw, stride, images2[i]);
    fh /= 2;
    fw /= 2;
#if SAVE_IMGS == 1 
    int *reducedimg = (int *)calloc(fh*fw,sizeof(int));
    int *reducedimg2 = (int *)calloc(fh*fw,sizeof(int));
    for(int y =0; y < fh; y++){
      for(int x = 0; x < fw; x++){
        reducedimg[y*fw + x] = images1[i][y*fw+x];
        reducedimg2[y*fw + x] = images2[i][y*fw+x];
      }
    }
    save_img(reducedimg, fh, fw, i, "reducedimg1_");
    save_img(reducedimg2, fh, fw, i, "reducedimg2_");
#endif
    YV12_BUFFER_CONFIG a = {.y_buffer = images1[i], .y_crop_width = fw,
                  .y_crop_height = fh, .y_stride=fw};
    YV12_BUFFER_CONFIG b = {.y_buffer = images2[i], .y_crop_width = fw,
                  .y_crop_height = fh, .y_stride=fw};
    buffers1[i] = a;
    buffers2[i] = b;
  }
  //Compute corners for specific frame
  /*int maxcorners = frame_to_filter->y_crop_width*frame_to_filter->y_crop_height; //MAX_CORNERS
  int *ref_corners = malloc(maxcorners*2*sizeof(int));
  int increment = 7;//pow(2,level);
  int num_ref_corners = all_points_as_corners(ref_corners, maxcorners,
                       frame_to_filter->y_crop_height,frame_to_filter->y_crop_width, increment);
  */
  int ref_corners[2 * MAX_CORNERS];
  int num_ref_corners = av1_fast_corner_detect(
      frame_to_filter->y_buffer, frame_to_filter->y_width,
      frame_to_filter->y_height, frame_to_filter->y_stride, ref_corners,
      MAX_CORNERS);
  
  //num_ref_corners = corners_from_file(ref_corners, maxcorners,
  //                   frame_height,frame_width, increment);
  printf("CORNERS COUNT %d\n", num_ref_corners); 
  for(int i = 0; i < frame_width*frame_height; i++){
      MV mv = mvs[i];
      LOCALMV localmv = {.row = mv.row/8.0, .col = mv.col/8.0};
      localmvs[i] = localmv;
  }
  for(int i = levels-1; i >= 0; i--){
    lucas_kanade(&buffers1[i], &buffers2[i], frame_idx, ref_frame_idx, localmvs,cpi, i, levels, num_ref_corners, ref_corners,
                 buffers1[0].y_crop_width, buffers1[0].y_crop_height);
  }
  for(int i = 0; i < frame_to_filter->y_crop_height*frame_to_filter->y_crop_width; i++){
    MV mv = {.row = round(8*localmvs[i].row), .col = round(8*localmvs[i].col)}; 
    mvs[i] = mv;
  }
 // free(ref_corners);
}
#endif
void initialize_mvs(AV1_COMP *cpi, YV12_BUFFER_CONFIG **frames, const int num_frames,
                    const int filter_frame_idx, const int mb_rows, const int mb_cols, 
                    const int mb_height, const int mb_width, 
                    const int mi_h, const int mi_w, const BLOCK_SIZE block_size,
                    MV *mvs){
  const YV12_BUFFER_CONFIG *const frame_to_filter = frames[filter_frame_idx];
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
  MACROBLOCK *const mb = &cpi->td.mb;
  MACROBLOCKD *const mbd = &mb->e_mbd;
  for (int mb_row = 0; mb_row < mb_rows; mb_row++) {
    av1_set_mv_row_limits(&cpi->common.mi_params, &mb->mv_limits,
                          (mb_row << mi_h), (mb_height >> MI_SIZE_LOG2),
                          cpi->oxcf.border_in_pixels);
    for (int mb_col = 0; mb_col < mb_cols; mb_col++) {
      av1_set_mv_col_limits(&cpi->common.mi_params, &mb->mv_limits,
                            (mb_col << mi_w), (mb_width >> MI_SIZE_LOG2),
                            cpi->oxcf.border_in_pixels);
      MV ref_mv = kZeroMv;  // Reference motion vector passed down along frames.
      // Perform temporal filtering frame by frame.
      for (int frame = 0; frame < num_frames; frame++) {
        if (frames[frame] == NULL) continue;

        // Motion search.
        MV subblock_mvs[4] = { kZeroMv, kZeroMv, kZeroMv, kZeroMv };
        int subblock_mses[4] = { INT_MAX, INT_MAX, INT_MAX, INT_MAX };
        if (frame == filter_frame_idx) {  // Frame to be filtered.
          // Change ref_mv sign for following frames.
          ref_mv.row *= -1;
          ref_mv.col *= -1;
        } else {  // Other reference frames.
          tf_motion_search(cpi, frame_to_filter, frames[frame], block_size,
                           mb_row, mb_col, &ref_mv, subblock_mvs,
                           subblock_mses);
        }
        const int y_height = mb_height >> mbd->plane[0].subsampling_y;
        const int y_width = mb_width >> mbd->plane[0].subsampling_x;
        const int source_y_stride = frame_to_filter->y_stride;
        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 2; j++) {
            int subblock_idx = i * 2 + j;
            const MV mv = subblock_mvs[subblock_idx];
            for(int yy = 0; yy < 0.5*mb_height; yy++){
              for(int xx = 0; xx < 0.5*mb_width; xx++){
                const int img_idx =
                    mb_row * mb_height * frame_width + mb_col * mb_width +
                    (i * 1.0 / 2 * mb_height + yy)* frame_width + (j * 1.0 / 2 * mb_width + xx);
                if(img_idx < frame_width*frame_height)
                  mvs[img_idx + frame * frame_height * frame_width] = mv;
              }
            }
          }
        }
      }
    }
  }
}

// Does temporal filter for a particular frame.
// Inputs:
//   cpi: Pointer to the composed information of input video.
//   frames: Frame buffers used for temporal filtering.
//   num_frames: Number of frames in the frame buffer.
//   filter_frame_idx: Index of the frame to be filtered.
//   is_key_frame: Whether the to-filter is a key frame.
//   block_size: Block size used for temporal filtering.
//   scale: Scaling factor.
//   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
//                 with each plane (in Y, U, V order).
// Returns:
//   Difference between filtered frame and the original frame.
static FRAME_DIFF tf_do_filtering(AV1_COMP *cpi, YV12_BUFFER_CONFIG **frames,
                                  const int num_frames,
                                  const int filter_frame_idx,
                                  const int is_key_frame,
                                  const BLOCK_SIZE block_size,
                                  const struct scale_factors *scale,
                                  const double *noise_levels) {
#if PRINT_STATEMENTS == 1
  printf("\nfiltering frame # %d with %d frames\n", filter_frame_idx, num_frames);
#endif
  // Basic information.
  const YV12_BUFFER_CONFIG *const frame_to_filter = frames[filter_frame_idx];
  const int frame_height = frame_to_filter->y_crop_height;
  const int frame_width = frame_to_filter->y_crop_width;
  const int mb_height = block_size_high[block_size];
  const int mb_width = block_size_wide[block_size];
  const int mb_pels = mb_height * mb_width;
  const int mb_rows = get_num_blocks(frame_height, mb_height);
  const int mb_cols = get_num_blocks(frame_width, mb_width);
  const int num_planes = av1_num_planes(&cpi->common);
  const int mi_h = mi_size_high_log2[block_size];
  const int mi_w = mi_size_wide_log2[block_size];
  assert(num_planes >= 1 && num_planes <= MAX_MB_PLANE);
  const int is_high_bitdepth = is_frame_high_bitdepth(frame_to_filter);
  // Quantization factor used in temporal filtering.
  const int q_factor = get_q(cpi);
  // Factor to control the filering strength.
  const int filter_strength = cpi->oxcf.arnr_strength;

  // Save input state.
  MACROBLOCK *const mb = &cpi->td.mb;
  MACROBLOCKD *const mbd = &mb->e_mbd;
  uint8_t *input_buffer[MAX_MB_PLANE];
  for (int i = 0; i < num_planes; i++) {
    input_buffer[i] = mbd->plane[i].pre[0].buf;
  }
  MB_MODE_INFO **input_mb_mode_info = mbd->mi;

  // Determine whether the video is with `YUV 4:2:2` format, since the avx2/sse2
  // function only supports square block size. We will use C function instead
  // for videos with `YUV 4:2:2` format.
  int is_yuv422_format = 0;
  for (int plane = 1; plane < num_planes; ++plane) {
    if (mbd->plane[plane].subsampling_x != mbd->plane[plane].subsampling_y) {
      is_yuv422_format = 1;
      break;
    }
  }

  // Setup.
  mbd->block_ref_scale_factors[0] = scale;
  mbd->block_ref_scale_factors[1] = scale;
  // A temporary block info used to store state in temporal filtering process.
  MB_MODE_INFO *tmp_mb_mode_info = (MB_MODE_INFO *)malloc(sizeof(MB_MODE_INFO));
  memset(tmp_mb_mode_info, 0, sizeof(MB_MODE_INFO));
  mbd->mi = &tmp_mb_mode_info;
  mbd->mi[0]->motion_mode = SIMPLE_TRANSLATION;
  // Allocate memory for predictor, accumulator and count.
  uint8_t *pred8 = aom_memalign(32, num_planes * mb_pels * sizeof(uint8_t));
  uint16_t *pred16 = aom_memalign(32, num_planes * mb_pels * sizeof(uint16_t));
  uint32_t *accum = aom_memalign(16, num_planes * mb_pels * sizeof(uint32_t));
  uint16_t *count = aom_memalign(16, num_planes * mb_pels * sizeof(uint16_t));
  memset(pred8, 0, num_planes * mb_pels * sizeof(pred8[0]));
  memset(pred16, 0, num_planes * mb_pels * sizeof(pred16[0]));
  uint8_t *const pred = is_high_bitdepth ? CONVERT_TO_BYTEPTR(pred16) : pred8;
  // Do filtering.
  FRAME_DIFF diff = { 0, 0 };
#if TF_DEBUG_MODE > 0
  int *img = calloc(frame_height * frame_width, sizeof(int));
  int *imgrow = calloc(frame_height * frame_width, sizeof(int));
  int *imgcol = calloc(frame_height * frame_width, sizeof(int));
  int countnewmvs = 0;
  printf("\nPixel mvs for frame #%d with %d frames to compare\n", filter_frame_idx, num_frames);
  MV *allframes_pixelmvs =
      malloc(num_frames * frame_height * frame_width * sizeof(MV));
  initialize_mvs(cpi,frames, num_frames,filter_frame_idx, mb_rows, mb_cols,
                    mb_height, mb_width, mi_h, mi_w, block_size,allframes_pixelmvs);
  printf("size of allframes %d and last idx %d\n", num_frames*frame_height*frame_width, 
         1*frame_height*frame_width+(frame_width*frame_height-1)); 
  for(int i = 0; i < frame_width*frame_height && num_frames >= 2; i++){
      MV mv = allframes_pixelmvs[1*frame_height*frame_width + i];
      const double distance = sqrt(pow(mv.row, 2) + pow(mv.col, 2));
      img[i] = distance;
      imgrow[i] = abs(mv.row);
      imgcol[i] = abs(mv.col);
  }
#endif
#if TF_DEBUG_MODE > 1
  for (int i = 0; i < num_frames; i++){
    optflow(cpi, frames[filter_frame_idx], frames[i], filter_frame_idx, i,
            (allframes_pixelmvs + i * frame_height * frame_width));
  }
  for(int i = 0; i < frame_width*frame_height && num_frames >= 2; i++){
      MV mv = allframes_pixelmvs[1*frame_height*frame_width + i];
      if(imgrow[i] != abs(mv.row) || imgcol[i] != abs(mv.col))
        countnewmvs++;
  }
  printf("NEW MVS %d\n", countnewmvs);
#else
  printf("\nBlock mvs for frame #%d with %d frames to compare\n", filter_frame_idx, num_frames);
#endif
  // Perform temporal filtering block by block.
  for (int mb_row = 0; mb_row < mb_rows; mb_row++) {
    av1_set_mv_row_limits(&cpi->common.mi_params, &mb->mv_limits,
                          (mb_row << mi_h), (mb_height >> MI_SIZE_LOG2),
                          cpi->oxcf.border_in_pixels);
    for (int mb_col = 0; mb_col < mb_cols; mb_col++) {
      av1_set_mv_col_limits(&cpi->common.mi_params, &mb->mv_limits,
                            (mb_col << mi_w), (mb_width >> MI_SIZE_LOG2),
                            cpi->oxcf.border_in_pixels);
      memset(accum, 0, num_planes * mb_pels * sizeof(accum[0]));
      memset(count, 0, num_planes * mb_pels * sizeof(count[0]));
      MV ref_mv = kZeroMv;  // Reference motion vector passed down along frames.
      // Perform temporal filtering frame by frame.
      for (int frame = 0; frame < num_frames; frame++) {
        if (frames[frame] == NULL) continue;

        // Motion search.
        MV subblock_mvs[4] = { kZeroMv, kZeroMv, kZeroMv, kZeroMv };
        int subblock_mses[4] = { INT_MAX, INT_MAX, INT_MAX, INT_MAX };
        if (frame == filter_frame_idx) {  // Frame to be filtered.
          // Change ref_mv sign for following frames.
          ref_mv.row *= -1;
          ref_mv.col *= -1;
        } else {  // Other reference frames.
          tf_motion_search(cpi, frame_to_filter, frames[frame], block_size,
                           mb_row, mb_col, &ref_mv, subblock_mvs,
                           subblock_mses);
        }
#if TF_DEBUG_MODE == 2
        MV *pixelmvs;
        pixelmvs = allframes_pixelmvs + frame * frame_height * frame_width;
        MV *blockpixel_mvs = malloc(mb_height * mb_width * sizeof(MV));
        for (int yy = 0; yy < mb_height; yy++) {
          for (int xx = 0; xx < mb_width; xx++) {
            MV sblock_orig;
            if (yy < mb_height / 2 && xx < mb_height / 2)
              sblock_orig = subblock_mvs[0];
            else if (yy < mb_height / 2)
              sblock_orig = subblock_mvs[1];
            else if (xx < mb_height / 2)
              sblock_orig = subblock_mvs[2];
            else
              sblock_orig = subblock_mvs[3];
            const int pixelidx = mb_row * mb_height * frame_width +
                                 mb_col * mb_width + yy * frame_width + xx;
            MV pixelmv = pixelmvs[pixelidx];
            double magn = sqrt(pow(pixelmv.row, 2) + pow(pixelmv.col, 2));
            blockpixel_mvs[yy * mb_width + xx] =
                fabs(magn) > 0. ? pixelmv : sblock_orig;
          }
        }
        tf_build_predictor(frames[frame], mbd, block_size, mb_row, mb_col,
                           num_planes, scale, blockpixel_mvs, pred,
                           frame_to_filter,subblock_mses);
#elif TF_DEBUG_MODE == 4
        MV *pixelmvs;
        pixelmvs = allframes_pixelmvs + frame * frame_height * frame_width;
        MV *blockpixel_mvs = malloc(mb_height * mb_width * sizeof(MV));
        for (int yy = 0; yy < mb_height; yy++) {
          for (int xx = 0; xx < mb_width; xx++) {
            const int pixelidx = mb_row * mb_height * frame_width +
                                 mb_col * mb_width + yy * frame_width + xx;
            MV pixelmv = pixelmvs[pixelidx];
            blockpixel_mvs[yy * mb_width + xx] = pixelmv;
            double magn = sqrt(pow(pixelmv.row, 2) + pow(pixelmv.col, 2));
          }
        }
        int *pixel_mses = malloc(mb_height*mb_width*sizeof(int));
        tf_build_predictor(frames[frame], mbd, block_size, mb_row, mb_col,
                           num_planes, scale, blockpixel_mvs, pred,
                           frame_to_filter,pixel_mses);
#else
        tf_build_predictor(frames[frame], mbd, block_size, mb_row, mb_col,
                           num_planes, scale, subblock_mvs, pred,
                           frame_to_filter,subblock_mses);
#endif
        // Perform weighted averaging.
        if (frame == filter_frame_idx) {  // Frame to be filtered.
          tf_apply_temporal_filter_self(mbd, block_size, num_planes, pred,
                                        accum, count);
        } else {  // Other reference frames.
          // TODO(any): avx2/sse2 version should be changed to align with C
          // function before using. In particular, current avx2/sse2 function
          // only supports 32x32 block size, 5x5 filtering window, 8-bit
          // encoding, and the case when the video is not with `YUV 4:2:2`
          // format.
          if (TF_BLOCK_SIZE == BLOCK_32X32 && TF_WINDOW_LENGTH == 5 &&
              !is_frame_high_bitdepth(frame_to_filter) && !is_yuv422_format) {
#if TF_DEBUG_MODE != 4
            av1_apply_temporal_filter(frame_to_filter, mbd, block_size, mb_row,
                                      mb_col, num_planes, noise_levels,
                                      subblock_mvs, subblock_mses, q_factor,
                                      filter_strength, pred, accum, count);
#else
            av1_apply_temporal_filter(frame_to_filter, mbd, block_size, mb_row,
                                      mb_col, num_planes, noise_levels,
                                      blockpixel_mvs, pixel_mses, q_factor,
                                      filter_strength, pred, accum, count);
#endif
          } else {
            av1_apply_temporal_filter_c(
                frame_to_filter, mbd, block_size, mb_row, mb_col, num_planes,
                noise_levels, subblock_mvs, subblock_mses, q_factor,
                filter_strength, pred, accum, count);
          }
        }
      }

      tf_normalize_filtered_frame(mbd, block_size, mb_row, mb_col, num_planes,
                                  accum, count, &cpi->alt_ref_buffer);

      if (!is_key_frame && cpi->sf.hl_sf.adaptive_overlay_encoding) {
        const int y_height = mb_height >> mbd->plane[0].subsampling_y;
        const int y_width = mb_width >> mbd->plane[0].subsampling_x;
        const int source_y_stride = frame_to_filter->y_stride;
        const int filter_y_stride = cpi->alt_ref_buffer.y_stride;
        const int source_offset =
            mb_row * y_height * source_y_stride + mb_col * y_width;
        const int filter_offset =
            mb_row * y_height * filter_y_stride + mb_col * y_width;
        unsigned int sse = 0;
        cpi->fn_ptr[block_size].vf(frame_to_filter->y_buffer + source_offset,
                                   source_y_stride,
                                   cpi->alt_ref_buffer.y_buffer + filter_offset,
                                   filter_y_stride, &sse);
        diff.sum += sse;
        diff.sse += sse * sse;
      }
    }
  }
#if TF_DEBUG_MODE > 0
#if SAVE_IMGS == 1
    if(filter_frame_idx == 0){
      save_img(img, frame_height, frame_width, filter_frame_idx, "subblockmvs");
      save_img(imgrow, frame_height, frame_width, filter_frame_idx,"subblockmvs_x"); 
      save_img(imgcol, frame_height, frame_width, filter_frame_idx, "subblockmvs_y");

#if TF_DEBUG_MODE == 1 
      //int frame_to_compare = 1;
      //printf("Trying out opt flow\n");
      //MV *pixelmvs = malloc(429*417*sizeof(MV));
      //MV *pixelmvs = malloc(frame_height*frame_width*sizeof(MV));
      //optflow_test(cpi, frames[filter_frame_idx], frames[frame_to_compare],filter_frame_idx, frame_to_compare, pixelmvs);
#endif
      int *tf_diff = calloc(frame_width*frame_height, sizeof(int));
      for (int yy = 0; yy < frame_height; yy++){
        for(int xx = 0; xx < frame_width; xx++){
          tf_diff[yy*frame_width + xx] =
            fabs(frame_to_filter->y_buffer[yy*frame_to_filter->y_stride + xx] -
                  cpi->alt_ref_buffer.y_buffer[yy*cpi->alt_ref_buffer.y_stride + xx]);
        }
      }
      save_img(tf_diff, frame_height, frame_width, filter_frame_idx, "tf_diff");
      free(tf_diff);
      dump_frame(cpi, frame_to_filter, filter_frame_idx, "_beforetf");
      dump_frame(cpi, &(cpi->alt_ref_buffer), filter_frame_idx, "_aftertf");
      dump_frame(cpi, frames[1], 1, "_beforetf");
    }
#endif
  free(img);
  free(imgrow);
  free(imgcol);
#endif

  // Restore input state
  for (int i = 0; i < num_planes; i++) {
    mbd->plane[i].pre[0].buf = input_buffer[i];
  }
  mbd->mi = input_mb_mode_info;

  free(tmp_mb_mode_info);
  aom_free(pred8);
  aom_free(pred16);
  aom_free(accum);
  aom_free(count);

  return diff;
}

// Setups the frame buffer for temporal filtering. Basically, this fuction
// determines how many frames will be used for temporal filtering and then
// groups them into a buffer. This function will also estimate the noise level
// of the to-filter frame.
// Inputs:
//   cpi: Pointer to the composed information of input video.
//   filter_frame_lookahead_idx: The index of the to-filter frame in the
//                               lookahead buffer `cpi->lookahead`.
//   is_second_arf: Whether the to-filter frame is the second ARF. This field
//                  will affect the number of frames used for filtering.
//   frames: Pointer to the frame buffer to setup.
//   num_frames_for_filtering: Number of frames used for filtering.
//   filter_frame_idx: Index of the to-filter frame in the setup frame buffer.
//   noise_levels: Pointer to the noise levels of the to-filter frame, estimated
//                 with each plane (in Y, U, V order).
// Returns:
//   Nothing will be returned. But the frame buffer `frames`, number of frames
//   in the buffer `num_frames_for_filtering`, and the index of the to-filter
//   frame in the buffer `filter_frame_idx` will be updated in this function.
//   Estimated noise levels for YUV planes will be saved in `noise_levels`.
static void tf_setup_filtering_buffer(const AV1_COMP *cpi,
                                      const int filter_frame_lookahead_idx,
                                      const int is_second_arf,
                                      YV12_BUFFER_CONFIG **frames,
                                      int *num_frames_for_filtering,
                                      int *filter_frame_idx,
                                      double *noise_levels) {
  // Number of frames used for filtering. Set `arnr_max_frames` as 1 to disable
  // temporal filtering.
  int num_frames = AOMMAX(cpi->oxcf.arnr_max_frames, 1);
  int num_before = 0;  // Number of filtering frames before the to-filter frame.
  int num_after = 0;   // Number of filtering frames after the to-filer frame.
  const int lookahead_depth =
      av1_lookahead_depth(cpi->lookahead, cpi->compressor_stage);
  // Number of buffered frames before the to-filter frame.
  const int max_before = filter_frame_lookahead_idx < -1
                             ? -filter_frame_lookahead_idx + 1
                             : filter_frame_lookahead_idx + 1;
  // Number of buffered frames after the to-filter frame.
  const int max_after = lookahead_depth - max_before;

  const int filter_frame_offset = filter_frame_lookahead_idx < -1
                                      ? -filter_frame_lookahead_idx
                                      : filter_frame_lookahead_idx;

  // Estimate noises for each plane.
  const struct lookahead_entry *to_filter_buf = av1_lookahead_peek(
      cpi->lookahead, filter_frame_offset, cpi->compressor_stage);
  assert(to_filter_buf != NULL);
  const YV12_BUFFER_CONFIG *to_filter_frame = &to_filter_buf->img;
  const int num_planes = av1_num_planes(&cpi->common);
  for (int plane = 0; plane < num_planes; ++plane) {
    noise_levels[plane] = av1_estimate_noise_from_single_plane(
        to_filter_frame, plane, cpi->common.seq_params.bit_depth);
  }
  // Get quantization factor.
  const int q = get_q(cpi);

  // Adjust number of filtering frames based on noise and quantization factor.
  // Basically, we would like to use more frames to filter low-noise frame such
  // that the filtered frame can provide better predictions for more frames.
  // Also, when the quantization factor is small enough (lossless compression),
  // we will not change the number of frames for key frame filtering, which is
  // to avoid visual quality drop.
  int adjust_num = 0;
  if (num_frames == 1) {  // `arnr_max_frames = 1` is used to disable filtering.
    adjust_num = 0;
  } else if (filter_frame_lookahead_idx < 0 && q <= 10) {
    adjust_num = 0;
  } else if (noise_levels[0] < 0.5) {
    adjust_num = 6;
  } else if (noise_levels[0] < 1.0) {
    adjust_num = 4;
  } else if (noise_levels[0] < 2.0) {
    adjust_num = 2;
  }
  num_frames = AOMMIN(num_frames + adjust_num, lookahead_depth + 1);

  if (filter_frame_lookahead_idx == -1) {  // Key frame.
    num_before = 0;
    num_after = AOMMIN(num_frames - 1, max_after);
  } else if (filter_frame_lookahead_idx < -1) {  // Key frame in one-pass mode.
    num_before = AOMMIN(num_frames - 1, max_before);
    num_after = 0;
  } else {
    num_frames = AOMMIN(num_frames, cpi->rc.gfu_boost / 150);
    num_frames += !(num_frames & 1);  // Make the number odd.
    // Only use 2 neighbours for the second ARF.
    if (is_second_arf) num_frames = AOMMIN(num_frames, 3);
    num_before = AOMMIN(num_frames >> 1, max_before);
    num_after = AOMMIN(num_frames >> 1, max_after);
  }
  num_frames = num_before + 1 + num_after;

  // Setup the frame buffer.
  for (int frame = 0; frame < num_frames; ++frame) {
    const int lookahead_idx = frame - num_before + filter_frame_offset;
    struct lookahead_entry *buf = av1_lookahead_peek(
        cpi->lookahead, lookahead_idx, cpi->compressor_stage);
    assert(buf != NULL);
    frames[frame] = &buf->img;
  }
  *num_frames_for_filtering = num_frames;
  *filter_frame_idx = num_before;
  assert(frames[*filter_frame_idx] == to_filter_frame);
}

// A constant number, sqrt(pi / 2),  used for noise estimation.
static const double SQRT_PI_BY_2 = 1.25331413732;

double av1_estimate_noise_from_single_plane(const YV12_BUFFER_CONFIG *frame,
                                            const int plane,
                                            const int bit_depth) {
  const int is_y_plane = (plane == 0);
  const int height = frame->crop_heights[is_y_plane ? 0 : 1];
  const int width = frame->crop_widths[is_y_plane ? 0 : 1];
  const int stride = frame->strides[is_y_plane ? 0 : 1];
  const uint8_t *src = frame->buffers[plane];
  const uint16_t *src16 = CONVERT_TO_SHORTPTR(src);
  const int is_high_bitdepth = is_frame_high_bitdepth(frame);

  int64_t accum = 0;
  int count = 0;
  for (int i = 1; i < height - 1; ++i) {
    for (int j = 1; j < width - 1; ++j) {
      // Setup a small 3x3 matrix.
      const int center_idx = i * stride + j;
      int mat[3][3];
      for (int ii = -1; ii <= 1; ++ii) {
        for (int jj = -1; jj <= 1; ++jj) {
          const int idx = center_idx + ii * stride + jj;
          mat[ii + 1][jj + 1] = is_high_bitdepth ? src16[idx] : src[idx];
        }
      }
      // Compute sobel gradients.
      const int Gx = (mat[0][0] - mat[0][2]) + (mat[2][0] - mat[2][2]) +
                     2 * (mat[1][0] - mat[1][2]);
      const int Gy = (mat[0][0] - mat[2][0]) + (mat[0][2] - mat[2][2]) +
                     2 * (mat[0][1] - mat[2][1]);
      const int Ga = ROUND_POWER_OF_TWO(abs(Gx) + abs(Gy), bit_depth - 8);
      // Accumulate Laplacian.
      if (Ga < NOISE_ESTIMATION_EDGE_THRESHOLD) {  // Only count smooth pixels.
        const int v = 4 * mat[1][1] -
                      2 * (mat[0][1] + mat[2][1] + mat[1][0] + mat[1][2]) +
                      (mat[0][0] + mat[0][2] + mat[2][0] + mat[2][2]);
        accum += ROUND_POWER_OF_TWO(abs(v), bit_depth - 8);
        ++count;
      }
    }
  }

  // Return -1.0 (unreliable estimation) if there are too few smooth pixels.
  return (count < 16) ? -1.0 : (double)accum / (6 * count) * SQRT_PI_BY_2;
}

int av1_temporal_filter(AV1_COMP *cpi, const int filter_frame_lookahead_idx,
                        int *show_existing_arf) {
  // Basic informaton of the current frame.
  const GF_GROUP *const gf_group = &cpi->gf_group;
  const uint8_t group_idx = gf_group->index;
  const FRAME_UPDATE_TYPE update_type = gf_group->update_type[group_idx];
  COUNT = COUNT + 1;
  // Filter one more ARF if the lookahead index is leq 7 (w.r.t. 9-th frame).
  // This frame is ALWAYS a show existing frame.
  const int is_second_arf = (update_type == INTNL_ARF_UPDATE) &&
                            (filter_frame_lookahead_idx >= 7) &&
                            cpi->sf.hl_sf.second_alt_ref_filtering;
  // TODO(anyone): Currently, we enforce the filtering strength on internal
  // ARFs except the second ARF to be zero. We should investigate in which case
  // it is more beneficial to use non-zero strength filtering.
  if (update_type == INTNL_ARF_UPDATE && !is_second_arf) {
    return 0;
  }

  // TODO(yunqing): For INTNL_ARF_UPDATE type, the following me initialization
  // is used somewhere unexpectedly. Should be resolved later.
  // Initialize errorperbit and sadperbit
  const int rdmult = av1_compute_rd_mult_based_on_qindex(cpi, TF_QINDEX);
  MvCosts *mv_costs = &cpi->td.mb.mv_costs;
  av1_set_error_per_bit(mv_costs, rdmult);
  av1_set_sad_per_bit(cpi, mv_costs, TF_QINDEX);
  av1_fill_mv_costs(cpi->common.fc,
                    cpi->common.features.cur_frame_force_integer_mv,
                    cpi->common.features.allow_high_precision_mv, mv_costs);

  // Setup frame buffer for filtering.
  YV12_BUFFER_CONFIG *frames[MAX_LAG_BUFFERS] = { NULL };
  int num_frames_for_filtering = 0;
  int filter_frame_idx = -1;
  double noise_levels[MAX_MB_PLANE] = { 0 };
  tf_setup_filtering_buffer(cpi, filter_frame_lookahead_idx, is_second_arf,
                            frames, &num_frames_for_filtering,
                            &filter_frame_idx, noise_levels);
  assert(num_frames_for_filtering > 0);
  assert(filter_frame_idx < num_frames_for_filtering);

  // Set showable frame.
  if (filter_frame_lookahead_idx >= 0) {
    cpi->common.showable_frame =
        num_frames_for_filtering == 1 || is_second_arf ||
        (cpi->oxcf.enable_overlay == 0 || cpi->sf.hl_sf.disable_overlay_frames);
  }

  // Do filtering.
  const int is_key_frame = (filter_frame_lookahead_idx < 0);
  // Setup scaling factors. Scaling on each of the arnr frames is not
  // supported.
  // ARF is produced at the native frame size and resized when coded.
  struct scale_factors sf;
  av1_setup_scale_factors_for_frame(
      &sf, frames[0]->y_crop_width, frames[0]->y_crop_height,
      frames[0]->y_crop_width, frames[0]->y_crop_height);
  const FRAME_DIFF diff =
      tf_do_filtering(cpi, frames, num_frames_for_filtering, filter_frame_idx,
                      is_key_frame, TF_BLOCK_SIZE, &sf, noise_levels);

  if (is_key_frame) {  // Key frame should always be filtered.
    return 1;
  }

  if ((show_existing_arf != NULL && cpi->sf.hl_sf.adaptive_overlay_encoding) ||
      is_second_arf) {
    const int frame_height = frames[filter_frame_idx]->y_crop_height;
    const int frame_width = frames[filter_frame_idx]->y_crop_width;
    const int block_height = block_size_high[TF_BLOCK_SIZE];
    const int block_width = block_size_wide[TF_BLOCK_SIZE];
    const int mb_rows = get_num_blocks(frame_height, block_height);
    const int mb_cols = get_num_blocks(frame_width, block_width);
    const int num_mbs = AOMMAX(1, mb_rows * mb_cols);
    const float mean = (float)diff.sum / num_mbs;
    const float std = (float)sqrt((float)diff.sse / num_mbs - mean * mean);

    aom_clear_system_state();
    // TODO(yunqing): This can be combined with TPL q calculation later.
    cpi->rc.base_frame_target = gf_group->bit_allocation[group_idx];
    av1_set_target_rate(cpi, cpi->common.width, cpi->common.height);
    int top_index = 0;
    int bottom_index = 0;
    const int q = av1_rc_pick_q_and_bounds(cpi, &cpi->rc, cpi->oxcf.width,
                                           cpi->oxcf.height, group_idx,
                                           &bottom_index, &top_index);
    const int ac_q = av1_ac_quant_QTX(q, 0, cpi->common.seq_params.bit_depth);
    const float threshold = 0.7f * ac_q * ac_q;

    if (!is_second_arf) {
      *show_existing_arf = 0;
      if (mean < threshold && std < mean * 1.2) {
        *show_existing_arf = 1;
      }
      cpi->common.showable_frame |= *show_existing_arf;
    } else {
      // Use source frame if the filtered frame becomes very different.
      if (!(mean < threshold && std < mean * 1.2)) {
        return 0;
      }
    }
  }

  return 1;
}
