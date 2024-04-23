/*
 *  Copyright (c) 2015 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <limits.h>
#include <math.h>

#include "av1/common/blockd.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/skin_detection.h"

static const uint8_t b_width_log2_lookup[BLOCK_SIZES] = { 0, 0, 1, 1, 1, 2,
                                                          2, 2, 3, 3, 3, 4,
                                                          4, 4, 5, 5 };
static const uint8_t b_height_log2_lookup[BLOCK_SIZES] = { 0, 1, 0, 1, 2, 1,
                                                           2, 3, 2, 3, 4, 3,
                                                           4, 5, 4, 5 };

int av1_compute_skin_block(const uint8_t *y, const uint8_t *u, const uint8_t *v,
                           int stride, int strideuv, int bsize,
                           int consec_zeromv, int low_motion) {
  // No skin if block has been zero/small motion for long consecutive time.
  if (consec_zeromv > 30 && low_motion) {
    return 0;
  } else {
    int motion = 1;
    // Take center pixel in block to determine is_skin.
    const int y_width_shift = (4 << b_width_log2_lookup[bsize]) >> 1;
    const int y_height_shift = (4 << b_height_log2_lookup[bsize]) >> 1;
    const int uv_width_shift = y_width_shift >> 1;
    const int uv_height_shift = y_height_shift >> 1;
    const uint8_t ysource = y[y_height_shift * stride + y_width_shift];
    const uint8_t usource = u[uv_height_shift * strideuv + uv_width_shift];
    const uint8_t vsource = v[uv_height_shift * strideuv + uv_width_shift];

    if (consec_zeromv > 25 && low_motion) motion = 0;
    return aom_skin_pixel(ysource, usource, vsource, motion);
  }
}

void av1_compute_skin_sb(AV1_COMP *const cpi, BLOCK_SIZE bsize, int mi_row,
                         int mi_col, int low_motion) {
  assert(bsize == BLOCK_8X8 || bsize == BLOCK_16X16);
  int i, j, i2, j2, num_bl;
  const uint8_t *src_y = cpi->source->y_buffer;
  const uint8_t *src_u = cpi->source->u_buffer;
  const uint8_t *src_v = cpi->source->v_buffer;
  const int src_ystride = cpi->source->y_stride;
  const int src_uvstride = cpi->source->uv_stride;
  const int y_bsize = 4 << b_width_log2_lookup[bsize];
  const int uv_bsize = y_bsize >> 1;
  const int shy = (y_bsize == 8) ? 3 : 4;
  const int shuv = shy - 1;
  const int fac = y_bsize / 4;
  const int y_shift = src_ystride * (mi_row << 2) + (mi_col << 2);
  const int uv_shift = src_uvstride * (mi_row << 1) + (mi_col << 1);
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  const int mi_row_limit = AOMMIN(mi_row + 16, mi_params->mi_rows - 2);
  const int mi_col_limit = AOMMIN(mi_col + 16, mi_params->mi_cols - 2);
  src_y += y_shift;
  src_u += uv_shift;
  src_v += uv_shift;

  for (i = mi_row; i < mi_row_limit; i += fac) {
    num_bl = 0;
    i2 = i >> 1;
    for (j = mi_col; j < mi_col_limit; j += fac) {
      j2 = j >> 1;
      const int bl_index = i * mi_params->mi_cols + j;
      const int bl_index2 = i2 * (mi_params->mi_cols >> 1) + j2;
      int consec_zeromv = 0;
      if (cpi->consec_zero_mv != NULL)
        consec_zeromv = cpi->consec_zero_mv[bl_index2];
      cpi->skin_map[bl_index] =
          av1_compute_skin_block(src_y, src_u, src_v, src_ystride, src_uvstride,
                                 bsize, consec_zeromv, low_motion);
      num_bl++;
      src_y += y_bsize;
      src_u += uv_bsize;
      src_v += uv_bsize;
    }
    src_y += (src_ystride << shy) - (num_bl << shy);
    src_u += (src_uvstride << shuv) - (num_bl << shuv);
    src_v += (src_uvstride << shuv) - (num_bl << shuv);
  }

  // Remove isolated skin blocks (none of its neighbors are skin) and isolated
  // non-skin blocks (all of its neighbors are skin).
  // Skip 4 corner blocks which have only 3 neighbors to remove isolated skin
  // blocks. Skip superblock borders to remove isolated non-skin blocks.
  for (i = mi_row; i < mi_row_limit; i += fac) {
    for (j = mi_col; j < mi_col_limit; j += fac) {
      int bl_index = i * mi_params->mi_cols + j;
      int num_neighbor = 0;
      int mi, mj;
      int non_skin_threshold = 8;
      // Skip 4 corners.
      if ((i == mi_row && (j == mi_col || j == mi_col_limit - fac)) ||
          (i == mi_row_limit - fac && (j == mi_col || j == mi_col_limit - fac)))
        continue;
      // There are only 5 neighbors for non-skin blocks on the border.
      if (i == mi_row || i == mi_row_limit - fac || j == mi_col ||
          j == mi_col_limit - fac)
        non_skin_threshold = 5;

      for (mi = -fac; mi <= fac; mi += fac) {
        for (mj = -fac; mj <= fac; mj += fac) {
          if (i + mi >= mi_row && i + mi < mi_row_limit && j + mj >= mi_col &&
              j + mj < mi_col_limit) {
            int bl_neighbor_index = (i + mi) * mi_params->mi_cols + j + mj;
            if (cpi->skin_map[bl_neighbor_index]) num_neighbor++;
          }
        }
      }

      if (cpi->skin_map[bl_index] && num_neighbor < 2)
        cpi->skin_map[bl_index] = 0;
      if (!cpi->skin_map[bl_index] && num_neighbor == non_skin_threshold)
        cpi->skin_map[bl_index] = 1;
    }
  }
}

#ifdef OUTPUT_YUV_SKINMAP
// For viewing skin map on input source.
void aom_write_one_yuv_frame(AV1_COMMON *cm, YV12_BUFFER_CONFIG *s,
                             FILE *map_file) {
  uint8_t *src = s->y_buffer;
  int h = cm->height;
  if (map_file == NULL) return;
  do {
    fwrite(src, s->y_width, 1, map_file);
    src += s->y_stride;
  } while (--h);

  src = s->u_buffer;
  h = s->uv_height;

  do {
    fwrite(src, s->uv_width, 1, map_file);
    src += s->uv_stride;
  } while (--h);

  src = s->v_buffer;
  h = s->uv_height;

  do {
    fwrite(src, s->uv_width, 1, map_file);
    src += s->uv_stride;
  } while (--h);

  fflush(map_file);
}

void av1_output_skin_map(AV1_COMP *const cpi, FILE *yuv_skinmap_file) {
  int i, j, mi_row, mi_col, num_bl;
  AV1_COMMON *const cm = &cpi->common;
  const CommonModeInfoParams *const mi_params = &cpi->common.mi_params;
  uint8_t *y;
  const uint8_t *src_y = cpi->source->y_buffer;
  const int src_ystride = cpi->source->y_stride;
  const int y_bsize = 8;  // Use 8x8 or 16x16.
  const int shy = (y_bsize == 8) ? 3 : 4;
  const int fac = y_bsize / 4;
  YV12_BUFFER_CONFIG skinmap;
  memset(&skinmap, 0, sizeof(YV12_BUFFER_CONFIG));
  const SequenceHeader *seq_params = cm->seq_params;
  if (aom_alloc_frame_buffer(
          &skinmap, cm->width, cm->height, seq_params->subsampling_x,
          seq_params->subsampling_y, seq_params->use_highbitdepth,
          cpi->oxcf.border_in_pixels, cm->features.byte_alignment, false, 0)) {
    aom_internal_error(cm->error, AOM_CODEC_MEM_ERROR,
                       "Failed to allocate buffer for skin map");
    return;
  }
  memset(skinmap.buffer_alloc, 128, skinmap.frame_size);
  y = skinmap.y_buffer;
  // Loop through blocks and set skin map based on center pixel of block.
  // Set y to white for skin block, otherwise set to source with gray scale.
  // Ignore rightmost/bottom boundary blocks.
  for (mi_row = 0; mi_row < mi_params->mi_rows - 1; mi_row += fac) {
    num_bl = 0;
    for (mi_col = 0; mi_col < mi_params->mi_cols - 1; mi_col += fac) {
      const int block_index = mi_row * mi_params->mi_cols + mi_col;
      const int is_skin = cpi->skin_map[block_index];
      for (i = 0; i < y_bsize; i++) {
        for (j = 0; j < y_bsize; j++) {
          y[i * skinmap.y_stride + j] =
              is_skin ? 255 : src_y[i * src_ystride + j];
        }
      }
      num_bl++;
      y += y_bsize;
      src_y += y_bsize;
    }
    y += (skinmap.y_stride << shy) - (num_bl << shy);
    src_y += (src_ystride << shy) - (num_bl << shy);
  }
  aom_write_one_yuv_frame(cm, &skinmap, yuv_skinmap_file);
  aom_free_frame_buffer(&skinmap);
}
#endif  // OUTPUT_YUV_SKINMAP
