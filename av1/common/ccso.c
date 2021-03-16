#include <assert.h>
#include <math.h>
#include <string.h>

#include "config/aom_scale_rtcd.h"

#include "aom/aom_integer.h"
#include "av1/common/ccso.h"
#include "av1/common/reconinter.h"

void extend_ccso_border(uint16_t **buf, int d, AV1_COMMON *cm,
                        MACROBLOCKD *xd) {
  int s = xd->plane[0].dst.width;
  for (int pli = 0; pli < av1_num_planes(cm); pli++) {
    uint16_t *p = &buf[pli][d * s + d];
    int h = xd->plane[pli].dst.height;
    int w = xd->plane[pli].dst.width;
    for (int y = 0; y < h; y++) {
      for (int x = 0; x < d; x++) {
        *(p - d + x) = p[0];
        p[w + x] = p[w - 1];
      }
      p += s;
    }
    p -= (s + d);
    for (int y = 0; y < d; y++) {
      memcpy(p + (y + 1) * s, p, sizeof(uint16_t) * (w + (d << 1)));
    }
    p -= ((h - 1) * s);
    for (int y = 0; y < d; y++) {
      memcpy(p - (y + 1) * s, p, sizeof(uint16_t) * (w + (d << 1)));
    }
  }
}

void cal_filter_support(int rec_luma_idx[2], const uint16_t *rec_y,
                        uint8_t quant_step_size, int inv_quant_step,
                        int rec_idx[2]) {
  for (int i = 0; i < 2; i++) {
    int d = rec_y[rec_idx[i]] - rec_y[0];
    if (d > quant_step_size)
      rec_luma_idx[i] = 2;
    else if (d < inv_quant_step)
      rec_luma_idx[i] = 0;
    else
      rec_luma_idx[i] = 1;
  }
}

void apply_ccso_filter(AV1_COMMON *cm, MACROBLOCKD *xd, const int plane,
                       uint16_t *temp_rec_y_buf, uint8_t *dstYuv8,
                       const int dst_stride,
                       int filter_offset[CCSO_LUT_EXT_SIZE],
                       uint8_t quant_step_size, uint8_t ext_filter_support) {
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int ccso_stride = xd->plane[0].dst.width;
  const int pic_height_c = xd->plane[1].dst.height;
  const int pic_width_c = xd->plane[1].dst.width;
  int rec_luma_idx[2];
  int inv_quant_step = quant_step_size * -1;
  int rec_idx[2];
#if CCSO_FILTER_SUPPORT_EXT
  if (ext_filter_support == 0) {
    rec_idx[0] = -1 * ccso_stride;
    rec_idx[1] = 1 * ccso_stride;
  } else if (ext_filter_support == 1) {
    rec_idx[0] = -1 * ccso_stride - 1;
    rec_idx[1] = 1 * ccso_stride + 1;
  } else if (ext_filter_support == 2) {
    rec_idx[0] = 0 * ccso_stride - 1;
    rec_idx[1] = 0 * ccso_stride + 1;
  } else if (ext_filter_support == 3) {
    rec_idx[0] = 1 * ccso_stride - 1;
    rec_idx[1] = -1 * ccso_stride + 1;
  } else if (ext_filter_support == 4) {
    rec_idx[0] = -4 * ccso_stride;
    rec_idx[1] = 4 * ccso_stride;
  } else if (ext_filter_support == 5) {
    rec_idx[0] = 0 * ccso_stride - 4;
    rec_idx[1] = 0 * ccso_stride + 4;
  } else if (ext_filter_support == 6) {
    rec_idx[0] = -7 * ccso_stride - 7;
    rec_idx[1] = 7 * ccso_stride + 7;
  } else {  // if(ext_filter_support == 7) {
    rec_idx[0] = 7 * ccso_stride - 7;
    rec_idx[1] = -7 * ccso_stride + 7;
  }
#endif
  int *offset_buf;
  if (plane > 0) {
    offset_buf = cm->ccso_info.filter_offset[plane - 1];
  } else {
    offset_buf = filter_offset;
  }
  int ccso_stride_idx[1 << CCSO_BLK_SIZE];
  int dst_stride_idx[1 << CCSO_BLK_SIZE];
  for (int i = 0; i < (1 << CCSO_BLK_SIZE); i++) {
    ccso_stride_idx[i] = ccso_stride * i;
    dst_stride_idx[i] = dst_stride * i;
  }
  int pad_stride = CCSO_PADDING_SIZE * ccso_stride + CCSO_PADDING_SIZE;
  int y_uv_scale = xd->plane[1].subsampling_x;
  for (int y = 0; y < pic_height_c; y += (1 << CCSO_BLK_SIZE)) {
    for (int x = 0; x < pic_width_c; x += (1 << CCSO_BLK_SIZE)) {
      if (plane > 0) {
        bool use_ccso =
            (plane == 1)
                ? mi_params
                      ->mi_grid_base
                          [(1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_y)) *
                               (y >> CCSO_BLK_SIZE) * mi_params->mi_stride +
                           (1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_x)) *
                               (x >> CCSO_BLK_SIZE)]
                      ->ccso_blk_u
                : mi_params
                      ->mi_grid_base
                          [(1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_y)) *
                               (y >> CCSO_BLK_SIZE) * mi_params->mi_stride +
                           (1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_x)) *
                               (x >> CCSO_BLK_SIZE)]
                      ->ccso_blk_v;
        if (!use_ccso) continue;
      }
      int y_offset;
      int x_offset;
      if (y + (1 << CCSO_BLK_SIZE) >= pic_height_c)
        y_offset = pic_height_c - y;
      else
        y_offset = (1 << CCSO_BLK_SIZE);

      if (x + (1 << CCSO_BLK_SIZE) >= pic_width_c)
        x_offset = pic_width_c - x;
      else
        x_offset = (1 << CCSO_BLK_SIZE);

      for (int yOff = 0; yOff < y_offset; yOff++) {
        for (int xOff = 0; xOff < x_offset; xOff++) {
          cal_filter_support(rec_luma_idx,
                             &temp_rec_y_buf[((ccso_stride_idx[yOff] + x + xOff)
                                              << y_uv_scale) +
                                             pad_stride],
                             quant_step_size, inv_quant_step, rec_idx);
          int offset_val = offset_buf[(rec_luma_idx[0] << 2) + rec_luma_idx[1]];
          dstYuv8[dst_stride_idx[yOff] + x + xOff] =
              clamp(offset_val + dstYuv8[dst_stride_idx[yOff] + x + xOff], 0,
                    (1 << cm->seq_params.bit_depth) - 1);
        }
      }
    }
    temp_rec_y_buf += (ccso_stride << (CCSO_BLK_SIZE + 1));
    dstYuv8 += (dst_stride << CCSO_BLK_SIZE);
  }
}

void apply_ccso_filter_hbd(AV1_COMMON *cm, MACROBLOCKD *xd, const int plane,
                           uint16_t *temp_rec_y_buf, uint16_t *rec_yuv_16,
                           const int dst_stride,
                           int filter_offset[CCSO_LUT_EXT_SIZE],
                           uint8_t quant_step_size,
                           uint8_t ext_filter_support) {
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  const int ccso_stride = xd->plane[0].dst.width;
  const int pic_height_c = xd->plane[1].dst.height;
  const int pic_width_c = xd->plane[1].dst.width;
  int rec_luma_idx[2];
  int inv_quant_step = quant_step_size * -1;
  int rec_idx[2];
#if CCSO_FILTER_SUPPORT_EXT
  if (ext_filter_support == 0) {
    rec_idx[0] = -1 * ccso_stride;
    rec_idx[1] = 1 * ccso_stride;
  } else if (ext_filter_support == 1) {
    rec_idx[0] = -1 * ccso_stride - 1;
    rec_idx[1] = 1 * ccso_stride + 1;
  } else if (ext_filter_support == 2) {
    rec_idx[0] = 0 * ccso_stride - 1;
    rec_idx[1] = 0 * ccso_stride + 1;
  } else if (ext_filter_support == 3) {
    rec_idx[0] = 1 * ccso_stride - 1;
    rec_idx[1] = -1 * ccso_stride + 1;
  } else if (ext_filter_support == 4) {
    rec_idx[0] = -4 * ccso_stride;
    rec_idx[1] = 4 * ccso_stride;
  } else if (ext_filter_support == 5) {
    rec_idx[0] = 0 * ccso_stride - 4;
    rec_idx[1] = 0 * ccso_stride + 4;
  } else if (ext_filter_support == 6) {
    rec_idx[0] = -7 * ccso_stride - 7;
    rec_idx[1] = 7 * ccso_stride + 7;
  } else {  // if(ext_filter_support == 7) {
    rec_idx[0] = 7 * ccso_stride - 7;
    rec_idx[1] = -7 * ccso_stride + 7;
  }
#endif
  int *offset_buf;
  if (plane > 0) {
    offset_buf = cm->ccso_info.filter_offset[plane - 1];
  } else {
    offset_buf = filter_offset;
  }
  int ccso_stride_idx[1 << CCSO_BLK_SIZE];
  int dst_stride_idx[1 << CCSO_BLK_SIZE];
  for (int i = 0; i < (1 << CCSO_BLK_SIZE); i++) {
    ccso_stride_idx[i] = ccso_stride * i;
    dst_stride_idx[i] = dst_stride * i;
  }
  int pad_stride = CCSO_PADDING_SIZE * ccso_stride + CCSO_PADDING_SIZE;
  int y_uv_scale = xd->plane[1].subsampling_x;
  for (int y = 0; y < pic_height_c; y += (1 << CCSO_BLK_SIZE)) {
    for (int x = 0; x < pic_width_c; x += (1 << CCSO_BLK_SIZE)) {
      if (plane > 0) {
        bool use_ccso =
            (plane == 1)
                ? mi_params
                      ->mi_grid_base
                          [(1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_y)) *
                               (y >> CCSO_BLK_SIZE) * mi_params->mi_stride +
                           (1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_x)) *
                               (x >> CCSO_BLK_SIZE)]
                      ->ccso_blk_u
                : mi_params
                      ->mi_grid_base
                          [(1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_y)) *
                               (y >> CCSO_BLK_SIZE) * mi_params->mi_stride +
                           (1 << CCSO_BLK_SIZE >>
                            (MI_SIZE_LOG2 - xd->plane[plane].subsampling_x)) *
                               (x >> CCSO_BLK_SIZE)]
                      ->ccso_blk_v;
        if (!use_ccso) continue;
      }
      int y_offset;
      int x_offset;
      if (y + (1 << CCSO_BLK_SIZE) >= pic_height_c)
        y_offset = pic_height_c - y;
      else
        y_offset = (1 << CCSO_BLK_SIZE);

      if (x + (1 << CCSO_BLK_SIZE) >= pic_width_c)
        x_offset = pic_width_c - x;
      else
        x_offset = (1 << CCSO_BLK_SIZE);

      for (int yOff = 0; yOff < y_offset; yOff++) {
        for (int xOff = 0; xOff < x_offset; xOff++) {
          cal_filter_support(rec_luma_idx,
                             &temp_rec_y_buf[((ccso_stride_idx[yOff] + x + xOff)
                                              << y_uv_scale) +
                                             pad_stride],
                             quant_step_size, inv_quant_step, rec_idx);
          int offset_val = offset_buf[(rec_luma_idx[0] << 2) + rec_luma_idx[1]];
          rec_yuv_16[dst_stride_idx[yOff] + x + xOff] =
              clamp(offset_val + rec_yuv_16[dst_stride_idx[yOff] + x + xOff], 0,
                    (1 << cm->seq_params.bit_depth) - 1);
        }
      }
    }
    temp_rec_y_buf += (ccso_stride << (CCSO_BLK_SIZE + 1));
    rec_yuv_16 += (dst_stride << CCSO_BLK_SIZE);
  }
}

void ccso_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm, MACROBLOCKD *xd) {
  uint16_t *ext_rec_yuv[3];
  const int num_planes = av1_num_planes(cm);
  av1_setup_dst_planes(xd->plane, cm->seq_params.sb_size, frame, 0, 0, 0,
                       num_planes);
  const int ccso_stride = xd->plane[0].dst.width;
  for (int pli = 0; pli < num_planes; pli++) {
    ext_rec_yuv[pli] =
        aom_malloc(sizeof(*ext_rec_yuv) *
                   (xd->plane[0].dst.height + CCSO_PADDING_SIZE * 2) *
                   (ccso_stride + CCSO_PADDING_SIZE * 2));
    int pic_height = xd->plane[pli].dst.height;
    int pic_width = xd->plane[pli].dst.width;
    const int dst_stride = xd->plane[pli].dst.stride;
    for (int r = 0; r < pic_height; ++r) {
      for (int c = 0; c < pic_width; ++c) {
        if (cm->seq_params.use_highbitdepth) {
          ext_rec_yuv[pli][(r + CCSO_PADDING_SIZE) * ccso_stride + c +
                           CCSO_PADDING_SIZE] =
              CONVERT_TO_SHORTPTR(xd->plane[pli].dst.buf)[r * dst_stride + c];
        } else {
          ext_rec_yuv[pli][(r + CCSO_PADDING_SIZE) * ccso_stride + c +
                           CCSO_PADDING_SIZE] =
              xd->plane[pli].dst.buf[r * dst_stride + c];
        }
      }
    }
  }
  extend_ccso_border(ext_rec_yuv, CCSO_PADDING_SIZE, cm, xd);
  uint8_t quant_sz[4] = { 16, 8, 32, 64 };
  for (int plane = 1; plane < 3; plane++) {
    const int dst_stride = xd->plane[plane].dst.stride;
    uint8_t quant_step_size =
        CCSO_QUANT_STEP ? quant_sz[cm->ccso_info.quant_idx[plane - 1]] : 16;
    if (cm->ccso_info.ccso_enable[plane - 1]) {
      if (cm->seq_params.use_highbitdepth) {
        apply_ccso_filter_hbd(cm, xd, plane, ext_rec_yuv[0],
                              &CONVERT_TO_SHORTPTR(xd->plane[plane].dst.buf)[0],
                              dst_stride, NULL, quant_step_size,
                              cm->ccso_info.ext_filter_support[plane - 1]);
      } else {
        apply_ccso_filter(cm, xd, plane, ext_rec_yuv[0],
                          &xd->plane[plane].dst.buf[0], dst_stride, NULL,
                          quant_step_size,
                          cm->ccso_info.ext_filter_support[plane - 1]);
      }
    }
  }
  for (int pli = 0; pli < num_planes; pli++) {
    aom_free(ext_rec_yuv[pli]);
  }
}