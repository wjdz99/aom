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

#include <stdio.h>
#include <av1/encoder/encoder.h>

#include "av1/common/av1_common_int.h"
#include "av1/common/blockd.h"
#include "av1/common/enums.h"
#include "av1/common/debugmodes.h"

static void log_frame_info(AV1_COMMON *cm, const char *str, FILE *f) {
  fprintf(f, "%s", str);
  fprintf(f, "(Frame %d, Show:%d, Q:%d): \n", cm->current_frame.frame_number,
          cm->show_frame, cm->quant_params.base_qindex);
}

void log_rd_info(const RD_STATS *rd, const char *str, FILE *f) {
  fprintf(f, "[%s] rate %d, dist %ld, rdcost %ld, sse %ld\n", str, rd->rate,
          rd->dist, rd->rdcost, rd->sse);
}

void log_mi_info(const AV1_COMMON *cm, BLOCK_SIZE bsize,
                 PARTITION_TYPE partition_type, int mi_row, int mi_col,
                 DSPL_TYPE dspl_type, int skip_txfm, const char *str,
                 int indent, FILE *f) {
  fprintf(f,
          "\n%*s[%s] frame_no=%d, qindex=%d, block_size=%s, partition_type=%s, "
          "mi_row:=%d, "
          "mi_col=%d, dspl_type=%s, skip_txfm=%d\n",
          indent, "", str, cm->current_frame.frame_number,
          cm->quant_params.base_qindex, block_size_to_str(bsize),
          partition_type_to_str(partition_type), mi_row, mi_col,
          dspl_type_to_str(dspl_type), skip_txfm);
}

const char *block_size_to_str(BLOCK_SIZE bsize) {
  const char *block_size_to_str_[] = {
    "BLOCK_4X4",   "BLOCK_4X8",    "BLOCK_8X4",      "BLOCK_8X8",
    "BLOCK_8X16",  "BLOCK_16X8",   "BLOCK_16X16",    "BLOCK_16X32",
    "BLOCK_32X16", "BLOCK_32X32",  "BLOCK_32X64",    "BLOCK_64X32",
    "BLOCK_64X64", "BLOCK_64X128", "BLOCK_128X64",   "BLOCK_128X128",
    "BLOCK_4X16",  "BLOCK_16X4",   "BLOCK_8X32",     "BLOCK_32X8",
    "BLOCK_16X64", "BLOCK_64X16",  "BLOCK_SIZES_ALL"
  };

  return block_size_to_str_[bsize];
}

const char *partition_type_to_str(PARTITION_TYPE partition_type) {
  const char *partition_type_to_str_[] = {
    "PARTITION_NONE",   "PARTITION_HORZ",   "PARTITION_VERT",
    "PARTITION_SPLIT",  "PARTITION_HORZ_A", "PARTITION_HORZ_B",
    "PARTITION_VERT_A", "PARTITION_VERT_B", "PARTITION_HORZ_4",
    "PARTITION_VERT_4"
  };
  return partition_type_to_str_[partition_type];
}

#define FN_PRINT_ARRAY2D_DEFN(TYPE)                                     \
  void print_array2d_##TYPE(const TYPE *base, int w, int h, int stride, \
                            const char *str, FILE *f) {                 \
    fprintf(f, "\n[%s, %dx%d]\n", str, w, h);                           \
    for (int i = 0; i < h; ++i) {                                       \
      for (int j = 0; j < w; ++j)                                       \
        fprintf(f, "%*d", 5, *(base + i * stride + j));                 \
      fprintf(f, "\n");                                                 \
    }                                                                   \
  }

FN_PRINT_ARRAY2D_DEFN(int8_t)
FN_PRINT_ARRAY2D_DEFN(uint8_t)
FN_PRINT_ARRAY2D_DEFN(int16_t)
FN_PRINT_ARRAY2D_DEFN(uint16_t)
FN_PRINT_ARRAY2D_DEFN(int32_t)

const char *dspl_type_to_str(DSPL_TYPE dspl_type) {
  const char *dspl_type_to_str_[] = { "", "DSPL_BAD", "DSPL_NO_TXFM",
                                      "DSPL_TXFM", "DSPL_END" };
  return dspl_type_to_str_[dspl_type];
}

/* This function dereferences a pointer to the mbmi structure
 * and uses the passed in member offset to print out the value of an integer
 * for each mbmi member value in the mi structure.
 */
static void print_mi_data(AV1_COMMON *cm, FILE *file, const char *descriptor,
                          size_t member_offset) {
  const CommonModeInfoParams *const mi_params = &cm->mi_params;
  MB_MODE_INFO **mi = mi_params->mi_grid_base;
  int rows = mi_params->mi_rows;
  int cols = mi_params->mi_cols;
  char prefix = descriptor[0];

  log_frame_info(cm, descriptor, file);
  for (int mi_row = 0; mi_row < rows; mi_row++) {
    fprintf(file, "%c ", prefix);
    for (int mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(file, "%2d ", *((char *)((char *)(mi[0]) + member_offset)));
      mi++;
    }
    fprintf(file, "\n");
    mi += mi_params->mi_stride - cols;
  }
  fprintf(file, "\n");
}

void av1_print_modes_and_motion_vectors(AV1_COMMON *cm, const char *file) {
  CommonModeInfoParams *mi_params = &cm->mi_params;
  FILE *mvs = fopen(file, "a");
  MB_MODE_INFO **mi = mi_params->mi_grid_base;
  const int rows = mi_params->mi_rows;
  const int cols = mi_params->mi_cols;

  print_mi_data(cm, mvs, "Partitions:", offsetof(MB_MODE_INFO, sb_type));
  print_mi_data(cm, mvs, "Modes:", offsetof(MB_MODE_INFO, mode));
  print_mi_data(cm, mvs, "Ref frame:", offsetof(MB_MODE_INFO, ref_frame[0]));
  print_mi_data(cm, mvs, "Transform:", offsetof(MB_MODE_INFO, tx_size));
  print_mi_data(cm, mvs, "UV Modes:", offsetof(MB_MODE_INFO, uv_mode));

  // output skip infomation.
  log_frame_info(cm, "Skips:", mvs);
  for (int mi_row = 0; mi_row < rows; mi_row++) {
    fprintf(mvs, "S ");
    for (int mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%2d ", mi[0]->skip_txfm);
      mi++;
    }
    fprintf(mvs, "\n");
    mi += mi_params->mi_stride - cols;
  }
  fprintf(mvs, "\n");

  // output motion vectors.
  log_frame_info(cm, "Vectors ", mvs);
  mi = mi_params->mi_grid_base;
  for (int mi_row = 0; mi_row < rows; mi_row++) {
    fprintf(mvs, "V ");
    for (int mi_col = 0; mi_col < cols; mi_col++) {
      fprintf(mvs, "%4d:%4d ", mi[0]->mv[0].as_mv.row, mi[0]->mv[0].as_mv.col);
      mi++;
    }
    fprintf(mvs, "\n");
    mi += mi_params->mi_stride - cols;
  }
  fprintf(mvs, "\n");

  fclose(mvs);
}

void av1_print_uncompressed_frame_header(const uint8_t *data, int size,
                                         const char *filename) {
  FILE *hdrFile = fopen(filename, "w");
  fwrite(data, size, sizeof(uint8_t), hdrFile);

  // Reset order hints(7bit + a previous bit) to 0, so that all camera frame
  // headers are identical in large scale coding.
  uint8_t zero = 0;
  fseek(hdrFile, 1, SEEK_SET);
  // Reset second byte.
  fwrite(&zero, 1, sizeof(uint8_t), hdrFile);
  fclose(hdrFile);
}

void av1_print_frame_contexts(const FRAME_CONTEXT *fc, const char *filename) {
  FILE *fcFile = fopen(filename, "w");
  const uint16_t *fcp = (uint16_t *)fc;
  const unsigned int n_contexts = sizeof(FRAME_CONTEXT) / sizeof(uint16_t);
  unsigned int i;

  for (i = 0; i < n_contexts; ++i) fprintf(fcFile, "%d ", *fcp++);
  fclose(fcFile);
}
