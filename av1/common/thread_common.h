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

#ifndef AOM_AV1_COMMON_THREAD_COMMON_H_
#define AOM_AV1_COMMON_THREAD_COMMON_H_

#include "config/aom_config.h"

#include "av1/common/av1_loopfilter.h"
#include "av1/common/cdef.h"
#include "aom_util/aom_thread.h"

#ifdef __cplusplus
extern "C" {
#endif

struct AV1Common;

typedef struct AV1LfMTInfo {
  int mi_row;
  int plane;
  int dir;
  int lpf_opt_level;
} AV1LfMTInfo;

// Loopfilter row synchronization
typedef struct AV1LfSyncData {
#if CONFIG_MULTITHREAD
  pthread_mutex_t *mutex_[MAX_MB_PLANE];
  pthread_cond_t *cond_[MAX_MB_PLANE];
#endif
  // Allocate memory to store the loop-filtered superblock index in each row.
  int *cur_sb_col[MAX_MB_PLANE];
  // The optimal sync_range for different resolution and platform should be
  // determined by testing. Currently, it is chosen to be a power-of-2 number.
  int sync_range;
  int rows;

  // Row-based parallel loopfilter data
  LFWorkerData *lfdata;
  int num_workers;

#if CONFIG_MULTITHREAD
  pthread_mutex_t *job_mutex;
#endif
  AV1LfMTInfo *job_queue;
  int jobs_enqueued;
  int jobs_dequeued;
} AV1LfSync;

typedef struct AV1LrMTInfo {
  int v_start;
  int v_end;
  int lr_unit_row;
  int plane;
  int sync_mode;
  int v_copy_start;
  int v_copy_end;
} AV1LrMTInfo;

typedef struct LoopRestorationWorkerData {
  int32_t *rst_tmpbuf;
  void *rlbs;
  void *lr_ctxt;
  int do_extend_border;
} LRWorkerData;

// Looprestoration row synchronization
typedef struct AV1LrSyncData {
#if CONFIG_MULTITHREAD
  pthread_mutex_t *mutex_[MAX_MB_PLANE];
  pthread_cond_t *cond_[MAX_MB_PLANE];
#endif
  // Allocate memory to store the loop-restoration block index in each row.
  int *cur_sb_col[MAX_MB_PLANE];
  // The optimal sync_range for different resolution and platform should be
  // determined by testing. Currently, it is chosen to be a power-of-2 number.
  int sync_range;
  int rows;
  int num_planes;

  int num_workers;

#if CONFIG_MULTITHREAD
  pthread_mutex_t *job_mutex;
#endif
  // Row-based parallel loopfilter data
  LRWorkerData *lrworkerdata;

  AV1LrMTInfo *job_queue;
  int jobs_enqueued;
  int jobs_dequeued;
} AV1LrSync;

typedef struct AV1CdefWorker {
  AV1_COMMON *cm;
  MACROBLOCKD *xd;
  uint16_t *colbuf[MAX_MB_PLANE];
  uint16_t *srcbuf;
  uint16_t *linebuf[MAX_MB_PLANE];
  cdef_init_fb_row_t cdef_init_fb_row_fn;
  int do_extend_border;
} AV1CdefWorkerData;

typedef struct AV1CdefRowSync {
#if CONFIG_MULTITHREAD
  pthread_mutex_t *row_mutex_;
  pthread_cond_t *row_cond_;
#endif  // CONFIG_MULTITHREAD
  int is_row_done;
} AV1CdefRowSync;

// Data related to CDEF search multi-thread synchronization.
typedef struct AV1CdefSyncData {
#if CONFIG_MULTITHREAD
  // Mutex lock used while dispatching jobs.
  pthread_mutex_t *mutex_;
#endif  // CONFIG_MULTITHREAD
  // Data related to CDEF row mt sync information
  AV1CdefRowSync *cdef_row_mt;
  // Flag to indicate all blocks are processed and end of frame is reached
  int end_of_frame;
  // Row index in units of 64x64 block
  int fbr;
  // Column index in units of 64x64 block
  int fbc;
} AV1CdefSync;

void av1_cdef_frame_mt(AV1_COMMON *const cm, MACROBLOCKD *const xd,
                       AV1CdefWorkerData *const cdef_worker,
                       AVxWorker *const workers, AV1CdefSync *const cdef_sync,
                       int num_workers, cdef_init_fb_row_t cdef_init_fb_row_fn,
                       int do_extend_border);
void av1_cdef_init_fb_row_mt(const AV1_COMMON *const cm,
                             const MACROBLOCKD *const xd,
                             CdefBlockInfo *const fb_info,
                             uint16_t **const linebuf, uint16_t *const src,
                             struct AV1CdefSyncData *const cdef_sync, int fbr);
void av1_cdef_copy_sb8_16(const AV1_COMMON *const cm, uint16_t *const dst,
                          int dstride, const uint8_t *src, int src_voffset,
                          int src_hoffset, int sstride, int vsize, int hsize);
void av1_alloc_cdef_sync(AV1_COMMON *const cm, AV1CdefSync *cdef_sync,
                         int num_workers);
void av1_free_cdef_sync(AV1CdefSync *cdef_sync);

// Deallocate loopfilter synchronization related mutex and data.
void av1_loop_filter_dealloc(AV1LfSync *lf_sync);
void av1_loop_filter_alloc(AV1LfSync *lf_sync, AV1_COMMON *cm, int rows,
                           int width, int num_workers);

void av1_loop_filter_frame_mt(YV12_BUFFER_CONFIG *frame, struct AV1Common *cm,
                              struct macroblockd *xd, int plane_start,
                              int plane_end, int partial_frame,
                              AVxWorker *workers, int num_workers,
                              AV1LfSync *lf_sync, int lpf_opt_level);

void av1_loop_restoration_filter_frame_mt(YV12_BUFFER_CONFIG *frame,
                                          struct AV1Common *cm,
                                          int optimized_lr, AVxWorker *workers,
                                          int num_workers, AV1LrSync *lr_sync,
                                          void *lr_ctxt, int do_extend_border);
void av1_loop_restoration_dealloc(AV1LrSync *lr_sync, int num_workers);
void av1_loop_restoration_alloc(AV1LrSync *lr_sync, AV1_COMMON *cm,
                                int num_workers, int num_rows_lr,
                                int num_planes, int width);
int av1_get_intrabc_extra_top_right_sb_delay(const AV1_COMMON *cm);

void av1_loop_filter_data_reset(LFWorkerData *lf_data,
                                YV12_BUFFER_CONFIG *frame_buffer,
                                struct AV1Common *cm, MACROBLOCKD *xd);

void av1_lf_mt_init(AV1_COMMON *cm, AV1LfSync *lf_sync, int *planes_to_lf,
                    int num_workers, int lpf_opt_level);

AV1LfMTInfo *av1_get_lf_job_info(AV1LfSync *lf_sync);

void av1_thread_loop_filter_rows(
    const YV12_BUFFER_CONFIG *const frame_buffer, AV1_COMMON *const cm,
    struct macroblockd_plane *planes, MACROBLOCKD *xd, int mi_row, int plane,
    int dir, int lpf_opt_level, AV1LfSync *const lf_sync,
    AV1_DEBLOCKING_PARAMETERS *params_buf, TX_SIZE *tx_buf, int mib_size_log2);

static AOM_INLINE int planes_to_loop_filter(struct loopfilter *lf,
                                            int *planes_to_lf, int plane_start,
                                            int plane_end) {
  // For each luma and chroma plane, whether to filter it or not.
  planes_to_lf[0] = (lf->filter_level[0] || lf->filter_level[1]) &&
                    plane_start <= 0 && 0 < plane_end;
  planes_to_lf[1] = lf->filter_level_u && plane_start <= 1 && 1 < plane_end;
  planes_to_lf[2] = lf->filter_level_v && plane_start <= 2 && 2 < plane_end;
  // If the luma plane is purposely not filtered, neither are the chroma
  // planes.
  if (!planes_to_lf[0] && plane_start <= 0 && 0 < plane_end) return 0;
  // Early exit.
  if (!planes_to_lf[0] && !planes_to_lf[1] && !planes_to_lf[2]) return 0;
  return 1;
}

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_COMMON_THREAD_COMMON_H_
