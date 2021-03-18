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
#ifndef AOM_AV1_ENCODER_PICKCDEF_H_
#define AOM_AV1_ENCODER_PICKCDEF_H_

#include "av1/common/cdef.h"
#include "av1/encoder/speed_features.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!\brief AV1 CDEF parameter search
 *
 * \ingroup in_loop_cdef
 *
 * Searches for optimal CDEF parameters for frame
 *
 * \param[in]      frame        Compressed frame buffer
 * \param[in]      ref          Source frame buffer
 * \param[in,out]  cm           Pointer to top level common structure
 * \param[in]      xd           Pointer to common current coding block structure
 * \param[in]      pick_method  The method used to select params
 * \param[in]      rdmult       rd multiplier to use in making param choices
 *
 * \return Nothing is returned. Instead, optimal CDEF parameters are stored
 * in the \c cdef_info structure of type \ref CdefInfo inside \c cm:
 * \arg \c cdef_bits: Bits of strength parameters
 * \arg \c nb_cdef_strengths: Number of strength parameters
 * \arg \c cdef_strengths: list of \c nb_cdef_strengths strength parameters
 * for the luma plane.
 * \arg \c uv_cdef_strengths: list of \c nb_cdef_strengths strength parameters
 * for the chroma planes.
 * \arg \c damping_factor: CDEF damping factor.
 *
 */
void av1_cdef_search(const YV12_BUFFER_CONFIG *frame,
                     const YV12_BUFFER_CONFIG *ref, AV1_COMMON *cm,
                     MACROBLOCKD *xd, CDEF_PICK_METHOD pick_method, int rdmult);

#if CONFIG_CC_CDEF
typedef struct {
  double cross_correlation[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];
  double auto_correlation[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                         [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF];
} cccdefStats;

void av1_cc_cdef_search(const YV12_BUFFER_CONFIG *frame,
                        const YV12_BUFFER_CONFIG *ref, AV1_COMMON *cm,
                        MACROBLOCKD *xd, CDEF_PICK_METHOD pick_method,
                        int rdmult, int key_freq_max, int key_freq_min);

void av1_get_filter_statstics_64x64(
    uint16_t *dst16, int dstride, uint16_t *in, int instride, uint8_t *org8,
    uint16_t *org, int orgstride, int xdec, int ydec,
    int dir[CDEF_NBLOCKS][CDEF_NBLOCKS], cdef_list *dlist, int cdef_count,
    bool is_rdo, cccdefStats filter_block_stats[MAX_NUMBER_OF_DIRECTIONS]);

void cccdef_filter_statistics_8x8_block(
    uint8_t *dst8, uint16_t *dst16, int dstride, const uint16_t *in,
    const int instride, uint8_t *org8, uint16_t *org16, int orgstride, int dir,
    int bsize, int xdec, int ydec,
    cccdefStats filter_block_stats[MAX_NUMBER_OF_DIRECTIONS]);
void get_frame_filter_stats(
    AV1_COMMON *cm, cccdefStats filter_block_stats[][MAX_NUMBER_OF_DIRECTIONS],
    cccdefStats frame_stats[MAX_NUMBER_OF_DIRECTIONS], int direction_merge);
void derive_filter_coeffs_from_stats(
    cccdefStats frame_stats[MAX_NUMBER_OF_DIRECTIONS],
    short f[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF * MAX_NUMBER_OF_DIRECTIONS],
    int direction_merge);

int matrix_decomposition(double inpMatr[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                                       [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF],
                         double outMatr[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                                       [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF],
                         int numEq);

int derive_ccdef_filter_coefficients(
    double auto_correlation[MAX_NUMBER_OF_CCCDEF_FILTER_COEFF]
                           [MAX_NUMBER_OF_CCCDEF_FILTER_COEFF],
    double *cross_correlation, double *filter_coefficients,
    int number_of_coeff);

uint64_t joint_strength_search_ccdef(int *best_lev, int nb_strengths,
                                     uint64_t mse[][TOTAL_CCCDEF_STRENGTHS],
                                     int sb_count);
uint64_t search_one_cccdef(int *lev, int nb_strengths,
                           uint64_t mse[][TOTAL_CCCDEF_STRENGTHS],
                           int sb_count);
int get_ccdef_coeff_golomb_length(short coeff);
int get_ccdef_coeff_total_bit_length(const short *const coeff,
                                     const int direction_merge_mode);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AOM_AV1_ENCODER_PICKCDEF_H_
