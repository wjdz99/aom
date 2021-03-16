#ifndef AOM_AV1_COMMON_CCF_H_
#define AOM_AV1_COMMON_CCF_H_

#define min_ccf(X, Y) (((X) < (Y)) ? (X) : (Y))
#define max_ccf(X, Y) (((X) > (Y)) ? (X) : (Y))

#include <float.h>
#include "config/aom_config.h"
#include "aom/aom_integer.h"
#include "aom_ports/mem.h"
#include "av1/common/av1_common_int.h"


#ifdef __cplusplus
extern "C" {
#endif

void extend_ccso_border(uint16_t **buf, int d, AV1_COMMON *cm, MACROBLOCKD *xd);

void cal_filter_support(int rec_luma_idx[2], const uint16_t *rec_y,
                        uint8_t quant_step_size, int inv_quant_step,
                        int rec_idx[2]);

void apply_ccso_filter(AV1_COMMON *cm, MACROBLOCKD *xd, const int plane,
                       uint16_t *temp_rec_y_buf, uint8_t *dstYuv8,
                       const int dst_stride,
                       int filter_offset[CCSO_LUT_EXT_SIZE],
                       uint8_t quant_step_size, uint8_t ext_filter_support);

void apply_ccso_filter_hbd(AV1_COMMON *cm, MACROBLOCKD *xd, const int plane,
                           uint16_t *temp_rec_y_buf, uint16_t *rec_yuv_16,
                           const int dst_stride,
                           int filter_offset[CCSO_LUT_EXT_SIZE],
                           uint8_t quant_step_size, uint8_t ext_filter_support);

void ccso_frame(YV12_BUFFER_CONFIG *frame, AV1_COMMON *cm, MACROBLOCKD *xd);

#ifdef __cplusplus
}  // extern "C"
#endif
#endif  // AOM_AV1_COMMON_CDEF_H_
