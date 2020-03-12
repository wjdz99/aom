#include "av1/common/reconinter.h"
#include "av1/encoder/reconinter_enc.h"
#include "av1/encoder/save_data.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define SOURCE_IMAGE_SIZE (16 * 16)
#define INTRAPRED_LSHAPE_SIZE (20 * 4 + 4 * 16)
#define INTERPRED_SIZE (20 * 20)

typedef struct DataInfo {
  uint8_t source_image[SOURCE_IMAGE_SIZE];
  uint8_t intrapred_lshape[INTRAPRED_LSHAPE_SIZE];
  uint8_t interpred[INTERPRED_SIZE];
  int32_t lambda;
  uint8_t ref_q;
  uint8_t base_q;
} DataInfo;

static FILE *FP = NULL;
static DataInfo *INFO = NULL;

static void EK_CHECK(bool b, const char *msg) {
  if (!b) {
    printf("Check failed: %s\n", msg);
    exit(2);
  }
}

void open_output(const char *fname) {
  EK_CHECK(FP == NULL, "open_output: file pointer is null");
  FP = fopen(fname, "wb");
  EK_CHECK(FP != NULL, "open_output: file pointer is not null");
}

void close_output() {
  EK_CHECK(FP != NULL, "close_output: file pointer is not null");
  EK_CHECK(fclose(FP) == 0, "close_output: fclose success");
  FP = NULL;
}

static void copy_source_image(MACROBLOCK *const x) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_x = xd->mi_col * MI_SIZE;
  const int mi_y = xd->mi_row * MI_SIZE;
  const struct buf_2d *ref = &x->plane[0].src;
  const int stride = ref->stride;
  uint8_t *buf = ref->buf0 + mi_y * stride + mi_x;
  for (int j = 0; j < 16; ++j) {
    for (int i = 0; i < 16; ++i) {
      INFO->source_image[j * 16 + i] = buf[j * stride + i];
    }
  }
}

static void copy_intrapred_lshape(MACROBLOCK *const x) {
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_x = xd->mi_col * MI_SIZE;
  const int mi_y = xd->mi_row * MI_SIZE;
  const struct buf_2d *ref = &xd->plane[0].dst;
  const int stride = ref->stride;
  uint8_t *buf = xd->plane[0].dst.buf0 + mi_y * stride + mi_x;
  // Point back at the start of the L-region that is 4 pixels wide.
  buf -= (4 * stride + 4);

  // Copy over the top part.
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 20; ++i) {
      INFO->intrapred_lshape[j * 20 + i] = buf[j * stride + i];
    }
  }

  // Copy over the side pixels.
  const int offset = 4 * 20;
  for (int j = 4; j < 16; ++j) {
    for (int i = 0; i < 4; ++i) {
      INFO->intrapred_lshape[offset + (j - 4) * 4 + i] = buf[j * stride + i];
    }
  }
}

static void copy_interpred(const AV1_COMP *const cpi, MACROBLOCK *const x) {
  const AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  const int mi_x = xd->mi_col * MI_SIZE;
  const int mi_y = xd->mi_row * MI_SIZE;
  DECLARE_ALIGNED(32, uint8_t, dst[128 * 128]);
  InterPredExt ext = { .border_top = 8, .border_left = 8,
                       .border_right = 0, .border_bottom = 0 };
  av1_build_inter_predictors(cm, xd, 0, xd->mi[0], false, 16, 16,
                             mi_x, mi_y, enc_calc_subpel_params, NULL,
                             dst, 128, &ext);
  // Copy over the predicted region. Since border is in multiples of 8,
  // but we only want 4 pixels, skip the first 4 rows and offset by 4.
  uint8_t *buf = dst + 128 * 4 + 4;
  for (int j = 0; j < 20; ++j) {
    for (int i = 0; i < 20; ++i) {
      INFO->interpred[j * 20 + i] = buf[j * 128 + i];
    }
  }
}

static void copy_lambda(const AV1_COMP *const cpi, MACROBLOCK *const x) {
  const AV1_COMMON *const cm = &cpi->common;
  MACROBLOCKD *const xd = &x->e_mbd;
  RefCntBuffer *refbuf = get_ref_frame_buf(cm, xd->mi[0]->ref_frame[0]);
  assert(refbuf != NULL);
  INFO->ref_q = refbuf->base_qindex;
  INFO->base_q = cm->base_qindex;
  int lambda = av1_compute_rd_mult_based_on_qindex(cpi, cm->base_qindex);
  EK_CHECK(lambda > 0, "lambda is positive");
  EK_CHECK(lambda < 1000000000, "lambda less than 1 billion");
  INFO->lambda = av1_compute_rd_mult_based_on_qindex(cpi, cm->base_qindex);
}

void process_block(const AV1_COMP *const cpi, MACROBLOCK *const x,
                   BLOCK_SIZE bsize) {
  EK_CHECK(INFO == NULL, "process_block: info is null");
  MACROBLOCKD *const xd = &x->e_mbd;
  MB_MODE_INFO *mbmi = xd->mi[0];
  const int is_compound = has_second_ref(mbmi);
  const int is_inter = is_inter_block(mbmi);
  if (!is_inter || is_compound || is_intrabc_block(mbmi) ||
      bsize != BLOCK_16X16) {
    return;
  }

  INFO = malloc(sizeof(DataInfo));
  EK_CHECK(INFO != NULL, "process_block: malloc succeeded");

  copy_source_image(x);
  copy_intrapred_lshape(x);
  copy_interpred(cpi, x);
  copy_lambda(cpi, x);
}

void write_network_order(int32_t value) {
  EK_CHECK(value >= 0, "write_network_order: value is non-negative");
  for (size_t i = 0; i < sizeof(value); ++i) {
    const long shift = 8 * (sizeof(value) - i - 1);
    uint8_t byte = 0xff & (value >> shift);
    EK_CHECK(1 == fwrite(&byte, 1, sizeof(byte), FP),
             "write_network_order: byte written successfully");
  }
}

void process_block_eobs(int eobs) {
  if (INFO == NULL) {
    return;
  }
  if (eobs > 0) {
    EK_CHECK(SOURCE_IMAGE_SIZE ==
             fwrite(INFO->source_image, 1, SOURCE_IMAGE_SIZE, FP),
             "process_block_eobs: source image written");
    EK_CHECK(INTRAPRED_LSHAPE_SIZE ==
             fwrite(INFO->intrapred_lshape, 1, INTRAPRED_LSHAPE_SIZE, FP),
             "process_block_eobs: intrapred-lshape written");
    EK_CHECK(INTERPRED_SIZE == fwrite(INFO->interpred, 1, INTERPRED_SIZE, FP),
             "process_block_eobs: interpred written");
    write_network_order(INFO->lambda);
    EK_CHECK(1 == fwrite(&(INFO->ref_q), 1, 1, FP),
             "process_block_eobs: ref_q written");
    EK_CHECK(1 == fwrite(&(INFO->base_q), 1, 1, FP),
             "process_block_eobs: base_q written");
  }
  printf("ref_q: %d\n", INFO->ref_q);
  printf("base_q: %d\n", INFO->base_q);
  printf("Lambda: %d\n", INFO->lambda);
  free(INFO);
  INFO = NULL;
}
