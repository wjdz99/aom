#ifndef AOM_INSPECTION_H_
#define AOM_INSPECTION_H_

#if CONFIG_ACCOUNTING
#include "av1/common/accounting.h"
#endif

typedef struct insp_frame_data insp_frame_data;

#include "av1/decoder/decoder.h"

typedef struct insp_mi_data insp_mi_data;

struct insp_mi_data {
  int8_t mode;
  int8_t skip;
  int8_t tx_type;
  int8_t tx_size;
};

struct insp_frame_data {
#if CONFIG_ACCOUNTING
  Accounting *accounting;
#endif
  insp_mi_data *mi_grid;
  int show_frame;
  int frame_type;
  int mi_rows;
  int mi_cols;
};

void init_frame_data(insp_frame_data *frame_data, int frame_width,
 int frame_height);
void inspect_frame(insp_frame_data *frame_data, struct AV1Decoder *pbi);

#endif
