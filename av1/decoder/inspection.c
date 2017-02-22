#include "av1/decoder/decoder.h"
#include "av1/decoder/inspection.h"
#include "av1/common/enums.h"

void init_frame_data(insp_frame_data *frame_data, int frame_width,
 int frame_height) {
  frame_data->mi_cols =
   ALIGN_POWER_OF_TWO(frame_width, MI_SIZE_LOG2) >> MI_SIZE_LOG2;
  frame_data->mi_rows =
   ALIGN_POWER_OF_TWO(frame_height, MI_SIZE_LOG2) >> MI_SIZE_LOG2;
  frame_data->mi_grid = (insp_mi_data *)aom_malloc(sizeof(insp_mi_data)
   *frame_data->mi_rows*frame_data->mi_cols);
}

void inspect_frame(insp_frame_data *frame_data, struct AV1Decoder *pbi) {
  AV1_COMMON *const cm = &pbi->common;
  frame_data->show_frame = cm->show_frame;
  frame_data->frame_type = cm->frame_type;
}
