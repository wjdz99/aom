#if CONFIG_DSPL_RESIDUAL && CONFIG_DSPL_DEBUG
#ifndef AV1_COMMON_DEBUGMODES_H_
#define AV1_COMMON_DEBUGMODES_H_

#include "av1/common/av1_common_int.h"
#include "av1/encoder/ratectrl.h"

const char *block_size_to_str(BLOCK_SIZE bsize);
const char *update_type_to_str(FRAME_UPDATE_TYPE update_type);
const char *dspl_type_to_str(DSPL_TYPE dspl_type);
const char *partition_type_to_str(PARTITION_TYPE partition_type);

void log_rd_info(const RD_STATS *rd, const char *str, FILE *f);
void log_mi_info(const AV1_COMMON *cm, BLOCK_SIZE bsize,
                 PARTITION_TYPE partition_type, int mi_row, int mi_col,
                 DSPL_TYPE dspl_type, int skip_txfm, const char *str,
                 int indent, FILE *f);

#define FN_PRINT_ARRAY2D_DECL(TYPE)                                     \
  void print_array2d_##TYPE(const TYPE *base, int w, int h, int stride, \
                            const char *str, FILE *f);
FN_PRINT_ARRAY2D_DECL(int8_t)
FN_PRINT_ARRAY2D_DECL(uint8_t)
FN_PRINT_ARRAY2D_DECL(int16_t)
FN_PRINT_ARRAY2D_DECL(uint16_t)
FN_PRINT_ARRAY2D_DECL(int32_t)

#endif  // AV1_COMMON_DEBUGMODES_H_
#endif
