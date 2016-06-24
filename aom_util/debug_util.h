#include <stdio.h>
#include "./aom_config.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CONFIG_BITSTREAM_DEBUG
void bitstream_queue_record();
void bitstream_queue_reset();
void bitstream_queue_pop(int* result, int* prob);
void bitstream_queue_push(int result, int prob);
#endif

#ifdef __cplusplus
}  // extern "C"
#endif
