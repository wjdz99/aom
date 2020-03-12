#ifndef AOM_AV1_SAVE_DATA_H_
#define AOM_AV1_SAVE_DATA_H_

#include <stdint.h>
#include "av1/encoder/encoder.h"

#ifdef __cplusplus
extern "C" {
#endif

// Open the output file for recording.
void open_output(const char *fname);

// Record information about the current block. Note that this must be followed
// by a call to process_block_eob -- this information is computed later in the
// pipeline, thus must be made separately.
void process_block(const AV1_COMP *const cpi, MACROBLOCK *const x,
                   BLOCK_SIZE bsize);

// Called after the process_block call. With the completed information,
// the block may be saved. Note that extra calls to process_block_eob
// without a corresponding call to process_block are ignored.
void process_block_eobs(int eobs);

// Close / finalize the output file.
void close_output();

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_SAVE_DATA_H_
