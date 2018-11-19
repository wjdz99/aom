#ifndef AOM_AV1_ENCODER_NEW_ENCODE_STRUCTURES_H_
#define AOM_AV1_ENCODER_NEW_ENCODE_STRUCTURES_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stddef.h>

#define ENCODE_PACK_SIZE 65536

struct AV1_COMP;

typedef struct {
  unsigned int frame_number;
  int reduced_tx_set_used;
  int monochrome;
  int show_frame;
  int order_hint;
} EncodeInputFrame;

typedef struct {
} EncodedFrame;

typedef struct {
  uint8_t *pack;
  size_t pack_length;
  size_t pack_buffer_length;
} PackedFrame;

EncodeInputFrame *create_EncodeInputFrame();
void free_EncodeInputFrame(EncodeInputFrame *eif);
void av1_comp_to_EncodeInputFrame(struct AV1_COMP *cpi, EncodeInputFrame *eif);
int EncodeInputFrame_get_num_planes(EncodeInputFrame *eif);

EncodedFrame *create_EncodedFrame();
void free_EncodedFrame(EncodedFrame *ef);

PackedFrame *create_PackedFrame();
void PackedFrame_allocate_pack_buffer(PackedFrame *pf, size_t sz);
void free_PackedFrame(PackedFrame *pf);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_AV1_ENCODER_BITSTREAM_H_
