#include "aom_mem/aom_mem.h"
#include "av1/encoder/encoder.h"
#include "new_encode_structures.h"

EncodeInputFrame *create_EncodeInputFrame() {
  EncodeInputFrame *eif;
  eif = (EncodeInputFrame *)aom_malloc(sizeof(EncodeInputFrame));
  return eif;
}

void free_EncodeInputFrame(EncodeInputFrame *eif) { aom_free(eif); }

void av1_comp_to_EncodeInputFrame(AV1_COMP *cpi, EncodeInputFrame *eif) {
  // This should take all the data out of AV1_COMP needed to set up
  // EncodeInputFrame and transfer it into the supplied object. No shared
  // buffers should be reused and there should be no pointers between the two
  // structures.
}

EncodedFrame *create_EncodedFrame() {
  EncodedFrame *ef;
  ef = (EncodedFrame *)aom_malloc(sizeof(EncodedFrame));
  return ef;
}

void free_EncodedFrame(EncodedFrame *ef) { aom_free(ef); }

PackedFrame *create_PackedFrame() {
  PackedFrame *pf;
  pf = (PackedFrame *)aom_malloc(sizeof(PackedFrame));
  pf->pack = NULL;
  pf->pack_length = 0;
  pf->pack_buffer_length = 0;
  return pf;
}

void PackedFrame_allocate_pack_buffer(PackedFrame *pf, size_t sz) {
  if (pf->pack) {
    aom_free(pf->pack);
    pf->pack = NULL;
    pf->pack_buffer_length = 0;
  }
  pf->pack = (uint8_t *)aom_malloc(sizeof(uint8_t) * sz);
  pf->pack_length = 0;
  pf->pack_buffer_length = sz;
}

void free_PackedFrame(PackedFrame *pf) {
  if (pf->pack) {
    aom_free(pf->pack);
    pf->pack = NULL;
    pf->pack_length = 0;
    pf->pack_buffer_length = 0;
  }
  aom_free(pf);
}
