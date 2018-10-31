#include "aom_mem/aom_mem.h"
#include "av1/encoder/encoder.h"
#include "av1/encoder/new_encode_structures.h"

EncodeInputFrame *create_EncodeInputFrame() {
  EncodeInputFrame *eif;
  eif = (EncodeInputFrame *)aom_malloc(sizeof(EncodeInputFrame));
  eif->reduced_tx_set_used = 0;
  return eif;
}

void free_EncodeInputFrame(EncodeInputFrame *eif) { aom_free(eif); }

void av1_comp_to_EncodeInputFrame(AV1_COMP *cpi, EncodeInputFrame *eif) {
  // This should take all the data out of AV1_COMP needed to set up
  // EncodeInputFrame and transfer it into the supplied object. No shared
  // buffers should be reused and there should be no pointers between the two
  // structures.
  eif->frame_number = cpi->common.current_frame.frame_number;
  eif->monochrome = cpi->common.seq_params.monochrome;
  eif->show_frame = cpi->common.show_frame;
}

int EncodeInputFrame_get_num_planes(EncodeInputFrame *eif) {
  return eif->monochrome ? 1 : MAX_MB_PLANE;
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
