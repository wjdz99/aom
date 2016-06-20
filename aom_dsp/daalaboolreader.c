#include "daalaboolreader.h"

int aom_daala_reader_init(daala_reader *r, const uint8_t *buffer, int size) {
  if (size && !buffer) {
    return 1;
  }
  r->buffer_end = buffer + size;
  r->buffer = buffer;
  od_ec_dec_init(&r->ec, buffer, size - 1);
  return 0;
}

const uint8_t *aom_daala_reader_find_end(daala_reader *r) {
  return r->buffer_end;
}
