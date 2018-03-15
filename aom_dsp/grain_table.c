#include <string.h>
#include <stdio.h>
#include "aom_dsp/aom_dsp_common.h"
#include "aom_dsp/grain_table.h"
#include "aom_mem/aom_mem.h"

static const char kFileMagic[8] = "filmgrn1";

void aom_film_grain_table_append(aom_film_grain_table_t *t, int64_t time_stamp,
                                 int64_t end_time,
                                 const aom_film_grain_t *grain) {
  if (!t->tail || memcmp(grain, &t->tail->params, sizeof(*grain))) {
    aom_film_grain_table_entry_t *new_tail = aom_malloc(sizeof(*new_tail));
    memset(new_tail, 0, sizeof(*new_tail));
    if (t->tail) t->tail->next = new_tail;
    if (!t->head) t->head = new_tail;
    t->tail = new_tail;

    new_tail->start_time = time_stamp;
    new_tail->end_time = end_time;
    new_tail->params = *grain;
  } else {
    t->tail->end_time = AOMMAX(t->tail->end_time, end_time);
    t->tail->start_time = AOMMIN(t->tail->start_time, time_stamp);
  }
}

int aom_film_grain_table_lookup(aom_film_grain_table_t *t, int64_t time_stamp,
                                int64_t end_time, int erase,
                                aom_film_grain_t *grain) {
  aom_film_grain_table_entry_t *entry = t->head;
  aom_film_grain_table_entry_t *prev_entry = 0;
  int16_t random_seed = grain ? grain->random_seed : 0;
  if (grain) memset(grain, 0, sizeof(*grain));

  while (entry) {
    aom_film_grain_table_entry_t *next = entry->next;
    if (time_stamp >= entry->start_time && time_stamp < entry->end_time) {
      if (grain) {
        *grain = entry->params;
        if (time_stamp != 0) grain->random_seed = random_seed;
      }
      if (!erase) return 1;

      const int64_t entry_end_time = entry->end_time;
      if (time_stamp <= entry->start_time && end_time >= entry->end_time) {
        if (t->tail == entry) t->tail = prev_entry;
        if (prev_entry) {
          prev_entry->next = entry->next;
        } else {
          t->head = entry->next;
        }
        aom_free(entry);
      } else if (time_stamp <= entry->start_time &&
                 end_time < entry->end_time) {
        entry->start_time = end_time;
      } else if (time_stamp > entry->start_time &&
                 end_time >= entry->end_time) {
        entry->end_time = time_stamp;
      } else {
        aom_film_grain_table_entry_t *new_entry =
            aom_malloc(sizeof(*new_entry));
        new_entry->next = entry->next;
        new_entry->start_time = end_time;
        new_entry->end_time = entry->end_time;
        new_entry->params = entry->params;
        entry->next = new_entry;
        entry->end_time = time_stamp;
        if (t->tail == entry) t->tail = new_entry;
      }
      // If segments aren't aligned, delete from the beggining of subsequent
      // segments
      if (end_time > entry_end_time) {
        aom_film_grain_table_lookup(t, entry->end_time, end_time, 1, 0);
      }
      return 1;
    }
    prev_entry = entry;
    entry = next;
  }
  memset(grain, 0, sizeof(*grain));
  return 0;
}

aom_codec_err_t aom_film_grain_table_read(
    aom_film_grain_table_t *t, const char *filename,
    struct aom_internal_error_info *error_info) {
  FILE *file = fopen(filename, "rb");
  if (!file) {
    aom_internal_error(error_info, AOM_CODEC_ERROR, "Unable to open %s",
                       filename);
    return error_info->error_code;
  }
  error_info->error_code = AOM_CODEC_OK;

  char magic[8];
  if (!fread(magic, 8, 1, file) || memcmp(magic, kFileMagic, 8)) {
    aom_internal_error(error_info, AOM_CODEC_ERROR,
                       "Unable to read (or invalid) file magic");
    fclose(file);
    return error_info->error_code;
  }

  uint8_t *buffer = aom_malloc(65536);
  if (!buffer) {
    aom_internal_error(error_info, AOM_CODEC_MEM_ERROR,
                       "Unable to create read buffer");
    return error_info->error_code;
  }
  aom_film_grain_table_entry_t *prev_entry = 0;
  while (!feof(file)) {
    uint16_t grain_size = 0;
    int64_t timestamps[2];
    const int num_timestamps =
        fread(timestamps, sizeof(timestamps[0]), 2, file);
    if (num_timestamps != 2) {
      if (num_timestamps == 0 && feof(file)) break;
      aom_internal_error(error_info, AOM_CODEC_ERROR,
                         "Unable to read timestamps");
      break;
    }
    if (1 != fread(&grain_size, sizeof(grain_size), 1, file)) {
      aom_internal_error(error_info, AOM_CODEC_ERROR,
                         "Unable to read grain size");
      break;
    }
    if (grain_size != fread(buffer, 1, grain_size, file)) {
      aom_internal_error(error_info, AOM_CODEC_ERROR,
                         "Unable to read grain parameters");
      break;
    }

    struct aom_read_bit_buffer rb = { 0 };
    rb.bit_buffer = buffer;
    rb.bit_buffer_end = buffer + grain_size;

    aom_film_grain_table_entry_t *entry = aom_malloc(sizeof(*entry));
    entry->start_time = timestamps[0];
    entry->end_time = timestamps[1];
    entry->params.apply_grain = aom_rb_read_bit(&rb);
    entry->params.random_seed = aom_rb_read_literal(&rb, 16);
    entry->params.update_parameters = aom_rb_read_bit(&rb);
    entry->next = 0;
    av1_film_grain_read_updated(&entry->params, 0, &rb, error_info);
    if (error_info->error_code != AOM_CODEC_OK) break;

    if (prev_entry) prev_entry->next = entry;
    if (!t->head) t->head = entry;
    t->tail = entry;
    prev_entry = entry;
  }

  aom_free(buffer);
  fclose(file);
  return error_info->error_code;
}

aom_codec_err_t aom_film_grain_table_write(
    const aom_film_grain_table_t *t, const char *filename,
    struct aom_internal_error_info *error_info) {
  error_info->error_code = AOM_CODEC_OK;

  FILE *file = fopen(filename, "wb");
  if (!file) {
    aom_internal_error(error_info, AOM_CODEC_ERROR, "Unable to open file %s",
                       filename);
    return error_info->error_code;
  }

  if (!fwrite(kFileMagic, 8, 1, file)) {
    aom_internal_error(error_info, AOM_CODEC_ERROR,
                       "Unable to write file magic");
    fclose(file);
    return error_info->error_code;
  }

  uint8_t *buffer = aom_malloc(65536);
  if (!buffer) {
    aom_internal_error(error_info, AOM_CODEC_MEM_ERROR,
                       "Unable to create write buffer");
    return error_info->error_code;
  }
  aom_film_grain_table_entry_t *entry = t->head;
  while (entry) {
    struct aom_write_bit_buffer wb = { buffer, 0 };
    aom_wb_write_bit(&wb, entry->params.apply_grain);
    aom_wb_write_literal(&wb, entry->params.random_seed, 16);
    aom_wb_write_bit(&wb, entry->params.update_parameters);
    av1_film_grain_write_updated(&entry->params, 0, &wb);

    const int64_t timestamps[2] = { entry->start_time, entry->end_time };
    if (2 != fwrite(timestamps, sizeof(timestamps[0]), 2, file)) {
      aom_internal_error(error_info, AOM_CODEC_ERROR,
                         "Unable to write timestamps");
      break;
    }
    const uint16_t grain_size = (wb.bit_offset + 7) / 8;
    if (1 != fwrite(&grain_size, sizeof(grain_size), 1, file)) {
      aom_internal_error(error_info, AOM_CODEC_ERROR,
                         "Unable to write grain_size");
      break;
    }
    if (grain_size != fwrite(buffer, 1, grain_size, file)) {
      aom_internal_error(error_info, AOM_CODEC_ERROR,
                         "Unable to write grain data");
      break;
    }
    entry = entry->next;
  }
  aom_free(buffer);
  fclose(file);
  return error_info->error_code;
}

void aom_film_grain_table_free(aom_film_grain_table_t *t) {
  aom_film_grain_table_entry_t *entry = t->head;
  while (entry) {
    aom_film_grain_table_entry_t *next = entry->next;
    aom_free(entry);
    entry = next;
  }
  memset(t, 0, sizeof(*t));
}
