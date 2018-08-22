/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#include <stdio.h>
#include <string.h>

#include "aom/aom_image.h"
#include "aom/aom_integer.h"
#include "aom_dsp/binary_codes_reader.h"
#include "aom_dsp/bitreader_buffer.h"
#include "av1/decoder/obu.h"
#include "common/av1c.h"
#include "config/aom_config.h"

// Helper macros to reduce verbosity required to check for read errors.
#define AV1C_READ_BIT_OR_FAIL(field, reader)                                   \
  int field = 0;                                                               \
  do {                                                                         \
    field = aom_rb_read_bit(reader);                                           \
    if (result == -1) {                                                        \
      fprintf(stderr,                                                          \
              "av1c: Error reading bit for " #field ", value=%d result=%d.\n", \
              field, result);                                                  \
      return -1;                                                               \
    }                                                                          \
  } while (0)

#define AV1C_READ_BITS_OR_FAIL(field, reader, length) \
  int field = 0;                                      \
  do {                                                \
    field = aom_rb_read_literal(reader, length);      \
    if (result == -1) {                               \
      fprintf(stderr,                                 \
              "av1c: Could not read bits for " #field \
              ", value=%d result=%d.\n",              \
              field, result);                         \
      return -1;                                      \
    }                                                 \
  } while (0)

// Helper macros for setting/restoring the error handler data in
// aom_read_bit_buffer.
#define AV1C_PUSH_ERROR_HANDLER_DATA(reader, new_data)        \
  void *original_error_handler_data = NULL;                   \
  do {                                                        \
    original_error_handler_data = reader->error_handler_data; \
    reader->error_handler_data = &new_data;                   \
  } while (0)

#define AV1C_POP_ERROR_HANDLER_DATA(reader)                   \
  do {                                                        \
    reader->error_handler_data = original_error_handler_data; \
  } while (0)

static const size_t kAv1cSize = 4;

static void bitreader_error_handler(void *data) {
  int *error_val = (int *)data;
  *error_val = 1;
}

// Parse the AV1 timing_info() structure:
// timing_info( ) {
//   num_units_in_display_tick       f(32)
//   time_scale                      f(32)
//   equal_picture_interval          f(1)
//   if (equal_picture_interval)
//     num_ticks_per_picture_minus_1 uvlc()
//   }
static int parse_timing_info(struct aom_read_bit_buffer *reader) {
  int result = 0;
  AV1C_PUSH_ERROR_HANDLER_DATA(reader, result);

  AV1C_READ_BITS_OR_FAIL(num_units_in_display_tick, reader, 32);
  AV1C_READ_BITS_OR_FAIL(time_scale, reader, 32);

  AV1C_READ_BIT_OR_FAIL(equal_picture_interval, reader);
  if (equal_picture_interval) {
    uint32_t num_ticks_per_picture_minus_1 = aom_rb_read_uvlc(reader);
    if (result == -1) {
      fprintf(stderr,
              "av1c: Could not read bits for "
              "num_ticks_per_picture_minus_1, value=%u.\n",
              num_ticks_per_picture_minus_1);
      return result;
    }
  }

  AV1C_POP_ERROR_HANDLER_DATA(reader);
  return result;
}

// Parse the AV1 decoder_model_info() structure:
// decoder_model_info( ) {
//   buffer_delay_length_minus_1            f(5)
//   num_units_in_decoding_tick             f(32)
//   buffer_removal_time_length_minus_1     f(5)
//   frame_presentation_time_length_minus_1 f(5)
// }
//
// Returns -1 upon failure, or the value of buffer_delay_length_minus_1.
static int parse_decoder_model_info(struct aom_read_bit_buffer *reader) {
  int result = 0;
  AV1C_PUSH_ERROR_HANDLER_DATA(reader, result);

  AV1C_READ_BITS_OR_FAIL(buffer_delay_length_minus_1, reader, 5);
  AV1C_READ_BITS_OR_FAIL(num_units_in_decoding_tick, reader, 32);
  AV1C_READ_BITS_OR_FAIL(buffer_removal_time_length_minus_1, reader, 5);
  AV1C_READ_BITS_OR_FAIL(frame_presentation_time_length_minus_1, reader, 5);

  AV1C_POP_ERROR_HANDLER_DATA(reader);
  return buffer_delay_length_minus_1;
}

// Parse the AV1 operating_parameters_info() structure:
// operating_parameters_info( op ) {
//   n = buffer_delay_length_minus_1 + 1
//   decoder_buffer_delay[ op ] f(n)
//   encoder_buffer_delay[ op ] f(n)
//   low_delay_mode_flag[ op ] f(1)
// }
static int parse_operating_parameters_info(struct aom_read_bit_buffer *reader,
                                           int buffer_delay_length_minus_1) {
  int result = 0;
  AV1C_PUSH_ERROR_HANDLER_DATA(reader, result);

  const int buffer_delay_length = buffer_delay_length_minus_1 + 1;
  AV1C_READ_BITS_OR_FAIL(decoder_buffer_delay, reader, buffer_delay_length);
  AV1C_READ_BITS_OR_FAIL(encoder_buffer_delay, reader, buffer_delay_length);
  AV1C_READ_BIT_OR_FAIL(low_delay_mode_flag, reader);

  AV1C_POP_ERROR_HANDLER_DATA(reader);
  return result;
}

// Parse the AV1 color_config() structure..See:
// https://aomediacodec.github.io/av1-spec/av1-spec.pdf#page=44
static int parse_color_config(struct aom_read_bit_buffer *reader,
                              av1c *config) {
  int result = 0;
  AV1C_PUSH_ERROR_HANDLER_DATA(reader, result);

  AV1C_READ_BIT_OR_FAIL(high_bitdepth, reader);
  config->high_bitdepth = high_bitdepth;

  int bit_depth = 0;
  if (config->seq_profile == 2 && config->high_bitdepth) {
    AV1C_READ_BIT_OR_FAIL(twelve_bit, reader);
    config->twelve_bit = twelve_bit;
    bit_depth = config->twelve_bit ? 12 : 10;
  } else {
    bit_depth = config->high_bitdepth ? 10 : 8;
  }

  if (config->seq_profile != 1) {
    AV1C_READ_BIT_OR_FAIL(mono_chrome, reader);
    config->monochrome = mono_chrome;
  }

  int color_primaries = AOM_CICP_CP_UNSPECIFIED;
  int transfer_characteristics = AOM_CICP_TC_UNSPECIFIED;
  int matrix_coefficients = AOM_CICP_MC_UNSPECIFIED;

  AV1C_READ_BIT_OR_FAIL(color_description_present_flag, reader);
  if (color_description_present_flag) {
    AV1C_READ_BITS_OR_FAIL(color_primaries_val, reader, 8);
    color_primaries = color_primaries_val;
    AV1C_READ_BITS_OR_FAIL(transfer_characteristics_val, reader, 8);
    transfer_characteristics = transfer_characteristics_val;
    AV1C_READ_BITS_OR_FAIL(matrix_coefficients_val, reader, 8);
    matrix_coefficients = matrix_coefficients_val;
  }

  if (config->monochrome) {
    AV1C_READ_BIT_OR_FAIL(color_range, reader);
    config->chroma_subsampling_x = 1;
    config->chroma_subsampling_y = 1;
  } else if (color_primaries == AOM_CICP_CP_BT_709 &&
             transfer_characteristics == AOM_CICP_TC_SRGB &&
             matrix_coefficients == AOM_CICP_MC_IDENTITY) {
    config->chroma_subsampling_x = 0;
    config->chroma_subsampling_y = 0;
  } else {
    AV1C_READ_BIT_OR_FAIL(color_range, reader);
    if (config->seq_profile == 0) {
      config->chroma_subsampling_x = 1;
      config->chroma_subsampling_y = 1;
    } else if (config->seq_profile == 1) {
      config->chroma_subsampling_x = 0;
      config->chroma_subsampling_y = 0;
    } else {
      if (bit_depth == 12) {
        AV1C_READ_BIT_OR_FAIL(subsampling_x, reader);
        config->chroma_subsampling_x = subsampling_x;
        if (subsampling_x) {
          AV1C_READ_BIT_OR_FAIL(subsampling_y, reader);
          config->chroma_subsampling_y = subsampling_y;
        } else {
          config->chroma_subsampling_y = 0;
        }
      } else {
        config->chroma_subsampling_x = 1;
        config->chroma_subsampling_y = 0;
      }
    }

    if (config->chroma_subsampling_x && config->chroma_subsampling_y) {
      AV1C_READ_BITS_OR_FAIL(chroma_sample_position, reader, 2);
      config->chroma_sample_position = chroma_sample_position;
    }
  }

  if (!config->monochrome) {
    AV1C_READ_BIT_OR_FAIL(separate_uv_delta_q, reader);
  }

  AV1C_POP_ERROR_HANDLER_DATA(reader);
  return result;
}

// Parse AV1 Sequence Header OBU. See:
// https://aomediacodec.github.io/av1-spec/av1-spec.pdf#page=41
static int parse_sequence_header(const uint8_t *const buffer, size_t length,
                                 av1c *config) {
  int result = 0;
  struct aom_read_bit_buffer reader_instance = {
    buffer, buffer + length, 0, &result, bitreader_error_handler
  };
  struct aom_read_bit_buffer *reader = &reader_instance;

  AV1C_READ_BITS_OR_FAIL(seq_profile, reader, 3);
  config->seq_profile = seq_profile;
  AV1C_READ_BIT_OR_FAIL(still_picture, reader);
  AV1C_READ_BIT_OR_FAIL(reduced_still_picture_header, reader);
  if (reduced_still_picture_header) {
    config->initial_presentation_delay_present = 0;
    AV1C_READ_BITS_OR_FAIL(seq_level_idx_0, reader, 5);
    config->seq_level_idx_0 = seq_level_idx_0;
    config->seq_tier_0 = 0;
  } else {
    int has_decoder_model = 0;
    int buffer_delay_length_minus_1 = 0;

    AV1C_READ_BIT_OR_FAIL(timing_info_present_flag, reader);
    if (timing_info_present_flag) {
      if (parse_timing_info(reader) != 0) return -1;

      AV1C_READ_BIT_OR_FAIL(decoder_model_info_present_flag, reader);
      if (decoder_model_info_present_flag &&
          (buffer_delay_length_minus_1 = parse_decoder_model_info(reader)) ==
              -1) {
        return -1;
      }
      has_decoder_model = 1;
    }

    AV1C_READ_BIT_OR_FAIL(initial_presentation_delay_present, reader);
    config->initial_presentation_delay_present =
        initial_presentation_delay_present;

    AV1C_READ_BITS_OR_FAIL(operating_points_cnt_minus_1, reader, 5);
    const int num_operating_points = operating_points_cnt_minus_1 + 1;

    for (int op_index = 0; op_index < num_operating_points; ++op_index) {
      AV1C_READ_BITS_OR_FAIL(operating_point_idc, reader, 12);
      AV1C_READ_BITS_OR_FAIL(seq_level_idx, reader, 5);

      int seq_tier = 0;
      if (seq_level_idx > 7) {
        AV1C_READ_BIT_OR_FAIL(seq_tier_this_op, reader);
        seq_tier = seq_tier_this_op;
      }

      if (has_decoder_model) {
        AV1C_READ_BIT_OR_FAIL(decoder_model_present_for_op, reader);
        if (decoder_model_present_for_op) {
          if (parse_operating_parameters_info(
                  reader, buffer_delay_length_minus_1) == -1)
            return -1;
        }
      }

      if (config->initial_presentation_delay_present) {
        // Skip the initial presentation delay bits if present since this
        // function has no access to the data required to properly set the
        // field.
        AV1C_READ_BIT_OR_FAIL(initial_presentation_delay_present_for_this_op,
                              reader);
        if (initial_presentation_delay_present_for_this_op) {
          AV1C_READ_BITS_OR_FAIL(initial_presentation_delay_minus_1, reader, 4);
        }
      }

      if (op_index == 0) {
        // av1c needs only the values from the first operating point.
        config->seq_level_idx_0 = seq_level_idx;
        config->seq_tier_0 = seq_tier;
        config->initial_presentation_delay_present = 0;
        config->initial_presentation_delay_minus_one = 0;
      }
    }
  }

  AV1C_READ_BITS_OR_FAIL(frame_width_bits_minus_1, reader, 4);
  AV1C_READ_BITS_OR_FAIL(frame_height_bits_minus_1, reader, 4);
  AV1C_READ_BITS_OR_FAIL(max_frame_width_minus_1, reader,
                         frame_width_bits_minus_1 + 1);
  AV1C_READ_BITS_OR_FAIL(max_frame_height_minus_1, reader,
                         frame_height_bits_minus_1 + 1);

  int frame_id_numbers_present = 0;
  if (!reduced_still_picture_header) {
    AV1C_READ_BIT_OR_FAIL(frame_id_numbers_present_flag, reader);
    frame_id_numbers_present = frame_id_numbers_present_flag;
  }

  if (frame_id_numbers_present) {
    AV1C_READ_BITS_OR_FAIL(delta_frame_id_length_minus_2, reader, 4);
    AV1C_READ_BITS_OR_FAIL(additional_frame_id_length_minus_1, reader, 3);
  }

  AV1C_READ_BIT_OR_FAIL(use_128x128_superblock, reader);
  AV1C_READ_BIT_OR_FAIL(enable_filter_intra, reader);
  AV1C_READ_BIT_OR_FAIL(enable_intra_edge_filter, reader);

  if (!reduced_still_picture_header) {
    AV1C_READ_BIT_OR_FAIL(enable_interintra_compound, reader);
    AV1C_READ_BIT_OR_FAIL(enable_masked_compound, reader);
    AV1C_READ_BIT_OR_FAIL(enable_warped_motion, reader);
    AV1C_READ_BIT_OR_FAIL(enable_dual_filter, reader);

    AV1C_READ_BIT_OR_FAIL(enable_order_hint, reader);
    if (enable_order_hint) {
      AV1C_READ_BIT_OR_FAIL(enable_jnt_comp, reader);
      AV1C_READ_BIT_OR_FAIL(enable_ref_frame_mvs, reader);
    }

    const int SELECT_SCREEN_CONTENT_TOOLS = 2;
    int seq_force_screen_content_tools = SELECT_SCREEN_CONTENT_TOOLS;
    AV1C_READ_BIT_OR_FAIL(seq_choose_screen_content_tools, reader);
    if (!seq_choose_screen_content_tools) {
      AV1C_READ_BIT_OR_FAIL(seq_force_screen_content_tools_val, reader);
      seq_force_screen_content_tools = seq_force_screen_content_tools_val;
    }

    if (seq_force_screen_content_tools > 0) {
      AV1C_READ_BIT_OR_FAIL(seq_choose_integer_mv, reader);

      const int SELECT_INTEGER_MV = 2;
      int seq_force_integer_mv = SELECT_INTEGER_MV;
      if (!seq_choose_integer_mv) {
        AV1C_READ_BIT_OR_FAIL(seq_force_integer_mv_val, reader);
        seq_force_integer_mv = seq_force_integer_mv_val;
      }
    }

    if (enable_order_hint) {
      AV1C_READ_BITS_OR_FAIL(order_hint_bits_minus_1, reader, 3);
    }
  }

  AV1C_READ_BIT_OR_FAIL(enable_superres, reader);
  AV1C_READ_BIT_OR_FAIL(enable_cdef, reader);
  AV1C_READ_BIT_OR_FAIL(enable_restoration, reader);

  return parse_color_config(reader, config);
}

int get_av1c_from_obu(const uint8_t *const buffer, size_t length, int is_annexb,
                      av1c *config) {
  if (!buffer || length == 0 || !config) {
    return -1;
  }

  ObuHeader obu_header;
  memset(&obu_header, 0, sizeof(obu_header));
  size_t consumed = 0;
  size_t sequence_header_length = 0;

  if (is_annexb) {
    // In an Annex-B stream there are 2 LEB128 size fields before the Sequence
    // Header size: temporal_unit_size and frame_unit_size. Consume them.
    size_t size_length = 0;
    uint64_t leb128_value = 0;

    for (int i = 0; i < 2; ++i) {
      if (aom_uleb_decode(buffer, length, &leb128_value, &size_length) != 0) {
        return -1;
      }
      consumed += size_length;
    }
  }

  // TODO(tomfinegan): Support consuming a Temporal Delimiter OBU.

  size_t obu_header_length = 0;
  if (aom_read_obu_header_and_size(
          buffer + consumed, length - consumed, is_annexb, &obu_header,
          &sequence_header_length, &obu_header_length) != AOM_CODEC_OK ||
      obu_header.type != OBU_SEQUENCE_HEADER ||
      sequence_header_length + consumed > length) {
    return -1;
  }
  consumed += obu_header_length;

  memset(config, 0, sizeof(*config));
  config->marker = 1;
  config->version = 1;
  return parse_sequence_header(buffer + consumed, sequence_header_length,
                               config);
}

int read_av1c(const uint8_t *const buffer, size_t buffer_length,
              size_t *bytes_read, av1c *config) {
  if (!buffer || buffer_length < kAv1cSize || !bytes_read || !config) return -1;

  *bytes_read = 0;

  int result = 0;
  struct aom_read_bit_buffer reader_instance = {
    buffer, buffer + buffer_length, 0, &result, bitreader_error_handler
  };
  struct aom_read_bit_buffer *reader = &reader_instance;

  AV1C_READ_BIT_OR_FAIL(marker, reader);
  config->marker = marker;

  AV1C_READ_BITS_OR_FAIL(version, reader, 7);
  config->version = version;

  AV1C_READ_BITS_OR_FAIL(seq_profile, reader, 3);
  config->seq_profile = seq_profile;

  AV1C_READ_BITS_OR_FAIL(seq_level_idx_0, reader, 5);
  config->seq_level_idx_0 = seq_level_idx_0;

  AV1C_READ_BIT_OR_FAIL(seq_tier_0, reader);
  config->seq_tier_0 = seq_tier_0;

  AV1C_READ_BIT_OR_FAIL(high_bitdepth, reader);
  config->high_bitdepth = high_bitdepth;

  AV1C_READ_BIT_OR_FAIL(twelve_bit, reader);
  config->twelve_bit = twelve_bit;

  AV1C_READ_BIT_OR_FAIL(monochrome, reader);
  config->monochrome = monochrome;

  AV1C_READ_BIT_OR_FAIL(chroma_subsampling_x, reader);
  config->chroma_subsampling_x = chroma_subsampling_x;

  AV1C_READ_BIT_OR_FAIL(chroma_subsampling_y, reader);
  config->chroma_subsampling_y = chroma_subsampling_y;

  AV1C_READ_BITS_OR_FAIL(chroma_sample_position, reader, 2);
  config->chroma_sample_position = chroma_sample_position;

  AV1C_READ_BITS_OR_FAIL(reserved, reader, 3);

  AV1C_READ_BIT_OR_FAIL(initial_presentation_delay_present, reader);
  config->initial_presentation_delay_present =
      initial_presentation_delay_present;

  AV1C_READ_BITS_OR_FAIL(initial_presentation_delay_minus_one, reader, 2);
  config->initial_presentation_delay_minus_one =
      initial_presentation_delay_minus_one;

  *bytes_read = reader->bit_offset >> 3;

  return 0;
}

int write_av1c(const av1c *const config, size_t capacity, size_t *bytes_written,
               uint8_t *buffer) {
  if (!buffer || capacity < kAv1cSize || !bytes_written) return -1;

  *bytes_written = 0;
  memset(buffer, 0, kAv1cSize);

  buffer[0] |= 1 << 7;
  buffer[0] |= config->version & 0x7F;
  buffer[1] |= (config->seq_profile & 0x7) << 5;
  buffer[1] |= (config->seq_level_idx_0 & 0x1F);
  buffer[2] |= config->seq_tier_0 << 7;
  buffer[2] |= (config->high_bitdepth << 6) & 0x40;
  buffer[2] |= (config->twelve_bit << 5) & 0x20;
  buffer[2] |= (config->monochrome << 4) & 0x10;
  buffer[2] |= (config->chroma_subsampling_x << 3) & 0x8;
  buffer[2] |= (config->chroma_subsampling_y << 2) & 0x4;
  buffer[2] |= config->chroma_sample_position & 0x3;
  buffer[3] |= (config->initial_presentation_delay_present << 4) & 0x10;

  if (config->initial_presentation_delay_present) {
    buffer[3] |= config->initial_presentation_delay_minus_one & 0xF;
  }

  *bytes_written = 4;
  return 0;
}
