/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom/aom_codec.h"
#include "aom/aomdx.h"
#include "aom_mem/aom_mem.h"
#include "av1/encoder/firstpass.h"
#include "av1/encoder/thirdpass.h"
#include "av1/common/blockd.h"
#include "common/ivfdec.h"

static int setup_two_pass_stream_input(struct AvxInputContext **input_ctx_ptr,
                                       const char *input_file_name) {
  FILE *infile;
  infile = fopen(input_file_name, "rb");
  if (!infile) {
    fprintf(stderr, "Failed to open input file '%s'.\n", input_file_name);
    return -1;
  }
  struct AvxInputContext *aom_input_ctx = aom_malloc(sizeof(*aom_input_ctx));
  if (!aom_input_ctx) return -1;
  memset(aom_input_ctx, 0, sizeof(*aom_input_ctx));
  aom_input_ctx->filename = input_file_name;
  aom_input_ctx->file = infile;

  if (file_is_ivf(aom_input_ctx)) {
    aom_input_ctx->file_type = FILE_TYPE_IVF;
  } else {
    fprintf(stderr, "Unrecognized input file type.\n");
    return -1;
  }
  *input_ctx_ptr = aom_input_ctx;
  return 0;
}

static int init_third_pass(THIRD_PASS_DEC_CTX *ctx) {
  if (!ctx->input_ctx) {
    if (ctx->file == NULL) {
      fprintf(stderr, "No third pass input specified.\n");
      return -1;
    }
    if (setup_two_pass_stream_input(&ctx->input_ctx, ctx->file) != 0) {
      return -1;
    }
  }

#if CONFIG_AV1_DECODER
  if (!ctx->codec.iface) {
    aom_codec_iface_t *decoder = aom_codec_av1_dx();
    if (!decoder) {
      fprintf(stderr, "Unknown input codec.\n");
      return -1;
    }
    if (aom_codec_dec_init(&ctx->codec, decoder, NULL, 0)) {
      fprintf(stderr, "Failed to initialize decoder.\n");
      return -1;
    }
  }
#else
  fprintf(stderr,
          "To utilize three-pass encoding, libaom must be built with "
          "CONFIG_AV1_DECODER=1.\n");
  return -1;
#endif
  return 0;
}

static int read_frame(THIRD_PASS_DEC_CTX *ctx) {
  if (!ctx->input_ctx || !ctx->codec.iface) {
    if (init_third_pass(ctx) < 0) return -1;
  }
  if (!ctx->have_frame) {
    if (ivf_read_frame(ctx->input_ctx->file, &ctx->buf, &ctx->bytes_in_buffer,
                       &ctx->buffer_size, NULL) < 0) {
      if (feof(ctx->input_ctx->file)) {
        return 1;
      } else {
        return -1;
      }
    }
    ctx->frame = ctx->buf;
    ctx->end_frame = ctx->frame + ctx->bytes_in_buffer;
    ctx->have_frame = 1;
  }
#if !CONFIG_INSPECTION
  // TODO(bohanli): this is because currently only with CONFIG_INSPECTION can we
  // decode and stop at no-show frames. It could be better if we expose such a
  // function via the decoder API and use that instead, so we do not rely on the
  // CONFIG_INSPECTION flag.
  fprintf(stderr,
          "To utilize three-pass encoding, libaom must be built with "
          "CONFIG_INSPECTION=1.\n");
  return -1;
#endif
  Av1DecodeReturn adr;
  if (aom_codec_decode(&ctx->codec, ctx->frame,
                       (unsigned int)ctx->bytes_in_buffer,
                       &adr) != AOM_CODEC_OK) {
    fprintf(stderr, "Failed to decode frame for third pass.\n");
    return -1;
  }
  ctx->frame = adr.buf;
  ctx->bytes_in_buffer = ctx->end_frame - ctx->frame;
  if (ctx->frame == ctx->end_frame) ctx->have_frame = 0;
  return 0;
}

// This function gets the information needed from the recently decoded frame,
// via various decoder APIs, and saves the info into ctx->frame_info.
static int get_frame_info(THIRD_PASS_DEC_CTX *ctx) {
  int ret = read_frame(ctx);
  if (ret != 0) return ret;
  int cur = ctx->num_frame_info_left;
  int frame_type_flags = 0;
  if (aom_codec_control(&ctx->codec, AOMD_GET_FRAME_FLAGS, &frame_type_flags) !=
      AOM_CODEC_OK) {
    fprintf(stderr, "Failed to read frame flags.\n");
    return -1;
  }
  if (frame_type_flags & AOM_FRAME_IS_KEY) {
    ctx->frame_info[cur].frame_type = KEY_FRAME;
  } else if (frame_type_flags & AOM_FRAME_IS_INTRAONLY) {
    ctx->frame_info[cur].frame_type = INTRA_ONLY_FRAME;
  } else if (frame_type_flags & AOM_FRAME_IS_SWITCH) {
    ctx->frame_info[cur].frame_type = S_FRAME;
  } else {
    ctx->frame_info[cur].frame_type = INTER_FRAME;
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_BASE_Q_IDX,
                        &ctx->frame_info[cur].base_q_idx) != AOM_CODEC_OK) {
    fprintf(stderr, "Failed to read base q index.\n");
    return -1;
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_SHOW_EXISTING_FRAME_FLAG,
                        &ctx->frame_info[cur].is_show_existing) !=
      AOM_CODEC_OK) {
    fprintf(stderr, "Failed to read show existing frame flag.\n");
    return -1;
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_SHOW_FRAME_FLAG,
                        &ctx->frame_info[cur].is_show_frame) != AOM_CODEC_OK) {
    fprintf(stderr, "Failed to read show frame flag.\n");
    return -1;
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_ORDER_HINT,
                        &ctx->frame_info[cur].order_hint) != AOM_CODEC_OK) {
    fprintf(stderr, "Failed to read order hint.\n");
    return -1;
  }
  ctx->num_frame_info_left++;
  return 0;
}

// Parse the frames in the gop and determine the last frame of the current GOP.
// Decode more frames if necessary. The variable max_num is the maximum static
// GOP length if we detect an IPPP structure, and it is expected that max_mum >=
// MAX_GF_INTERVAL.
static int get_current_gop_end(THIRD_PASS_DEC_CTX *ctx, int max_num,
                               int *last_idx) {
  assert(max_num >= MAX_GF_INTERVAL);
  *last_idx = 0;
  int cur_idx = 0;
  int arf_order_hint = -1;
  int num_show_frames = 0;
  while (num_show_frames < max_num) {
    assert(cur_idx < MAX_THIRD_PASS_BUF);
    // Read in from bitstream if needed.
    if (cur_idx >= ctx->num_frame_info_left) {
      int ret = get_frame_info(ctx);
      if (ret < 0) return -1;
      if (ret == 1) {
        // At the end of the file, GOP ends in the prev frame.
        *last_idx = cur_idx - 1;
        return 0;
      }
    }

    // TODO(bohanli): verify that fwd_kf works here.
    if (ctx->frame_info[cur_idx].frame_type == KEY_FRAME &&
        ctx->frame_info[cur_idx].is_show_frame == 1) {
      if (ctx->frame_info[cur_idx].order_hint != 0) {
        // If this is a key frame and is not the first kf in this kf group, we
        // have reached the next key frame. Stop here.
        *last_idx = cur_idx - 1;
        return 0;
      }
    } else if (ctx->frame_info[cur_idx].is_show_frame == 0 &&
               arf_order_hint == -1) {
      // If this is an arf (the first no show)
      if (num_show_frames <= 1) {
        // This is an arf and we should end the GOP with its overlay.
        arf_order_hint = ctx->frame_info[cur_idx].order_hint;
      } else {
        // There are multiple show frames before the this arf, so we treat the
        // frames previous to this arf as a GOP.
        *last_idx = cur_idx - 1;
        return 0;
      }
    } else if (arf_order_hint >= 0 &&
               ctx->frame_info[cur_idx].order_hint == arf_order_hint) {
      // If this is the overlay/show existing of the arf
      assert(ctx->frame_info[cur_idx].is_show_frame);
      *last_idx = cur_idx;
      return 0;
    } else {
      // This frame is part of the GOP.
      if (ctx->frame_info[cur_idx].is_show_frame) num_show_frames++;
    }
    cur_idx++;
  }
  // This is a long IPPP GOP and we will use a length of max_num here.
  assert(arf_order_hint < 0);
  *last_idx = max_num - 1;
  return 0;
}

int av1_set_gop_third_pass(THIRD_PASS_DEC_CTX *ctx, GF_GROUP *gf_group,
                           int order_hint_bits, int *gf_len) {
  // Read in future frames and find the last frame in the current GOP.
  int last_idx;
  if (get_current_gop_end(ctx, MAX_GF_INTERVAL, &last_idx) < 0) return -1;

  // Determine the GOP length.
  // TODO(bohanli): Define and set the GOP structure here. Then we also
  // dont't need to store prev_gop_end here.
  (void)gf_group;
  *gf_len = ctx->frame_info[last_idx].order_hint - ctx->prev_gop_end;
  *gf_len = (*gf_len + (1 << order_hint_bits)) % (1 << order_hint_bits);

  ctx->prev_gop_end = ctx->frame_info[last_idx].order_hint;
  return 0;
}

int av1_pop_third_pass_info(THIRD_PASS_DEC_CTX *ctx) {
  if (ctx->num_frame_info_left == 0) return -1;
  ctx->num_frame_info_left--;
  for (int i = 0; i < ctx->num_frame_info_left; i++) {
    ctx->frame_info[i] = ctx->frame_info[i + 1];
  }
  return 0;
}

int av1_init_thirdpass_ctx(THIRD_PASS_DEC_CTX **ctx, const char *file) {
  if (*ctx == NULL) {
    *ctx = aom_calloc(1, sizeof(**ctx));
  }
  if (!ctx) return -1;
  THIRD_PASS_DEC_CTX *ctx_ptr = *ctx;
  ctx_ptr->file = file;
  ctx_ptr->input_ctx = NULL;
  ctx_ptr->num_frame_info_left = 0;
  ctx_ptr->have_frame = 0;
  ctx_ptr->buf = NULL;
  ctx_ptr->bytes_in_buffer = 0;
  ctx_ptr->buffer_size = 0;
  ctx_ptr->prev_gop_end = -1;
  return 0;
}

void av1_free_thirdpass_ctx(THIRD_PASS_DEC_CTX *ctx) {
  if (ctx == NULL) return;
  if (ctx->codec.iface) {
    aom_codec_destroy(&ctx->codec);
  }
  if (ctx->input_ctx && ctx->input_ctx->file) fclose(ctx->input_ctx->file);
  aom_free(ctx->input_ctx);
  if (ctx->buf) free(ctx->buf);
  aom_free(ctx);
}
