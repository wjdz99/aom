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
#include "av1/encoder/firstpass.h"
#include "av1/encoder/thirdpass.h"
#include "av1/common/blockd.h"
#include "common/video_reader.h"

static void codec_die(aom_codec_ctx_t *ctx, const char *s) {
  const char *detail = aom_codec_error_detail(ctx);

  printf("%s: %s\n", s, aom_codec_error(ctx));
  if (detail) printf("    %s\n", detail);
  exit(EXIT_FAILURE);
}

static int open_file(THIRD_PASS_DEC_CTX *ctx) {
  if (ctx->file == NULL) {
    codec_die(&ctx->codec, "No third pass input specified. ");
  }
  ctx->reader = aom_video_reader_open(ctx->file);
  if (!ctx->reader)
    codec_die(&ctx->codec, "Failed to open file for third pass.");
#if CONFIG_AV1_DECODER
  aom_codec_iface_t *decoder = aom_codec_av1_dx();
  if (!decoder) codec_die(&ctx->codec, "Unknown input codec.");
  if (aom_codec_dec_init(&ctx->codec, decoder, NULL, 0))
    codec_die(&ctx->codec, "Failed to initialize decoder.");
#else
  codec_die(&ctx->codec,
            "To utilize three-pass encoding, libaom must be built with "
            "CONFIG_AV1_DECODER=1.");
#endif
  return 0;
}

static int read_frame(THIRD_PASS_DEC_CTX *ctx) {
  if (ctx->reader == NULL) {
    open_file(ctx);
  }
  if (!ctx->have_frame) {
    if (!aom_video_reader_read_frame(ctx->reader)) return EXIT_FAILURE;
    ctx->frame = aom_video_reader_get_frame(ctx->reader, &ctx->frame_size);
    ctx->end_frame = ctx->frame + ctx->frame_size;
    ctx->have_frame = 1;
  }
  Av1DecodeReturn adr;
  if (aom_codec_decode(&ctx->codec, ctx->frame, (unsigned int)ctx->frame_size,
                       &adr) != AOM_CODEC_OK) {
    codec_die(&ctx->codec, "Failed to decode frame for third pass.");
  }
  ctx->frame = adr.buf;
  ctx->frame_size = ctx->end_frame - ctx->frame;
  if (ctx->frame == ctx->end_frame) ctx->have_frame = 0;
  return 0;
}

static int get_frame_gop_info(THIRD_PASS_DEC_CTX *ctx) {
  int cur = ctx->num_gop_info_left;
  int frame_type_flags = 0;
  if (aom_codec_control(&ctx->codec, AOMD_GET_FRAME_FLAGS, &frame_type_flags) !=
      AOM_CODEC_OK) {
    codec_die(&ctx->codec, "failed to read frame flags.");
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
    codec_die(&ctx->codec, "failed to read base q index.");
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_SHOW_EXISTING_FRAME_FLAG,
                        &ctx->frame_info[cur].is_show_existing) !=
      AOM_CODEC_OK) {
    codec_die(&ctx->codec, "failed to read show existing frame flag.");
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_SHOW_FRAME_FLAG,
                        &ctx->frame_info[cur].is_show_frame) != AOM_CODEC_OK) {
    codec_die(&ctx->codec, "failed to read show frame flag.");
  }

  if (aom_codec_control(&ctx->codec, AOMD_GET_ORDER_HINT,
                        &ctx->frame_info[cur].order_hint) != AOM_CODEC_OK) {
    codec_die(&ctx->codec, "failed to read order hint.");
  }
  ctx->num_gop_info_left++;
  return 0;
}

// TODO(bohanli): define an enum for this return
// return :
// -1: failed
// 1: no arf
// 2: with arf
static int read_gop_frames(THIRD_PASS_DEC_CTX *ctx, int max_num,
                           int *last_idx) {
  *last_idx = 0;
  int cur_idx = 0;
  int arf_order_hint = -1;
  // int kf_offset = 0;
  // int has_kf = 0;
  int num_show_frames = 0;
  while (num_show_frames < max_num) {
    // read in from bitstream if needed
    if (cur_idx >= ctx->num_gop_info_left) {
      read_frame(ctx);
      get_frame_gop_info(ctx);
    }

    // TODO(bohanli): verify that fwd_kf works here.
    if (ctx->frame_info[cur_idx].frame_type == KEY_FRAME &&
        ctx->frame_info[cur_idx].is_show_frame == 1) {
      // if this is a key frame,
      if (cur_idx != 0) {
        // we have reached the next key frame. Stop here.
        *last_idx = cur_idx - 1;
        return 1;
      }
    } else if (ctx->frame_info[cur_idx].is_show_frame == 0 &&
               arf_order_hint == -1) {
      // if this is an arf (the first no show)
      if (num_show_frames <= 1) {
        // this is an arf and we should end the GOP with its overlay
        arf_order_hint = ctx->frame_info[cur_idx].order_hint;
      } else {
        // we treat the frames previous to this arf as a gop.
        *last_idx = cur_idx - 1;
        return 1;
      }
    } else if (arf_order_hint >= 0 &&
               ctx->frame_info[cur_idx].order_hint == arf_order_hint) {
      // if this is the overlay/show existing of the arf
      assert(ctx->frame_info[cur_idx].is_show_frame);
      *last_idx = cur_idx;
      return 2;
    } else {
      // this frame is part of the GOP.
      if (ctx->frame_info[cur_idx].is_show_frame) num_show_frames++;
    }
    cur_idx++;
  }
  assert(arf_order_hint < 0);
  *last_idx = max_num - 1;
  return 1;
}

int av1_set_gop_third_pass(THIRD_PASS_DEC_CTX *ctx, GF_GROUP *gf_group,
                           int *gf_len) {
  (void)gf_group;
  // Read in future frames and find the last frame in the GOP
  int last_idx;
  read_gop_frames(ctx, 32, &last_idx);

  // TODO(bohanli): Define and set the GOP structure here too
  *gf_len =
      (ctx->frame_info[last_idx].order_hint - ctx->prev_gop_end) % (1 << 7);

  ctx->prev_gop_end = ctx->frame_info[last_idx].order_hint;

  // Reset the frame info ptr to the next frame.
  ctx->num_gop_info_left -= (last_idx + 1);
  for (int i = 0; i < ctx->num_gop_info_left; i++) {
    ctx->frame_info[i] = ctx->frame_info[i + last_idx + 1];
  }
  return 0;
}

void av1_init_thirdpass_ctx(THIRD_PASS_DEC_CTX **ctx, const char *file) {
  if (*ctx == NULL) {
    *ctx = aom_malloc(sizeof(**ctx));
  }
  THIRD_PASS_DEC_CTX *ctx_ptr = *ctx;
  ctx_ptr->file = file;
  ctx_ptr->reader = NULL;
  ctx_ptr->num_gop_info_left = 0;
  ctx_ptr->prev_gop_end = -1;
  ctx_ptr->have_frame = 0;
}

void av1_free_thirdpass_ctx(THIRD_PASS_DEC_CTX *ctx) {
  if (ctx == NULL) return;
  if (ctx->codec.iface && aom_codec_destroy(&ctx->codec)) {
    codec_die(&ctx->codec, "Failded to destroy decoder codec for thirdpass.");
  }
  aom_video_reader_close(ctx->reader);
  ctx->reader = NULL;
  ctx->num_gop_info_left = 0;
  aom_free(ctx);
}
