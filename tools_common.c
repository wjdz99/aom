/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "./tools_common.h"

#if CONFIG_AV1_ENCODER
#include "aom/aomcx.h"
#endif

#if CONFIG_AV1_DECODER
#include "aom/aomdx.h"
#endif

#if defined(_WIN32) || defined(__OS2__)
#include <io.h>
#include <fcntl.h>

#ifdef __OS2__
#define _setmode setmode
#define _fileno fileno
#define _O_BINARY O_BINARY
#endif
#endif

#define LOG_ERROR(label)               \
  do {                                 \
    const char *l = label;             \
    va_list ap;                        \
    va_start(ap, fmt);                 \
    if (l) fprintf(stderr, "%s: ", l); \
    vfprintf(stderr, fmt, ap);         \
    fprintf(stderr, "\n");             \
    va_end(ap);                        \
  } while (0)

FILE *set_binary_mode(FILE *stream) {
  (void)stream;
#if defined(_WIN32) || defined(__OS2__)
  _setmode(_fileno(stream), _O_BINARY);
#endif
  return stream;
}

void die(const char *fmt, ...) {
  LOG_ERROR(NULL);
  usage_exit();
}

void fatal(const char *fmt, ...) {
  LOG_ERROR("Fatal");
  exit(EXIT_FAILURE);
}

void warn(const char *fmt, ...) { LOG_ERROR("Warning"); }

void die_codec(aom_codec_ctx_t *ctx, const char *s) {
  const char *detail = aom_codec_error_detail(ctx);

  printf("%s: %s\n", s, aom_codec_error(ctx));
  if (detail) printf("    %s\n", detail);
  exit(EXIT_FAILURE);
}

int read_yuv_frame(struct AvxInputContext *input_ctx, aom_image_t *yuv_frame) {
  FILE *f = input_ctx->file;
  struct FileTypeDetectionBuffer *detect = &input_ctx->detect;
  int plane = 0;
  int shortread = 0;
  const int bytespp = (yuv_frame->fmt & AOM_IMG_FMT_HIGHBITDEPTH) ? 2 : 1;

  for (plane = 0; plane < 3; ++plane) {
    uint8_t *ptr;
    const int w = aom_img_plane_width(yuv_frame, plane);
    const int h = aom_img_plane_height(yuv_frame, plane);
    int r;

    /* Determine the correct plane based on the image format. The for-loop
     * always counts in Y,U,V order, but this may not match the order of
     * the data on disk.
     */
    switch (plane) {
      case 1:
        ptr =
            yuv_frame->planes[yuv_frame->fmt == AOM_IMG_FMT_YV12 ? AOM_PLANE_V
                                                                 : AOM_PLANE_U];
        break;
      case 2:
        ptr =
            yuv_frame->planes[yuv_frame->fmt == AOM_IMG_FMT_YV12 ? AOM_PLANE_U
                                                                 : AOM_PLANE_V];
        break;
      default: ptr = yuv_frame->planes[plane];
    }

    for (r = 0; r < h; ++r) {
      size_t needed = w * bytespp;
      size_t buf_position = 0;
      const size_t left = detect->buf_read - detect->position;
      if (left > 0) {
        const size_t more = (left < needed) ? left : needed;
        memcpy(ptr, detect->buf + detect->position, more);
        buf_position = more;
        needed -= more;
        detect->position += more;
      }
      if (needed > 0) {
        shortread |= (fread(ptr + buf_position, 1, needed, f) < needed);
      }

      ptr += yuv_frame->stride[plane];
    }
  }

  return shortread;
}

#if CONFIG_ENCODERS

static const AvxInterface aom_encoders[] = {
#if CONFIG_AV1_ENCODER
  { "av1", AV1_FOURCC, &aom_codec_av1_cx },
#endif
};

int get_aom_encoder_count(void) {
  return sizeof(aom_encoders) / sizeof(aom_encoders[0]);
}

const AvxInterface *get_aom_encoder_by_index(int i) { return &aom_encoders[i]; }

const AvxInterface *get_aom_encoder_by_name(const char *name) {
  int i;

  for (i = 0; i < get_aom_encoder_count(); ++i) {
    const AvxInterface *encoder = get_aom_encoder_by_index(i);
    if (strcmp(encoder->name, name) == 0) return encoder;
  }

  return NULL;
}

#endif  // CONFIG_ENCODERS

#if CONFIG_DECODERS

static const AvxInterface aom_decoders[] = {
#if CONFIG_AV1_DECODER
  { "av1", AV1_FOURCC, &aom_codec_av1_dx },
#endif
};

int get_aom_decoder_count(void) {
  return sizeof(aom_decoders) / sizeof(aom_decoders[0]);
}

const AvxInterface *get_aom_decoder_by_index(int i) { return &aom_decoders[i]; }

const AvxInterface *get_aom_decoder_by_name(const char *name) {
  int i;

  for (i = 0; i < get_aom_decoder_count(); ++i) {
    const AvxInterface *const decoder = get_aom_decoder_by_index(i);
    if (strcmp(decoder->name, name) == 0) return decoder;
  }

  return NULL;
}

const AvxInterface *get_aom_decoder_by_fourcc(uint32_t fourcc) {
  int i;

  for (i = 0; i < get_aom_decoder_count(); ++i) {
    const AvxInterface *const decoder = get_aom_decoder_by_index(i);
    if (decoder->fourcc == fourcc) return decoder;
  }

  return NULL;
}

#endif  // CONFIG_DECODERS

// TODO(dkovalev): move this function to aom_image.{c, h}, so it will be part
// of aom_image_t support
int aom_img_plane_width(const aom_image_t *img, int plane) {
  if (plane > 0 && img->x_chroma_shift > 0)
    return (img->d_w + 1) >> img->x_chroma_shift;
  else
    return img->d_w;
}

int aom_img_plane_height(const aom_image_t *img, int plane) {
  if (plane > 0 && img->y_chroma_shift > 0)
    return (img->d_h + 1) >> img->y_chroma_shift;
  else
    return img->d_h;
}

void aom_img_write(const aom_image_t *img, FILE *file) {
  int plane;

  for (plane = 0; plane < 3; ++plane) {
    const unsigned char *buf = img->planes[plane];
    const int stride = img->stride[plane];
    const int w = aom_img_plane_width(img, plane) *
                  ((img->fmt & AOM_IMG_FMT_HIGHBITDEPTH) ? 2 : 1);
    const int h = aom_img_plane_height(img, plane);
    int y;

    for (y = 0; y < h; ++y) {
      fwrite(buf, 1, w, file);
      buf += stride;
    }
  }
}

int aom_img_read(aom_image_t *img, FILE *file) {
  int plane;

  for (plane = 0; plane < 3; ++plane) {
    unsigned char *buf = img->planes[plane];
    const int stride = img->stride[plane];
    const int w = aom_img_plane_width(img, plane) *
                  ((img->fmt & AOM_IMG_FMT_HIGHBITDEPTH) ? 2 : 1);
    const int h = aom_img_plane_height(img, plane);
    int y;

    for (y = 0; y < h; ++y) {
      if (fread(buf, 1, w, file) != (size_t)w) return 0;
      buf += stride;
    }
  }

  return 1;
}

#define mmin(a, b) ((a) < (b) ? (a) : (b))

#if CONFIG_HIGHBITDEPTH
// TODO(urvang): Refactor the highbd and lowbd version.
void aom_find_mismatch_high(const aom_image_t *const img1,
                            const aom_image_t *const img2, int yloc[4],
                            int uloc[4], int vloc[4]) {
  uint16_t *plane1, *plane2;
  uint32_t stride1, stride2;
  const uint32_t bsize = 64;
  const uint32_t bsizey = bsize >> img1->y_chroma_shift;
  const uint32_t bsizex = bsize >> img1->x_chroma_shift;
  const uint32_t c_w =
      (img1->d_w + img1->x_chroma_shift) >> img1->x_chroma_shift;
  const uint32_t c_h =
      (img1->d_h + img1->y_chroma_shift) >> img1->y_chroma_shift;
  int match = 1;
  uint32_t i, j;
  yloc[0] = yloc[1] = yloc[2] = yloc[3] = -1;
  plane1 = (uint16_t *)img1->planes[AOM_PLANE_Y];
  plane2 = (uint16_t *)img2->planes[AOM_PLANE_Y];
  stride1 = img1->stride[AOM_PLANE_Y] / 2;
  stride2 = img2->stride[AOM_PLANE_Y] / 2;
  for (i = 0, match = 1; match && i < img1->d_h; i += bsize) {
    for (j = 0; match && j < img1->d_w; j += bsize) {
      int k, l;
      const int si = mmin(i + bsize, img1->d_h) - i;
      const int sj = mmin(j + bsize, img1->d_w) - j;
      for (k = 0; match && k < si; ++k) {
        for (l = 0; match && l < sj; ++l) {
          if (*(plane1 + (i + k) * stride1 + j + l) !=
              *(plane2 + (i + k) * stride2 + j + l)) {
            yloc[0] = i + k;
            yloc[1] = j + l;
            yloc[2] = *(plane1 + (i + k) * stride1 + j + l);
            yloc[3] = *(plane2 + (i + k) * stride2 + j + l);
            match = 0;
            break;
          }
        }
      }
    }
  }

  uloc[0] = uloc[1] = uloc[2] = uloc[3] = -1;
  plane1 = (uint16_t *)img1->planes[AOM_PLANE_U];
  plane2 = (uint16_t *)img2->planes[AOM_PLANE_U];
  stride1 = img1->stride[AOM_PLANE_U] / 2;
  stride2 = img2->stride[AOM_PLANE_U] / 2;
  for (i = 0, match = 1; match && i < c_h; i += bsizey) {
    for (j = 0; match && j < c_w; j += bsizex) {
      int k, l;
      const int si = mmin(i + bsizey, c_h - i);
      const int sj = mmin(j + bsizex, c_w - j);
      for (k = 0; match && k < si; ++k) {
        for (l = 0; match && l < sj; ++l) {
          if (*(plane1 + (i + k) * stride1 + j + l) !=
              *(plane2 + (i + k) * stride2 + j + l)) {
            uloc[0] = i + k;
            uloc[1] = j + l;
            uloc[2] = *(plane1 + (i + k) * stride1 + j + l);
            uloc[3] = *(plane2 + (i + k) * stride2 + j + l);
            match = 0;
            break;
          }
        }
      }
    }
  }

  vloc[0] = vloc[1] = vloc[2] = vloc[3] = -1;
  plane1 = (uint16_t *)img1->planes[AOM_PLANE_V];
  plane2 = (uint16_t *)img2->planes[AOM_PLANE_V];
  stride1 = img1->stride[AOM_PLANE_V] / 2;
  stride2 = img2->stride[AOM_PLANE_V] / 2;
  for (i = 0, match = 1; match && i < c_h; i += bsizey) {
    for (j = 0; match && j < c_w; j += bsizex) {
      int k, l;
      const int si = mmin(i + bsizey, c_h - i);
      const int sj = mmin(j + bsizex, c_w - j);
      for (k = 0; match && k < si; ++k) {
        for (l = 0; match && l < sj; ++l) {
          if (*(plane1 + (i + k) * stride1 + j + l) !=
              *(plane2 + (i + k) * stride2 + j + l)) {
            vloc[0] = i + k;
            vloc[1] = j + l;
            vloc[2] = *(plane1 + (i + k) * stride1 + j + l);
            vloc[3] = *(plane2 + (i + k) * stride2 + j + l);
            match = 0;
            break;
          }
        }
      }
    }
  }
}
#endif

void aom_find_mismatch(const aom_image_t *const img1,
                       const aom_image_t *const img2, int yloc[4], int uloc[4],
                       int vloc[4]) {
  const uint32_t bsize = 64;
  const uint32_t bsizey = bsize >> img1->y_chroma_shift;
  const uint32_t bsizex = bsize >> img1->x_chroma_shift;
  const uint32_t c_w =
      (img1->d_w + img1->x_chroma_shift) >> img1->x_chroma_shift;
  const uint32_t c_h =
      (img1->d_h + img1->y_chroma_shift) >> img1->y_chroma_shift;
  int match = 1;
  uint32_t i, j;
  yloc[0] = yloc[1] = yloc[2] = yloc[3] = -1;
  for (i = 0, match = 1; match && i < img1->d_h; i += bsize) {
    for (j = 0; match && j < img1->d_w; j += bsize) {
      int k, l;
      const int si = mmin(i + bsize, img1->d_h) - i;
      const int sj = mmin(j + bsize, img1->d_w) - j;
      for (k = 0; match && k < si; ++k) {
        for (l = 0; match && l < sj; ++l) {
          if (*(img1->planes[AOM_PLANE_Y] +
                (i + k) * img1->stride[AOM_PLANE_Y] + j + l) !=
              *(img2->planes[AOM_PLANE_Y] +
                (i + k) * img2->stride[AOM_PLANE_Y] + j + l)) {
            yloc[0] = i + k;
            yloc[1] = j + l;
            yloc[2] = *(img1->planes[AOM_PLANE_Y] +
                        (i + k) * img1->stride[AOM_PLANE_Y] + j + l);
            yloc[3] = *(img2->planes[AOM_PLANE_Y] +
                        (i + k) * img2->stride[AOM_PLANE_Y] + j + l);
            match = 0;
            break;
          }
        }
      }
    }
  }

  uloc[0] = uloc[1] = uloc[2] = uloc[3] = -1;
  for (i = 0, match = 1; match && i < c_h; i += bsizey) {
    for (j = 0; match && j < c_w; j += bsizex) {
      int k, l;
      const int si = mmin(i + bsizey, c_h - i);
      const int sj = mmin(j + bsizex, c_w - j);
      for (k = 0; match && k < si; ++k) {
        for (l = 0; match && l < sj; ++l) {
          if (*(img1->planes[AOM_PLANE_U] +
                (i + k) * img1->stride[AOM_PLANE_U] + j + l) !=
              *(img2->planes[AOM_PLANE_U] +
                (i + k) * img2->stride[AOM_PLANE_U] + j + l)) {
            uloc[0] = i + k;
            uloc[1] = j + l;
            uloc[2] = *(img1->planes[AOM_PLANE_U] +
                        (i + k) * img1->stride[AOM_PLANE_U] + j + l);
            uloc[3] = *(img2->planes[AOM_PLANE_U] +
                        (i + k) * img2->stride[AOM_PLANE_U] + j + l);
            match = 0;
            break;
          }
        }
      }
    }
  }
  vloc[0] = vloc[1] = vloc[2] = vloc[3] = -1;
  for (i = 0, match = 1; match && i < c_h; i += bsizey) {
    for (j = 0; match && j < c_w; j += bsizex) {
      int k, l;
      const int si = mmin(i + bsizey, c_h - i);
      const int sj = mmin(j + bsizex, c_w - j);
      for (k = 0; match && k < si; ++k) {
        for (l = 0; match && l < sj; ++l) {
          if (*(img1->planes[AOM_PLANE_V] +
                (i + k) * img1->stride[AOM_PLANE_V] + j + l) !=
              *(img2->planes[AOM_PLANE_V] +
                (i + k) * img2->stride[AOM_PLANE_V] + j + l)) {
            vloc[0] = i + k;
            vloc[1] = j + l;
            vloc[2] = *(img1->planes[AOM_PLANE_V] +
                        (i + k) * img1->stride[AOM_PLANE_V] + j + l);
            vloc[3] = *(img2->planes[AOM_PLANE_V] +
                        (i + k) * img2->stride[AOM_PLANE_V] + j + l);
            match = 0;
            break;
          }
        }
      }
    }
  }
}

int aom_compare_img(const aom_image_t *const img1,
                    const aom_image_t *const img2) {
  uint32_t l_w = img1->d_w;
  uint32_t c_w = (img1->d_w + img1->x_chroma_shift) >> img1->x_chroma_shift;
  const uint32_t c_h =
      (img1->d_h + img1->y_chroma_shift) >> img1->y_chroma_shift;
  uint32_t i;
  int match = 1;

  match &= (img1->fmt == img2->fmt);
  match &= (img1->d_w == img2->d_w);
  match &= (img1->d_h == img2->d_h);
#if CONFIG_HIGHBITDEPTH
  if (img1->fmt & AOM_IMG_FMT_HIGHBITDEPTH) {
    l_w *= 2;
    c_w *= 2;
  }
#endif

  for (i = 0; i < img1->d_h; ++i)
    match &= (memcmp(img1->planes[AOM_PLANE_Y] + i * img1->stride[AOM_PLANE_Y],
                     img2->planes[AOM_PLANE_Y] + i * img2->stride[AOM_PLANE_Y],
                     l_w) == 0);

  for (i = 0; i < c_h; ++i)
    match &= (memcmp(img1->planes[AOM_PLANE_U] + i * img1->stride[AOM_PLANE_U],
                     img2->planes[AOM_PLANE_U] + i * img2->stride[AOM_PLANE_U],
                     c_w) == 0);

  for (i = 0; i < c_h; ++i)
    match &= (memcmp(img1->planes[AOM_PLANE_V] + i * img1->stride[AOM_PLANE_V],
                     img2->planes[AOM_PLANE_V] + i * img2->stride[AOM_PLANE_V],
                     c_w) == 0);

  return match;
}

// TODO(dkovalev) change sse_to_psnr signature: double -> int64_t
double sse_to_psnr(double samples, double peak, double sse) {
  static const double kMaxPSNR = 100.0;

  if (sse > 0.0) {
    const double psnr = 10.0 * log10(samples * peak * peak / sse);
    return psnr > kMaxPSNR ? kMaxPSNR : psnr;
  } else {
    return kMaxPSNR;
  }
}

// TODO(debargha): Consolidate the functions below into a separate file.
#if CONFIG_HIGHBITDEPTH
static void highbd_img_upshift(aom_image_t *dst, aom_image_t *src,
                               int input_shift) {
  // Note the offset is 1 less than half.
  const int offset = input_shift > 0 ? (1 << (input_shift - 1)) - 1 : 0;
  int plane;
  if (dst->d_w != src->d_w || dst->d_h != src->d_h ||
      dst->x_chroma_shift != src->x_chroma_shift ||
      dst->y_chroma_shift != src->y_chroma_shift || dst->fmt != src->fmt ||
      input_shift < 0) {
    fatal("Unsupported image conversion");
  }
  switch (src->fmt) {
    case AOM_IMG_FMT_I42016:
    case AOM_IMG_FMT_I42216:
    case AOM_IMG_FMT_I44416:
    case AOM_IMG_FMT_I44016: break;
    default: fatal("Unsupported image conversion"); break;
  }
  for (plane = 0; plane < 3; plane++) {
    int w = src->d_w;
    int h = src->d_h;
    int x, y;
    if (plane) {
      w = (w + src->x_chroma_shift) >> src->x_chroma_shift;
      h = (h + src->y_chroma_shift) >> src->y_chroma_shift;
    }
    for (y = 0; y < h; y++) {
      uint16_t *p_src =
          (uint16_t *)(src->planes[plane] + y * src->stride[plane]);
      uint16_t *p_dst =
          (uint16_t *)(dst->planes[plane] + y * dst->stride[plane]);
      for (x = 0; x < w; x++) *p_dst++ = (*p_src++ << input_shift) + offset;
    }
  }
}

static void lowbd_img_upshift(aom_image_t *dst, aom_image_t *src,
                              int input_shift) {
  // Note the offset is 1 less than half.
  const int offset = input_shift > 0 ? (1 << (input_shift - 1)) - 1 : 0;
  int plane;
  if (dst->d_w != src->d_w || dst->d_h != src->d_h ||
      dst->x_chroma_shift != src->x_chroma_shift ||
      dst->y_chroma_shift != src->y_chroma_shift ||
      dst->fmt != src->fmt + AOM_IMG_FMT_HIGHBITDEPTH || input_shift < 0) {
    fatal("Unsupported image conversion");
  }
  switch (src->fmt) {
    case AOM_IMG_FMT_I420:
    case AOM_IMG_FMT_I422:
    case AOM_IMG_FMT_I444:
    case AOM_IMG_FMT_I440: break;
    default: fatal("Unsupported image conversion"); break;
  }
  for (plane = 0; plane < 3; plane++) {
    int w = src->d_w;
    int h = src->d_h;
    int x, y;
    if (plane) {
      w = (w + src->x_chroma_shift) >> src->x_chroma_shift;
      h = (h + src->y_chroma_shift) >> src->y_chroma_shift;
    }
    for (y = 0; y < h; y++) {
      uint8_t *p_src = src->planes[plane] + y * src->stride[plane];
      uint16_t *p_dst =
          (uint16_t *)(dst->planes[plane] + y * dst->stride[plane]);
      for (x = 0; x < w; x++) {
        *p_dst++ = (*p_src++ << input_shift) + offset;
      }
    }
  }
}

void aom_img_upshift(aom_image_t *dst, aom_image_t *src, int input_shift) {
  if (src->fmt & AOM_IMG_FMT_HIGHBITDEPTH) {
    highbd_img_upshift(dst, src, input_shift);
  } else {
    lowbd_img_upshift(dst, src, input_shift);
  }
}

void aom_img_truncate_16_to_8(aom_image_t *dst, aom_image_t *src) {
  int plane;
  if (dst->fmt + AOM_IMG_FMT_HIGHBITDEPTH != src->fmt || dst->d_w != src->d_w ||
      dst->d_h != src->d_h || dst->x_chroma_shift != src->x_chroma_shift ||
      dst->y_chroma_shift != src->y_chroma_shift) {
    fatal("Unsupported image conversion");
  }
  switch (dst->fmt) {
    case AOM_IMG_FMT_I420:
    case AOM_IMG_FMT_I422:
    case AOM_IMG_FMT_I444:
    case AOM_IMG_FMT_I440: break;
    default: fatal("Unsupported image conversion"); break;
  }
  for (plane = 0; plane < 3; plane++) {
    int w = src->d_w;
    int h = src->d_h;
    int x, y;
    if (plane) {
      w = (w + src->x_chroma_shift) >> src->x_chroma_shift;
      h = (h + src->y_chroma_shift) >> src->y_chroma_shift;
    }
    for (y = 0; y < h; y++) {
      uint16_t *p_src =
          (uint16_t *)(src->planes[plane] + y * src->stride[plane]);
      uint8_t *p_dst = dst->planes[plane] + y * dst->stride[plane];
      for (x = 0; x < w; x++) {
        *p_dst++ = (uint8_t)(*p_src++);
      }
    }
  }
}

static void highbd_img_downshift(aom_image_t *dst, aom_image_t *src,
                                 int down_shift) {
  int plane;
  if (dst->d_w != src->d_w || dst->d_h != src->d_h ||
      dst->x_chroma_shift != src->x_chroma_shift ||
      dst->y_chroma_shift != src->y_chroma_shift || dst->fmt != src->fmt ||
      down_shift < 0) {
    fatal("Unsupported image conversion");
  }
  switch (src->fmt) {
    case AOM_IMG_FMT_I42016:
    case AOM_IMG_FMT_I42216:
    case AOM_IMG_FMT_I44416:
    case AOM_IMG_FMT_I44016: break;
    default: fatal("Unsupported image conversion"); break;
  }
  for (plane = 0; plane < 3; plane++) {
    int w = src->d_w;
    int h = src->d_h;
    int x, y;
    if (plane) {
      w = (w + src->x_chroma_shift) >> src->x_chroma_shift;
      h = (h + src->y_chroma_shift) >> src->y_chroma_shift;
    }
    for (y = 0; y < h; y++) {
      uint16_t *p_src =
          (uint16_t *)(src->planes[plane] + y * src->stride[plane]);
      uint16_t *p_dst =
          (uint16_t *)(dst->planes[plane] + y * dst->stride[plane]);
      for (x = 0; x < w; x++) *p_dst++ = *p_src++ >> down_shift;
    }
  }
}

static void lowbd_img_downshift(aom_image_t *dst, aom_image_t *src,
                                int down_shift) {
  int plane;
  if (dst->d_w != src->d_w || dst->d_h != src->d_h ||
      dst->x_chroma_shift != src->x_chroma_shift ||
      dst->y_chroma_shift != src->y_chroma_shift ||
      src->fmt != dst->fmt + AOM_IMG_FMT_HIGHBITDEPTH || down_shift < 0) {
    fatal("Unsupported image conversion");
  }
  switch (dst->fmt) {
    case AOM_IMG_FMT_I420:
    case AOM_IMG_FMT_I422:
    case AOM_IMG_FMT_I444:
    case AOM_IMG_FMT_I440: break;
    default: fatal("Unsupported image conversion"); break;
  }
  for (plane = 0; plane < 3; plane++) {
    int w = src->d_w;
    int h = src->d_h;
    int x, y;
    if (plane) {
      w = (w + src->x_chroma_shift) >> src->x_chroma_shift;
      h = (h + src->y_chroma_shift) >> src->y_chroma_shift;
    }
    for (y = 0; y < h; y++) {
      uint16_t *p_src =
          (uint16_t *)(src->planes[plane] + y * src->stride[plane]);
      uint8_t *p_dst = dst->planes[plane] + y * dst->stride[plane];
      for (x = 0; x < w; x++) {
        *p_dst++ = *p_src++ >> down_shift;
      }
    }
  }
}

void aom_img_downshift(aom_image_t *dst, aom_image_t *src, int down_shift) {
  if (dst->fmt & AOM_IMG_FMT_HIGHBITDEPTH) {
    highbd_img_downshift(dst, src, down_shift);
  } else {
    lowbd_img_downshift(dst, src, down_shift);
  }
}
#endif  // CONFIG_HIGHBITDEPTH
