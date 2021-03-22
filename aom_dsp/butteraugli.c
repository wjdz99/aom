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

#include <assert.h>
#include <jxl/butteraugli.h>
#include <math.h>

#include "aom_dsp/butteraugli.h"
#include "aom_mem/aom_mem.h"
#include "third_party/libyuv/include/libyuv/convert_argb.h"

// BT.709 color space matrices. Limited range.
static const float RGBtoYUVMatrix[9] = { .183f,  .614f,  .062f,  //
                                         -.101f, -.339f, .439f,  //
                                         .439f,  -.399f, -.040f };
static const float RGBtoYUVMatrixAdd[3] = { 16.f, 128.f, 128.f };
static const float YUVtoRGBMatrix[9] = {
  1.164f, 0.000f,  1.793f,   //
  1.164f, -0.213f, -0.533f,  //
  1.164f, 2.112f,  0.000f,
};

static uint8_t ClampRound8(float x) {
  int rounded = (int)(x + 0.5f);
  if (rounded < 0) {
    rounded = 0;
  } else if (rounded > 255) {
    rounded = 255;
  }
  return rounded;
}

static void RGBPixelToYUV(float r, float g, float b, uint8_t *y, uint8_t *u,
                          uint8_t *v) {
  *y = ClampRound8(RGBtoYUVMatrixAdd[0] + RGBtoYUVMatrix[0] * r +
                   RGBtoYUVMatrix[1] * g + RGBtoYUVMatrix[2] * b);
  *u = ClampRound8(RGBtoYUVMatrixAdd[1] + RGBtoYUVMatrix[3] * r +
                   RGBtoYUVMatrix[4] * g + RGBtoYUVMatrix[5] * b);
  *v = ClampRound8(RGBtoYUVMatrixAdd[2] + RGBtoYUVMatrix[6] * r +
                   RGBtoYUVMatrix[7] * g + RGBtoYUVMatrix[8] * b);
}

static void YUVPixelToRGB(uint8_t yv, uint8_t uv, uint8_t vv, float *r,
                          float *g, float *b) {
  const float y = (float)yv - RGBtoYUVMatrixAdd[0];
  const float u = (float)uv - RGBtoYUVMatrixAdd[1];
  const float v = (float)vv - RGBtoYUVMatrixAdd[2];
  *r = YUVtoRGBMatrix[0] * y + YUVtoRGBMatrix[1] * u + YUVtoRGBMatrix[2] * v;
  *g = YUVtoRGBMatrix[3] * y + YUVtoRGBMatrix[4] * u + YUVtoRGBMatrix[5] * v;
  *b = YUVtoRGBMatrix[6] * y + YUVtoRGBMatrix[7] * u + YUVtoRGBMatrix[8] * v;
}

static void RGBtoYUV(const float *rgb, uint8_t *y, int y_stride, uint8_t *u,
                     int u_stride, uint8_t *v, int v_stride, int width,
                     int height) {
  for (int r = 0; r < height; ++r) {
    for (int c = 0; c < width; ++c) {
      RGBPixelToYUV(rgb[0], rgb[1], rgb[2], &y[c], &u[c], &v[c]);
      rgb += 3;
    }
    y += y_stride;
    u += u_stride;
    v += v_stride;
  }
}

// TODO(sdeng): We should be able to use libyuv for this conversion.
static void YUVtoRGB(const uint8_t *y, int y_stride, const uint8_t *u,
                     int u_stride, const uint8_t *v, int v_stride, float *rgb,
                     int width, int height) {
  for (int r = 0; r < height; ++r) {
    for (int c = 0; c < width; ++c) {
      YUVPixelToRGB(y[c], u[c], v[c], &rgb[0], &rgb[1], &rgb[2]);
      rgb += 3;
    }
    y += y_stride;
    u += u_stride;
    v += v_stride;
  }
}

static const float kOpsinMatrix[3][3] = {
  { 0.300000011920929f, 0.621999979019165f, 0.078000001609325f },
  { 0.230000004172325f, 0.691999971866608f, 0.078000001609325f },
  { 0.243422687053680f, 0.204767450690269f, 0.551809847354889f }
};

static const float kOpsinMatrixInverse[3][3] = {
  { 11.031567f, -9.866942f, -0.16462275f },
  { -3.2541473f, 4.41877f, -0.16462313f },
  { -3.6588511f, 2.7129228f, 1.9459282f }
};

static const float kOpsinBias = 0.003793073119596f;
static const float kCbrtOpsinBias = 0.155954198689738f;  // cbrtf(kOpsinBias)

// Converts a sRGB (0-255) array of size 3 to the XYB colorspace.
static void SrgbToJxlXyb(const float srgb[3], float xyb[3]) {
  float srgb01[3] = { srgb[0] / 255.0f, srgb[1] / 255.0f, srgb[2] / 255.0f };
  // See https://en.wikipedia.org/wiki/SRGB#Specification_of_the_transformation
  // for the formula
  float lrgb[3];
  for (int i = 0; i < 3; ++i) {
    lrgb[i] = srgb01[i] < 0.04045f ? srgb01[i] / 12.92f
                                   : powf((srgb01[i] + 0.055f) / 1.055f, 2.4f);
  }
  float opsin[3];
  for (int i = 0; i < 3; ++i) {
    opsin[i] = kOpsinBias;
    for (int j = 0; j < 3; ++j) {
      opsin[i] += lrgb[j] * kOpsinMatrix[i][j];  // kOpsinMatrix transpose
    }
  }
  float opsin_cbrt[3];
  for (int i = 0; i < 3; ++i) {
    opsin_cbrt[i] = cbrtf(opsin[i] < 0.0f ? 0.0f : opsin[i]) - kCbrtOpsinBias;
  }
  xyb[0] = 0.5f * (opsin_cbrt[0] - opsin_cbrt[1]);
  xyb[1] = 0.5f * (opsin_cbrt[0] + opsin_cbrt[1]);
  xyb[2] = opsin_cbrt[2];
}

// Converts an XYB array of size 3 to a sRGB 0-255 array.
static void JxlXybToSrgb(const float xyb[3], float srgb[3]) {
  float opsin_cbrt[3] = { xyb[0] + xyb[1], xyb[1] - xyb[0], xyb[2] };
  float opsin[3];
  for (int i = 0; i < 3; ++i) {
    opsin[i] = opsin_cbrt[i] + kCbrtOpsinBias;
    opsin[i] *= opsin[i] * opsin[i];
  }
  float lrgb[3];
  for (int i = 0; i < 3; ++i) {
    lrgb[i] = 0.0f;
    for (int j = 0; j < 3; ++j) {
      lrgb[i] += (opsin[j] - kOpsinBias) *
                 kOpsinMatrixInverse[i][j];  // kOpsinMatrixInverse transpose
    }
  }
  // See https://en.wikipedia.org/wiki/SRGB#Specification_of_the_transformation
  // for the formula
  for (int i = 0; i < 3; ++i) {
    float srgb01 = lrgb[i] < 0.0031308f
                       ? 12.92f * lrgb[i]
                       : 1.055f * powf(lrgb[i], 1 / 2.4f) - 0.055f;
    srgb[i] = srgb01 * 255.0f;
  }
}

static void RGBToXYB(const float *rgb, int stride_rgb, uint8_t *x, int x_stride,
                     uint8_t *y, int y_stride, uint8_t *b, int b_stride,
                     int width, int height) {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      const float *srgb = &rgb[i * 3 + (j * stride_rgb)];
      float xyb[3];
      SrgbToJxlXyb(srgb, xyb);
      const float xPrime = (xyb[0] + 0.015386134f) * 22.995788804f;
      const float yPrime = xyb[1] * 1.183000077f;
      const float bPrime = ((xyb[2] - xyb[1]) + 0.27770459f) * 1.502141333f;
      y[i + (j * y_stride)] = ClampRound8(yPrime * 255.0f);
      b[i + (j * b_stride)] = ClampRound8(bPrime * 255.0f);
      x[i + (j * x_stride)] = ClampRound8(xPrime * 255.0f);
    }
  }
}

static void XYBToRGB(const uint8_t *x, int x_stride, const uint8_t *y,
                     int y_stride, const uint8_t *b, int b_stride, float *rgb,
                     int stride_argb, int width, int height) {
  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      const float xPrime = (float)x[i + j * x_stride] / 255.0f;
      const float yPrime = (float)y[i + j * y_stride] / 255.0f;
      const float bPrime = (float)b[i + j * b_stride] / 255.0f;
      float xyb[3];
      xyb[0] = (xPrime / 22.995788804f) - 0.015386134f;
      xyb[1] = yPrime / 1.183000077f;
      xyb[2] = (bPrime / 1.502141333f) - 0.27770459f + xyb[1];
      float *srgb = &rgb[i * 3 + (j * stride_argb)];
      JxlXybToSrgb(xyb, srgb);
    }
  }
}

void aom_yuv_to_xyb(const YV12_BUFFER_CONFIG *yuv, YV12_BUFFER_CONFIG *xyb) {
  const int width = yuv->y_crop_width;
  const int height = yuv->y_crop_height;
  const int ss_x = yuv->subsampling_x;
  const int ss_y = yuv->subsampling_y;
  (void)ss_x;
  (void)ss_y;
  assert(ss_x == 0 && ss_y == 0);

  const size_t stride_rgb = width * 3;
  const size_t buffer_size = height * stride_rgb;
  float *rgb = (float *)aom_malloc(buffer_size * sizeof(*rgb));

  YUVtoRGB(yuv->y_buffer, yuv->y_stride, yuv->u_buffer, yuv->uv_stride,
           yuv->v_buffer, yuv->uv_stride, rgb, width, height);

  RGBToXYB(rgb, stride_rgb, xyb->v_buffer, xyb->uv_stride, xyb->y_buffer,
           xyb->y_stride, xyb->u_buffer, xyb->uv_stride, width, height);

  aom_free(rgb);
}

void aom_xyb_to_yuv(const YV12_BUFFER_CONFIG *xyb, YV12_BUFFER_CONFIG *yuv) {
  const int width = xyb->y_crop_width;
  const int height = xyb->y_crop_height;
  const int ss_x = xyb->subsampling_x;
  const int ss_y = xyb->subsampling_y;
  (void)ss_x;
  (void)ss_y;
  assert(ss_x == 0 && ss_y == 0);

  const size_t stride_rgb = width * 3;
  const size_t buffer_size = height * stride_rgb;
  float *rgb = (float *)aom_malloc(buffer_size * sizeof(*rgb));

  XYBToRGB(xyb->v_buffer, xyb->uv_stride, xyb->y_buffer, xyb->y_stride,
           xyb->u_buffer, xyb->uv_stride, rgb, stride_rgb, width, height);

  RGBtoYUV(rgb, yuv->y_buffer, yuv->y_stride, yuv->u_buffer, yuv->uv_stride,
           yuv->v_buffer, yuv->uv_stride, width, height);

  aom_free(rgb);
}

int aom_calc_butteraugli(const YV12_BUFFER_CONFIG *source,
                         const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                         float *dist_map) {
  (void)bit_depth;
  assert(bit_depth == 8);
  const int width = source->y_crop_width;
  const int height = source->y_crop_height;
  const int ss_x = source->subsampling_x;
  const int ss_y = source->subsampling_y;

  const size_t stride_argb = width * 4;
  const size_t buffer_size = height * stride_argb;
  uint8_t *src_argb = (uint8_t *)aom_malloc(buffer_size);
  uint8_t *distorted_argb = (uint8_t *)aom_malloc(buffer_size);
  if (!src_argb || !distorted_argb) {
    aom_free(src_argb);
    aom_free(distorted_argb);
    return 0;
  }

  if (ss_x == 1 && ss_y == 1) {
    I420ToARGBMatrix(source->y_buffer, source->y_stride, source->u_buffer,
                     source->uv_stride, source->v_buffer, source->uv_stride,
                     src_argb, stride_argb, &kYuvH709Constants, width, height);
    I420ToARGBMatrix(distorted->y_buffer, distorted->y_stride,
                     distorted->u_buffer, distorted->uv_stride,
                     distorted->v_buffer, distorted->uv_stride, distorted_argb,
                     stride_argb, &kYuvH709Constants, width, height);
  } else if (ss_x == 1 && ss_y == 0) {
    I422ToARGBMatrix(source->y_buffer, source->y_stride, source->u_buffer,
                     source->uv_stride, source->v_buffer, source->uv_stride,
                     src_argb, stride_argb, &kYuvH709Constants, width, height);
    I422ToARGBMatrix(distorted->y_buffer, distorted->y_stride,
                     distorted->u_buffer, distorted->uv_stride,
                     distorted->v_buffer, distorted->uv_stride, distorted_argb,
                     stride_argb, &kYuvH709Constants, width, height);
  } else if (ss_x == 0 && ss_y == 0) {
    I444ToARGBMatrix(source->y_buffer, source->y_stride, source->u_buffer,
                     source->uv_stride, source->v_buffer, source->uv_stride,
                     src_argb, stride_argb, &kYuvH709Constants, width, height);
    I444ToARGBMatrix(distorted->y_buffer, distorted->y_stride,
                     distorted->u_buffer, distorted->uv_stride,
                     distorted->v_buffer, distorted->uv_stride, distorted_argb,
                     stride_argb, &kYuvH709Constants, width, height);
  } else {
    aom_free(src_argb);
    aom_free(distorted_argb);
    return 0;
  }

  JxlPixelFormat pixel_format = { 4, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0 };
  JxlButteraugliApi *api = JxlButteraugliApiCreate(NULL);
  JxlButteraugliApiSetHFAsymmetry(api, 0.8f);

  JxlButteraugliResult *result = JxlButteraugliCompute(
      api, width, height, &pixel_format, src_argb, buffer_size, &pixel_format,
      distorted_argb, buffer_size);

  const float *distmap = NULL;
  uint32_t row_stride;
  JxlButteraugliResultGetDistmap(result, &distmap, &row_stride);
  if (distmap == NULL) {
    JxlButteraugliApiDestroy(api);
    JxlButteraugliResultDestroy(result);
    aom_free(src_argb);
    aom_free(distorted_argb);
    return 0;
  }

  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      dist_map[j * width + i] = distmap[j * row_stride + i];
    }
  }

  JxlButteraugliApiDestroy(api);
  JxlButteraugliResultDestroy(result);
  aom_free(src_argb);
  aom_free(distorted_argb);
  return 1;
}
