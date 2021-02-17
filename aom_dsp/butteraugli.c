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

// TODO(sdeng): update the jxl api.
#include <assert.h>
#include <jxl/butteraugli.h>

#include "aom_dsp/butteraugli.h"
#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"

// BT.709 color space matrices.
static const double RGBtoYUVMatrixAdd[3] = { 16 / 256., 128 / 256.,
                                             128 / 256. };

static const double YUVtoRGBMatrix[9] = {
  1.164 * 255.,  0.000 * 255., 1.793 * 255., 1.164 * 255., -0.213 * 255.,
  -0.533 * 255., 1.164 * 255., 2.112 * 255., 0.000 * 255.,
};

// Clamp between 0 and m, and rounding to a uint16.
uint16_t Clamp(double v, uint16_t m) {
  if (v < 0) {
    return 0;
  }
  if (v > m) {
    return m;
  }
  return (uint16_t)(v + 0.5);
}

static void YUVPixelToRGB(uint16_t yv, uint16_t uv, uint16_t vv, uint8_t *r,
                          uint8_t *g, uint8_t *b, int bits) {
  const double maxv = (1 << bits) - 1;
  const double y = yv - RGBtoYUVMatrixAdd[0] * maxv;
  const double u = uv - RGBtoYUVMatrixAdd[1] * maxv;
  const double v = vv - RGBtoYUVMatrixAdd[2] * maxv;
  *r = Clamp(
      (YUVtoRGBMatrix[0] * y + YUVtoRGBMatrix[1] * u + YUVtoRGBMatrix[2] * v) /
          maxv,
      255);
  *g = Clamp(
      (YUVtoRGBMatrix[3] * y + YUVtoRGBMatrix[4] * u + YUVtoRGBMatrix[5] * v) /
          maxv,
      255);
  *b = Clamp(
      (YUVtoRGBMatrix[6] * y + YUVtoRGBMatrix[7] * u + YUVtoRGBMatrix[8] * v) /
          maxv,
      255);
}

// Converts YUV input buffers into RGB image. Width and height parameters are
// the width and height of the output rgb buffer.
static void ConvertYuv420ToRGB(const uint8_t *buffer_y, const uint8_t *buffer_u,
                               const uint8_t *buffer_v, int height, int width,
                               int row_stride_y, int row_stride_uv,
                               uint8_t *buffer_rgb) {
  assert(buffer_y);
  assert(buffer_u);
  assert(buffer_v);
  assert(buffer_rgb);

  int y_index = 0;
  int uv_index = 0;
  int rgb_index = 0;
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      const int16_t value_y = buffer_y[y_index];
      const int16_t value_u = buffer_u[uv_index];
      const int16_t value_v = buffer_v[uv_index];

      uint8_t value_r, value_g, value_b;
      YUVPixelToRGB(value_y, value_u, value_v, &value_r, &value_g, &value_b, 8);

      buffer_rgb[rgb_index++] = value_r;
      buffer_rgb[rgb_index++] = value_g;
      buffer_rgb[rgb_index++] = value_b;

      // Only progress by UV pixel stride if we've moved across 2 Y pixels.
      if (y_index % 2 == 1) {
        uv_index++;
      }
      y_index++;
    }
    y_index = i * row_stride_y;
    // Only progress by UV row stride if we've moved down 2 Y rows.
    uv_index = i / 2 * row_stride_uv;
  }
}

void aom_calc_butteraugli(const YV12_BUFFER_CONFIG *source,
                          const YV12_BUFFER_CONFIG *distorted, int bit_depth,
                          float *dist_map) {
  assert(bit_depth == 8);
  assert(source->y_width == source->uv_width * 2);
  uint8_t *src_y = source->y_buffer;
  uint8_t *src_u = source->u_buffer;
  uint8_t *src_v = source->v_buffer;
  uint8_t *distorted_y = distorted->y_buffer;
  uint8_t *distorted_u = distorted->u_buffer;
  uint8_t *distorted_v = distorted->v_buffer;
  const int width = source->y_width;
  const int height = source->y_height;

  size_t buffer_size = width * height * 3;
  uint8_t *src_rgb = (uint8_t *)aom_malloc(buffer_size);
  uint8_t *distorted_rgb = (uint8_t *)aom_malloc(buffer_size);
  ConvertYuv420ToRGB(src_y, src_u, src_v, height, width, source->y_stride,
                     source->uv_stride, src_rgb);
  ConvertYuv420ToRGB(distorted_y, distorted_u, distorted_v, height, width,
                     distorted->y_stride, distorted->uv_stride, distorted_rgb);

  JxlPixelFormat pixel_format = { 3, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0 };
  JxlButteraugliApi *api = JxlButteraugliApiCreate(NULL);
  JxlButteraugliApiSetHFAsymmetry(api, 0.8f);

  JxlButteraugliResult *result = JxlButteraugliCompute(
      api, width, height, &pixel_format, src_rgb, buffer_size, &pixel_format,
      distorted_rgb, buffer_size);

  JxlButteraugliResultGetDistance(result, /*pnorm=*/6.0);

  const float *distmap;
  uint32_t row_stride;
  JxlButteraugliResultGetDistmap(result, &distmap, &row_stride);

  for (int j = 0; j < height; ++j) {
    for (int i = 0; i < width; ++i) {
      dist_map[j * width + i] = distmap[j * row_stride + i];
    }
  }

  JxlButteraugliResultDestroy(result);
  aom_free(src_rgb);
  aom_free(distorted_rgb);
  (void)bit_depth;
}
