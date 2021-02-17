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

// Converts YUV input buffers into RGB image. Width and height parameters are
// the width and height of the output rgb buffer.
void static ConvertYuv420ToRGB(const uint8_t* buffer_y, const uint8_t* buffer_u,
                        const uint8_t* buffer_v, int height, int width,
                        int row_stride_y, int row_stride_uv,
                        uint8_t* buffer_rgb) {
  JXL_CHECK(buffer_y);
  JXL_CHECK(buffer_u);
  JXL_CHECK(buffer_v);
  JXL_CHECK(buffer_rgb);

  int y_index = 0;
  int uv_index = 0;
  int rgb_index = 0;
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      // Converts from YUV420 to RGB888
      // There are four times more Y pixels than U or V pixels.
      // See:
      // https://en.wikipedia.org/wiki/YUV#Y.E2.80.B2UV444_to_RGB888_conversion

      const int32_t value_y = buffer_y[y_index] - 16;
      const int32_t value_u = buffer_u[uv_index] - 128;
      const int32_t value_v = buffer_v[uv_index] - 128;

      const int32_t y_298_p_128 = value_y * 298 + 128;
      const int32_t value_r = (y_298_p_128 + value_v * 409) >> 8;
      const int32_t value_g = (y_298_p_128 - value_u * 100 - value_v * 208) >> 8;
      const int32_t value_b = (y_298_p_128 + value_u * 516) >> 8;

      buffer_rgb[rgb_index++] = Clamp(value_r);
      buffer_rgb[rgb_index++] = Clamp(value_g);
      buffer_rgb[rgb_index++] = Clamp(value_b);

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



  JxlPixelFormat pixel_format = {3, JXL_TYPE_UINT8, JXL_NATIVE_ENDIAN, 0};
  JxlButteraugliApiStruct* api = JxlButteraugliApiCreate(nullptr);
  JxlButteraugliApiSetHFAsymmetry(api, 0.8f);

  JxlButteraugliResult* result = JxlButteraugliCompute(
      api, width, height, &pixel_format, src_y,
      buffer_size, &pixel_format, distorted_y, buffer_size);

  double distance_p = JxlButteraugliResultGetDistance(result.get(),
                                                      /*pnorm=*/6);

}
