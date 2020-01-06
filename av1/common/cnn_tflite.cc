/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <vector>
#include "cnn_tflite.h"

#include "av1/common/onyxc_int.h"
#include "av1/tflite_models/intra_frame_model/qp32.h"
#include "av1/tflite_models/intra_frame_model/qp32_op_registrations.h"
#include "third_party/tensorflow/tensorflow/lite/interpreter.h"
#include "third_party/tensorflow/tensorflow/lite/kernels/kernel_util.h"
#include "third_party/tensorflow/tensorflow/lite/model.h"

// Builds and returns the TFlite interpreter.
static std::unique_ptr<tflite::Interpreter> get_tflite_interpreter(int qindex,
                                                                   int width,
                                                                   int height) {
  assert(av1_use_cnn_tflite(qindex));
  // TODO(urvang): Switch TF-lite model and op registrations based on qindex.
  (void)qindex;
  auto model = tflite::GetModel(qp32_model_tflite_data);
  tflite::MutableOpResolver resolver;
  RegisterSelectedOpsQp32(&resolver);
  tflite::InterpreterBuilder builder(model, resolver);
  std::unique_ptr<tflite::Interpreter> interpreter;
  builder(&interpreter);

  // TODO(urvang): Set to available number of threads in 'cm'.
  interpreter->SetNumThreads(1);

  tflite::ErrorReporter *reporter = tflite::DefaultErrorReporter();

  // Dimension order: batch_size, height, width, num_channels.
  // Note: height comes before width here!
  const std::vector<int> in_out_dims = { 1, height, width, 1 };
  // We only need to resize the input tensor. All other tensors (including
  // output tensor) will be resized automatically.
  if (interpreter->ResizeInputTensor(interpreter->inputs()[0], in_out_dims) !=
      kTfLiteOk) {
    reporter->Report("Failed at input tensor resize");
    return nullptr;
  }

  if (interpreter->AllocateTensors() != kTfLiteOk) {
    reporter->Report("Failed at tensor allocation");
    return nullptr;
  }
  return interpreter;
}

extern "C" int av1_restore_cnn_img_tflite(int qindex, const uint8_t *dgd,
                                          int width, int height, int dgd_stride,
                                          uint8_t *rst, int rst_stride) {
  std::unique_ptr<tflite::Interpreter> interpreter =
      get_tflite_interpreter(qindex, width, height);

  // Prepare input.
  const float max_val = 255.0f;
  const int in_stride = width;
  auto input = interpreter->typed_input_tensor<float>(0);
  for (int r = 0; r < height; ++r) {
    for (int c = 0; c < width; ++c) {
      input[r * in_stride + c] =
          static_cast<float>(dgd[r * dgd_stride + c]) / max_val;
      assert(input[r * in_stride + c] >= 0.0f);
      assert(input[r * in_stride + c] <= 1.0f);
    }
  }

  // Invoke TFlite inference.
  tflite::ErrorReporter *reporter = tflite::DefaultErrorReporter();
  auto status = interpreter->Invoke();
  if (status != kTfLiteOk) {
    reporter->Report("Failed at interpreter invocation");
    return 0;
  }

  // Use the output to restore 'dgd' and store in 'rst'.
  const auto output = interpreter->typed_output_tensor<float>(0);
  const int out_stride = width;
  // TODO(urvang): Try using the model which directly gives clip-by-value and
  // multiply.
  for (int r = 0; r < height; ++r) {
    for (int c = 0; c < width; ++c) {
      const int residue =
          static_cast<int>(output[r * out_stride + c] * max_val + 0.5);
      rst[r * rst_stride + c] = clip_pixel(dgd[r * dgd_stride + c] + residue);
    }
  }
  return 1;
}

extern "C" int av1_restore_cnn_img_tflite_highbd(int qindex,
                                                 const uint16_t *dgd, int width,
                                                 int height, int dgd_stride,
                                                 uint16_t *rst, int rst_stride,
                                                 int bit_depth) {
  assert(av1_use_cnn_tflite(qindex));
  std::unique_ptr<tflite::Interpreter> interpreter =
      get_tflite_interpreter(qindex, width, height);

  // Prepare input.
  const auto max_val = static_cast<float>((1 << bit_depth) - 1);
  const int in_stride = width;
  auto input = interpreter->typed_input_tensor<float>(0);
  for (int r = 0; r < height; ++r) {
    for (int c = 0; c < width; ++c) {
      input[r * in_stride + c] =
          static_cast<float>(dgd[r * dgd_stride + c]) / max_val;
      assert(input[r * in_stride + c] >= 0.0f);
      assert(input[r * in_stride + c] <= 1.0f);
    }
  }

  // Invoke TFlite inference.
  tflite::ErrorReporter *reporter = tflite::DefaultErrorReporter();
  auto status = interpreter->Invoke();
  if (status != kTfLiteOk) {
    reporter->Report("Failed at interpreter invocation");
    return 0;
  }

  // Use the output to restore 'dgd' and store in 'rst'.
  const auto output = interpreter->typed_output_tensor<float>(0);
  const int out_stride = width;
  for (int r = 0; r < height; ++r) {
    for (int c = 0; c < width; ++c) {
      const int residue =
          static_cast<int>(output[r * out_stride + c] * max_val + 0.5);
      rst[r * rst_stride + c] =
          clip_pixel_highbd(dgd[r * dgd_stride + c] + residue, bit_depth);
    }
  }
  return 1;
}

/*
extern "C" void av1_restore_cnn_plane_tflite(const AV1_COMMON *cm, int plane,
                                  int is_residue) {
  YV12_BUFFER_CONFIG *buf = &cm->cur_frame->buf;
  if (cm->seq_params.use_highbitdepth) {
    switch (plane) {
      case AOM_PLANE_Y:
        av1_restore_cnn_img_tflite_highbd(
            cm->base_qindex, CONVERT_TO_SHORTPTR(buf->y_buffer),
            buf->y_crop_width, buf->y_crop_height, buf->y_stride,
            cm->seq_params.bit_depth, is_residue,
            CONVERT_TO_SHORTPTR(buf->y_buffer), buf->y_stride);
        break;
      case AOM_PLANE_U:
        av1_restore_cnn_img_tflite_highbd(
            cm->base_qindex, CONVERT_TO_SHORTPTR(buf->u_buffer),
            buf->uv_crop_width, buf->uv_crop_height, buf->uv_stride,
            cm->seq_params.bit_depth, is_residue,
            CONVERT_TO_SHORTPTR(buf->u_buffer), buf->uv_stride);
        break;
      case AOM_PLANE_V:
        av1_restore_cnn_img_tflite_highbd(
            cm->base_qindex, CONVERT_TO_SHORTPTR(buf->v_buffer),
            buf->uv_crop_width, buf->uv_crop_height, buf->uv_stride,
            cm->seq_params.bit_depth, is_residue,
            CONVERT_TO_SHORTPTR(buf->u_buffer), buf->uv_stride);
        break;
      default: assert(0 && "Invalid plane index");
    }
  } else {
    assert(cm->seq_params.bit_depth == 8);
    switch (plane) {
      case AOM_PLANE_Y:
        av1_restore_cnn_img_tflite(cm->base_qindex, buf->y_buffer,
                                   buf->y_crop_width, buf->y_crop_height,
                                   buf->y_stride, is_residue, buf->y_buffer,
                                   buf->y_stride);
        break;
      case AOM_PLANE_U:
        av1_restore_cnn_img_tflite(cm->base_qindex, buf->u_buffer,
                                   buf->uv_crop_width, buf->uv_crop_height,
                                   buf->uv_stride, is_residue, buf->u_buffer,
                                   buf->uv_stride);
        break;
      case AOM_PLANE_V:
        av1_restore_cnn_img_tflite(cm->base_qindex, buf->v_buffer,
                                   buf->uv_crop_width, buf->uv_crop_height,
                                   buf->uv_stride, is_residue, buf->v_buffer,
                                   buf->uv_stride);
        break;
      default: assert(0 && "Invalid plane index");
    }
  }
}
*/