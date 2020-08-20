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

#include <memory>

#include "av1/common/interintra_ml.h"
#include "av1/common/interintra_ml_model.h"
#include "av1/common/reconintra.h"
#include "common/tf_lite_includes.h"

namespace {

// Assumes entire program is single-threaded.
tflite::Interpreter *interpreter_ = nullptr;
tflite::ErrorReporter *reporter_ = nullptr;

void add_resolver_builtins(::tflite::MutableOpResolver *resolver) {
  resolver->AddBuiltin(::tflite::BuiltinOperator_ADD,
                       ::tflite::ops::builtin::Register_ADD());
  resolver->AddBuiltin(::tflite::BuiltinOperator_CAST,
                       ::tflite::ops::builtin::Register_CAST());
  resolver->AddBuiltin(::tflite::BuiltinOperator_CONCATENATION,
                       ::tflite::ops::builtin::Register_CONCATENATION());
  resolver->AddBuiltin(::tflite::BuiltinOperator_CONV_2D,
                       ::tflite::ops::builtin::Register_CONV_2D());
  resolver->AddBuiltin(::tflite::BuiltinOperator_EQUAL,
                       ::tflite::ops::builtin::Register_EQUAL());
  resolver->AddBuiltin(::tflite::BuiltinOperator_FILL,
                       ::tflite::ops::builtin::Register_FILL());
  resolver->AddBuiltin(::tflite::BuiltinOperator_GATHER,
                       ::tflite::ops::builtin::Register_GATHER());
  resolver->AddBuiltin(::tflite::BuiltinOperator_IF,
                       ::tflite::ops::builtin::Register_IF());
  resolver->AddBuiltin(::tflite::BuiltinOperator_LEAKY_RELU,
                       ::tflite::ops::builtin::Register_LEAKY_RELU());
  resolver->AddBuiltin(::tflite::BuiltinOperator_LESS,
                       ::tflite::ops::builtin::Register_LESS());
  resolver->AddBuiltin(::tflite::BuiltinOperator_LOGICAL_AND,
                       ::tflite::ops::builtin::Register_LOGICAL_AND());
  resolver->AddBuiltin(::tflite::BuiltinOperator_RESHAPE,
                       ::tflite::ops::builtin::Register_RESHAPE());
  resolver->AddBuiltin(::tflite::BuiltinOperator_SHAPE,
                       ::tflite::ops::builtin::Register_SHAPE());
  resolver->AddBuiltin(::tflite::BuiltinOperator_SLICE,
                       ::tflite::ops::builtin::Register_SLICE());
  resolver->AddBuiltin(::tflite::BuiltinOperator_STRIDED_SLICE,
                       ::tflite::ops::builtin::Register_STRIDED_SLICE());
  resolver->AddBuiltin(::tflite::BuiltinOperator_TRANSPOSE,
                       ::tflite::ops::builtin::Register_TRANSPOSE());
  resolver->AddBuiltin(::tflite::BuiltinOperator_UNPACK,
                       ::tflite::ops::builtin::Register_UNPACK(), 3, 3);
  resolver->AddBuiltin(::tflite::BuiltinOperator_WHILE,
                       ::tflite::ops::builtin::Register_WHILE());
}

// If the error report has not been initialized, initialize it. Otherwise,
// return it.
tflite::ErrorReporter *get_or_init_reporter() {
  if (reporter_ != nullptr) {
    return reporter_;
  }
  reporter_ = tflite::DefaultErrorReporter();
  return reporter_;
}

// If the interpreter has not been initialized yet, initializes it and returns
// it. Otherwise, returns the initialized interpreter.
tflite::Interpreter *get_or_init_interpreter() {
  if (interpreter_ != nullptr) {
    return interpreter_;
  }

  auto model = tflite::GetModel(decode_13759197_5_tflite_data);
  tflite::MutableOpResolver resolver;
  add_resolver_builtins(&resolver);
  tflite::InterpreterBuilder builder(model, resolver);
  std::unique_ptr<tflite::Interpreter> interpreter;
  builder(&interpreter);
  tflite::ErrorReporter *reporter = get_or_init_reporter();

  if (interpreter->AllocateTensors() != kTfLiteOk) {
    reporter->Report("Failed");
    assert(false);
  }
  interpreter_ = interpreter.release();
  return interpreter_;
}

// Copy a blank square into the region. Needed as default behavior if
// the interintra ML model does not support a particular use case.
void copy_blank_square(uint8_t *dst, int stride, BLOCK_SIZE bsize,
                       bool is_hbd) {
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  for (int j = 0; j < bh; ++j) {
    av1_bd_memset(dst + j * stride, 0, bw, is_hbd);
  }
}

// Load the inputs (inter-predictor + border, intra-predictor border)
// into the interpreter.
void load_inputs(tflite::Interpreter *interpreter, INTERINTRA_MODE mode,
                 BLOCK_SIZE bsize, const uint8_t *inter_pred, int inter_stride,
                 const uint8_t *intra_pred, int intra_stride) {
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];

  // Load the inter-predictor and border.
  int *mode_input = interpreter->typed_input_tensor<int>(0);
  *mode_input = mode - II_ML_PRED0;  // Normalize so 0 is the first mode.
  int *inter_input = interpreter->typed_input_tensor<int>(1);
  inter_pred -= INTERINTRA_ML_BORDER * (1 + inter_stride);
  for (int j = 0; j < bh + INTERINTRA_ML_BORDER; ++j) {
    for (int i = 0; i < bw + INTERINTRA_ML_BORDER; ++i) {
      inter_input[i + j * (bw + INTERINTRA_ML_BORDER)] =
          inter_pred[i + j * inter_stride];
    }
  }

  // Load the top-part of the intra-predictor border.
  int *intra_input = interpreter->typed_input_tensor<int>(2);
  int interpreter_offset = 0;
  intra_pred -= INTERINTRA_ML_BORDER * (1 + intra_stride);
  for (int j = 0; j < INTERINTRA_ML_BORDER; ++j) {
    for (int i = 0; i < bw + INTERINTRA_ML_BORDER; ++i) {
      intra_input[interpreter_offset] = intra_pred[i + j * intra_stride];
      ++interpreter_offset;
    }
  }

  // Load the left columns of the intra-predictor border.
  intra_pred += INTERINTRA_ML_BORDER * intra_stride;
  for (int j = 0; j < bh; ++j) {
    for (int i = 0; i < INTERINTRA_ML_BORDER; ++i) {
      intra_input[interpreter_offset] = intra_pred[i + j * intra_stride];
      ++interpreter_offset;
    }
  }
}

// Copy the output of the interpreter into the destination buffer.
void copy_to_output(tflite::Interpreter *interpreter, BLOCK_SIZE bsize,
                    uint8_t *comp_pred, int comp_stride) {
  const int bw = block_size_wide[bsize];
  const int bh = block_size_high[bsize];
  for (int j = 0; j < bh; ++j) {
    for (int i = 0; i < bw; ++i) {
      comp_pred[i + j * comp_stride] =
          interpreter->typed_output_tensor<int>(0)[i + j * bw];
    }
  }
}

}  // namespace

void av1_combine_interintra_ml(INTERINTRA_MODE mode, BLOCK_SIZE plane_bsize,
                               uint8_t *comp_pred, int comp_stride,
                               const uint8_t *inter_pred, int inter_stride,
                               const uint8_t *intra_pred, int intra_stride,
                               int border) {
  (void)border;
  assert(border >= INTERINTRA_ML_BORDER);
  if (plane_bsize != BLOCK_16X16) {
    // Not yet implemented. Just copy a blank square into the predictor.
    copy_blank_square(comp_pred, comp_stride, plane_bsize, false);
    return;
  }
  tflite::Interpreter *interpreter = get_or_init_interpreter();
  load_inputs(interpreter, mode, plane_bsize, inter_pred, inter_stride,
              intra_pred, intra_stride);
  auto status = interpreter->Invoke();
  if (status != kTfLiteOk) {
    tflite::ErrorReporter *reporter = get_or_init_reporter();
    reporter->Report("Failed");
    assert(false);
  }
  copy_to_output(interpreter, plane_bsize, comp_pred, comp_stride);
}

void av1_combine_interintra_ml_highbd(
    INTERINTRA_MODE mode, BLOCK_SIZE plane_bsize, uint8_t *comp_pred8,
    int comp_stride, const uint8_t *inter_pred8, int inter_stride,
    const uint8_t *intra_pred8, int intra_stride, int bd, int border) {
  (void)mode;
  (void)inter_pred8;
  (void)inter_stride;
  (void)intra_pred8;
  (void)intra_stride;
  (void)bd;
  (void)border;
  assert(border >= INTERNTRA_ML_BORDER);
  // Not yet implemented. Just copy a blank square into the predictor.
  copy_blank_square(comp_pred8, comp_stride, plane_bsize, true);
}
