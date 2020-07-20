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

#ifndef AOM_COMMON_TF_LITE_INCLUDES_H_
#define AOM_COMMON_TF_LITE_INCLUDES_H_

// TensorFlow Lite has several unused parameters that are
// exposed as part of the API. In the AOM build process, this
// will cause failures when -Wunused-parameter is set.
// Since TF Lite is external code, instruct the compiler to
// ignore this warning when including it.
// Note that Clang supports this GCC pragma.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "third_party/tensorflow/tensorflow/lite/interpreter.h"
#include "third_party/tensorflow/tensorflow/lite/model.h"

#pragma GCC diagnostic pop

#endif  // AOM_COMMON_TF_LITE_INCLUDES_H_
