/*
 * Copyright (c) 2019, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

/*!\file
 * \brief This is an sample binary showing how to load a TF-lite model.
 *
 * 1. Build your model like normal.
 * 2. Follow the steps at https://www.tensorflow.org/lite/convert/python_api
 *    to convert into a TF-lite model.
 * 3. Run `xxd -i model.tflite > model.c` to make it a CC file.
 * 4. Change the declaration to be const.
 * 4. Create a .h file that exposes the array and length.
 * 5. Add appropriate copyright headers.
 */

#include <cstdio>
#include <cstdlib>
#include <memory>

#include "examples/tf_lite_model_data.h"
#include "third_party/tensorflow/tensorflow/lite/interpreter.h"
#include "third_party/tensorflow/tensorflow/lite/kernels/register.h"
#include "third_party/tensorflow/tensorflow/lite/model.h"

int main(int argc, char *argv[]) {
  // Model that takes in 64 floats and returns 1 float.
  auto model = tflite::GetModel(tf_lite_model_data);

  // Build the interpreter.
  tflite::ops::builtin::BuiltinOpResolver resolver;
  tflite::InterpreterBuilder builder(model, resolver);
  std::unique_ptr<tflite::Interpreter> interpreter;
  builder(&interpreter);

  std::unique_ptr<tflite::ErrorReporter> reporter(
      tflite::DefaultErrorReporter());

  if (interpreter->AllocateTensors() != kTfLiteOk) {
    reporter->Report("Failed");
    return EXIT_FAILURE;
  }

  float *input = interpreter->typed_tensor<float>(0);
  for (int i = 0; i < 64; ++i) {
    input[i] = i;
  }
  auto status = interpreter->Invoke();
  if (status != kTfLiteOk) {
    reporter->Report("Failed");
    return EXIT_FAILURE;
  }
  float y = interpreter->typed_output_tensor<float>(0)[0];
  printf("Sample output: %f\n", y);

  return EXIT_SUCCESS;
}
