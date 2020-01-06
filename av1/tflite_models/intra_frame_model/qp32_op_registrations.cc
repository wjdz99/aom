#include "qp32_op_registrations.h"

#include "tensorflow/lite/kernels/builtin_op_kernels.h"
#include "tensorflow/lite/model.h"

void RegisterSelectedOpsQp32(::tflite::MutableOpResolver *resolver) {
  resolver->AddBuiltin(::tflite::BuiltinOperator_ADD,
                       ::tflite::ops::builtin::Register_ADD());
  resolver->AddBuiltin(::tflite::BuiltinOperator_CONV_2D,
                       ::tflite::ops::builtin::Register_CONV_2D());
  resolver->AddBuiltin(::tflite::BuiltinOperator_DEPTHWISE_CONV_2D,
                       ::tflite::ops::builtin::Register_DEPTHWISE_CONV_2D());
  resolver->AddBuiltin(::tflite::BuiltinOperator_MAXIMUM,
                       ::tflite::ops::builtin::Register_MAXIMUM());
  resolver->AddBuiltin(::tflite::BuiltinOperator_MINIMUM,
                       ::tflite::ops::builtin::Register_MINIMUM());
}
