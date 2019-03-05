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

#include <stdio.h>

#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/public/session.h"

int print_that_status() {
  tensorflow::Session *session;
  tensorflow::Status status =
      NewSession(tensorflow::SessionOptions(), &session);
  if (!status.ok()) {
    printf("Session failed to create: %s\n", status.ToString());
    return 1;
  }
  printf("Session successfully created.\n");
  return 0;
}
