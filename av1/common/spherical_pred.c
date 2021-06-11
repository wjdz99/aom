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
#include <math.h>

#include "av1/common/common.h"
#include "av1/common/spherical_pred.h"

void av1_sphere_to_plane_erp(double phi, double theta, int width, int height,
                             double *x, double *y) {
  // The 0.01 is for issues when comparing floating point numbers
  // TODO(yaoyaogoogle): Find a better way to solve the issue
  assert(phi >= -PI * 0.5 - 0.01 && phi <= PI * 0.5 + 0.01 &&
         theta >= -PI - 0.01 && theta <= PI + 0.01);

  // This should actually be a range related to 1/cos(phi) since x is distorted
  // TODO(yaoyaogoogle): Adjust the width of x according to the interpolation
  // mode
  *x = theta / PI * width * 0.5 + width * 0.5;

  // Loop x to the other side if out of range
  if (*x < 0) {
    *x = *x + width;
  } else if (*x > width) {
    *x = *x - width;
  }
  // No minus sign for y since we only use an imaginary upside-down globe
  *y = phi / (PI * 0.5) * height * 0.5 + height * 0.5;

  // Limit y to the border if out of range
  if (*y < 0) {
    *y = 0;
  } else if (*y > height) {
    *y = height;
  }
}

void av1_plane_to_sphere_erp(double x, double y, int width, int height,
                             double *phi, double *theta) {
  assert(x < width && x >= 0 && y < height && y >= 0);

  *theta = (x - width * 0.5) / width * 2 * PI;

  // Loop theta to the other side if out of range
  if (*theta < -PI) {
    *theta = *theta + 2 * PI;
  } else if (*theta > PI) {
    *theta = *theta - 2 * PI;
  }
  *phi = (y - height * 0.5) / height * PI;

  // Limit phi to the border if out of range
  if (*phi < -0.5 * PI) {
    *phi = -0.5 * PI;
  } else if (*phi > 0.5 * PI) {
    *phi = 0.5 * PI;
  }
}
