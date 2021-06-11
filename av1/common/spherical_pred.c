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
  double phi_mod = fmod(phi, 0.5 * PI);
  if (phi_mod == 0 && fmod(phi, PI) != 0) {
    if (fmod(phi, 2 * PI) == -1.5 * PI) {
      phi_mod = 0.5 * PI;
    } else if (fmod(phi, 2 * PI) == 1.5 * PI) {
      phi_mod = -0.5 * PI;
    } else {
      phi_mod = phi > 0 ? 0.5 * PI : -0.5 * PI;
    }
  }

  double theta_mod = fmod(theta, PI);
  if (theta_mod == 0 && fmod(theta, 2 * PI) != 0) {
    if (theta > 0) {
      theta_mod = PI;
    } else if (theta < 0) {
      theta_mod = -PI;
    }
  }
  // This should actually be a range related to 1/cos(phi) since x is distorted
  // TODO(yaoyaogoogle): Adjust the width of x according to the interpolation
  // mode
  *x = theta_mod / PI * (width - 1) * 0.5 + (width - 1) * 0.5;

  // No minus sign for y since we only use an imaginary upside-down globe
  *y = phi_mod / (PI * 0.5) * (height - 1) * 0.5 + (height - 1) * 0.5;
}

void av1_plane_to_sphere_erp(double x, double y, int width, int height,
                             double *phi, double *theta) {
  assert(y >= 0 && y <= height - 1);

  double x_mod = fmod(x, width);
  x_mod = x_mod >= 0 ? x_mod : x_mod + width;

  *theta = (x_mod - (width - 1) * 0.5) / (width - 1) * 2 * PI;
  *phi = (y - (height - 1) * 0.5) / (height - 1) * PI;
}
