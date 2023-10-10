/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_DSP_MATHUTILS_H_
#define AOM_AOM_DSP_MATHUTILS_H_

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "aom_dsp/aom_dsp_common.h"
#include "aom_mem/aom_mem.h"

static const double TINY_NEAR_ZERO = 1.0E-16;

// Solves the matrix equation AX = B, where X and B are matrices of size n x m
// and A is n x n. This can be used to solve multiple systems of linear
// equations with the same left-hand-side (A matrix) at the same time, with less
// computational cost than solving each problem independently.
//
// The argument pivot determines whether we use pivoting. This improves
// numerical stability by rearranging A and B at each step to maximize the
// diagonal elements of A, but at the cost of additional compute time.
// In some cases, this may not be necessary - in particular, least squares
// matrices already have A[i][i] >= A[i][j] for all i, j, by construction.
// This property is not guaranteed to be maintained throughout solving,
// but it means that least squares matrices may not need pivoting,
// saving some compute time.
static INLINE int linsolve_multi(int n, int m, double *A, int A_stride,
                                 double *B, int B_stride, double *X,
                                 int X_stride, bool pivot) {
  int i, j, k;
  double c;
  // Forward elimination
  for (k = 0; k < n - 1; k++) {
    if (pivot) {
      // Bring the largest magnitude to the diagonal position
      for (i = n - 1; i > k; i--) {
        if (fabs(A[(i - 1) * A_stride + k]) < fabs(A[i * A_stride + k])) {
          // Swap rows i and i-1
          for (j = 0; j < n; j++) {
            c = A[i * A_stride + j];
            A[i * A_stride + j] = A[(i - 1) * A_stride + j];
            A[(i - 1) * A_stride + j] = c;
          }
          for (j = 0; j < m; j++) {
            c = B[i * B_stride + j];
            B[i * B_stride + j] = B[(i - 1) * B_stride + j];
            B[(i - 1) * B_stride + j] = c;
          }
        }
      }
    }

    // Subtract this row from all subsequent rows
    for (i = k; i < n - 1; i++) {
      if (fabs(A[k * A_stride + k]) < TINY_NEAR_ZERO) return 0;
      c = A[(i + 1) * A_stride + k] / A[k * A_stride + k];
      for (j = 0; j < n; j++)
        A[(i + 1) * A_stride + j] -= c * A[k * A_stride + j];
      for (j = 0; j < m; j++)
        B[(i + 1) * B_stride + j] -= c * B[k * B_stride + j];
    }
  }

  // Backward substitution
  for (i = n - 1; i >= 0; i--) {
    if (fabs(A[i * A_stride + i]) < TINY_NEAR_ZERO) return 0;
    for (j = 0; j < m; j++) {
      c = 0;
      for (k = i + 1; k <= n - 1; k++)
        c += A[i * A_stride + k] * X[k * X_stride + j];
      X[i * X_stride + j] = (B[i * B_stride + j] - c) / A[i * A_stride + i];
    }
  }

  return 1;
}

static INLINE int linsolve(int n, double *A, int stride, double *b, double *x) {
  return linsolve_multi(n, 1, A, stride, b, 1, x, 1, true);
}

////////////////////////////////////////////////////////////////////////////////
// Least-squares
// Solves for n-dim x in a least squares sense to minimize |Ax - b|^2
// The solution is simply x = (A'A)^-1 A'b or simply the solution for
// the system: A'A x = A'b
//
// This process is split into three steps in order to avoid needing to
// explicitly allocate the A matrix, which may be very large if there
// are many equations to solve.
//
// We can also allow multiple right-hand-sides to be solved simultaneously,
// to reduce computational complexity
//
// The process for using this is (in pseudocode):
//
// Allocate mat (size n*n), y (size n*m), a (size n), x (size nxm)
// least_squares_init(mat, y, n, m)
// for each (set of) equation(s) a . x = b {
//    least_squares_accumulate(mat, y, a, b, n, m)
// }
// least_squares_solve(mat, y, x, n, m)
//
// where:
// * mat, y are accumulators for the values A'A and A'b respectively,
// * a, b are the coefficients of each individual equation,
// * x is the result vector
// * n is the problem size
// * m is the number of simultaneous problems
//
// The solution to the j'th problem is stored in the j'th column of x.
// In other words, the solution to the j'th problem is
//   { x[j], x[j+m], ... x[(n-1)*m + j] }
static INLINE void least_squares_init(double *mat, double *y, int n, int m) {
  memset(mat, 0, n * n * sizeof(double));
  memset(y, 0, n * m * sizeof(double));
}

// Round the given positive value to nearest integer
static AOM_FORCE_INLINE int iroundpf(float x) {
  assert(x >= 0.0);
  return (int)(x + 0.5f);
}

static INLINE void least_squares_accumulate(double *mat, double *y,
                                            const double *a, double *b, int n,
                                            int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      mat[i * n + j] += a[i] * a[j];
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      y[i * m + j] += a[i] * b[j];
    }
  }
}

static INLINE int least_squares_solve(double *mat, double *y, double *x, int n,
                                      int m) {
  return linsolve_multi(n, m, mat, n, y, m, x, m, true);
}

// Matrix multiply
static INLINE void multiply_mat(const double *m1, const double *m2, double *res,
                                const int m1_rows, const int inner_dim,
                                const int m2_cols) {
  double sum;

  int row, col, inner;
  for (row = 0; row < m1_rows; ++row) {
    for (col = 0; col < m2_cols; ++col) {
      sum = 0;
      for (inner = 0; inner < inner_dim; ++inner)
        sum += m1[row * inner_dim + inner] * m2[inner * m2_cols + col];
      *(res++) = sum;
    }
  }
}

static AOM_INLINE float approx_exp(float y) {
#define A ((1 << 23) / 0.69314718056f)  // (1 << 23) / ln(2)
#define B \
  127  // Offset for the exponent according to IEEE floating point standard.
#define C 60801  // Magic number controls the accuracy of approximation
  union {
    float as_float;
    int32_t as_int32;
  } container;
  container.as_int32 = ((int32_t)(y * A)) + ((B << 23) - C);
  return container.as_float;
#undef A
#undef B
#undef C
}
#endif  // AOM_AOM_DSP_MATHUTILS_H_
