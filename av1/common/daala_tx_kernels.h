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

/* clang-format off */

#ifndef AOM_DSP_DAALA_TX_KERNELS_H_
#define AOM_DSP_DAALA_TX_KERNELS_H_

#include "aom_dsp/aom_dsp_common.h"
#include "av1/common/odintrin.h"

#define AVG_BIAS (0)

static INLINE od_coeff od_add(od_coeff p0, od_coeff p1) {
  return p0 + p1;
}

static INLINE od_coeff od_sub(od_coeff p0, od_coeff p1) {
  return p0 - p1;
}

static INLINE od_coeff od_add_avg(od_coeff p0, od_coeff p1) {
  return (od_add(p0, p1) + AVG_BIAS) >> 1;
}

static INLINE od_coeff od_sub_avg(od_coeff p0, od_coeff p1) {
  return (od_sub(p0, p1) + AVG_BIAS) >> 1;
}

static INLINE od_coeff od_rshift1(od_coeff v) {
  return (v + (v < 0)) >> 1;
}

/* Fixed point multiply. */
static INLINE od_coeff od_mul(od_coeff n, int c, int q) {
  return (n*c + ((1 << q) >> 1)) >> q;
}

/* Two multiply rotation primative (used when rotating by Pi/4). */
static INLINE void od_rot2(od_coeff *p0, od_coeff *p1, od_coeff t, int c0,
 int q0, int c1, int q1) {
  *p1 = od_mul(*p0, c0, q0);
  *p0 = od_mul(t, c1, q1);
}

/* Three multiply rotation primative. */
static INLINE void od_rot3(od_coeff *p0, od_coeff *p1, od_coeff *t, od_coeff *u,
 int c0, int q0, int c1, int q1, int c2, int q2) {
  *u = od_mul(*p0, c0, q0);
  *p0 = od_mul(*p1, c1, q1);
  *t = od_mul(*t, c2, q2);
}

#define NONE (0)
#define AVG (!NONE)
#define SHIFT (!NONE)

#define ADD (0)
#define SUB (1)

/* Rotate by Pi/4 and add. */
static INLINE void od_rotate_pi4_kernel(od_coeff *p0, od_coeff *p1, int c0,
 int q0, int c1, int q1, int type, int avg) {
  od_coeff t;
  t = type == ADD ?
   avg ? od_add_avg(*p1, *p0) : od_add(*p1, *p0) :
   avg ? od_sub_avg(*p1, *p0) : od_sub(*p1, *p0);
  od_rot2(p0, p1, t, c0, q0, c1, q1);
  *p1 = type == ADD ? od_sub(*p1, *p0) : od_add(*p1, *p0);
}

#define od_rotate_pi4_add(p0, p1, c0, q0, c1, q1) \
 od_rotate_pi4_kernel(p0, p1, c0, q0, c1, q1, ADD, NONE)
#define od_rotate_pi4_sub(p0, p1, c0, q0, c1, q1) \
 od_rotate_pi4_kernel(p0, p1, c0, q0, c1, q1, SUB, NONE)

#define od_rotate_pi4_add_avg(p0, p1, c0, q0, c1, q1) \
 od_rotate_pi4_kernel(p0, p1, c0, q0, c1, q1, ADD, AVG)
#define od_rotate_pi4_sub_avg(p0, p1, c0, q0, c1, q1) \
 od_rotate_pi4_kernel(p0, p1, c0, q0, c1, q1, SUB, AVG)

/* Rotate and add. */
static INLINE void od_rotate_kernel(od_coeff *p0, od_coeff *p1, od_coeff v,
 int c0, int q0, int c1, int q1, int c2, int q2, int type, int avg, int shift) {
  od_coeff u;
  od_coeff t;
  t = type == ADD ?
   avg ? od_add_avg(*p1, v) : od_add(*p1, v) :
   avg ? od_sub_avg(*p1, v) : od_sub(*p1, v);
  od_rot3(p0, p1, &t, &u, c0, q0, c1, q1, c2, q2);
  *p0 = od_add(*p0, t);
  if (shift) t = od_rshift1(t);
  *p1 = type == ADD ? od_sub(u, t) : od_add(u, t);
}

#define od_rotate_add(p0, p1, c0, q0, c1, q1, c2, q2, shift) \
 od_rotate_kernel(p0, p1, *p0, c0, q0, c1, q1, c2, q2, ADD, NONE, shift)
#define od_rotate_sub(p0, p1, c0, q0, c1, q1, c2, q2, shift) \
 od_rotate_kernel(p0, p1, *p0, c0, q0, c1, q1, c2, q2, SUB, NONE, shift)

#define od_rotate_add_avg(p0, p1, c0, q0, c1, q1, c2, q2, shift) \
 od_rotate_kernel(p0, p1, *p0, c0, q0, c1, q1, c2, q2, ADD, AVG, shift)
#define od_rotate_sub_avg(p0, p1, c0, q0, c1, q1, c2, q2, shift) \
 od_rotate_kernel(p0, p1, *p0, c0, q0, c1, q1, c2, q2, SUB, AVG, shift)

#define od_rotate_add_half(p0, p1, v, c0, q0, c1, q1, c2, q2, shift) \
 od_rotate_kernel(p0, p1, v, c0, q0, c1, q1, c2, q2, ADD, NONE, shift)
#define od_rotate_sub_half(p0, p1, v, c0, q0, c1, q1, c2, q2, shift) \
 od_rotate_kernel(p0, p1, v, c0, q0, c1, q1, c2, q2, SUB, NONE, shift)

/* Rotate and subtract with negation. */
static INLINE void od_rotate_neg_kernel(od_coeff *p0, od_coeff *p1,
 int c0, int q0, int c1, int q1, int c2, int q2, int avg) {
  od_coeff u;
  od_coeff t;
  t = avg ? od_sub_avg(*p0, *p1) : od_sub(*p0, *p1);
  od_rot3(p0, p1, &t, &u, c0, q0, c1, q1, c2, q2);
  *p0 = od_sub(*p0, t);
  *p1 = od_sub(t, u);
}

#define od_rotate_neg(p0, p1, c0, q0, c1, q1, c2, q2) \
 od_rotate_neg_kernel(p0, p1, c0, q0, c1, q1, c2, q2, NONE)
#define od_rotate_neg_avg(p0, p1, c0, q0, c1, q1, c2, q2) \
 od_rotate_neg_kernel(p0, p1, c0, q0, c1, q1, c2, q2, AVG)

/* Computes the +/- addition butterfly (asymmetric output).
   The inverse to this function is od_butterfly_add_asym().

    p0 = p0 + p1;
    p1 = p1 - p0/2; */
static INLINE void od_butterfly_add(od_coeff *p0, od_coeff *p0h, od_coeff *p1) {
  od_coeff p0h_;
  *p0 = od_add(*p0, *p1);
  p0h_ = od_rshift1(*p0);
  *p1 = od_sub(*p1, p0h_);
  if (p0h != NULL) *p0h = p0h_;
}

/* Computes the asymmetric +/- addition butterfly (unscaled output).
   The inverse to this function is od_butterfly_add().

    p1 = p1 + p0/2;
    p0 = p0 - p1; */
static INLINE void od_butterfly_add_asym(od_coeff *p0, od_coeff p0h,
 od_coeff *p1) {
  *p1 = od_add(*p1, p0h);
  *p0 = od_sub(*p0, *p1);
}

/* Computes the +/- subtraction butterfly (asymmetric output).
   The inverse to this function is od_butterfly_sub_asym().

    p0 = p0 - p1;
    p1 = p1 + p0/2; */
static INLINE void od_butterfly_sub(od_coeff *p0, od_coeff *p0h, od_coeff *p1) {
  od_coeff p0h_;
  *p0 = od_sub(*p0, *p1);
  p0h_ = od_rshift1(*p0);
  *p1 = od_add(*p1, p0h_);
  if (p0h != NULL) *p0h = p0h_;
}

/* Computes the asymmetric +/- subtraction butterfly (unscaled output).
   The inverse to this function is od_butterfly_sub().

    p1 = p1 - p0/2;
    p0 = p0 + p1; */
static INLINE void od_butterfly_sub_asym(od_coeff *p0, od_coeff p0h,
 od_coeff *p1) {
  *p1 = od_sub(*p1, p0h);
  *p0 = od_add(*p0, *p1);
}

/* Computes the +/- subtract and negate butterfly (asymmetric output).
   The inverse to this function is od_butterfly_neg_asym().

    p1 = p1 - p0;
    p0 = p0 + p1/2;
    p1 = -p1; */
static INLINE void od_butterfly_neg(od_coeff *p0, od_coeff *p1, od_coeff *p1h) {
  *p1 = od_sub(*p0, *p1);
  *p1h = od_rshift1(*p1);
  *p0 = od_sub(*p0, *p1h);
}

/*  Computes the asymmetric +/- negate and subtract butterfly (unscaled output).
    The inverse to this function is od_butterfly_neg().

    p1 = -p1;
    p0 = p0 - p1/2;
    p1 = p1 + p0; */
static INLINE void od_butterfly_neg_asym(od_coeff *p0, od_coeff *p1,
 od_coeff p1h) {
  *p0 = od_add(*p0, p1h);
  *p1 = od_sub(*p0, *p1);
}

/* --- 2-point Transforms --- */

/**
 * 2-point orthonormal Type-II fDCT
 */
static INLINE void od_fdct_2(od_coeff *p0, od_coeff *p1) {
  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4]  = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]            = 1.4142135623730951 */
  od_rotate_pi4_sub_avg(p1, p0, 11585, 13, 11585, 13);
}

/**
 * 2-point orthonormal Type-II iDCT
 */
static INLINE void od_idct_2(od_coeff *p0, od_coeff *p1) {
  /*  11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/16384 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(p0, p1, 11585, 13, 11585, 14);
}

/**
 * 2-point asymmetric Type-II fDCT
 */
static INLINE void od_fdct_2_asym(od_coeff *p0, od_coeff *p1, od_coeff p1h) {
  od_butterfly_neg_asym(p0, p1, p1h);
}

/**
 * 2-point asymmetric Type-II iDCT
 */
static INLINE void od_idct_2_asym(od_coeff *p0, od_coeff *p1, od_coeff *p1h) {
  od_butterfly_neg(p0, p1, p1h);
}

/**
 * 2-point orthonormal Type-IV fDCT
 */
static INLINE void od_fdst_2(od_coeff *p0, od_coeff *p1) {

  /* Stage 0 */

  /* 10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8]  = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8]  = 0.5411961001461971 */
  /*  3135/4096 = 2*Cos[3*Pi/8]              = 0.7653668647301796 */
  od_rotate_add_avg(p0, p1, 10703, 13, 8867, 14, 3135, 12, NONE);
}

/**
 * 2-point orthonormal Type-IV iDCT
 */
static INLINE void od_idst_2(od_coeff *p0, od_coeff *p1) {
  od_fdst_2(p0, p1);
}

/**
 * 2-point asymmetric Type-IV fDCT
 */
static INLINE void od_fdst_2_asym(od_coeff *p0, od_coeff p0h, od_coeff *p1) {

  /* Stage 0 */

  /*   473/512 = (Sin[3*Pi/8] + Cos[3*Pi/8])/Sqrt[2] = 0.9238795325112867 */
  /* 3135/4096 = (Sin[3*Pi/8] - Cos[3*Pi/8])*Sqrt[2] = 0.7653668647301795 */
  /* 4433/8192 = Cos[3*Pi/8]*Sqrt[2]                 = 0.5411961001461971 */
  od_rotate_add_half(p0, p1, p0h, 473, 9, 3135, 12, 4433, 13, NONE);
}

/**
 * 2-point asymmetric Type-IV iDCT
 */
static INLINE void od_idst_2_asym(od_coeff *p0, od_coeff *p1) {

  /* Stage 0 */

  /*   473/512 = (Sin[3*Pi/8] + Cos[3*Pi/8])/Sqrt[2] = 0.9238795325112867 */
  /* 3135/4096 = (Sin[3*Pi/8] - Cos[3*Pi/8])*Sqrt[2] = 0.7653668647301795 */
  /* 8867/8192 = 2*Cos[3*Pi/8]*Sqrt[2]               = 1.0823922002923940 */
  od_rotate_add_avg(p0, p1, 473, 9, 3135, 12, 8867, 13, SHIFT);
}

/* --- 4-point Transforms --- */

/**
 * 4-point orthonormal Type-II fDCT
 */
static INLINE void od_fdct_4(od_coeff *q0, od_coeff *q1, od_coeff *q2,
 od_coeff *q3) {
  od_coeff q1h;
  od_coeff q3h;

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_neg(q0, q3, &q3h);
  od_butterfly_add(q1, &q1h, q2);

  /* Embedded 2-point transforms with asymmetric input. */
  od_fdct_2_asym(q0, q1, q1h);
  od_fdst_2_asym(q3, q3h, q2);
}

/**
 * 4-point orthonormal Type-II iDCT
 */
static INLINE void od_idct_4(od_coeff *q0, od_coeff *q2,
                             od_coeff *q1, od_coeff *q3)  {
  od_coeff q1h;

  /* Embedded 2-point transforms with asymmetric output. */
  od_idst_2_asym(q3, q2);
  od_idct_2_asym(q0, q1, &q1h);

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_add_asym(q1, q1h, q2);
  od_butterfly_neg_asym(q0, q3, od_rshift1(*q3));
}

/**
 * 4-point asymmetric Type-II fDCT
 */
static INLINE void od_fdct_4_asym(od_coeff *q0, od_coeff *q1, od_coeff q1h,
                                  od_coeff *q2, od_coeff *q3, od_coeff q3h) {

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_neg_asym(q0, q3, q3h);
  od_butterfly_sub_asym(q1, q1h, q2);

  /* Embedded 2-point orthonormal transforms. */
  od_fdct_2(q0, q1);
  od_fdst_2(q3, q2);
}

/**
 * 4-point asymmetric Type-II iDCT
 */
static INLINE void od_idct_4_asym(od_coeff *q0, od_coeff *q2,
                                  od_coeff *q1, od_coeff *q1h,
                                  od_coeff *q3, od_coeff *q3h)  {

  /* Embedded 2-point orthonormal transforms. */
  od_idst_2(q3, q2);
  od_idct_2(q0, q1);

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_sub(q1, q1h, q2);
  od_butterfly_neg(q0, q3, q3h);
}

/**
 * 4-point orthonormal Type-IV fDST
 */
static INLINE void od_fdst_4(od_coeff *q0, od_coeff *q1,
                             od_coeff *q2, od_coeff *q3) {

  /* Stage 0 */

  /* 13623/16384 = (Sin[7*Pi/16] + Cos[7*Pi/16])/Sqrt[2] = 0.831469612302545 */
  /*   4551/4096 = (Sin[7*Pi/16] - Cos[7*Pi/16])*Sqrt[2] = 1.111140466039204 */
  /*  9041/32768 = Cos[7*Pi/16]*Sqrt[2]                  = 0.275899379282943 */
  od_rotate_add(q0, q3, 13623, 14, 4551, 12, 565, 11, SHIFT);

  /* 16069/16384 = (Sin[5*Pi/16] + Cos[5*Pi/16])/Sqrt[2] = 0.9807852804032304 */
  /* 12785/32768 = (Sin[5*Pi/16] - Cos[5*Pi/16])*Sqrt[2] = 0.3901806440322566 */
  /*   1609/2048 = Cos[5*Pi/16]*Sqrt[2]                  = 0.7856949583871021 */
  od_rotate_sub(q2, q1, 16069, 14, 12785, 15, 1609, 11, SHIFT);

  /* Stage 1 */

  od_butterfly_sub_asym(q0, od_rshift1(*q0), q1);
  od_butterfly_sub_asym(q2, od_rshift1(*q2), q3);

  /* Stage 2 */

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(q2, q1, 5793, 12, 11585, 13);
}

/**
 * 4-point orthonormal Type-IV iDST
 */
static INLINE void od_idst_4(od_coeff *q0, od_coeff *q2,
                             od_coeff *q1, od_coeff *q3) {
  od_coeff q0h;
  od_coeff q2h;

  /* Stage 0 */

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(q2, q1, 5793, 12, 11585, 13);

  /* Stage 1 */

  od_butterfly_sub(q2, &q2h, q3);
  od_butterfly_sub(q0, &q0h, q1);

  /* Stage 2 */

  /*   8035/8192 = (Sin[5*Pi/16] + Cos[5*Pi/16])/Sqrt[2] = 0.9807852804032304 */
  /* 12785/32768 = (Sin[5*Pi/16] - Cos[5*Pi/16])*Sqrt[2] = 0.3901806440322566 */
  /* 12873/16384 = Cos[5*Pi/16]*Sqrt[2]                  = 0.7856949583871021 */
  od_rotate_sub_half(q2, q1, q2h, 8035, 13, 12785, 15, 12873, 14, NONE);

  /*   6811/8192 = (Sin[7*Pi/16] + Cos[7*Pi/16])/Sqrt[2] = 0.831469612302545 */
  /* 18205/16384 = (Sin[7*Pi/16] - Cos[7*Pi/16])*Sqrt[2] = 1.111140466039204 */
  /*  9041/32768 = Cos[7*Pi/16]*Sqrt[2]                  = 0.275899379282943 */
  od_rotate_add_half(q0, q3, q0h, 6811, 13, 18205, 14, 9041, 15, NONE);
}

/**
 * 4-point asymmetric Type-IV fDST
 */
static INLINE void od_fdst_4_asym(od_coeff *q0, od_coeff q0h, od_coeff *q1,
                                  od_coeff *q2, od_coeff q2h, od_coeff *q3) {

  /* Stage 0 */

  /*  9633/16384 = (Sin[7*Pi/16] + Cos[7*Pi/16])/2 = 0.5879378012096793 */
  /*  12873/8192 = (Sin[7*Pi/16] - Cos[7*Pi/16])*2 = 1.5713899167742045 */
  /* 12785/32768 = Cos[7*Pi/16]*2                  = 0.3901806440322565 */
  od_rotate_add_half(q0, q3, q0h, 9633, 14, 12873, 13, 12785, 15, SHIFT);

  /* 11363/16384 = (Sin[5*Pi/16] + Cos[5*Pi/16])/2 = 0.6935199226610738 */
  /* 18081/32768 = (Sin[5*Pi/16] - Cos[5*Pi/16])*2 = 0.5517987585658861 */
  /*  4551/4096 = Cos[5*Pi/16]*2                  = 1.1111404660392044 */
  od_rotate_sub_half(q2, q1, q2h, 11363, 14, 18081, 15, 4551, 12, SHIFT);

  /* Stage 1 */

  od_butterfly_sub_asym(q0, od_rshift1(*q0), q1);
  od_butterfly_sub_asym(q2, od_rshift1(*q2), q3);

  /* Stage 2 */

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(q2, q1, 11585, 13, 11585, 13);
}

/**
 * 4-point asymmetric Type-IV iDST
 */
static INLINE void od_idst_4_asym(od_coeff *q0, od_coeff *q2,
                                  od_coeff *q1, od_coeff *q3) {
  od_coeff q0h;
  od_coeff q2h;

  /* Stage 0 */

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(q2, q1, 11585, 13, 11585, 13);

  /* Stage 1 */

  od_butterfly_sub(q2, &q2h, q3);
  od_butterfly_sub(q0, &q0h, q1);

  /* Stage 2 */

  /* 11363/16384 = (Sin[5*Pi/16] + Cos[5*Pi/16])/2 = 0.6935199226610738 */
  /* 18081/32768 = (Sin[5*Pi/16] - Cos[5*Pi/16])*2 = 0.5517987585658861 */
  /* 18205/16384 = Cos[5*Pi/16]*2                  = 1.1111404660392044 */
  od_rotate_sub_half(q2, q1, q2h, 11363, 14, 18081, 15, 18205, 14, SHIFT);

  /*  9633/16384 = (Sin[7*Pi/16] + Cos[7*Pi/16])/2 = 0.5879378012096793 */
  /*  12873/8192 = (Sin[7*Pi/16] - Cos[7*Pi/16])*2 = 1.5713899167742045 */
  /* 12785/32768 = Cos[7*Pi/16]*2                  = 0.3901806440322565 */
  od_rotate_add_half(q0, q3, q0h, 9633, 14, 12873, 13, 12785, 15, SHIFT);
}

/* --- 8-point Transforms --- */

/**
 * 8-point orthonormal Type-II fDCT
 */
static INLINE void od_fdct_8(od_coeff *r0, od_coeff *r1,
                             od_coeff *r2, od_coeff *r3,
                             od_coeff *r4, od_coeff *r5,
                             od_coeff *r6, od_coeff *r7) {
  od_coeff r1h;
  od_coeff r3h;
  od_coeff r5h;
  od_coeff r7h;

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_neg(r0, r7, &r7h);
  od_butterfly_add(r1, &r1h, r6);
  od_butterfly_neg(r2, r5, &r5h);
  od_butterfly_add(r3, &r3h, r4);

  /* Embedded 4-point forward transforms with asymmetric input. */
  od_fdct_4_asym(r0, r1, r1h, r2, r3, r3h);
  od_fdst_4_asym(r7, r7h, r6, r5, r5h, r4);
}

/**
 * 8-point orthonormal Type-II iDCT
 */
static INLINE void od_idct_8(od_coeff *r0, od_coeff *r4,
                             od_coeff *r2, od_coeff *r6,
                             od_coeff *r1, od_coeff *r5,
                             od_coeff *r3, od_coeff *r7) {
  od_coeff r1h;
  od_coeff r3h;

  /* Embedded 4-point inverse transforms with asymmetric output. */
  od_idst_4_asym(r7, r5, r6, r4);
  od_idct_4_asym(r0, r2, r1, &r1h, r3, &r3h);

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_add_asym(r3, r3h, r4);
  od_butterfly_neg_asym(r2, r5, od_rshift1(*r5));
  od_butterfly_add_asym(r1, r1h, r6);
  od_butterfly_neg_asym(r0, r7, od_rshift1(*r7));
}

/**
 * 8-point asymmetric Type-II fDCT
 */
static INLINE void od_fdct_8_asym(od_coeff *r0, od_coeff *r1, od_coeff r1h,
                                  od_coeff *r2, od_coeff *r3, od_coeff r3h,
                                  od_coeff *r4, od_coeff *r5, od_coeff r5h,
                                  od_coeff *r6, od_coeff *r7, od_coeff r7h) {

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_neg_asym(r0, r7, r7h);
  od_butterfly_sub_asym(r1, r1h, r6);
  od_butterfly_neg_asym(r2, r5, r5h);
  od_butterfly_sub_asym(r3, r3h, r4);

  /* Embedded 4-point orthonormal transforms. */
  od_fdct_4(r0, r1, r2, r3);
  od_fdst_4(r7, r6, r5, r4);
}

/**
 * 8-point asymmetric Type-II iDCT
 */
static INLINE void od_idct_8_asym(od_coeff *r0, od_coeff *r4,
                                  od_coeff *r2, od_coeff *r6,
                                  od_coeff *r1, od_coeff *r1h,
                                  od_coeff *r5, od_coeff *r5h,
                                  od_coeff *r3, od_coeff *r3h,
                                  od_coeff *r7, od_coeff *r7h)  {

  /* Embedded 4-point inverse orthonormal transforms. */
  od_idst_4(r7, r5, r6, r4);
  od_idct_4(r0, r2, r1, r3);

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_sub(r3, r3h, r4);
  od_butterfly_neg(r2, r5, r5h);
  od_butterfly_sub(r1, r1h, r6);
  od_butterfly_neg(r0, r7, r7h);
}

/**
 * 8-point orthonormal Type-IV fDST
 */
static INLINE void od_fdst_8(od_coeff *r0, od_coeff *r1,
                             od_coeff *r2, od_coeff *r3,
                             od_coeff *r4, od_coeff *r5,
                             od_coeff *r6, od_coeff *r7) {
  od_coeff r0h;
  od_coeff r2h;
  od_coeff r5h;
  od_coeff r7h;

  /* Stage 0 */

  /* 17911/16384 = Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* 14699/16384 = Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /*    803/8192 = Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_add(r0, r7, 17911, 14, 14699, 14, 803, 13, NONE);

  /* 20435/16384 = Sin[13*Pi/32] + Cos[13*Pi/32] = 1.24722501298667123 */
  /* 21845/32768 = Sin[13*Pi/32] - Cos[13*Pi/32] = 0.66665565847774650 */
  /*   1189/4096 = Cos[13*Pi/32]                 = 0.29028467725446233 */
  od_rotate_sub(r6, r1, 20435, 14, 21845, 15, 1189, 12, NONE);

  /* 22173/16384 = Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /*   3363/8192 = Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* 15447/32768 = Cos[11*Pi/32]                 = 0.47139673682599764 */
  od_rotate_add(r2, r5, 22173, 14, 3363, 13, 15447, 15, NONE);

  /* 23059/16384 = Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /*  2271/16384 = Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /*   5197/8192 = Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_sub(r4, r3, 23059, 14, 2271, 14, 5197, 13, NONE);

  /* Stage 1 */

  od_butterfly_add(r0, &r0h, r3);
  od_butterfly_sub(r2, &r2h, r1);
  od_butterfly_add(r5, &r5h, r6);
  od_butterfly_sub(r7, &r7h, r4);

  /* Stage 2 */

  od_butterfly_add_asym(r7, r7h, r6);
  od_butterfly_add_asym(r5, r5h, r3);
  od_butterfly_add_asym(r2, r2h, r4);
  od_butterfly_sub_asym(r0, r0h, r1);

  /* Stage 3 */

  /* 10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_sub_avg(r3, r4, 10703, 13, 8867, 14, 3135, 12, NONE);

  /* 10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_neg_avg(r2, r5, 10703, 13, 8867, 14, 3135, 12);

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_sub_avg(r1, r6, 11585, 13, 11585, 13);
}

/**
 * 8-point orthonormal Type-IV iDST
 */
static INLINE void od_idst_8(od_coeff *r0, od_coeff *r4,
                             od_coeff *r2, od_coeff *r6,
                             od_coeff *r1, od_coeff *r5,
                             od_coeff *r3, od_coeff *r7) {
  od_coeff r0h;
  od_coeff r2h;
  od_coeff r5h;
  od_coeff r7h;

  /* Stage 3 */

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(r6, r1, 11585, 13, 11585, 13);

  /* 10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_neg_avg(r5, r2, 10703, 13, 8867, 14, 3135, 12);

  /* 10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_add_avg(r4, r3, 10703, 13, 8867, 14, 3135, 12, NONE);

  /* Stage 2 */

  od_butterfly_sub(r0, &r0h, r1);
  od_butterfly_add(r2, &r2h, r4);
  od_butterfly_add(r5, &r5h, r3);
  od_butterfly_add(r7, &r7h, r6);

  /* Stage 1 */

  od_butterfly_sub_asym(r7, r7h, r4);
  od_butterfly_add_asym(r5, r5h, r6);
  od_butterfly_sub_asym(r2, r2h, r1);
  od_butterfly_add_asym(r0, r0h, r3);

  /* Stage 0 */

  /* 23059/16384 = Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /*  2271/16384 = Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /*   5197/8192 = Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_sub(r4, r3, 23059, 14, 2271, 14, 5197, 13, NONE);

  /* 22173/16384 = Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /*   3363/8192 = Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* 15447/32768 = Cos[11*Pi/32]                 = 0.47139673682599764 */
  od_rotate_add(r2, r5, 22173, 14, 3363, 13, 15447, 15, NONE);

  /* 20435/16384 = Sin[13*Pi/32] + Cos[13*Pi/32] = 1.24722501298667123 */
  /* 21845/32768 = Sin[13*Pi/32] - Cos[13*Pi/32] = 0.66665565847774650 */
  /*   1189/4096 = Cos[13*Pi/32]                 = 0.29028467725446233 */
  od_rotate_sub(r6, r1, 20435, 14, 21845, 15, 1189, 12, NONE);

  /* 17911/16384 = Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* 14699/16384 = Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /*    803/8192 = Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_add(r0, r7, 17911, 14, 14699, 14, 803, 13, NONE);
}

/**
 * 8-point asymmetric Type-IV fDST
 */
static INLINE void od_fdst_8_asym(od_coeff *r0, od_coeff r0h, od_coeff *r1,
                                  od_coeff *r2, od_coeff r2h, od_coeff *r3,
                                  od_coeff *r4, od_coeff r4h, od_coeff *r5,
                                  od_coeff *r6, od_coeff r6h, od_coeff *r7) {
  od_coeff r5h;
  od_coeff r7h;

  /* Stage 0 */

  /* 12665/16384 = (Sin[15*Pi/32] + Cos[15*Pi/32])/Sqrt[2] = 0.77301045336274 */
  /*   5197/4096 = (Sin[15*Pi/32] - Cos[15*Pi/32])*Sqrt[2] = 1.26878656832729 */
  /*  2271/16384 = Cos[15*Pi/32]*Sqrt[2]                   = 0.13861716919909 */
  od_rotate_add_half(r0, r7, r0h, 12665, 14, 5197, 12, 2271, 14, NONE);

  /* 14449/16384 = Sin[13*Pi/32] + Cos[13*Pi/32])/Sqrt[2] = 0.881921264348355 */
  /* 30893/32768 = Sin[13*Pi/32] - Cos[13*Pi/32])*Sqrt[2] = 0.942793473651995 */
  /*   3363/8192 = Cos[13*Pi/32]*Sqrt[2]                  = 0.410524527522357 */
  od_rotate_sub_half(r6, r1, r6h, 14449, 14, 30893, 15, 3363, 13, NONE);

  /* 15679/16384 = Sin[11*Pi/32] + Cos[11*Pi/32])/Sqrt[2] = 0.956940335732209 */
  /*   1189/2048 = Sin[11*Pi/32] - Cos[11*Pi/32])*Sqrt[2] = 0.580569354508925 */
  /*   5461/8192 = Cos[11*Pi/32]*Sqrt[2]                  = 0.666655658477747 */
  od_rotate_add_half(r2, r5, r2h, 15679, 14, 1189, 11, 5461, 13, NONE);

  /* 16305/16384 = (Sin[9*Pi/32] + Cos[9*Pi/32])/Sqrt[2] = 0.9951847266721969 */
  /*    803/4096 = (Sin[9*Pi/32] - Cos[9*Pi/32])*Sqrt[2] = 0.1960342806591213 */
  /* 14699/16384 = Cos[9*Pi/32]*Sqrt[2]                  = 0.8971675863426364 */
  od_rotate_sub_half(r4, r3, r4h, 16305, 14, 803, 12, 14699, 14, NONE);

  /* Stage 1 */

  od_butterfly_add(r0, &r0h, r3);
  od_butterfly_sub(r2, &r2h, r1);
  od_butterfly_add(r5, &r5h, r6);
  od_butterfly_sub(r7, &r7h, r4);

  /* Stage 2 */

  od_butterfly_add_asym(r7, r7h, r6);
  od_butterfly_add_asym(r5, r5h, r3);
  od_butterfly_add_asym(r2, r2h, r4);
  od_butterfly_sub_asym(r0, r0h, r1);

  /* Stage 3 */

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_sub_avg(r3, r4, 669, 9, 8867, 14, 3135, 12, NONE);

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_neg_avg(r2, r5, 669, 9, 8867, 14, 3135, 12);

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_sub_avg(r1, r6, 5793, 12, 11585, 13);
}

/**
 * 8-point asymmetric Type-IV iDST
 */
static INLINE void od_idst_8_asym(od_coeff *r0, od_coeff *r4,
                                  od_coeff *r2, od_coeff *r6,
                                  od_coeff *r1, od_coeff *r5,
                                  od_coeff *r3, od_coeff *r7) {
  od_coeff r0h;
  od_coeff r2h;
  od_coeff r5h;
  od_coeff r7h;

  /* Stage 3 */

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(r6, r1, 5793, 12, 11585, 13);

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_neg_avg(r5, r2, 669, 9, 8867, 14, 3135, 12);

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_add_avg(r4, r3, 669, 9, 8867, 14, 3135, 12, NONE);

  /* Stage 2 */

  od_butterfly_sub(r0, &r0h, r1);
  od_butterfly_add(r2, &r2h, r4);
  od_butterfly_add(r5, &r5h, r3);
  od_butterfly_add(r7, &r7h, r6);

  /* Stage 1 */

  od_butterfly_sub_asym(r7, r7h, r4);
  od_butterfly_add_asym(r5, r5h, r6);
  od_butterfly_sub_asym(r2, r2h, r1);
  od_butterfly_add_asym(r0, r0h, r3);

  /* Stage 0 */

  /* 16305/16384 = (Sin[9*Pi/32] + Cos[9*Pi/32])/Sqrt[2] = 0.9951847266721969 */
  /*    803/4096 = (Sin[9*Pi/32] - Cos[9*Pi/32])*Sqrt[2] = 0.1960342806591213 */
  /* 14699/16384 = Cos[9*Pi/32]*Sqrt[2]                  = 0.8971675863426364 */
  od_rotate_sub(r4, r3, 16305, 14, 803, 12, 14699, 14, SHIFT);

  /* 15679/16384 = Sin[11*Pi/32] + Cos[11*Pi/32])/Sqrt[2] = 0.956940335732209 */
  /*   1189/2048 = Sin[11*Pi/32] - Cos[11*Pi/32])*Sqrt[2] = 0.580569354508925 */
  /*   5461/8192 = Cos[11*Pi/32]*Sqrt[2]                  = 0.666655658477747 */
  od_rotate_add(r2, r5, 15679, 14, 1189, 11, 5461, 13, SHIFT);

  /* 14449/16384 = Sin[13*Pi/32] + Cos[13*Pi/32])/Sqrt[2] = 0.881921264348355 */
  /* 15447/16384 = Sin[13*Pi/32] - Cos[13*Pi/32])*Sqrt[2] = 0.942793473651995 */
  /*   3363/8192 = Cos[13*Pi/32]*Sqrt[2]                  = 0.410524527522357 */
  od_rotate_sub(r6, r1, 14449, 14, 15447, 14, 3363, 13, SHIFT);

  /* 12665/16384 = (Sin[15*Pi/32] + Cos[15*Pi/32])/Sqrt[2] = 0.77301045336274 */
  /*   5197/4096 = (Sin[15*Pi/32] - Cos[15*Pi/32])*Sqrt[2] = 1.26878656832729 */
  /*  2271/16384 = Cos[15*Pi/32]*Sqrt[2]                   = 0.13861716919909 */
  od_rotate_add(r0, r7, 12665, 14, 5197, 12, 2271, 14, SHIFT);
}

/* --- 16-point Transforms --- */

/**
 * 16-point orthonormal Type-II fDCT
 */
static INLINE void od_fdct_16(od_coeff *s0, od_coeff *s1,
                              od_coeff *s2, od_coeff *s3,
                              od_coeff *s4, od_coeff *s5,
                              od_coeff *s6, od_coeff *s7,
                              od_coeff *s8, od_coeff *s9,
                              od_coeff *sa, od_coeff *sb,
                              od_coeff *sc, od_coeff *sd,
                              od_coeff *se, od_coeff *sf) {
  od_coeff s1h;
  od_coeff s3h;
  od_coeff s5h;
  od_coeff s7h;
  od_coeff s9h;
  od_coeff sbh;
  od_coeff sdh;
  od_coeff sfh;

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_neg(s0, sf, &sfh);
  od_butterfly_add(s1, &s1h, se);
  od_butterfly_neg(s2, sd, &sdh);
  od_butterfly_add(s3, &s3h, sc);
  od_butterfly_neg(s4, sb, &sbh);
  od_butterfly_add(s5, &s5h, sa);
  od_butterfly_neg(s6, s9, &s9h);
  od_butterfly_add(s7, &s7h, s8);

  /* Embedded 8-point transforms with asymmetric input. */
  od_fdct_8_asym(s0, s1, s1h, s2, s3, s3h, s4, s5, s5h, s6, s7, s7h);
  od_fdst_8_asym(sf, sfh, se, sd, sdh, sc, sb, sbh, sa, s9, s9h, s8);
}

/**
 * 16-point orthonormal Type-II iDCT
 */
static INLINE void od_idct_16(od_coeff *s0, od_coeff *s8,
                              od_coeff *s4, od_coeff *sc,
                              od_coeff *s2, od_coeff *sa,
                              od_coeff *s6, od_coeff *se,
                              od_coeff *s1, od_coeff *s9,
                              od_coeff *s5, od_coeff *sd,
                              od_coeff *s3, od_coeff *sb,
                              od_coeff *s7, od_coeff *sf) {
  od_coeff s1h;
  od_coeff s3h;
  od_coeff s5h;
  od_coeff s7h;

  /* Embedded 8-point transforms with asymmetric output. */
  od_idst_8_asym(sf, sb, sd, s9, se, sa, sc, s8);
  od_idct_8_asym(s0, s4, s2, s6, s1, &s1h, s5, &s5h, s3, &s3h, s7, &s7h);

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_add_asym(s7, s7h, s8);
  od_butterfly_neg_asym(s6, s9, od_rshift1(*s9));
  od_butterfly_add_asym(s5, s5h, sa);
  od_butterfly_neg_asym(s4, sb, od_rshift1(*sb));
  od_butterfly_add_asym(s3, s3h, sc);
  od_butterfly_neg_asym(s2, sd, od_rshift1(*sd));
  od_butterfly_add_asym(s1, s1h, se);
  od_butterfly_neg_asym(s0, sf, od_rshift1(*sf));
}

/**
 * 16-point asymmetric Type-II fDCT
 */
static INLINE void od_fdct_16_asym(od_coeff *s0, od_coeff *s1, od_coeff s1h,
                                   od_coeff *s2, od_coeff *s3, od_coeff s3h,
                                   od_coeff *s4, od_coeff *s5, od_coeff s5h,
                                   od_coeff *s6, od_coeff *s7, od_coeff s7h,
                                   od_coeff *s8, od_coeff *s9, od_coeff s9h,
                                   od_coeff *sa, od_coeff *sb, od_coeff sbh,
                                   od_coeff *sc, od_coeff *sd, od_coeff sdh,
                                   od_coeff *se, od_coeff *sf, od_coeff sfh) {

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_neg_asym(s0, sf, sfh);
  od_butterfly_sub_asym(s1, s1h, se);
  od_butterfly_neg_asym(s2, sd, sdh);
  od_butterfly_sub_asym(s3, s3h, sc);
  od_butterfly_neg_asym(s4, sb, sbh);
  od_butterfly_sub_asym(s5, s5h, sa);
  od_butterfly_neg_asym(s6, s9, s9h);
  od_butterfly_sub_asym(s7, s7h, s8);

  /* Embedded 8-point orthonormal transforms. */
  od_fdct_8(s0, s1, s2, s3, s4, s5, s6, s7);
  od_fdst_8(sf, se, sd, sc, sb, sa, s9, s8);
}

/**
 * 16-point asymmetric Type-II iDCT
 */
static INLINE void od_idct_16_asym(od_coeff *s0, od_coeff *s8,
                                   od_coeff *s4, od_coeff *sc,
                                   od_coeff *s2, od_coeff *sa,
                                   od_coeff *s6, od_coeff *se,
                                   od_coeff *s1, od_coeff *s1h,
                                   od_coeff *s9, od_coeff *s9h,
                                   od_coeff *s5, od_coeff *s5h,
                                   od_coeff *sd, od_coeff *sdh,
                                   od_coeff *s3, od_coeff *s3h,
                                   od_coeff *sb, od_coeff *sbh,
                                   od_coeff *s7, od_coeff *s7h,
                                   od_coeff *sf, od_coeff *sfh) {

  /* Embedded 8-point orthonormal transforms. */
  od_idst_8(sf, sb, sd, s9, se, sa, sc, s8);
  od_idct_8(s0, s4, s2, s6, s1, s5, s3, s7);

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_sub(s7, s7h, s8);
  od_butterfly_neg(s6, s9, s9h);
  od_butterfly_sub(s5, s5h, sa);
  od_butterfly_neg(s4, sb, sbh);
  od_butterfly_sub(s3, s3h, sc);
  od_butterfly_neg(s2, sd, sdh);
  od_butterfly_sub(s1, s1h, se);
  od_butterfly_neg(s0, sf, sfh);
}

/**
 * 16-point orthonormal Type-IV fDST
 */
static INLINE void od_fdst_16(od_coeff *s0, od_coeff *s1,
                              od_coeff *s2, od_coeff *s3,
                              od_coeff *s4, od_coeff *s5,
                              od_coeff *s6, od_coeff *s7,
                              od_coeff *s8, od_coeff *s9,
                              od_coeff *sa, od_coeff *sb,
                              od_coeff *sc, od_coeff *sd,
                              od_coeff *se, od_coeff *sf) {
  od_coeff s0h;
  od_coeff s2h;
  od_coeff sdh;
  od_coeff sfh;

  /* Stage 0 */

  /* 24279/32768 = (Sin[31*Pi/64] + Cos[31*Pi/64])/Sqrt[2] = 0.74095112535496 */
  /*  11003/8192 = (Sin[31*Pi/64] - Cos[31*Pi/64])*Sqrt[2] = 1.34311790969404 */
  /*  1137/16384 = Cos[31*Pi/64]*Sqrt[2]                   = 0.06939217050794 */
  od_rotate_add(s0, sf, 24279, 15, 11003, 13, 1137, 14, SHIFT);

  /* 1645/2048 = (Sin[29*Pi/64] + Cos[29*Pi/64])/Sqrt[2] = 0.8032075314806449 */
  /*   305/256 = (Sin[29*Pi/64] - Cos[29*Pi/64])*Sqrt[2] = 1.1913986089848667 */
  /*  425/2048 = Cos[29*Pi/64]*Sqrt[2]                   = 0.2075082269882116 */
  od_rotate_sub(se, s1, 1645, 11, 305, 8, 425, 11, SHIFT);

  /* 14053/32768 = (Sin[27*Pi/64] + Cos[27*Pi/64])/Sqrt[2] = 0.85772861000027 */
  /*   8423/8192 = (Sin[27*Pi/64] - Cos[27*Pi/64])*Sqrt[2] = 1.02820548838644 */
  /*   2815/8192 = Cos[27*Pi/64]*Sqrt[2]                   = 0.34362586580705 */
  od_rotate_add(s2, sd, 14053, 14, 8423, 13, 2815, 13, SHIFT);

  /* 14811/16384 = (Sin[25*Pi/64] + Cos[25*Pi/64])/Sqrt[2] = 0.90398929312344 */
  /*   7005/8192 = (Sin[25*Pi/64] - Cos[25*Pi/64])*Sqrt[2] = 0.85511018686056 */
  /*   3903/8192 = Cos[25*Pi/64]*Sqrt[2]                   = 0.47643419969316 */
  od_rotate_sub(sc, s3, 14811, 14, 7005, 13, 3903, 13, SHIFT);

  /* 30853/32768 = (Sin[23*Pi/64] + Cos[23*Pi/64])/Sqrt[2] = 0.94154406518302 */
  /* 11039/16384 = (Sin[23*Pi/64] - Cos[23*Pi/64])*Sqrt[2] = 0.67377970678444 */
  /*  9907/16384 = Cos[23*Pi/64]*Sqrt[2]                   = 0.60465421179080 */
  od_rotate_add(s4, sb, 30853, 15, 11039, 14, 9907, 14, SHIFT);

  /* 15893/16384 = (Sin[21*Pi/64] + Cos[21*Pi/64])/Sqrt[2] = 0.97003125319454 */
  /*   3981/8192 = (Sin[21*Pi/64] - Cos[21*Pi/64])*Sqrt[2] = 0.89716758634264 */
  /*   1489/2048 = Cos[21*Pi/64]*Sqrt[2]                   = 0.72705107329128 */
  od_rotate_sub(sa, s5, 15893, 14, 3981, 13, 1489, 11, SHIFT);

  /* 32413/32768 = (Sin[19*Pi/64] + Cos[19*Pi/64])/Sqrt[2] = 0.98917650996478 */
  /*    601/2048 = (Sin[19*Pi/64] - Cos[19*Pi/64])*Sqrt[2] = 0.29346094891072 */
  /* 13803/16384 = Cos[19*Pi/64]*Sqrt[2]                   = 0.84244603550942 */
  od_rotate_add(s6, s9, 32413, 15, 601, 11, 13803, 14, SHIFT);

  /* 32729/32768 = (Sin[17*Pi/64] + Cos[17*Pi/64])/Sqrt[2] = 0.99879545620517 */
  /*    201/2048 = (Sin[17*Pi/64] - Cos[17*Pi/64])*Sqrt[2] = 0.09813534865484 */
  /*   1945/2048 = Cos[17*Pi/64]*Sqrt[2]                   = 0.94972778187775 */
  od_rotate_sub(s8, s7, 32729, 15, 201, 11, 1945, 11, SHIFT);

  /* Stage 1 */

  od_butterfly_sub_asym(s0, od_rshift1(*s0), s7);
  od_butterfly_sub_asym(s8, od_rshift1(*s8), sf);
  od_butterfly_add_asym(s4, od_rshift1(*s4), s3);
  od_butterfly_add_asym(sc, od_rshift1(*sc), sb);
  od_butterfly_sub_asym(s2, od_rshift1(*s2), s5);
  od_butterfly_sub_asym(sa, od_rshift1(*sa), sd);
  od_butterfly_add_asym(s6, od_rshift1(*s6), s1);
  od_butterfly_add_asym(se, od_rshift1(*se), s9);

  /* Stage 2 */

  od_butterfly_add(s8, NULL, s4);
  od_butterfly_add(s7, NULL, sb);
  od_butterfly_sub(sa, NULL, s6);
  od_butterfly_sub(s5, NULL, s9);
  od_butterfly_add(s0, &s0h, s3);
  od_butterfly_add(sd, &sdh, se);
  od_butterfly_sub(s2, &s2h, s1);
  od_butterfly_sub(sf, &sfh, sc);

  /* Stage 3 */

  /*     301/256 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /*   1609/2048 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* 12785/32768 = 2*Cos[7*Pi/16]              = 0.3901806440322565 */
  od_rotate_add_avg(s8, s7, 301, 8, 1609, 11, 12785, 15, NONE);

  /* 11363/8192 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* 9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*  4551/8192 = Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(s9, s6, 11363, 13, 9041, 15, 4551, 13, NONE);

  /*  5681/4096 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* 9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*  4551/4096 = 2*Cos[5*Pi/16]              = 1.1111404660392044 */
  od_rotate_neg_avg(s5, sa, 5681, 12, 9041, 15, 4551, 12);

  /*   9633/8192 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* 12873/16384 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /*  6393/32768 = Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(s4, sb, 9633, 13, 12873, 14, 6393, 15);

  /* Stage 4 */

  od_butterfly_add_asym(s2, s2h, sc);
  od_butterfly_sub_asym(s0, s0h, s1);
  od_butterfly_add_asym(sf, sfh, se);
  od_butterfly_add_asym(sd, sdh, s3);
  od_butterfly_add_asym(s7, od_rshift1(*s7), s6);
  od_butterfly_sub_asym(s8, od_rshift1(*s8), s9);
  od_butterfly_sub_asym(sa, od_rshift1(*sa), sb);
  od_butterfly_add_asym(s5, od_rshift1(*s5), s4);

  /* Stage 5 */

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[7*Pi/8]             = 0.7653668647301796 */
  od_rotate_add_avg(sc, s3, 669, 9, 8867, 14, 3135, 12, NONE);

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3870398453221475 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_neg_avg(s2, sd, 669, 9, 8867, 14, 3135, 12);

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(sa, s5, 5793, 12, 11585, 13);

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(s6, s9, 5793, 12, 11585, 13);

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(se, s1, 5793, 12, 11585, 13);
}

/**
 * 16-point orthonormal Type-IV iDST
 */
static INLINE void od_idst_16(od_coeff *s0, od_coeff *s8,
                              od_coeff *s4, od_coeff *sc,
                              od_coeff *s2, od_coeff *sa,
                              od_coeff *s6, od_coeff *se,
                              od_coeff *s1, od_coeff *s9,
                              od_coeff *s5, od_coeff *sd,
                              od_coeff *s3, od_coeff *sb,
                              od_coeff *s7, od_coeff *sf) {
  od_coeff s0h;
  od_coeff s2h;
  od_coeff s4h;
  od_coeff s6h;
  od_coeff s8h;
  od_coeff sah;
  od_coeff sch;
  od_coeff sdh;
  od_coeff seh;
  od_coeff sfh;

  /* Stage 5 */

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(s6, s9, 5793, 12, 11585, 13);

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(sa, s5, 5793, 12, 11585, 13);

  /*  5793/4096 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* 11585/8192 = 2*Cos[Pi/4]           = 1.4142135623730951 */
  od_rotate_pi4_add_avg(se, s1, 5793, 12, 11585, 13);

  /*     669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /*  8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*   3135/4096 = 2*Cos[7*Pi/8]             = 0.7653668647301796 */
  od_rotate_add_avg(sc, s3, 669, 9, 8867, 14, 3135, 12, NONE);

  /*    669/512 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3870398453221475 */
  /* 8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*  3135/4096 = 2*Cos[3*Pi/8]             = 0.7653668647301796 */
  od_rotate_neg_avg(sd, s2, 669, 9, 8867, 14, 3135, 12);

  /* Stage 4 */

  od_butterfly_add(s5, NULL, s4);
  od_butterfly_sub(sa, NULL, sb);
  od_butterfly_sub(s8, NULL, s9);
  od_butterfly_add(s7, NULL, s6);
  od_butterfly_add(sd, &sdh, s3);
  od_butterfly_add(sf, &sfh, se);
  od_butterfly_sub(s0, &s0h, s1);
  od_butterfly_add(s2, &s2h, sc);

  /* Stage 3 */

  /*     301/256 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /*   1609/2048 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* 12785/32768 = 2*Cos[7*Pi/16]              = 0.3901806440322565 */
  od_rotate_add_avg(s8, s7, 301, 8, 1609, 11, 12785, 15, NONE);

  /* 11363/8192 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* 9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*  4551/8192 = Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(s9, s6, 11363, 13, 9041, 15, 4551, 13, NONE);

  /*  5681/4096 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* 9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*  4551/4096 = 2*Cos[5*Pi/16]              = 1.1111404660392044 */
  od_rotate_neg_avg(sa, s5, 5681, 12, 9041, 15, 4551, 12);

  /*   9633/8192 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* 12873/16384 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /*  6393/32768 = Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(sb, s4, 9633, 13, 12873, 14, 6393, 15);

  /* Stage 2 */

  od_butterfly_add_asym(s8, od_rshift1(*s8), s4);
  od_butterfly_add_asym(s7, od_rshift1(*s7), sb);
  od_butterfly_sub_asym(sa, od_rshift1(*sa), s6);
  od_butterfly_sub_asym(s5, od_rshift1(*s5), s9);
  od_butterfly_add_asym(s0, s0h, s3);
  od_butterfly_add_asym(sd, sdh, se);
  od_butterfly_sub_asym(s2, s2h, s1);
  od_butterfly_sub_asym(sf, sfh, sc);

  /* Stage 1 */

  od_butterfly_sub(s0, &s0h, s7);
  od_butterfly_sub(s8, &s8h, sf);
  od_butterfly_add(s4, &s4h, s3);
  od_butterfly_add(sc, &sch, sb);
  od_butterfly_sub(s2, &s2h, s5);
  od_butterfly_sub(sa, &sah, sd);
  od_butterfly_add(s6, &s6h, s1);
  od_butterfly_add(se, &seh, s9);

  /* Stage 0 */

  /*   4091/4096 = (Sin[17*Pi/64] + Cos[17*Pi/64])/Sqrt[2] = 0.99879545620517 */
  /*    201/2048 = (Sin[17*Pi/64] - Cos[17*Pi/64])*Sqrt[2] = 0.09813534865484 */
  /* 31121/32768 = Cos[17*Pi/64]*Sqrt[2]                   = 0.94972778187775 */
  od_rotate_sub_half(s8, s7, s8h, 4091, 12, 201, 11, 31121, 15, NONE);

  /* 16207/16384 = (Sin[19*Pi/64] + Cos[19*Pi/64])/Sqrt[2] = 0.98917650996478 */
  /*    601/2048 = (Sin[19*Pi/64] - Cos[19*Pi/64])*Sqrt[2] = 0.29346094891072 */
  /* 27605/32768 = Cos[19*Pi/64]*Sqrt[2]                   = 0.84244603550942 */
  od_rotate_add_half(s6, s9, s6h, 16207, 14, 601, 11, 27605, 15, NONE);

  /* 15893/16384 = (Sin[21*Pi/64] + Cos[21*Pi/64])/Sqrt[2] = 0.97003125319454 */
  /*   3981/8192 = (Sin[21*Pi/64] - Cos[21*Pi/64])*Sqrt[2] = 0.89716758634264 */
  /*   1489/2048 = Cos[21*Pi/64]*Sqrt[2]                   = 0.72705107329128 */
  od_rotate_sub_half(sa, s5, sah, 15893, 14, 3981, 13, 1489, 11, NONE);

  /*   7713/8192 = (Sin[23*Pi/64] + Cos[23*Pi/64])/Sqrt[2] = 0.94154406518302 */
  /* 11039/16384 = (Sin[23*Pi/64] - Cos[23*Pi/64])*Sqrt[2] = 0.67377970678444 */
  /* 19813/32768 = Cos[23*Pi/64]*Sqrt[2]                   = 0.60465421179080 */
  od_rotate_add_half(s4, sb, s4h, 7713, 13, 11039, 14, 19813, 15, NONE);

  /* 14811/16384 = (Sin[25*Pi/64] + Cos[25*Pi/64])/Sqrt[2] = 0.90398929312344 */
  /*   7005/8192 = (Sin[25*Pi/64] - Cos[25*Pi/64])*Sqrt[2] = 0.85511018686056 */
  /*   3903/8192 = Cos[25*Pi/64]*Sqrt[2]                   = 0.47643419969316 */
  od_rotate_sub_half(sc, s3, sch, 14811, 14, 7005, 13, 3903, 13, NONE);

  /* 14053/32768 = (Sin[27*Pi/64] + Cos[27*Pi/64])/Sqrt[2] = 0.85772861000027 */
  /*   8423/8192 = (Sin[27*Pi/64] - Cos[27*Pi/64])*Sqrt[2] = 1.02820548838644 */
  /*   2815/8192 = Cos[27*Pi/64]*Sqrt[2]                   = 0.34362586580705 */
  od_rotate_add_half(s2, sd, s2h, 14053, 14, 8423, 13, 2815, 13, NONE);

  /* 1645/2048 = (Sin[29*Pi/64] + Cos[29*Pi/64])/Sqrt[2] = 0.8032075314806449 */
  /*   305/256 = (Sin[29*Pi/64] - Cos[29*Pi/64])*Sqrt[2] = 1.1913986089848667 */
  /*  425/2048 = Cos[29*Pi/64]*Sqrt[2]                   = 0.2075082269882116 */
  od_rotate_sub_half(se, s1, seh, 1645, 11, 305, 8, 425, 11, NONE);

  /*  3035/32768 = (Sin[31*Pi/64] + Cos[31*Pi/64])/Sqrt[2] = 0.74095112535496 */
  /* 44011/32768 = (Sin[31*Pi/64] - Cos[31*Pi/64])*Sqrt[2] = 1.34311790969404 */
  /*  1137/16384 = Cos[31*Pi/64]*Sqrt[2]                   = 0.06939217050794 */
  od_rotate_add_half(s0, sf, s0h, 3035, 12, 44011, 15, 1137, 14, NONE);
}

/**
 * 16-point asymmetric Type-IV fDST
 */
static INLINE void od_fdst_16_asym(od_coeff *s0, od_coeff s0h, od_coeff *s1,
                                   od_coeff *s2, od_coeff s2h, od_coeff *s3,
                                   od_coeff *s4, od_coeff s4h, od_coeff *s5,
                                   od_coeff *s6, od_coeff s6h, od_coeff *s7,
                                   od_coeff *s8, od_coeff s8h, od_coeff *s9,
                                   od_coeff *sa, od_coeff sah, od_coeff *sb,
                                   od_coeff *sc, od_coeff sch, od_coeff *sd,
                                   od_coeff *se, od_coeff seh, od_coeff *sf) {
  od_coeff sdh;
  od_coeff sfh;

  /* Stage 0 */

  /*   1073/2048 = (Sin[31*Pi/64] + Cos[31*Pi/64])/2 = 0.5239315652662953 */
  /* 62241/32768 = (Sin[31*Pi/64] - Cos[31*Pi/64])*2 = 1.8994555637555088 */
  /*   201/16384 = Cos[31*Pi/64]*2                   = 0.0981353486548360 */
  od_rotate_add_half(s0, sf, s0h, 1073, 11, 62241, 15, 201, 11, SHIFT);

  /* 18611/32768 = (Sin[29*Pi/64] + Cos[29*Pi/64])/2 = 0.5679534922100714 */
  /* 55211/32768 = (Sin[29*Pi/64] - Cos[29*Pi/64])*2 = 1.6848920710188384 */
  /*    601/2048 = Cos[29*Pi/64]*2                   = 0.2934609489107235 */
  od_rotate_sub_half(se, s1, seh, 18611, 15, 55211, 15, 601, 11, SHIFT);

  /*  9937/16384 = (Sin[27*Pi/64] + Cos[27*Pi/64])/2 = 0.6065057165489039 */
  /*   1489/1024 = (Sin[27*Pi/64] - Cos[27*Pi/64])*2 = 1.4541021465825602 */
  /*   3981/8192 = Cos[27*Pi/64]*2                   = 0.4859603598065277 */
  od_rotate_add_half(s2, sd, s2h, 9937, 14, 1489, 10, 3981, 13, SHIFT);

  /* 10473/16384 = (Sin[25*Pi/64] + Cos[25*Pi/64])/2 = 0.6392169592876205 */
  /* 39627/32768 = (Sin[25*Pi/64] - Cos[25*Pi/64])*2 = 1.2093084235816014 */
  /* 11039/16384 = Cos[25*Pi/64]*2                   = 0.6737797067844401 */
  od_rotate_sub_half(sc, s3, sch, 10473, 14, 39627, 15, 11039, 14, SHIFT);

  /* 2727/4096 = (Sin[23*Pi/64] + Cos[23*Pi/64])/2 = 0.6657721932768628 */
  /* 3903/4096 = (Sin[23*Pi/64] - Cos[23*Pi/64])*2 = 0.9528683993863225 */
  /* 7005/8192 = Cos[23*Pi/64]*2                   = 0.8551101868605642 */
  od_rotate_add_half(s4, sb, s4h, 2727, 12, 3903, 12, 7005, 13, SHIFT);

  /* 5619/8192 = (Sin[21*Pi/64] + Cos[21*Pi/64])/2 = 0.6859156770967569 */
  /* 2815/4096 = (Sin[21*Pi/64] - Cos[21*Pi/64])*2 = 0.6872517316141069 */
  /* 8423/8192 = Cos[21*Pi/64]*2                   = 1.0282054883864433 */
  od_rotate_sub_half(sa, s5, sah, 5619, 13, 2815, 12, 8423, 13, SHIFT);

  /*   2865/4096 = (Sin[19*Pi/64] + Cos[19*Pi/64])/2 = 0.6994534179865391 */
  /* 13588/32768 = (Sin[19*Pi/64] - Cos[19*Pi/64])*2 = 0.4150164539764232 */
  /*     305/256 = Cos[19*Pi/64]*2                   = 1.1913986089848667 */
  od_rotate_add_half(s6, s9, s6h, 2865, 12, 13599, 15, 305, 8, SHIFT);

  /* 23143/32768 = (Sin[17*Pi/64] + Cos[17*Pi/64])/2 = 0.7062550401009887 */
  /*   1137/8192 = (Sin[17*Pi/64] - Cos[17*Pi/64])*2 = 0.1387843410158816 */
  /*  11003/8192 = Cos[17*Pi/64]*2                   = 1.3431179096940367 */
  od_rotate_sub_half(s8, s7, s8h, 23143, 15, 1137, 13, 11003, 13, SHIFT);

  /* Stage 1 */

  od_butterfly_sub_asym(s0, od_rshift1(*s0), s7);
  od_butterfly_sub_asym(s8, od_rshift1(*s8), sf);
  od_butterfly_add_asym(s4, od_rshift1(*s4), s3);
  od_butterfly_add_asym(sc, od_rshift1(*sc), sb);
  od_butterfly_sub_asym(s2, od_rshift1(*s2), s5);
  od_butterfly_sub_asym(sa, od_rshift1(*sa), sd);
  od_butterfly_add_asym(s6, od_rshift1(*s6), s1);
  od_butterfly_add_asym(se, od_rshift1(*se), s9);

  /* Stage 2 */

  od_butterfly_add(s8, NULL, s4);
  od_butterfly_add(s7, NULL, sb);
  od_butterfly_sub(sa, NULL, s6);
  od_butterfly_sub(s5, NULL, s9);
  od_butterfly_add(s0, &s0h, s3);
  od_butterfly_add(sd, &sdh, se);
  od_butterfly_sub(s2, &s2h, s1);
  od_butterfly_sub(sf, &sfh, sc);

  /* Stage 3 */

  /*   9633/8192 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* 12873/16384 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /*  6393/32768 = Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_add(s8, s7, 9633, 13, 12873, 14, 6393, 15, NONE);

  /* 45451/32768 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /*  9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*   4551/8192 = Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(s9, s6, 22725, 14, 9041, 15, 4551, 13, NONE);

  /*  11363/8192 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /*  9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*   4551/8192 = Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_neg(s5, sa, 11363, 13, 9041, 15, 4551, 13);

  /*  9633/32768 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* 12873/16384 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /*  6393/32768 = Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(s4, sb, 9633, 13, 12873, 14, 6393, 15);

  /* Stage 4 */

  od_butterfly_add_asym(s2, s2h, sc);
  od_butterfly_sub_asym(s0, s0h, s1);
  od_butterfly_add_asym(sf, sfh, se);
  od_butterfly_add_asym(sd, sdh, s3);
  od_butterfly_add_asym(s7, od_rshift1(*s7), s6);
  od_butterfly_sub_asym(s8, od_rshift1(*s8), s9);
  od_butterfly_sub_asym(sa, od_rshift1(*sa), sb);
  od_butterfly_add_asym(s5, od_rshift1(*s5), s4);

  /* Stage 5 */

  /*  10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /*  8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*   3135/8192 = Cos[7*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(sc, s3, 10703, 13, 8867, 14, 3135, 13, NONE);

  /*  10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3870398453221475 */
  /*  8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*   3135/8192 = Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(s2, sd, 10703, 13, 8867, 14, 3135, 13);

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /*  5793/8192 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(sa, s5, 11585, 13, 5793, 13);

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /*  5793/8192 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(s6, s9, 11585, 13, 5793, 13);

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /*  5793/8192 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(se, s1, 11585, 13, 5793, 13);
}

/**
 * 16-point asymmetric Type-IV iDST
 */
static INLINE void od_idst_16_asym(od_coeff *s0, od_coeff *s8,
                                   od_coeff *s4, od_coeff *sc,
                                   od_coeff *s2, od_coeff *sa,
                                   od_coeff *s6, od_coeff *se,
                                   od_coeff *s1, od_coeff *s9,
                                   od_coeff *s5, od_coeff *sd,
                                   od_coeff *s3, od_coeff *sb,
                                   od_coeff *s7, od_coeff *sf) {
  od_coeff s0h;
  od_coeff s2h;
  od_coeff s4h;
  od_coeff s6h;
  od_coeff s8h;
  od_coeff sah;
  od_coeff sch;
  od_coeff sdh;
  od_coeff seh;
  od_coeff sfh;

  /* Stage 5 */

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /*  5793/8192 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(s6, s9, 11585, 13, 5793, 13);

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /*  5793/8192 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(sa, s5, 11585, 13, 5793, 13);

  /* 11585/8192 = Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /*  5793/8192 = Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(se, s1, 11585, 13, 5793, 13);

  /*  10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /*  8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*   3135/8192 = Cos[7*Pi/8]               = 0.7653668647301796 */
  od_rotate_add(sc, s3, 10703, 13, 8867, 14, 3135, 13, NONE);

  /*  10703/8192 = Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3870398453221475 */
  /*  8867/16384 = Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /*   3135/8192 = Cos[3*Pi/8]               = 0.7653668647301796 */
  od_rotate_neg(sd, s2, 10703, 13, 8867, 14, 3135, 13);

  /* Stage 4 */

  od_butterfly_add(s5, NULL, s4);
  od_butterfly_sub(sa, NULL, sb);
  od_butterfly_sub(s8, NULL, s9);
  od_butterfly_add(s7, NULL, s6);
  od_butterfly_add(sd, &sdh, s3);
  od_butterfly_add(sf, &sfh, se);
  od_butterfly_sub(s0, &s0h, s1);
  od_butterfly_add(s2, &s2h, sc);

  /* Stage 3 */

  /*   9633/8192 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* 12873/16384 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /*  6393/32768 = Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(sb, s4, 9633, 13, 12873, 14, 6393, 15);

  /*  11363/8192 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /*  9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /*   4551/8192 = Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_neg(sa, s5, 11363, 13, 9041, 15, 4551, 13);

  /* 22725/16384 = Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /*  9041/32768 = Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* 18205/32768 = Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(s9, s6, 22725, 14, 9041, 15, 18205, 15, NONE);

  /*  9633/8192 = Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /*  1609/2048 = Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* 6393/32768 = Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_add(s8, s7, 9633, 13, 1609, 11, 6393, 15, NONE);

  /* Stage 2 */

  od_butterfly_add_asym(s8, od_rshift1(*s8), s4);
  od_butterfly_add_asym(s7, od_rshift1(*s7), sb);
  od_butterfly_sub_asym(sa, od_rshift1(*sa), s6);
  od_butterfly_sub_asym(s5, od_rshift1(*s5), s9);
  od_butterfly_add_asym(s0, s0h, s3);
  od_butterfly_add_asym(sd, sdh, se);
  od_butterfly_sub_asym(s2, s2h, s1);
  od_butterfly_sub_asym(sf, sfh, sc);

  /* Stage 1 */

  od_butterfly_sub(s0, &s0h, s7);
  od_butterfly_sub(s8, &s8h, sf);
  od_butterfly_add(s4, &s4h, s3);
  od_butterfly_add(sc, &sch, sb);
  od_butterfly_sub(s2, &s2h, s5);
  od_butterfly_sub(sa, &sah, sd);
  od_butterfly_add(s6, &s6h, s1);
  od_butterfly_add(se, &seh, s9);

  /* Stage 0 */

  /* 23143/32768 = (Sin[17*Pi/64] + Cos[17*Pi/64])/2 = 0.7062550401009887 */
  /*   1137/8192 = (Sin[17*Pi/64] - Cos[17*Pi/64])*2 = 0.1387843410158816 */
  /* 44011/32768 = Cos[17*Pi/64]*2                   = 1.3431179096940367 */
  od_rotate_sub_half(s8, s7, s8h, 23143, 15, 1137, 13, 44011, 15, SHIFT);

  /*   2865/4096 = (Sin[19*Pi/64] + Cos[19*Pi/64])/2 = 0.6994534179865391 */
  /* 13599/32768 = (Sin[19*Pi/64] - Cos[19*Pi/64])*2 = 0.4150164539764232 */
  /*     305/256 = Cos[19*Pi/64]*2                   = 1.1913986089848667 */
  od_rotate_add_half(s6, s9, s6h, 2865, 12, 13599, 15, 305, 8, SHIFT);

  /* 5619/8192 = (Sin[21*Pi/64] + Cos[21*Pi/64])/2 = 0.6859156770967569 */
  /* 2815/4096 = (Sin[21*Pi/64] - Cos[21*Pi/64])*2 = 0.6872517316141069 */
  /* 8423/8192 = Cos[21*Pi/64]*2                   = 1.0282054883864433 */
  od_rotate_sub_half(sa, s5, sah, 5619, 13, 2815, 12, 8423, 13, SHIFT);

  /* 2727/4096 = (Sin[23*Pi/64] + Cos[23*Pi/64])/2 = 0.6657721932768628 */
  /* 3903/4096 = (Sin[23*Pi/64] - Cos[23*Pi/64])*2 = 0.9528683993863225 */
  /* 7005/8192 = Cos[23*Pi/64]*2                   = 0.8551101868605642 */
  od_rotate_add_half(s4, sb, s4h, 2727, 12, 3903, 12, 7005, 13, SHIFT);

  /* 10473/16384 = (Sin[25*Pi/64] + Cos[25*Pi/64])/2 = 0.6392169592876205 */
  /* 39627/32768 = (Sin[25*Pi/64] - Cos[25*Pi/64])*2 = 1.2093084235816014 */
  /* 11039/16384 = Cos[25*Pi/64]*2                   = 0.6737797067844401 */
  od_rotate_sub_half(sc, s3, sch, 10473, 14, 39627, 15, 11039, 14, SHIFT);

  /*  9937/16384 = (Sin[27*Pi/64] + Cos[27*Pi/64])/2 = 0.6065057165489039 */
  /*   1489/1024 = (Sin[27*Pi/64] - Cos[27*Pi/64])*2 = 1.4541021465825602 */
  /*   3981/8192 = Cos[27*Pi/64]*2                   = 0.4859603598065277 */
  od_rotate_add_half(s2, sd, s2h, 9937, 14, 1489, 10, 3981, 13, SHIFT);

  /* 18611/32768 = (Sin[29*Pi/64] + Cos[29*Pi/64])/2 = 0.5679534922100714 */
  /* 55211/32768 = (Sin[29*Pi/64] - Cos[29*Pi/64])*2 = 1.6848920710188384 */
  /*    601/2048 = Cos[29*Pi/64]*2                   = 0.2934609489107235 */
  od_rotate_sub_half(se, s1, seh, 18611, 15, 55211, 15, 601, 11, SHIFT);

  /*   1073/2048 = (Sin[31*Pi/64] + Cos[31*Pi/64])/2 = 0.5239315652662953 */
  /* 62241/32768 = (Sin[31*Pi/64] - Cos[31*Pi/64])*2 = 1.8994555637555088 */
  /*    201/2048 = Cos[31*Pi/64]*2                   = 0.0981353486548360 */
  od_rotate_add_half(s0, sf, s0h, 1073, 11, 62241, 15, 201, 11, SHIFT);
}

/* --- 32-point Transforms --- */

/**
 * 32-point orthonormal Type-II fDCT
 */
static INLINE void od_fdct_32(od_coeff *t0, od_coeff *t1,
                              od_coeff *t2, od_coeff *t3,
                              od_coeff *t4, od_coeff *t5,
                              od_coeff *t6, od_coeff *t7,
                              od_coeff *t8, od_coeff *t9,
                              od_coeff *ta, od_coeff *tb,
                              od_coeff *tc, od_coeff *td,
                              od_coeff *te, od_coeff *tf,
                              od_coeff *tg, od_coeff *th,
                              od_coeff *ti, od_coeff *tj,
                              od_coeff *tk, od_coeff *tl,
                              od_coeff *tm, od_coeff *tn,
                              od_coeff *to, od_coeff *tp,
                              od_coeff *tq, od_coeff *tr,
                              od_coeff *ts, od_coeff *tt,
                              od_coeff *tu, od_coeff *tv) {
  od_coeff t1h;
  od_coeff t3h;
  od_coeff t5h;
  od_coeff t7h;
  od_coeff t9h;
  od_coeff tbh;
  od_coeff tdh;
  od_coeff tfh;
  od_coeff thh;
  od_coeff tjh;
  od_coeff tlh;
  od_coeff tnh;
  od_coeff tph;
  od_coeff trh;
  od_coeff tth;
  od_coeff tvh;

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_neg(t0, tv, &tvh);
  od_butterfly_add(t1, &t1h, tu);
  od_butterfly_neg(t2, tt, &tth);
  od_butterfly_add(t3, &t3h, ts);
  od_butterfly_neg(t4, tr, &trh);
  od_butterfly_add(t5, &t5h, tq);
  od_butterfly_neg(t6, tp, &tph);
  od_butterfly_add(t7, &t7h, to);
  od_butterfly_neg(t8, tn, &tnh);
  od_butterfly_add(t9, &t9h, tm);
  od_butterfly_neg(ta, tl, &tlh);
  od_butterfly_add(tb, &tbh, tk);
  od_butterfly_neg(tc, tj, &tjh);
  od_butterfly_add(td, &tdh, ti);
  od_butterfly_neg(te, th, &thh);
  od_butterfly_add(tf, &tfh, tg);

  /* Embedded 16-point transforms with asymmetric input. */
  od_fdct_16_asym(
   t0, t1, t1h, t2, t3, t3h, t4, t5, t5h, t6, t7, t7h,
   t8, t9, t9h, ta, tb, tbh, tc, td, tdh, te, tf, tfh);
  od_fdst_16_asym(
   tv, tvh, tu, tt, tth, ts, tr, trh, tq, tp, tph, to,
   tn, tnh, tm, tl, tlh, tk, tj, tjh, ti, th, thh, tg);
}

/**
 * 32-point orthonormal Type-II iDCT
 */
static INLINE void od_idct_32(od_coeff *t0, od_coeff *tg,
                              od_coeff *t8, od_coeff *to,
                              od_coeff *t4, od_coeff *tk,
                              od_coeff *tc, od_coeff *ts,
                              od_coeff *t2, od_coeff *ti,
                              od_coeff *ta, od_coeff *tq,
                              od_coeff *t6, od_coeff *tm,
                              od_coeff *te, od_coeff *tu,
                              od_coeff *t1, od_coeff *th,
                              od_coeff *t9, od_coeff *tp,
                              od_coeff *t5, od_coeff *tl,
                              od_coeff *td, od_coeff *tt,
                              od_coeff *t3, od_coeff *tj,
                              od_coeff *tb, od_coeff *tr,
                              od_coeff *t7, od_coeff *tn,
                              od_coeff *tf, od_coeff *tv) {
  od_coeff t1h;
  od_coeff t3h;
  od_coeff t5h;
  od_coeff t7h;
  od_coeff t9h;
  od_coeff tbh;
  od_coeff tdh;
  od_coeff tfh;

  /* Embedded 16-point transforms with asymmetric output. */
  od_idst_16_asym(
   tv, tn, tr, tj, tt, tl, tp, th, tu, tm, tq, ti, ts, tk, to, tg);
  od_idct_16_asym(
   t0, t8, t4, tc, t2, ta, t6, te,
   t1, &t1h, t9, &t9h, t5, &t5h, td, &tdh,
   t3, &t3h, tb, &tbh, t7, &t7h, tf, &tfh);

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_add_asym(tf, tfh, tg);
  od_butterfly_neg_asym(te, th, od_rshift1(*th));
  od_butterfly_add_asym(td, tdh, ti);
  od_butterfly_neg_asym(tc, tj, od_rshift1(*tj));
  od_butterfly_add_asym(tb, tbh, tk);
  od_butterfly_neg_asym(ta, tl, od_rshift1(*tl));
  od_butterfly_add_asym(t9, t9h, tm);
  od_butterfly_neg_asym(t8, tn, od_rshift1(*tn));
  od_butterfly_add_asym(t7, t7h, to);
  od_butterfly_neg_asym(t6, tp, od_rshift1(*tp));
  od_butterfly_add_asym(t5, t5h, tq);
  od_butterfly_neg_asym(t4, tr, od_rshift1(*tr));
  od_butterfly_add_asym(t3, t3h, ts);
  od_butterfly_neg_asym(t2, tt, od_rshift1(*tt));
  od_butterfly_add_asym(t1, t1h, tu);
  od_butterfly_neg_asym(t0, tv, od_rshift1(*tv));
}

/**
 * 32-point asymmetric Type-II fDCT
 */
static INLINE void od_fdct_32_asym(od_coeff *t0, od_coeff *t1, od_coeff t1h,
                                   od_coeff *t2, od_coeff *t3, od_coeff t3h,
                                   od_coeff *t4, od_coeff *t5, od_coeff t5h,
                                   od_coeff *t6, od_coeff *t7, od_coeff t7h,
                                   od_coeff *t8, od_coeff *t9, od_coeff t9h,
                                   od_coeff *ta, od_coeff *tb, od_coeff tbh,
                                   od_coeff *tc, od_coeff *td, od_coeff tdh,
                                   od_coeff *te, od_coeff *tf, od_coeff tfh,
                                   od_coeff *tg, od_coeff *th, od_coeff thh,
                                   od_coeff *ti, od_coeff *tj, od_coeff tjh,
                                   od_coeff *tk, od_coeff *tl, od_coeff tlh,
                                   od_coeff *tm, od_coeff *tn, od_coeff tnh,
                                   od_coeff *to, od_coeff *tp, od_coeff tph,
                                   od_coeff *tq, od_coeff *tr, od_coeff trh,
                                   od_coeff *ts, od_coeff *tt, od_coeff tth,
                                   od_coeff *tu, od_coeff *tv, od_coeff tvh) {

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_neg_asym(t0, tv, tvh);
  od_butterfly_sub_asym(t1, t1h, tu);
  od_butterfly_neg_asym(t2, tt, tth);
  od_butterfly_sub_asym(t3, t3h, ts);
  od_butterfly_neg_asym(t4, tr, trh);
  od_butterfly_sub_asym(t5, t5h, tq);
  od_butterfly_neg_asym(t6, tp, tph);
  od_butterfly_sub_asym(t7, t7h, to);
  od_butterfly_neg_asym(t8, tn, tnh);
  od_butterfly_sub_asym(t9, t9h, tm);
  od_butterfly_neg_asym(ta, tl, tlh);
  od_butterfly_sub_asym(tb, tbh, tk);
  od_butterfly_neg_asym(tc, tj, tjh);
  od_butterfly_sub_asym(td, tdh, ti);
  od_butterfly_neg_asym(te, th, thh);
  od_butterfly_sub_asym(tf, tfh, tg);

  /* Embedded 16-point orthonormal transforms. */
  od_fdct_16(t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, ta, tb, tc, td, te, tf);
  od_fdst_16(tv, tu, tt, ts, tr, tq, tp, to, tn, tm, tl, tk, tj, ti, th, tg);
}

/**
 * 32-point asymmetric Type-II iDCT
 */
static INLINE void od_idct_32_asym(od_coeff *t0, od_coeff *tg,
                                   od_coeff *t8, od_coeff *to,
                                   od_coeff *t4, od_coeff *tk,
                                   od_coeff *tc, od_coeff *ts,
                                   od_coeff *t2, od_coeff *ti,
                                   od_coeff *ta, od_coeff *tq,
                                   od_coeff *t6, od_coeff *tm,
                                   od_coeff *te, od_coeff *tu,
                                   od_coeff *t1, od_coeff *t1h,
                                   od_coeff *th, od_coeff *thh,
                                   od_coeff *t9, od_coeff *t9h,
                                   od_coeff *tp, od_coeff *tph,
                                   od_coeff *t5, od_coeff *t5h,
                                   od_coeff *tl, od_coeff *tlh,
                                   od_coeff *td, od_coeff *tdh,
                                   od_coeff *tt, od_coeff *tth,
                                   od_coeff *t3, od_coeff *t3h,
                                   od_coeff *tj, od_coeff *tjh,
                                   od_coeff *tb, od_coeff *tbh,
                                   od_coeff *tr, od_coeff *trh,
                                   od_coeff *t7, od_coeff *t7h,
                                   od_coeff *tn, od_coeff *tnh,
                                   od_coeff *tf, od_coeff *tfh,
                                   od_coeff *tv, od_coeff *tvh) {

  /* Embedded 16-point orthonormal transforms. */
  od_idst_16(tv, tn, tr, tj, tt, tl, tp, th, tu, tm, tq, ti, ts, tk, to, tg);
  od_idct_16(t0, t8, t4, tc, t2, ta, t6, te, t1, t9, t5, td, t3, tb, t7, tf);

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_sub(tf, tfh, tg);
  od_butterfly_neg(te, th, thh);
  od_butterfly_sub(td, tdh, ti);
  od_butterfly_neg(tc, tj, tjh);
  od_butterfly_sub(tb, tbh, tk);
  od_butterfly_neg(ta, tl, tlh);
  od_butterfly_sub(t9, t9h, tm);
  od_butterfly_neg(t8, tn, tnh);
  od_butterfly_sub(t7, t7h, to);
  od_butterfly_neg(t6, tp, tph);
  od_butterfly_sub(t5, t5h, tq);
  od_butterfly_neg(t4, tr, trh);
  od_butterfly_sub(t3, t3h, ts);
  od_butterfly_neg(t2, tt, tth);
  od_butterfly_sub(t1, t1h, tu);
  od_butterfly_neg(t0, tv, tvh);
}

/**
 * 32-point orthonormal Type-IV fDST
 */
static INLINE void od_fdst_32(od_coeff *t0, od_coeff *t1,
                              od_coeff *t2, od_coeff *t3,
                              od_coeff *t4, od_coeff *t5,
                              od_coeff *t6, od_coeff *t7,
                              od_coeff *t8, od_coeff *t9,
                              od_coeff *ta, od_coeff *tb,
                              od_coeff *tc, od_coeff *td,
                              od_coeff *te, od_coeff *tf,
                              od_coeff *tg, od_coeff *th,
                              od_coeff *ti, od_coeff *tj,
                              od_coeff *tk, od_coeff *tl,
                              od_coeff *tm, od_coeff *tn,
                              od_coeff *to, od_coeff *tp,
                              od_coeff *tq, od_coeff *tr,
                              od_coeff *ts, od_coeff *tt,
                              od_coeff *tu, od_coeff *tv) {
  od_coeff t0h;
  od_coeff t1h;
  od_coeff t2h;
  od_coeff t3h;
  od_coeff t4h;
  od_coeff t6h;
  od_coeff t8h;
  od_coeff t9h;
  od_coeff tah;
  od_coeff tbh;
  od_coeff tch;
  od_coeff tdh;
  od_coeff teh;
  od_coeff tfh;
  od_coeff tgh;
  od_coeff thh;
  od_coeff tih;
  od_coeff tjh;
  od_coeff tkh;
  od_coeff tlh;
  od_coeff tmh;
  od_coeff tnh;
  od_coeff tph;
  od_coeff trh;
  od_coeff tsh;
  od_coeff tth;
  od_coeff tuh;
  od_coeff tvh;

  /* Stage 0 */

  /* Sin[63*Pi/128] + Cos[63*Pi/128] = 1.0242400472191164 */
  /* Sin[63*Pi/128] - Cos[63*Pi/128] = 0.9751575901732919 */
  /* Cos[63*Pi/128]                  = 0.0245412285229123 */
  od_rotate_add(t0, tv, 16781, 14, 15977, 14, 201, 13, NONE);

  /* Sin[61*Pi/128] + Cos[61*Pi/128] = 1.0708550202783576 */
  /* Sin[61*Pi/128] - Cos[61*Pi/128] = 0.9237258930790228 */
  /* Cos[61*Pi/128]                  = 0.0735645635996674 */
  od_rotate_sub(tu, t1, 17545, 14, 30269, 15, 2411, 15, NONE);

  /* Sin[59*Pi/128] + Cos[59*Pi/128] = 1.1148902097979262 */
  /* Sin[59*Pi/128] - Cos[59*Pi/128] = 0.8700688593994937 */
  /* Cos[59*Pi/128]                  = 0.1224106751992162 */
  od_rotate_add(t2, tt, 36533, 15, 14255, 14, 4011, 15, NONE);

  /* Sin[57*Pi/128] + Cos[57*Pi/128] = 1.1562395311492424 */
  /* Sin[57*Pi/128] - Cos[57*Pi/128] = 0.8143157536286401 */
  /* Cos[57*Pi/128]                  = 0.1709618887603012 */
  od_rotate_sub(ts, t3, 37, 5, 26683, 15, 2801, 14, NONE);

  /* Sin[55*Pi/128] + Cos[55*Pi/128] = 1.1948033701953984 */
  /* Sin[55*Pi/128] - Cos[55*Pi/128] = 0.7566008898816587 */
  /* Cos[55*Pi/128]                  = 0.2191012401568698 */
  od_rotate_add(t4, tr, 39151, 15, 3099, 12, 1795, 13, NONE);

  /* Sin[53*Pi/128] + Cos[53*Pi/128] = 1.2304888232703382 */
  /* Sin[53*Pi/128] - Cos[53*Pi/128] = 0.6970633083205415 */
  /* Cos[53*Pi/128]                  = 0.2667127574748984 */
  od_rotate_sub(tq, t5, 40321, 15, 22841, 15, 2185, 13, NONE);

  /* Sin[51*Pi/128] + Cos[51*Pi/128] = 1.2632099209919283 */
  /* Sin[51*Pi/128] - Cos[51*Pi/128] = 0.6358464401941452 */
  /* Cos[51*Pi/128]                  = 0.3136817403988915 */
  od_rotate_add(t6, tp, 41393, 15, 20835, 15, 10279, 15, NONE);

  /* Sin[49*Pi/128] + Cos[49*Pi/128] = 1.2928878353697270 */
  /* Sin[49*Pi/128] - Cos[49*Pi/128] = 0.5730977622997508 */
  /* Cos[49*Pi/128]                  = 0.3598950365349881 */
  od_rotate_sub(to, t7, 42365, 15, 18778, 15, 11793, 15, NONE);

  /* Sin[47*Pi/128] + Cos[47*Pi/128] = 1.3194510697085207 */
  /* Sin[47*Pi/128] - Cos[47*Pi/128] = 0.5089684416985408 */
  /* Cos[47*Pi/128]                  = 0.4052413140049899 */
  od_rotate_add(t8, tn, 10809, 13, 8339, 14, 13279, 15, NONE);

  /* Sin[45*Pi/128] + Cos[45*Pi/128] = 1.3428356308501219 */
  /* Sin[45*Pi/128] - Cos[45*Pi/128] = 0.4436129715409088 */
  /* Cos[45*Pi/128]                  = 0.4496113296546065 */
  od_rotate_sub(tm, t9, 22001, 14, 1817, 12, 14733, 15, NONE);

  /* Sin[43*Pi/128] + Cos[43*Pi/128] = 1.3629851833384956 */
  /* Sin[43*Pi/128] - Cos[43*Pi/128] = 0.3771887988789274 */
  /* Cos[43*Pi/128]                  = 0.4928981922297840 */
  od_rotate_add(ta, tl, 22331, 14, 1545, 12, 16151, 15, NONE);

  /* Sin[41*Pi/128] + Cos[41*Pi/128] = 1.3798511851368043 */
  /* Sin[41*Pi/128] - Cos[41*Pi/128] = 0.3098559453626100 */
  /* Cos[41*Pi/128]                  = 0.5349976198870972 */
  od_rotate_sub(tk, tb, 45215, 15, 10153, 15, 17531, 15, NONE);

  /* Sin[39*Pi/128] + Cos[39*Pi/128] = 1.3933930045694290 */
  /* Sin[39*Pi/128] - Cos[39*Pi/128] = 0.2417766217337384 */
  /* Cos[39*Pi/128]                  = 0.5758081914178453 */
  od_rotate_add(tc, tj, 45659, 15, 7923, 15, 4717, 13, NONE);

  /* Sin[37*Pi/128] + Cos[37*Pi/128] = 1.4035780182072331 */
  /* Sin[37*Pi/128] - Cos[37*Pi/128] = 0.1731148370459795 */
  /* Cos[37*Pi/128]                  = 0.6152315905806268 */
  od_rotate_sub(ti, td, 5749, 12, 5673, 15, 315, 9, NONE);

  /* Sin[35*Pi/128] + Cos[35*Pi/128] = 1.4103816894602614 */
  /* Sin[35*Pi/128] - Cos[35*Pi/128] = 0.1040360035527078 */
  /* Cos[35*Pi/128]                  = 0.6531728429537768 */
  od_rotate_add(te, th, 46215, 15, 3409, 15, 21403, 15, NONE);

  /* Sin[33*Pi/128] + Cos[33*Pi/128] = 1.4137876276885337 */
  /* Sin[33*Pi/128] - Cos[33*Pi/128] = 0.0347065382144002 */
  /* Cos[33*Pi/128]                  = 0.6895405447370668 */
  od_rotate_sub(tg, tf, 46327, 15, 1137, 15, 22595, 15, NONE);

  /* Stage 1 */

  od_butterfly_add(t0, &t0h, tf);
  od_butterfly_sub(tv, &tvh, tg);
  od_butterfly_add(th, &thh, tu);
  od_butterfly_sub(te, &teh, t1);

  od_butterfly_add(t2, &t2h, td);
  od_butterfly_sub(tt, &tth, ti);
  od_butterfly_add(tj, &tjh, ts);
  od_butterfly_sub(tc, &tch, t3);

  od_butterfly_add(t4, &t4h, tb);
  od_butterfly_sub(tr, &trh, tk);
  od_butterfly_add(tl, &tlh, tq);
  od_butterfly_sub(ta, &tah, t5);

  od_butterfly_add(t6, &t6h, t9);
  od_butterfly_sub(tp, &tph, tm);
  od_butterfly_add(tn, &tnh, to);
  od_butterfly_sub(t8, &t8h, t7);

  /* Stage 2 */

  od_butterfly_sub_asym(t0, t0h, t7);
  od_butterfly_add_asym(tv, tvh, to);
  od_butterfly_sub_asym(tp, tph, tu);
  od_butterfly_add_asym(t6, t6h, t1);

  od_butterfly_sub_asym(t2, t2h, t5);
  od_butterfly_add_asym(tt, tth, tq);
  od_butterfly_sub_asym(tr, trh, ts);
  od_butterfly_add_asym(t4, t4h, t3);

  od_butterfly_add_asym(t8, t8h, tg);
  od_butterfly_sub_asym(te, teh, tm);
  od_butterfly_add_asym(tn, tnh, tf);
  od_butterfly_sub_asym(th, thh, t9);

  od_butterfly_add_asym(ta, tah, ti);
  od_butterfly_sub_asym(tc, tch, tk);
  od_butterfly_add_asym(tl, tlh, td);
  od_butterfly_sub_asym(tj, tjh, tb);

  /* Stage 3 */

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_sub(tf, tg, 17911, 14, 14699, 14, 803, 13, NONE);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_add(th, te, 20435, 14, 21845, 15, 1189, 12, NONE);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_add(ti, td, 22173, 14, 3363, 13, 15447, 15, NONE);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_sub(tc, tj, 23059, 14, 2271, 14, 5197, 13, NONE);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_neg(tb, tk, 23059, 14, 2271, 14, 5197, 13);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_neg(ta, tl, 22173, 14, 3363, 13, 15447, 15);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_neg(t9, tm, 20435, 14, 21845, 15, 1189, 12);

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_neg(t8, tn, 17911, 14, 14699, 14, 803, 13);

  /* Stage 4 */

  od_butterfly_sub(t3, &t3h, t0);
  od_butterfly_add(ts, &tsh, tv);
  od_butterfly_sub(tu, &tuh, tt);
  od_butterfly_add(t1, &t1h, t2);

  od_butterfly_add(to, NULL, t4);
  od_butterfly_sub(tq, NULL, t6);
  od_butterfly_add(t7, NULL, tr);
  od_butterfly_sub(t5, NULL, tp);

  od_butterfly_sub(tb, &tbh, t8);
  od_butterfly_add(tk, &tkh, tn);
  od_butterfly_sub(tm, &tmh, tl);
  od_butterfly_add(t9, &t9h, ta);

  od_butterfly_sub(tf, &tfh, tc);
  od_butterfly_add(tg, &tgh, tj);
  od_butterfly_sub(ti, &tih, th);
  od_butterfly_add(td, &tdh, te);

  /* Stage 5 */

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_add(to, t7, 9633, 13, 12873, 14, 6393, 15, NONE);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(tp, t6, 22725, 14, 9041, 15, 4551, 13, NONE);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_neg(t5, tq, 11363, 13, 9041, 15, 4551, 13);

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(t4, tr, 9633, 13, 12873, 14, 6393, 15);

  /* Stage 6 */

  od_butterfly_add_asym(t1, t1h, t0);
  od_butterfly_sub_asym(tu, tuh, tv);
  od_butterfly_sub_asym(ts, tsh, t2);
  od_butterfly_sub_asym(t3, t3h, tt);

  od_butterfly_add_asym(t5, od_rshift1(*t5), t4);
  od_butterfly_sub_asym(tq, od_rshift1(*tq), tr);
  od_butterfly_add_asym(t7, od_rshift1(*t7), t6);
  od_butterfly_sub_asym(to, od_rshift1(*to), tp);

  od_butterfly_add_asym(t9, t9h, t8);
  od_butterfly_sub_asym(tm, tmh, tn);
  od_butterfly_sub_asym(tk, tkh, ta);
  od_butterfly_sub_asym(tb, tbh, tl);

  od_butterfly_add_asym(ti, tih, tc);
  od_butterfly_add_asym(td, tdh, tj);
  od_butterfly_add_asym(tf, tfh, te);
  od_butterfly_sub_asym(tg, tgh, th);

  /* Stage 7 */

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(t2, tt, 10703, 13, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(ts, t3, 10703, 13, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(ta, tl, 10703, 13, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tk, tb, 10703, 13, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tc, tj, 10703, 13, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(ti, td, 10703, 13, 8867, 14, 3135, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tu, t1, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tq, t5, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_sub(tp, t6, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tm, t9, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(te, th, 11585, 13, 5793, 13);
}

/**
 * 32-point orthonormal Type-IV fDST
 */
static INLINE void od_idst_32(od_coeff *t0, od_coeff *tg,
                              od_coeff *t8, od_coeff *to,
                              od_coeff *t4, od_coeff *tk,
                              od_coeff *tc, od_coeff *ts,
                              od_coeff *t2, od_coeff *ti,
                              od_coeff *ta, od_coeff *tq,
                              od_coeff *t6, od_coeff *tm,
                              od_coeff *te, od_coeff *tu,
                              od_coeff *t1, od_coeff *th,
                              od_coeff *t9, od_coeff *tp,
                              od_coeff *t5, od_coeff *tl,
                              od_coeff *td, od_coeff *tt,
                              od_coeff *t3, od_coeff *tj,
                              od_coeff *tb, od_coeff *tr,
                              od_coeff *t7, od_coeff *tn,
                              od_coeff *tf, od_coeff *tv) {
  od_coeff t0h;
  od_coeff t1h;
  od_coeff t2h;
  od_coeff t3h;
  od_coeff t4h;
  od_coeff t6h;
  od_coeff t8h;
  od_coeff t9h;
  od_coeff tah;
  od_coeff tbh;
  od_coeff tch;
  od_coeff tdh;
  od_coeff teh;
  od_coeff tfh;
  od_coeff tgh;
  od_coeff thh;
  od_coeff tih;
  od_coeff tjh;
  od_coeff tkh;
  od_coeff tlh;
  od_coeff tmh;
  od_coeff tnh;
  od_coeff tph;
  od_coeff trh;
  od_coeff tsh;
  od_coeff tth;
  od_coeff tuh;
  od_coeff tvh;

  /* Stage 7 */

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(te, th, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tm, t9, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_sub(tp, t6, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tq, t5, 11585, 13, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tu, t1, 11585, 13, 5793, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(td, ti, 10703, 13, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tc, tj, 10703, 13, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tk, tb, 10703, 13, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(tl, ta, 10703, 13, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(ts, t3, 10703, 13, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(tt, t2, 10703, 13, 8867, 14, 3135, 13);

  /* Stage 6 */

  od_butterfly_sub(tg, &tgh, th);
  od_butterfly_add(tf, &tfh, te);
  od_butterfly_add(td, &tdh, tj);
  od_butterfly_add(ti, &tih, tc);

  od_butterfly_sub(tb, &tbh, tl);
  od_butterfly_sub(tk, &tkh, ta);
  od_butterfly_sub(tm, &tmh, tn);
  od_butterfly_add(t9, &t9h, t8);

  od_butterfly_sub(to, NULL, tp);
  od_butterfly_add(t7, NULL, t6);
  od_butterfly_sub(tq, NULL, tr);
  od_butterfly_add(t5, NULL, t4);

  od_butterfly_sub(t3, &t3h, tt);
  od_butterfly_sub(ts, &tsh, t2);
  od_butterfly_sub(tu, &tuh, tv);
  od_butterfly_add(t1, &t1h, t0);

  /* Stage 5 */

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(tr, t4, 9633, 13, 12873, 14, 6393, 15);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_neg(tq, t5, 11363, 13, 9041, 15, 4551, 13);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(tp, t6, 22725, 14, 9041, 15, 4551, 13, NONE);

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_add(to, t7, 9633, 13, 12873, 14, 6393, 15, NONE);

  /* Stage 4 */

  od_butterfly_add_asym(td, tdh, te);
  od_butterfly_sub_asym(ti, tih, th);
  od_butterfly_add_asym(tg, tgh, tj);
  od_butterfly_sub_asym(tf, tfh, tc);

  od_butterfly_add_asym(t9, t9h, ta);
  od_butterfly_sub_asym(tm, tmh, tl);
  od_butterfly_add_asym(tk, tkh, tn);
  od_butterfly_sub_asym(tb, tbh, t8);

  od_butterfly_sub_asym(t5, od_rshift1(*t5), tp);
  od_butterfly_add_asym(t7, od_rshift1(*t7), tr);
  od_butterfly_sub_asym(tq, od_rshift1(*tq), t6);
  od_butterfly_add_asym(to, od_rshift1(*to), t4);

  od_butterfly_add_asym(t1, t1h, t2);
  od_butterfly_sub_asym(tu, tuh, tt);
  od_butterfly_add_asym(ts, tsh, tv);
  od_butterfly_sub_asym(t3, t3h, t0);

  /* Stage 3 */

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_neg(tn, t8, 17911, 14, 14699, 14, 803, 13);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_neg(tm, t9, 20435, 14, 21845, 15, 1189, 12);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_neg(tl, ta, 22173, 14, 3363, 13, 15447, 15);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_neg(tk, tb, 23059, 14, 2271, 14, 5197, 13);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_sub(tc, tj, 23059, 14, 2271, 14, 5197, 13, NONE);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_add(ti, td, 22173, 14, 3363, 13, 15447, 15, NONE);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_add(th, te, 20435, 14, 21845, 15, 1189, 12, NONE);

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_sub(tf, tg, 17911, 14, 14699, 14, 803, 13, NONE);

  /* Stage 2 */

  od_butterfly_sub(tj, &tjh, tb);
  od_butterfly_add(tl, &tlh, td);
  od_butterfly_sub(tc, &tch, tk);
  od_butterfly_add(ta, &tah, ti);

  od_butterfly_sub(th, &thh, t9);
  od_butterfly_add(tn, &tnh, tf);
  od_butterfly_sub(te, &teh, tm);
  od_butterfly_add(t8, &t8h, tg);

  od_butterfly_add(t4, &t4h, t3);
  od_butterfly_sub(tr, &trh, ts);
  od_butterfly_add(tt, &tth, tq);
  od_butterfly_sub(t2, &t2h, t5);

  od_butterfly_add(t6, &t6h, t1);
  od_butterfly_sub(tp, &tph, tu);
  od_butterfly_add(tv, &tvh, to);
  od_butterfly_sub(t0, &t0h, t7);

  /* Stage 1 */

  od_butterfly_sub_asym(t8, t8h, t7);
  od_butterfly_add_asym(tn, tnh, to);
  od_butterfly_sub_asym(tp, tph, tm);
  od_butterfly_add_asym(t6, t6h, t9);

  od_butterfly_sub_asym(ta, tah, t5);
  od_butterfly_add_asym(tl, tlh, tq);
  od_butterfly_sub_asym(tr, trh, tk);
  od_butterfly_add_asym(t4, t4h, tb);

  od_butterfly_sub_asym(tc, tch, t3);
  od_butterfly_add_asym(tj, tjh, ts);
  od_butterfly_sub_asym(tt, tth, ti);
  od_butterfly_add_asym(t2, t2h, td);

  od_butterfly_sub_asym(te, teh, t1);
  od_butterfly_add_asym(th, thh, tu);
  od_butterfly_sub_asym(tv, tvh, tg);
  od_butterfly_add_asym(t0, t0h, tf);

  /* Stage 0 */

  /* Sin[33*Pi/128] + Cos[33*Pi/128] = 1.4137876276885337 */
  /* Sin[33*Pi/128] - Cos[33*Pi/128] = 0.0347065382144002 */
  /* Cos[33*Pi/128]                  = 0.6895405447370668 */
  od_rotate_sub(tg, tf, 46327, 15, 1137, 15, 22595, 15, NONE);

  /* Sin[35*Pi/128] + Cos[35*Pi/128] = 1.4103816894602614 */
  /* Sin[35*Pi/128] - Cos[35*Pi/128] = 0.1040360035527078 */
  /* Cos[35*Pi/128]                  = 0.6531728429537768 */
  od_rotate_add(te, th, 46215, 15, 3409, 15, 21403, 15, NONE);

  /* Sin[37*Pi/128] + Cos[37*Pi/128] = 1.4035780182072331 */
  /* Sin[37*Pi/128] - Cos[37*Pi/128] = 0.1731148370459795 */
  /* Cos[37*Pi/128]                  = 0.6152315905806268 */
  od_rotate_sub(ti, td, 5749, 12, 5673, 15, 315, 9, NONE);

  /* Sin[39*Pi/128] + Cos[39*Pi/128] = 1.3933930045694290 */
  /* Sin[39*Pi/128] - Cos[39*Pi/128] = 0.2417766217337384 */
  /* Cos[39*Pi/128]                  = 0.5758081914178453 */
  od_rotate_add(tc, tj, 45659, 15, 7923, 15, 4717, 13, NONE);

  /* Sin[41*Pi/128] + Cos[41*Pi/128] = 1.3798511851368043 */
  /* Sin[41*Pi/128] - Cos[41*Pi/128] = 0.3098559453626100 */
  /* Cos[41*Pi/128]                  = 0.5349976198870972 */
  od_rotate_sub(tk, tb, 45215, 15, 10153, 15, 17531, 15, NONE);

  /* Sin[43*Pi/128] + Cos[43*Pi/128] = 1.3629851833384956 */
  /* Sin[43*Pi/128] - Cos[43*Pi/128] = 0.3771887988789274 */
  /* Cos[43*Pi/128]                  = 0.4928981922297840 */
  od_rotate_add(ta, tl, 22331, 14, 1545, 12, 16151, 15, NONE);

  /* Sin[45*Pi/128] + Cos[45*Pi/128] = 1.3428356308501219 */
  /* Sin[45*Pi/128] - Cos[45*Pi/128] = 0.4436129715409088 */
  /* Cos[45*Pi/128]                  = 0.4496113296546065 */
  od_rotate_sub(tm, t9, 22001, 14, 1817, 12, 14733, 15, NONE);

  /* Sin[47*Pi/128] + Cos[47*Pi/128] = 1.3194510697085207 */
  /* Sin[47*Pi/128] - Cos[47*Pi/128] = 0.5089684416985408 */
  /* Cos[47*Pi/128]                  = 0.4052413140049899 */
  od_rotate_add(t8, tn, 10809, 13, 8339, 14, 13279, 15, NONE);

  /* Sin[49*Pi/128] + Cos[49*Pi/128] = 1.2928878353697270 */
  /* Sin[49*Pi/128] - Cos[49*Pi/128] = 0.5730977622997508 */
  /* Cos[49*Pi/128]                  = 0.3598950365349881 */
  od_rotate_sub(to, t7, 42365, 15, 18779, 15, 11793, 15, NONE);

  /* Sin[51*Pi/128] + Cos[51*Pi/128] = 1.2632099209919283 */
  /* Sin[51*Pi/128] - Cos[51*Pi/128] = 0.6358464401941452 */
  /* Cos[51*Pi/128]                  = 0.3136817403988915 */
  od_rotate_add(t6, tp, 41393, 15, 20835, 15, 10279, 15, NONE);

  /* Sin[53*Pi/128] + Cos[53*Pi/128] = 1.2304888232703382 */
  /* Sin[53*Pi/128] - Cos[53*Pi/128] = 0.6970633083205415 */
  /* Cos[53*Pi/128]                  = 0.2667127574748984 */
  od_rotate_sub(tq, t5, 40321, 15, 22841, 15, 2185, 13, NONE);

  /* Sin[55*Pi/128] + Cos[55*Pi/128] = 1.1948033701953984 */
  /* Sin[55*Pi/128] - Cos[55*Pi/128] = 0.7566008898816587 */
  /* Cos[55*Pi/128]                  = 0.2191012401568698 */
  od_rotate_add(t4, tr, 39151, 15, 3099, 12, 1795, 13, NONE);

  /* Sin[57*Pi/128] + Cos[57*Pi/128] = 1.1562395311492424 */
  /* Sin[57*Pi/128] - Cos[57*Pi/128] = 0.8143157536286401 */
  /* Cos[57*Pi/128]                  = 0.1709618887603012 */
  od_rotate_sub(ts, t3, 37, 5, 26683, 15, 2801, 14, NONE);

  /* Sin[59*Pi/128] + Cos[59*Pi/128] = 1.1148902097979262 */
  /* Sin[59*Pi/128] - Cos[59*Pi/128] = 0.8700688593994937 */
  /* Cos[59*Pi/128]                  = 0.1224106751992162 */
  od_rotate_add(t2, tt, 36533, 15, 14255, 14, 4011, 15, NONE);

  /* Sin[61*Pi/128] + Cos[61*Pi/128] = 1.0708550202783576 */
  /* Sin[61*Pi/128] - Cos[61*Pi/128] = 0.9237258930790228 */
  /* Cos[61*Pi/128]                  = 0.0735645635996674 */
  od_rotate_sub(tu, t1, 17545, 14, 30269, 15, 2411, 15, NONE);

  /* Sin[63*Pi/128] + Cos[63*Pi/128] = 1.0242400472191164 */
  /* Sin[63*Pi/128] - Cos[63*Pi/128] = 0.9751575901732919 */
  /* Cos[63*Pi/128]                  = 0.0245412285229123 */
  od_rotate_add(t0, tv, 16781, 14, 15977, 14, 201, 13, NONE);
}

/**
 * 32-point orthonormal Type-IV fDST
 */
static INLINE void od_fdst_32_asym(od_coeff *t0, od_coeff t0h, od_coeff *t1,
                                   od_coeff *t2, od_coeff t2h, od_coeff *t3,
                                   od_coeff *t4, od_coeff t4h, od_coeff *t5,
                                   od_coeff *t6, od_coeff t6h, od_coeff *t7,
                                   od_coeff *t8, od_coeff t8h, od_coeff *t9,
                                   od_coeff *ta, od_coeff tah, od_coeff *tb,
                                   od_coeff *tc, od_coeff tch, od_coeff *td,
                                   od_coeff *te, od_coeff teh, od_coeff *tf,
                                   od_coeff *tg, od_coeff tgh, od_coeff *th,
                                   od_coeff *ti, od_coeff tih, od_coeff *tj,
                                   od_coeff *tk, od_coeff tkh, od_coeff *tl,
                                   od_coeff *tm, od_coeff tmh, od_coeff *tn,
                                   od_coeff *to, od_coeff toh, od_coeff *tp,
                                   od_coeff *tq, od_coeff tqh, od_coeff *tr,
                                   od_coeff *ts, od_coeff tsh, od_coeff *tt,
                                   od_coeff *tu, od_coeff tuh, od_coeff *tv) {
  od_coeff t1h;
  od_coeff t3h;
  od_coeff t9h;
  od_coeff tbh;
  od_coeff tdh;
  od_coeff tfh;
  od_coeff thh;
  od_coeff tjh;
  od_coeff tlh;
  od_coeff tnh;
  od_coeff tph;
  od_coeff trh;
  od_coeff tth;
  od_coeff tvh;

  /* Stage 0 */

  /* Sin[63*Pi/128] + Cos[63*Pi/128] = 1.0242400472191164 */
  /* Sin[63*Pi/128] - Cos[63*Pi/128] = 0.9751575901732919 */
  /* Cos[63*Pi/128]                  = 0.0245412285229123 */
  od_rotate_add_half(t0, tv, t0h, 5933, 13, 22595, 14, 1137, 15, NONE);

  /* Sin[61*Pi/128] + Cos[61*Pi/128] = 1.0708550202783576 */
  /* Sin[61*Pi/128] - Cos[61*Pi/128] = 0.9237258930790228 */
  /* Cos[61*Pi/128]                  = 0.0735645635996674 */
  od_rotate_sub_half(tu, t1, tuh, 6203, 13, 21403, 14, 3409, 15, NONE);

  /* Sin[59*Pi/128] + Cos[59*Pi/128] = 1.1148902097979262 */
  /* Sin[59*Pi/128] - Cos[59*Pi/128] = 0.8700688593994937 */
  /* Cos[59*Pi/128]                  = 0.1224106751992162 */
  od_rotate_add_half(t2, tt, t2h, 25833, 15, 315, 8, 5673, 15, NONE);

  /* Sin[57*Pi/128] + Cos[57*Pi/128] = 1.1562395311492424 */
  /* Sin[57*Pi/128] - Cos[57*Pi/128] = 0.8143157536286401 */
  /* Cos[57*Pi/128]                  = 0.1709618887603012 */
  od_rotate_sub_half(ts, t3, tsh, 26791, 15, 4717, 12, 7923, 15, NONE);

  /* Sin[55*Pi/128] + Cos[55*Pi/128] = 1.1948033701953984 */
  /* Sin[55*Pi/128] - Cos[55*Pi/128] = 0.7566008898816587 */
  /* Cos[55*Pi/128]                  = 0.2191012401568698 */
  od_rotate_add_half(t4, tr, t4h, 6921, 13, 17531, 14, 10153, 15, NONE);

  /* Sin[53*Pi/128] + Cos[53*Pi/128] = 1.2304888232703382 */
  /* Sin[53*Pi/128] - Cos[53*Pi/128] = 0.6970633083205415 */
  /* Cos[53*Pi/128]                  = 0.2667127574748984 */
  od_rotate_sub_half(tq, t5, tqh, 28511, 15, 32303, 15, 1545, 12, NONE);

  /* Sin[51*Pi/128] + Cos[51*Pi/128] = 1.2632099209919283 */
  /* Sin[51*Pi/128] - Cos[51*Pi/128] = 0.6358464401941452 */
  /* Cos[51*Pi/128]                  = 0.3136817403988915 */
  od_rotate_add_half(t6, tp, t6h, 29269, 15, 14733, 14, 1817, 12, NONE);

  /* Sin[49*Pi/128] + Cos[49*Pi/128] = 1.2928878353697270 */
  /* Sin[49*Pi/128] - Cos[49*Pi/128] = 0.5730977622997508 */
  /* Cos[49*Pi/128]                  = 0.3598950365349881 */
  od_rotate_sub_half(to, t7, toh, 29957, 15, 13279, 14, 8339, 14, NONE);

  /* Sin[47*Pi/128] + Cos[47*Pi/128] = 1.3194510697085207 */
  /* Sin[47*Pi/128] - Cos[47*Pi/128] = 0.5089684416985408 */
  /* Cos[47*Pi/128]                  = 0.4052413140049899 */
  od_rotate_add_half(t8, tn, t8h, 7643, 13, 11793, 14, 18779, 15, NONE);

  /* Sin[45*Pi/128] + Cos[45*Pi/128] = 1.3428356308501219 */
  /* Sin[45*Pi/128] - Cos[45*Pi/128] = 0.4436129715409088 */
  /* Cos[45*Pi/128]                  = 0.4496113296546065 */
  od_rotate_sub_half(tm, t9, tmh, 15557, 14, 20557, 15, 20835, 15, NONE);

  /* Sin[43*Pi/128] + Cos[43*Pi/128] = 1.3629851833384956 */
  /* Sin[43*Pi/128] - Cos[43*Pi/128] = 0.3771887988789274 */
  /* Cos[43*Pi/128]                  = 0.4928981922297840 */
  od_rotate_add_half(ta, tl, tah, 31581, 15, 17479, 15, 22841, 15, NONE);

  /* Sin[41*Pi/128] + Cos[41*Pi/128] = 1.3798511851368043 */
  /* Sin[41*Pi/128] - Cos[41*Pi/128] = 0.3098559453626100 */
  /* Cos[41*Pi/128]                  = 0.5349976198870972 */
  od_rotate_sub_half(tk, tb, tkh, 7993, 13, 14359, 15, 3099, 12, NONE);

  /* Sin[39*Pi/128] + Cos[39*Pi/128] = 1.3933930045694290 */
  /* Sin[39*Pi/128] - Cos[39*Pi/128] = 0.2417766217337384 */
  /* Cos[39*Pi/128]                  = 0.5758081914178453 */
  od_rotate_add_half(tc, tj, tch, 16143, 14, 2801, 13, 26683, 15, NONE);

  /* Sin[37*Pi/128] + Cos[37*Pi/128] = 1.4035780182072331 */
  /* Sin[37*Pi/128] - Cos[37*Pi/128] = 0.1731148370459795 */
  /* Cos[37*Pi/128]                  = 0.6152315905806268 */
  od_rotate_sub_half(ti, td, tih, 16261, 14, 4011, 14, 14255, 14, NONE);

  /* Sin[35*Pi/128] + Cos[35*Pi/128] = 1.4103816894602614 */
  /* Sin[35*Pi/128] - Cos[35*Pi/128] = 0.1040360035527078 */
  /* Cos[35*Pi/128]                  = 0.6531728429537768 */
  od_rotate_add_half(te, th, teh, 32679, 15, 4821, 15, 30269, 15, NONE);

  /* Sin[33*Pi/128] + Cos[33*Pi/128] = 1.4137876276885337 */
  /* Sin[33*Pi/128] - Cos[33*Pi/128] = 0.0347065382144002 */
  /* Cos[33*Pi/128]                  = 0.6895405447370668 */
  od_rotate_sub_half(tg, tf, tgh, 16379, 14, 201, 12, 15977, 14, NONE);

  /* Stage 1 */

  od_butterfly_add(t0, &t0h, tf);
  od_butterfly_sub(tv, &tvh, tg);
  od_butterfly_add(th, &thh, tu);
  od_butterfly_sub(te, &teh, t1);

  od_butterfly_add(t2, &t2h, td);
  od_butterfly_sub(tt, &tth, ti);
  od_butterfly_add(tj, &tjh, ts);
  od_butterfly_sub(tc, &tch, t3);

  od_butterfly_add(t4, &t4h, tb);
  od_butterfly_sub(tr, &trh, tk);
  od_butterfly_add(tl, &tlh, tq);
  od_butterfly_sub(ta, &tah, t5);

  od_butterfly_add(t6, &t6h, t9);
  od_butterfly_sub(tp, &tph, tm);
  od_butterfly_add(tn, &tnh, to);
  od_butterfly_sub(t8, &t8h, t7);

  /* Stage 2 */

  od_butterfly_sub_asym(t0, t0h, t7);
  od_butterfly_add_asym(tv, tvh, to);
  od_butterfly_sub_asym(tp, tph, tu);
  od_butterfly_add_asym(t6, t6h, t1);

  od_butterfly_sub_asym(t2, t2h, t5);
  od_butterfly_add_asym(tt, tth, tq);
  od_butterfly_sub_asym(tr, trh, ts);
  od_butterfly_add_asym(t4, t4h, t3);

  od_butterfly_add_asym(t8, t8h, tg);
  od_butterfly_sub_asym(te, teh, tm);
  od_butterfly_add_asym(tn, tnh, tf);
  od_butterfly_sub_asym(th, thh, t9);

  od_butterfly_add_asym(ta, tah, ti);
  od_butterfly_sub_asym(tc, tch, tk);
  od_butterfly_add_asym(tl, tlh, td);
  od_butterfly_sub_asym(tj, tjh, tb);

  /* Stage 3 */

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_sub(tf, tg, 17911, 14, 14699, 14, 803, 13, NONE);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_add(th, te, 10217, 13, 5461, 13, 1189, 12, NONE);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_add(ti, td, 5543, 12, 3363, 13, 7723, 14, NONE);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_sub(tc, tj, 11529, 13, 2271, 14, 5197, 13, NONE);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_neg(tb, tk, 11529, 13, 2271, 14, 5197, 13);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_neg(ta, tl, 5543, 12, 3363, 13, 7723, 14);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_neg(t9, tm, 10217, 13, 5461, 13, 1189, 12);

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_neg(t8, tn, 17911, 14, 14699, 14, 803, 13);

  /* Stage 4 */

  od_butterfly_sub(t3, &t3h, t0);
  od_butterfly_add(ts, &tsh, tv);
  od_butterfly_sub(tu, &tuh, tt);
  od_butterfly_add(t1, &t1h, t2);

  od_butterfly_add(to, NULL, t4);
  od_butterfly_sub(tq, NULL, t6);
  od_butterfly_add(t7, NULL, tr);
  od_butterfly_sub(t5, NULL, tp);

  od_butterfly_sub(tb, &tbh, t8);
  od_butterfly_add(tk, &tkh, tn);
  od_butterfly_sub(tm, &tmh, tl);
  od_butterfly_add(t9, &t9h, ta);

  od_butterfly_sub(tf, &tfh, tc);
  od_butterfly_add(tg, &tgh, tj);
  od_butterfly_sub(ti, &tih, th);
  od_butterfly_add(td, &tdh, te);

  /* Stage 5 */

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_add(to, t7, 301, 8, 1609, 11, 6393, 15, NONE);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(tp, t6, 11363, 13, 9041, 15, 4551, 13, NONE);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_neg(t5, tq, 5681, 12, 9041, 15, 4551, 13);

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(t4, tr, 9633, 13, 12873, 14, 6393, 15);

  /* Stage 6 */

  od_butterfly_add_asym(t1, t1h, t0);
  od_butterfly_sub_asym(tu, tuh, tv);
  od_butterfly_sub_asym(ts, tsh, t2);
  od_butterfly_sub_asym(t3, t3h, tt);

  od_butterfly_add_asym(t5, od_rshift1(*t5), t4);
  od_butterfly_sub_asym(tq, od_rshift1(*tq), tr);
  od_butterfly_add_asym(t7, od_rshift1(*t7), t6);
  od_butterfly_sub_asym(to, od_rshift1(*to), tp);

  od_butterfly_add_asym(t9, t9h, t8);
  od_butterfly_sub_asym(tm, tmh, tn);
  od_butterfly_sub_asym(tk, tkh, ta);
  od_butterfly_sub_asym(tb, tbh, tl);

  od_butterfly_add_asym(ti, tih, tc);
  od_butterfly_add_asym(td, tdh, tj);
  od_butterfly_add_asym(tf, tfh, te);
  od_butterfly_sub_asym(tg, tgh, th);

  /* Stage 7 */

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(t2, tt, 669, 9, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(ts, t3, 669, 9, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(ta, tl, 669, 9, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tk, tb, 669, 9, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tc, tj, 669, 9, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(ti, td, 669, 9, 8867, 14, 3135, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tu, t1, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tq, t5, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_sub(tp, t6, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tm, t9, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(te, th, 5793, 12, 5793, 13);
}

/**
 * 32-point asymmetric Type-IV iDST
 */
static INLINE void od_idst_32_asym(od_coeff *t0, od_coeff *tg,
                                   od_coeff *t8, od_coeff *to,
                                   od_coeff *t4, od_coeff *tk,
                                   od_coeff *tc, od_coeff *ts,
                                   od_coeff *t2, od_coeff *ti,
                                   od_coeff *ta, od_coeff *tq,
                                   od_coeff *t6, od_coeff *tm,
                                   od_coeff *te, od_coeff *tu,
                                   od_coeff *t1, od_coeff *th,
                                   od_coeff *t9, od_coeff *tp,
                                   od_coeff *t5, od_coeff *tl,
                                   od_coeff *td, od_coeff *tt,
                                   od_coeff *t3, od_coeff *tj,
                                   od_coeff *tb, od_coeff *tr,
                                   od_coeff *t7, od_coeff *tn,
                                   od_coeff *tf, od_coeff *tv) {
  od_coeff t0h;
  od_coeff t1h;
  od_coeff t2h;
  od_coeff t3h;
  od_coeff t4h;
  od_coeff t6h;
  od_coeff t8h;
  od_coeff t9h;
  od_coeff tah;
  od_coeff tbh;
  od_coeff tch;
  od_coeff tdh;
  od_coeff teh;
  od_coeff tfh;
  od_coeff tgh;
  od_coeff thh;
  od_coeff tih;
  od_coeff tjh;
  od_coeff tkh;
  od_coeff tlh;
  od_coeff tmh;
  od_coeff tnh;
  od_coeff tph;
  od_coeff trh;
  od_coeff tsh;
  od_coeff tth;
  od_coeff tuh;
  od_coeff tvh;

  /* Stage 7 */

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(te, th, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tm, t9, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_sub(tp, t6, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tq, t5, 5793, 12, 5793, 13);

  /* Sin[Pi/4] + Cos[Pi/4] = 1.4142135623730951 */
  /* Cos[Pi/4]             = 0.7071067811865475 */
  od_rotate_pi4_add(tu, t1, 5793, 12, 5793, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(td, ti, 669, 9, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tc, tj, 669, 9, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(tk, tb, 669, 9, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(tl, ta, 669, 9, 8867, 14, 3135, 13);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_add(ts, t3, 669, 9, 8867, 14, 3135, 13, NONE);

  /* Sin[3*Pi/8] + Cos[3*Pi/8] = 1.3065629648763766 */
  /* Sin[3*Pi/8] - Cos[3*Pi/8] = 0.5411961001461969 */
  /* Cos[3*Pi/8]               = 0.3826834323650898 */
  od_rotate_neg(tt, t2, 669, 9, 8867, 14, 3135, 13);

  /* Stage 6 */

  od_butterfly_sub(tg, &tgh, th);
  od_butterfly_add(tf, &tfh, te);
  od_butterfly_add(td, &tdh, tj);
  od_butterfly_add(ti, &tih, tc);

  od_butterfly_sub(tb, &tbh, tl);
  od_butterfly_sub(tk, &tkh, ta);
  od_butterfly_sub(tm, &tmh, tn);
  od_butterfly_add(t9, &t9h, t8);

  od_butterfly_sub(to, NULL, tp);
  od_butterfly_add(t7, NULL, t6);
  od_butterfly_sub(tq, NULL, tr);
  od_butterfly_add(t5, NULL, t4);

  od_butterfly_sub(t3, &t3h, tt);
  od_butterfly_sub(ts, &tsh, t2);
  od_butterfly_sub(tu, &tuh, tv);
  od_butterfly_add(t1, &t1h, t0);

  /* Stage 5 */

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_neg(tr, t4, 9633, 13, 12873, 14, 6393, 15);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_neg(tq, t5, 5681, 12, 9041, 15, 4551, 13);

  /* Sin[5*Pi/16] + Cos[5*Pi/16] = 1.3870398453221475 */
  /* Sin[5*Pi/16] - Cos[5*Pi/16] = 0.2758993792829431 */
  /* Cos[5*Pi/16]                = 0.5555702330196022 */
  od_rotate_add(tp, t6, 11363, 13, 9041, 15, 4551, 13, NONE);

  /* Sin[7*Pi/16] + Cos[7*Pi/16] = 1.1758756024193586 */
  /* Sin[7*Pi/16] - Cos[7*Pi/16] = 0.7856949583871022 */
  /* Cos[7*Pi/16]                = 0.1950903220161283 */
  od_rotate_add(to, t7, 301, 8, 1609, 11, 6393, 15, NONE);

  /* Stage 4 */

  od_butterfly_add_asym(td, tdh, te);
  od_butterfly_sub_asym(ti, tih, th);
  od_butterfly_add_asym(tg, tgh, tj);
  od_butterfly_sub_asym(tf, tfh, tc);

  od_butterfly_add_asym(t9, t9h, ta);
  od_butterfly_sub_asym(tm, tmh, tl);
  od_butterfly_add_asym(tk, tkh, tn);
  od_butterfly_sub_asym(tb, tbh, t8);

  od_butterfly_sub_asym(t5, od_rshift1(*t5), tp);
  od_butterfly_add_asym(t7, od_rshift1(*t7), tr);
  od_butterfly_sub_asym(tq, od_rshift1(*tq), t6);
  od_butterfly_add_asym(to, od_rshift1(*to), t4);

  od_butterfly_add_asym(t1, t1h, t2);
  od_butterfly_sub_asym(tu, tuh, tt);
  od_butterfly_add_asym(ts, tsh, tv);
  od_butterfly_sub_asym(t3, t3h, t0);

  /* Stage 3 */

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_neg(tn, t8, 17911, 14, 14699, 14, 803, 13);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_neg(tm, t9, 10217, 13, 5461, 13, 1189, 12);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_neg(tl, ta, 5543, 12, 3363, 13, 7723, 14);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_neg(tk, tb, 11529, 13, 2271, 14, 5197, 13);

  /* Sin[9*Pi/32] + Cos[9*Pi/32] = 1.4074037375263826 */
  /* Sin[9*Pi/32] - Cos[9*Pi/32] = 0.1386171691990915 */
  /* Cos[9*Pi/32]                = 0.6343932841636455 */
  od_rotate_sub(tc, tj, 11529, 13, 2271, 14, 5197, 13, NONE);

  /* Sin[11*Pi/32] + Cos[11*Pi/32] = 1.3533180011743526 */
  /* Sin[11*Pi/32] - Cos[11*Pi/32] = 0.4105245275223574 */
  /* Cos[11*Pi/32]                 = 0.4713967368259976 */
  od_rotate_add(ti, td, 5543, 12, 3363, 13, 7723, 14, NONE);

  /* Sin[13*Pi/32] + Cos[13*Pi/32] = 1.2472250129866712 */
  /* Sin[13*Pi/32] - Cos[13*Pi/32] = 0.6666556584777465 */
  /* Cos[13*Pi/32]                 = 0.2902846772544623 */
  od_rotate_add(th, te, 10217, 13, 5461, 13, 1189, 12, NONE);

  /* Sin[15*Pi/32] + Cos[15*Pi/32] = 1.0932018670017576 */
  /* Sin[15*Pi/32] - Cos[15*Pi/32] = 0.8971675863426363 */
  /* Cos[15*Pi/32]                 = 0.0980171403295606 */
  od_rotate_sub(tf, tg, 17911, 14, 14699, 14, 803, 13, NONE);

  /* Stage 2 */

  od_butterfly_sub(tj, &tjh, tb);
  od_butterfly_add(tl, &tlh, td);
  od_butterfly_sub(tc, &tch, tk);
  od_butterfly_add(ta, &tah, ti);

  od_butterfly_sub(th, &thh, t9);
  od_butterfly_add(tn, &tnh, tf);
  od_butterfly_sub(te, &teh, tm);
  od_butterfly_add(t8, &t8h, tg);

  od_butterfly_add(t4, &t4h, t3);
  od_butterfly_sub(tr, &trh, ts);
  od_butterfly_add(tt, &tth, tq);
  od_butterfly_sub(t2, &t2h, t5);

  od_butterfly_add(t6, &t6h, t1);
  od_butterfly_sub(tp, &tph, tu);
  od_butterfly_add(tv, &tvh, to);
  od_butterfly_sub(t0, &t0h, t7);

  /* Stage 1 */

  od_butterfly_sub_asym(t8, t8h, t7);
  od_butterfly_add_asym(tn, tnh, to);
  od_butterfly_sub_asym(tp, tph, tm);
  od_butterfly_add_asym(t6, t6h, t9);

  od_butterfly_sub_asym(ta, tah, t5);
  od_butterfly_add_asym(tl, tlh, tq);
  od_butterfly_sub_asym(tr, trh, tk);
  od_butterfly_add_asym(t4, t4h, tb);

  od_butterfly_sub_asym(tc, tch, t3);
  od_butterfly_add_asym(tj, tjh, ts);
  od_butterfly_sub_asym(tt, tth, ti);
  od_butterfly_add_asym(t2, t2h, td);

  od_butterfly_sub_asym(te, teh, t1);
  od_butterfly_add_asym(th, thh, tu);
  od_butterfly_sub_asym(tv, tvh, tg);
  od_butterfly_add_asym(t0, t0h, tf);

  /* Stage 0 */

  /* Sin[33*Pi/128] + Cos[33*Pi/128] = 1.4137876276885337 */
  /* Sin[33*Pi/128] - Cos[33*Pi/128] = 0.0347065382144002 */
  /* Cos[33*Pi/128]                  = 0.6895405447370668 */
  od_rotate_sub(tg, tf, 16379, 14, 201, 12, 15977, 14, SHIFT);

  /* Sin[35*Pi/128] + Cos[35*Pi/128] = 1.4103816894602614 */
  /* Sin[35*Pi/128] - Cos[35*Pi/128] = 0.1040360035527078 */
  /* Cos[35*Pi/128]                  = 0.6531728429537768 */
  od_rotate_add(te, th, 32679, 15, 4821, 15, 30269, 15, SHIFT);

  /* Sin[37*Pi/128] + Cos[37*Pi/128] = 1.4035780182072331 */
  /* Sin[37*Pi/128] - Cos[37*Pi/128] = 0.1731148370459795 */
  /* Cos[37*Pi/128]                  = 0.6152315905806268 */
  od_rotate_sub(ti, td, 16261, 14, 4011, 14, 14255, 14, SHIFT);

  /* Sin[39*Pi/128] + Cos[39*Pi/128] = 1.3933930045694290 */
  /* Sin[39*Pi/128] - Cos[39*Pi/128] = 0.2417766217337384 */
  /* Cos[39*Pi/128]                  = 0.5758081914178453 */
  od_rotate_add(tc, tj, 16143, 14, 2801, 13, 26683, 15, SHIFT);

  /* Sin[41*Pi/128] + Cos[41*Pi/128] = 1.3798511851368043 */
  /* Sin[41*Pi/128] - Cos[41*Pi/128] = 0.3098559453626100 */
  /* Cos[41*Pi/128]                  = 0.5349976198870972 */
  od_rotate_sub(tk, tb, 7993, 13, 14359, 15, 3099, 12, SHIFT);

  /* Sin[43*Pi/128] + Cos[43*Pi/128] = 1.3629851833384956 */
  /* Sin[43*Pi/128] - Cos[43*Pi/128] = 0.3771887988789274 */
  /* Cos[43*Pi/128]                  = 0.4928981922297840 */
  od_rotate_add(ta, tl, 31581, 15, 17479, 15, 22841, 15, SHIFT);

  /* Sin[45*Pi/128] + Cos[45*Pi/128] = 1.3428356308501219 */
  /* Sin[45*Pi/128] - Cos[45*Pi/128] = 0.4436129715409088 */
  /* Cos[45*Pi/128]                  = 0.4496113296546065 */
  od_rotate_sub(tm, t9, 15557, 14, 20557, 15, 20835, 15, SHIFT);

  /* Sin[47*Pi/128] + Cos[47*Pi/128] = 1.3194510697085207 */
  /* Sin[47*Pi/128] - Cos[47*Pi/128] = 0.5089684416985408 */
  /* Cos[47*Pi/128]                  = 0.4052413140049899 */
  od_rotate_add(t8, tn, 7643, 13, 11793, 14, 18779, 15, SHIFT);

  /* Sin[49*Pi/128] + Cos[49*Pi/128] = 1.2928878353697270 */
  /* Sin[49*Pi/128] - Cos[49*Pi/128] = 0.5730977622997508 */
  /* Cos[49*Pi/128]                  = 0.3598950365349881 */
  od_rotate_sub(to, t7, 29957, 15, 13279, 14, 8339, 14, SHIFT);

  /* Sin[51*Pi/128] + Cos[51*Pi/128] = 1.2632099209919283 */
  /* Sin[51*Pi/128] - Cos[51*Pi/128] = 0.6358464401941452 */
  /* Cos[51*Pi/128]                  = 0.3136817403988915 */
  od_rotate_add(t6, tp, 29269, 15, 14733, 14, 1817, 12, SHIFT);

  /* Sin[53*Pi/128] + Cos[53*Pi/128] = 1.2304888232703382 */
  /* Sin[53*Pi/128] - Cos[53*Pi/128] = 0.6970633083205415 */
  /* Cos[53*Pi/128]                  = 0.2667127574748984 */
  od_rotate_sub(tq, t5, 28511, 15, 32303, 15, 1545, 12, SHIFT);

  /* Sin[55*Pi/128] + Cos[55*Pi/128] = 1.1948033701953984 */
  /* Sin[55*Pi/128] - Cos[55*Pi/128] = 0.7566008898816587 */
  /* Cos[55*Pi/128]                  = 0.2191012401568698 */
  od_rotate_add(t4, tr, 6921, 13, 17531, 14, 10153, 15, SHIFT);

  /* Sin[57*Pi/128] + Cos[57*Pi/128] = 1.1562395311492424 */
  /* Sin[57*Pi/128] - Cos[57*Pi/128] = 0.8143157536286401 */
  /* Cos[57*Pi/128]                  = 0.1709618887603012 */
  od_rotate_sub(ts, t3, 26791, 15, 4717, 12, 7923, 15, SHIFT);

  /* Sin[59*Pi/128] + Cos[59*Pi/128] = 1.1148902097979262 */
  /* Sin[59*Pi/128] - Cos[59*Pi/128] = 0.8700688593994937 */
  /* Cos[59*Pi/128]                  = 0.1224106751992162 */
  od_rotate_add(t2, tt, 25833, 15, 315, 8, 5673, 15, SHIFT);

  /* Sin[61*Pi/128] + Cos[61*Pi/128] = 1.0708550202783576 */
  /* Sin[61*Pi/128] - Cos[61*Pi/128] = 0.9237258930790228 */
  /* Cos[61*Pi/128]                  = 0.0735645635996674 */
  od_rotate_sub(tu, t1, 6203, 13, 21403, 14, 3409, 15, SHIFT);

  /* Sin[63*Pi/128] + Cos[63*Pi/128] = 1.0242400472191164 */
  /* Sin[63*Pi/128] - Cos[63*Pi/128] = 0.9751575901732919 */
  /* Cos[63*Pi/128]                  = 0.0245412285229123 */
  od_rotate_add(t0, tv, 5933, 13, 22595, 14, 1137, 15, SHIFT);
}

/* --- 64-point Transforms --- */

/**
 * 64-point orthonormal Type-II fDCT
 */
static INLINE void od_fdct_64(od_coeff *u0, od_coeff *u1,
                              od_coeff *u2, od_coeff *u3,
                              od_coeff *u4, od_coeff *u5,
                              od_coeff *u6, od_coeff *u7,
                              od_coeff *u8, od_coeff *u9,
                              od_coeff *ua, od_coeff *ub,
                              od_coeff *uc, od_coeff *ud,
                              od_coeff *ue, od_coeff *uf,
                              od_coeff *ug, od_coeff *uh,
                              od_coeff *ui, od_coeff *uj,
                              od_coeff *uk, od_coeff *ul,
                              od_coeff *um, od_coeff *un,
                              od_coeff *uo, od_coeff *up,
                              od_coeff *uq, od_coeff *ur,
                              od_coeff *us, od_coeff *ut,
                              od_coeff *uu, od_coeff *uv,
                              od_coeff *uw, od_coeff *ux,
                              od_coeff *uy, od_coeff *uz,
                              od_coeff *uA, od_coeff *uB,
                              od_coeff *uC, od_coeff *uD,
                              od_coeff *uE, od_coeff *uF,
                              od_coeff *uG, od_coeff *uH,
                              od_coeff *uI, od_coeff *uJ,
                              od_coeff *uK, od_coeff *uL,
                              od_coeff *uM, od_coeff *uN,
                              od_coeff *uO, od_coeff *uP,
                              od_coeff *uQ, od_coeff *uR,
                              od_coeff *uS, od_coeff *uT,
                              od_coeff *uU, od_coeff *uV,
                              od_coeff *uW, od_coeff *uX,
                              od_coeff *uY, od_coeff *uZ,
                              od_coeff *u_, od_coeff *u) {
  od_coeff u1h;
  od_coeff u3h;
  od_coeff u5h;
  od_coeff u7h;
  od_coeff u9h;
  od_coeff ubh;
  od_coeff udh;
  od_coeff ufh;
  od_coeff uhh;
  od_coeff ujh;
  od_coeff ulh;
  od_coeff unh;
  od_coeff uph;
  od_coeff urh;
  od_coeff uth;
  od_coeff uvh;
  od_coeff uxh;
  od_coeff uzh;
  od_coeff uBh;
  od_coeff uDh;
  od_coeff uFh;
  od_coeff uHh;
  od_coeff uJh;
  od_coeff uLh;
  od_coeff uNh;
  od_coeff uPh;
  od_coeff uRh;
  od_coeff uTh;
  od_coeff uVh;
  od_coeff uXh;
  od_coeff uZh;
  od_coeff uh_;

  /* +/- Butterflies with asymmetric output. */
  od_butterfly_neg(u0, u , &uh_);
  od_butterfly_add(u1, &u1h, u_);
  od_butterfly_neg(u2, uZ, &uZh);
  od_butterfly_add(u3, &u3h, uY);
  od_butterfly_neg(u4, uX, &uXh);
  od_butterfly_add(u5, &u5h, uW);
  od_butterfly_neg(u6, uV, &uVh);
  od_butterfly_add(u7, &u7h, uU);
  od_butterfly_neg(u8, uT, &uTh);
  od_butterfly_add(u9, &u9h, uS);
  od_butterfly_neg(ua, uR, &uRh);
  od_butterfly_add(ub, &ubh, uQ);
  od_butterfly_neg(uc, uP, &uPh);
  od_butterfly_add(ud, &udh, uO);
  od_butterfly_neg(ue, uN, &uNh);
  od_butterfly_add(uf, &ufh, uM);
  od_butterfly_neg(ug, uL, &uLh);
  od_butterfly_add(uh, &uhh, uK);
  od_butterfly_neg(ui, uJ, &uJh);
  od_butterfly_add(uj, &ujh, uI);
  od_butterfly_neg(uk, uH, &uHh);
  od_butterfly_add(ul, &ulh, uG);
  od_butterfly_neg(um, uF, &uFh);
  od_butterfly_add(un, &unh, uE);
  od_butterfly_neg(uo, uD, &uDh);
  od_butterfly_add(up, &uph, uC);
  od_butterfly_neg(uq, uB, &uBh);
  od_butterfly_add(ur, &urh, uA);
  od_butterfly_neg(us, uz, &uzh);
  od_butterfly_add(ut, &uth, uy);
  od_butterfly_neg(uu, ux, &uxh);
  od_butterfly_add(uv, &uvh, uw);

  /* Embedded 16-point transforms with asymmetric input. */
  od_fdct_32_asym(
   u0, u1, u1h, u2, u3, u3h, u4, u5, u5h, u6, u7, u7h,
   u8, u9, u9h, ua, ub, ubh, uc, ud, udh, ue, uf, ufh,
   ug, uh, uhh, ui, uj, ujh, uk, ul, ulh, um, un, unh,
   uo, up, uph, uq, ur, urh, us, ut, uth, uu, uv, uvh);

  od_fdst_32_asym(
   u , uh_, u_, uZ, uZh, uY, uX, uXh, uW, uV, uVh, uU,
   uT, uTh, uS, uR, uRh, uQ, uP, uPh, uO, uN, uNh, uM,
   uL, uLh, uK, uJ, uJh, uI, uH, uHh, uG, uF, uFh, uE,
   uD, uDh, uC, uB, uBh, uA, uz, uzh, uy, ux, uxh, uw);
}

/**
 * 64-point orthonormal Type-II iDCT
 */
static INLINE void od_idct_64(od_coeff *u0, od_coeff *uw,
                              od_coeff *ug, od_coeff *uM,
                              od_coeff *u8, od_coeff *uE,
                              od_coeff *uo, od_coeff *uU,
                              od_coeff *u4, od_coeff *uA,
                              od_coeff *uk, od_coeff *uQ,
                              od_coeff *uc, od_coeff *uI,
                              od_coeff *us, od_coeff *uY,
                              od_coeff *u2, od_coeff *uy,
                              od_coeff *ui, od_coeff *uO,
                              od_coeff *ua, od_coeff *uG,
                              od_coeff *uq, od_coeff *uW,
                              od_coeff *u6, od_coeff *uC,
                              od_coeff *um, od_coeff *uS,
                              od_coeff *ue, od_coeff *uK,
                              od_coeff *uu, od_coeff *u_,
                              od_coeff *u1, od_coeff *ux,
                              od_coeff *uh, od_coeff *uN,
                              od_coeff *u9, od_coeff *uF,
                              od_coeff *up, od_coeff *uV,
                              od_coeff *u5, od_coeff *uB,
                              od_coeff *ul, od_coeff *uR,
                              od_coeff *ud, od_coeff *uJ,
                              od_coeff *ut, od_coeff *uZ,
                              od_coeff *u3, od_coeff *uz,
                              od_coeff *uj, od_coeff *uP,
                              od_coeff *ub, od_coeff *uH,
                              od_coeff *ur, od_coeff *uX,
                              od_coeff *u7, od_coeff *uD,
                              od_coeff *un, od_coeff *uT,
                              od_coeff *uf, od_coeff *uL,
                              od_coeff *uv, od_coeff *u ) {
  od_coeff u1h;
  od_coeff u3h;
  od_coeff u5h;
  od_coeff u7h;
  od_coeff u9h;
  od_coeff ubh;
  od_coeff udh;
  od_coeff ufh;
  od_coeff uhh;
  od_coeff ujh;
  od_coeff ulh;
  od_coeff unh;
  od_coeff uph;
  od_coeff urh;
  od_coeff uth;
  od_coeff uvh;

  /* Embedded 32-point transforms with asymmetric output. */
  od_idst_32_asym(
   u , uL, uT, uD, uX, uH, uP, uz, uZ, uJ, uR, uB, uV, uF, uN, ux,
   u_, uK, uS, uC, uW, uG, uO, uy, uY, uI, uQ, uA, uU, uE, uM, uw);
  od_idct_32_asym(
   u0, ug, u8, uo, u4, uk, uc, us, u2, ui, ua, uq, u6, um, ue, uu,
   u1, &u1h, uh, &uhh, u9, &u9h, up, &uph,
   u5, &u5h, ul, &ulh, ud, &udh, ut, &uth,
   u3, &u3h, uj, &ujh, ub, &ubh, ur, &urh,
   u7, &u7h, un, &unh, uf, &ufh, uv, &uvh);

  /* +/- Butterflies with asymmetric input. */
  od_butterfly_neg_asym(u0, u , od_rshift1(*u));
  od_butterfly_add_asym(u1, u1h, u_);
  od_butterfly_neg_asym(u2, uZ, od_rshift1(*uZ));
  od_butterfly_add_asym(u3, u3h, uY);
  od_butterfly_neg_asym(u4, uX, od_rshift1(*uX));
  od_butterfly_add_asym(u5, u5h, uW);
  od_butterfly_neg_asym(u6, uV, od_rshift1(*uV));
  od_butterfly_add_asym(u7, u7h, uU);
  od_butterfly_neg_asym(u8, uT, od_rshift1(*uT));
  od_butterfly_add_asym(u9, u9h, uS);
  od_butterfly_neg_asym(ua, uR, od_rshift1(*uR));
  od_butterfly_add_asym(ub, ubh, uQ);
  od_butterfly_neg_asym(uc, uP, od_rshift1(*uP));
  od_butterfly_add_asym(ud, udh, uO);
  od_butterfly_neg_asym(ue, uN, od_rshift1(*uN));
  od_butterfly_add_asym(uf, ufh, uM);
  od_butterfly_neg_asym(ug, uL, od_rshift1(*uL));
  od_butterfly_add_asym(uh, uhh, uK);
  od_butterfly_neg_asym(ui, uJ, od_rshift1(*uJ));
  od_butterfly_add_asym(uj, ujh, uI);
  od_butterfly_neg_asym(uk, uH, od_rshift1(*uH));
  od_butterfly_add_asym(ul, ulh, uG);
  od_butterfly_neg_asym(um, uF, od_rshift1(*uF));
  od_butterfly_add_asym(un, unh, uE);
  od_butterfly_neg_asym(uo, uD, od_rshift1(*uD));
  od_butterfly_add_asym(up, uph, uC);
  od_butterfly_neg_asym(uq, uB, od_rshift1(*uB));
  od_butterfly_add_asym(ur, urh, uA);
  od_butterfly_neg_asym(us, uz, od_rshift1(*uz));
  od_butterfly_add_asym(ut, uth, uy);
  od_butterfly_neg_asym(uu, ux, od_rshift1(*ux));
  od_butterfly_add_asym(uv, uvh, uw);
}

#endif
