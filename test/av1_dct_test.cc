/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
*/

#include <math.h>
#include <stdlib.h>
#include <new>

#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/acm_random.h"
#include "test/util.h"
#include "./aom_config.h"
#include "aom_ports/msvc.h"

#undef CONFIG_COEFFICIENT_RANGE_CHECKING
#define CONFIG_COEFFICIENT_RANGE_CHECKING 1
#include "av1/encoder/dct.c"

using libaom_test::ACMRandom;

namespace {
void reference_dct_1d(const double *in, double *out, int size) {
  const double kInvSqrt2 = 0.707106781186547524400844362104;
  for (int k = 0; k < size; ++k) {
    out[k] = 0;
    for (int n = 0; n < size; ++n) {
      out[k] += in[n] * cos(PI * (2 * n + 1) * k / (2 * size));
    }
    if (k == 0) out[k] = out[k] * kInvSqrt2;
  }
}

typedef void (*FdctFuncRef)(const double *in, double *out, int size);
typedef void (*IdctFuncRef)(const double *in, double *out, int size);
typedef void (*FdctFunc)(const tran_low_t *in, tran_low_t *out);
typedef void (*IdctFunc)(const tran_low_t *in, tran_low_t *out);

class TransTestBase {
 public:
  virtual ~TransTestBase() {}

 protected:
  void RunFwdAccuracyCheck() {
    tran_low_t *input = new tran_low_t[txfm_size_];
    tran_low_t *output = new tran_low_t[txfm_size_];
    double *ref_input = new double[txfm_size_];
    double *ref_output = new double[txfm_size_];

    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 5000;
    for (int ti = 0; ti < count_test_block; ++ti) {
      for (int ni = 0; ni < txfm_size_; ++ni) {
        input[ni] = rnd.Rand8() - rnd.Rand8();
        ref_input[ni] = static_cast<double>(input[ni]);
      }

      fwd_txfm_(input, output);
      fwd_txfm_ref_(ref_input, ref_output, txfm_size_);

      for (int ni = 0; ni < txfm_size_; ++ni) {
        EXPECT_LE(
            abs(output[ni] - static_cast<tran_low_t>(round(ref_output[ni]))),
            max_error_);
      }
    }

    delete[] input;
    delete[] output;
    delete[] ref_input;
    delete[] ref_output;
  }

  void RunFwdExtremeCheck() {
    tran_low_t *input = new tran_low_t[txfm_size_];
    tran_low_t *temp_out = new tran_low_t[txfm_size_];
    tran_low_t *output = new tran_low_t[txfm_size_];
    double *ref_input = new double[txfm_size_];
    double *tmp_output = new double[txfm_size_];
    double *ref_output = new double[txfm_size_];

    ACMRandom rnd(ACMRandom::DeterministicSeed());
    const int count_test_block = 5;
    const tran_low_t mask = 0xFF;
    for (int ti = 0; ti < count_test_block; ++ti) {
      for (int ni = 0; ni < txfm_size_; ++ni) {
        input[ni] = rnd.Rand8() % 2 ? mask : -mask;
        ref_input[ni] = static_cast<double>(input[ni]);
      }

      if (ti == 0) {
        for (int j = 0; j < txfm_size_; ++j) {
          input[j] = mask;
          ref_input[j] = static_cast<double>(input[j]);
        }
      } else if (ti == 1) {
        for (int j = 0; j < txfm_size_; ++j) {
          input[j] = -mask;
          ref_input[j] = static_cast<double>(input[j]);
        }
      }

      scale_input(input);
      fwd_txfm_(input, temp_out);
      simulate_col_txfm(temp_out);
      fwd_txfm_(temp_out, output);
      scale_output(output);

      scale_input_ref(ref_input);
      fwd_txfm_ref_(ref_input, tmp_output, txfm_size_);
      simulate_col_txfm_ref(tmp_output);
      fwd_txfm_ref_(tmp_output, ref_output, txfm_size_);
      scale_output_ref(ref_output);

      for (int ni = 0; ni < txfm_size_; ++ni) {
        EXPECT_LE(
            abs(output[ni] - static_cast<tran_low_t>(round(ref_output[ni]))),
            max_error_);
      }
    }

    delete[] input;
    delete[] temp_out;
    delete[] output;
    delete[] ref_input;
    delete[] tmp_output;
    delete[] ref_output;
  }

  double max_error_;
  int txfm_size_;
  FdctFunc fwd_txfm_;
  FdctFuncRef fwd_txfm_ref_;

  void simulate_col_txfm(tran_low_t *out) {
    tran_low_t value = ROUND_POWER_OF_TWO(out[0], 2);
    for (int i = 0; i < txfm_size_; i++) {
      out[i] = value;
    }
  }
  void simulate_col_txfm_ref(double *out) {
    double value = out[0] / 4.0;
    for (int i = 0; i < txfm_size_; i++) {
      out[i] = value;
    }
  }
  void scale_input(tran_low_t *in) {
    for (int i = 0; i < txfm_size_; i++) {
      in[i] *= 4;
    }
  }
  void scale_input_ref(double *in) {
    for (int i = 0; i < txfm_size_; i++) {
      in[i] *= 4;
    }
  }
  void scale_output(tran_low_t *out) {
    for (int i = 0; i < txfm_size_; i++) {
      out[i] /= 4;
    }
  }
  void scale_output_ref(double *out) {
    for (int i = 0; i < txfm_size_; i++) {
      out[i] /= 4;
    }
  }
};

typedef std::tr1::tuple<FdctFunc, FdctFuncRef, int, int> FdctParam;
class AV1FwdTxfm : public TransTestBase,
                   public ::testing::TestWithParam<FdctParam> {
 public:
  virtual void SetUp() {
    fwd_txfm_ = GET_PARAM(0);
    fwd_txfm_ref_ = GET_PARAM(1);
    txfm_size_ = GET_PARAM(2);
    max_error_ = GET_PARAM(3);
  }
  virtual void TearDown() {}
};

TEST_P(AV1FwdTxfm, FwdAccuracyCheck) { RunFwdAccuracyCheck(); }
TEST_P(AV1FwdTxfm, FwdExtremeCheck) { RunFwdExtremeCheck(); }

INSTANTIATE_TEST_CASE_P(
    C, AV1FwdTxfm,
    ::testing::Values(FdctParam(&fdct4, &reference_dct_1d, 4, 1),
                      FdctParam(&fdct8, &reference_dct_1d, 8, 1),
                      FdctParam(&fdct16, &reference_dct_1d, 16, 2),
                      FdctParam(&fdct32, &reference_dct_1d, 32, 3)));
}  // namespace
