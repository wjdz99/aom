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

#ifndef SIMD_TEST_H_
#define SIMD_TEST_H_

#include "third_party/googletest/src/include/gtest/gtest.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"

namespace simd_test {

template <typename param_signature, typename simd_signature,
          typename ref_signature>
class TestIntrinsic : public ::testing::TestWithParam<param_signature> {
 public:
  virtual ~TestIntrinsic() {}
  virtual void SetUp() {
    simd = std::tr1::get<0>(this->GetParam());
    ref_simd = std::tr1::get<1>(this->GetParam());
    mask = std::tr1::get<2>(this->GetParam());
    maskwidth = std::tr1::get<3>(this->GetParam());
    name = std::tr1::get<4>(this->GetParam());
  }

  virtual void TearDown() { libaom_test::ClearSystemState(); }

 protected:
  simd_signature simd;
  ref_signature ref_simd;
  uint32_t mask, maskwidth;
  const char *name;
};

typedef uint32_t c_uint32_t;
typedef uint64_t c_uint64_t;
typedef int64_t c_int64_t;

// Create one typedef for each function signature
#define TYPEDEF_SIMD1(name, r, a)                                             \
  typedef TestIntrinsic<std::tr1::tuple<r (*)(a), c_##r (*)(c_##a), uint32_t, \
                                        uint32_t, const char *>,              \
                        r (*)(a), c_##r (*)(c_##a)>                           \
      name

#define TYPEDEF_SIMD2(name, r, a, b)                                          \
  typedef TestIntrinsic<std::tr1::tuple<r (*)(a, b), c_##r (*)(c_##a, c_##b), \
                                        uint32_t, uint32_t, const char *>,    \
                        r (*)(a, b), c_##r (*)(c_##a, c_##b)>                 \
      name

TYPEDEF_SIMD1(V64_V64, v64, v64);
TYPEDEF_SIMD1(V128_V64, v128, v64);
TYPEDEF_SIMD1(S64_V64, int64_t, v64);
TYPEDEF_SIMD1(V128_V128, v128, v128);
TYPEDEF_SIMD1(U64_V64, uint64_t, v64);
TYPEDEF_SIMD1(U64_V128, uint64_t, v128);
TYPEDEF_SIMD2(V64_V64V64, v64, v64, v64);
TYPEDEF_SIMD2(V128_V64V64, v128, v64, v64);
TYPEDEF_SIMD2(V64_V64U32, v64, v64, uint32_t);
TYPEDEF_SIMD2(U32_V64V64, uint32_t, v64, v64);
TYPEDEF_SIMD2(V128_V128V128, v128, v128, v128);
TYPEDEF_SIMD2(V128_V128U32, v128, v128, uint32_t);
TYPEDEF_SIMD2(U32_V128V128, uint32_t, v128, v128);

// Google Test allows up to 50 tests per case, so split the largest
typedef V64_V64 V64_V64_Part2;
typedef V128_V128 V128_V128_Part2;
typedef V128_V128 V128_V128_Part3;
typedef V64_V64V64 V64_V64V64_Part2;
typedef V128_V128V128 V128_V128V128_Part2;
}
#endif
