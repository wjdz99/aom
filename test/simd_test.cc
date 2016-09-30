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

#include <string>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "aom_dsp/aom_simd.h"
#include "aom_dsp/simd/v128_intrinsics_c.h"
#include "./simd_test_defs.h"

using libaom_test::ACMRandom;

namespace simd_test {

std::string print(uint8_t *a, int size) {
  std::string text = "0x";
  for (int i = 0; i < size; i++) {
    char buf[3];
    snprintf(buf, sizeof(buf), "%02x",
             a[!CONFIG_BIG_ENDIAN ? size - 1 - i : i]);
    text += buf;
  }

  return text;
}

void setmask(uint8_t *s, int size, uint32_t mask, uint32_t maskwidth) {
  switch (maskwidth) {
    case 0: break;
    case 8:
      for (int i = 0; i < size; i++) s[i] &= mask;
      break;
    case 16: {
      uint16_t *t = (uint16_t *)s;
      for (int i = 0; i < size / 2; i++) t[i] &= mask;
      break;
    }
    case 32: {
      uint32_t *t = (uint32_t *)s;
      for (int i = 0; i < size / 4; i++) t[i] &= mask;
      break;
    }
    case 64: {
      uint64_t *t = (uint64_t *)s;
      for (int i = 0; i < size / 8; i++) t[i] &= mask;
      break;
    }
    default: std::cerr << "Unsupported mask width" << std::endl; abort();
  }
}

template <class T, class U>
struct is_same {
  enum { value = 0 };
};

template <class T>
struct is_same<T, T> {
  enum { value = 1 };
};

void u64_store_aligned(void *p, uint64_t a) {
  v64_store_aligned(p, v64_from_64(a));
}

void c_u64_store_aligned(void *p, uint64_t a) {
  c_v64_store_aligned(p, c_v64_from_64(a));
}

// Ideally, the loads and stores should be function templates, but
// v64 and v128 could be typedefed to the same type (which is the
// case on x86) and then we can't instantiate both v64 and v128.  So
// this workaround takes the combination of all types as template
// arguments and the functions as void pointers so we get no
// matching issues, and then cast the pointers into the proper
// function pointers.
template <typename ret, typename arg, typename c_ret, typename c_arg>
int compare_simd1(void *store, void *load, void *simd, void *d, void *c_store,
                  void *c_load, void *c_simd, void *ref_d, const void *a) {
  void (*my_store)(void *, ret) = (void (*)(void *, ret))store;
  arg (*my_load)(const void *) = (arg(*)(const void *))load;
  ret (*my_simd)(arg) = (ret(*)(arg))simd;
  void (*my_c_store)(void *, c_ret) = (void (*)(void *, c_ret))c_store;
  c_arg (*my_c_load)(const void *) = (c_arg(*)(const void *))c_load;
  c_ret (*my_c_simd)(c_arg) = (c_ret(*)(c_arg))c_simd;

  // Call reference and intrinsic
  my_c_store(ref_d, my_c_simd(my_c_load(a)));
  my_store(d, my_simd(my_load(a)));

  // Compare results
  for (unsigned int pos = 0; pos < sizeof(c_ret) / sizeof(uint32_t); pos++)
    if (((uint32_t *)ref_d)[pos] != ((uint32_t *)d)[pos]) return 1;
  return 0;
}

template <typename ret, typename arg1, typename arg2, typename c_ret,
          typename c_arg1, typename c_arg2>
int compare_simd2(void *store, void *load1, void *load2, void *simd, void *d,
                  void *c_store, void *c_load1, void *c_load2, void *c_simd,
                  void *ref_d, const void *a, const void *b) {
  void (*my_store)(void *, ret) = (void (*)(void *, ret))store;
  arg1 (*my_load1)(const void *) = (arg1(*)(const void *))load1;
  arg2 (*my_load2)(const void *) = (arg2(*)(const void *))load2;
  ret (*my_simd)(arg1, arg2) = (ret(*)(arg1, arg2))simd;
  void (*my_c_store)(void *, c_ret) = (void (*)(void *, c_ret))c_store;
  c_arg1 (*my_c_load1)(const void *) = (c_arg1(*)(const void *))c_load1;
  c_arg2 (*my_c_load2)(const void *) = (c_arg2(*)(const void *))c_load2;
  c_ret (*my_c_simd)(c_arg1, c_arg2) = (c_ret(*)(c_arg1, c_arg2))c_simd;

  // Call reference and intrinsic
  my_c_store(ref_d, my_c_simd(my_c_load1(a), my_c_load2(b)));
  my_store(d, my_simd(my_load1(a), my_load2(b)));

  // Compare results
  for (unsigned int pos = 0; pos < sizeof(c_ret) / sizeof(uint32_t); pos++)
    if (((uint32_t *)ref_d)[pos] != ((uint32_t *)d)[pos]) return 1;
  return 0;
}

template <typename ret, typename arg, typename c_ret, typename c_arg>
void test_simd1(uint32_t iterations, uint32_t mask, uint32_t maskwidth,
                const char *name, ret (*simd)(arg), c_ret (*ref_simd)(c_arg)) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED(sizeof(arg), uint16_t, s[sizeof(arg) / sizeof(uint16_t)]);
  DECLARE_ALIGNED(sizeof(ret), uint8_t, d[sizeof(ret)]);
  DECLARE_ALIGNED(sizeof(ret), uint8_t, ref_d[sizeof(ret)]);
  memset(ref_d, 0, sizeof(ref_d));
  memset(d, 0, sizeof(d));

  int error = 0;
  for (unsigned int count = 0; count < iterations && !error; count++) {
    for (unsigned int c = 0; c < sizeof(arg) / sizeof(uint16_t); c++)
      s[c] = rnd.Rand16();

    if (maskwidth) setmask((uint8_t *)s, sizeof(arg), mask, maskwidth);

    // V64_V64
    if (is_same<c_ret, c_v64>::value && is_same<c_arg, c_v64>::value) {
      error = compare_simd1<ret, arg, c_ret, c_arg>(
          (void *)v64_store_aligned, (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s);
      // U64_V64
    } else if (is_same<c_ret, uint64_t>::value &&
               is_same<c_arg, c_v64>::value) {
      error = compare_simd1<ret, arg, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s);
      // S64_V64
    } else if (is_same<c_ret, int64_t>::value && is_same<c_arg, c_v64>::value) {
      error = compare_simd1<ret, arg, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s);
      // U64_V128
    } else if (is_same<c_ret, uint64_t>::value &&
               is_same<c_arg, c_v128>::value) {
      error = compare_simd1<ret, arg, c_ret, c_arg>(
          (void *)u64_store_aligned, (void *)v128_load_aligned, (void *)simd, d,
          (void *)c_u64_store_aligned, (void *)c_v128_load_aligned,
          (void *)ref_simd, ref_d, s);
      // V128_V128
    } else if (is_same<c_ret, c_v128>::value && is_same<c_arg, c_v128>::value) {
      error = compare_simd1<ret, arg, c_ret, c_arg>(
          (void *)v128_store_aligned, (void *)v128_load_aligned, (void *)simd,
          d, (void *)c_v128_store_aligned, (void *)c_v128_load_aligned,
          (void *)ref_simd, ref_d, s);
      // V128_V64
    } else if (is_same<c_ret, c_v128>::value && is_same<c_arg, c_v64>::value) {
      error = compare_simd1<ret, arg, c_ret, c_arg>(
          (void *)v128_store_aligned, (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_v128_store_aligned, (void *)c_v64_load_aligned,
          (void *)ref_simd, ref_d, s);
    } else {
      std::cerr << "Internal errer: Unknown intrinsic function "
                << typeid(c_ret).name() << " f(" << typeid(c_arg).name() << ")"
                << std::endl;
      abort();
    }
  }

  EXPECT_EQ(0, error) << "Error: mismatch for " << name << "("
                      << print((uint8_t *)s, sizeof(s)) << ") -> "
                      << print(d, sizeof(d)) << " (simd), "
                      << print(ref_d, sizeof(ref_d)) << " (ref)" << std::endl;
}

template <typename ret, typename arg1, typename arg2, typename c_ret,
          typename c_arg1, typename c_arg2>
void test_simd2(uint32_t iterations, uint32_t mask, uint32_t maskwidth,
                const char *name, ret (*simd)(arg1, arg2),
                c_ret (*ref_simd)(c_arg1, c_arg2)) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  DECLARE_ALIGNED(sizeof(arg1), uint16_t, s1[sizeof(arg1) / sizeof(uint16_t)]);
  DECLARE_ALIGNED(sizeof(arg2), uint16_t, s2[sizeof(arg2) / sizeof(uint16_t)]);
  DECLARE_ALIGNED(sizeof(ret), uint8_t, d[sizeof(ret)]);
  DECLARE_ALIGNED(sizeof(ret), uint8_t, ref_d[sizeof(ret)]);
  memset(ref_d, 0, sizeof(ref_d));
  memset(d, 0, sizeof(d));

  int error = 0;
  for (unsigned int count = 0; count < iterations && !error; count++) {
    for (unsigned int c = 0; c < sizeof(arg1) / sizeof(uint16_t); c++)
      s1[c] = rnd.Rand16();

    for (unsigned int c = 0; c < sizeof(arg2) / sizeof(uint16_t); c++)
      s2[c] = rnd.Rand16();

    if (maskwidth) setmask((uint8_t *)s2, sizeof(arg2), mask, maskwidth);

    // V64_V64V64
    if (is_same<c_ret, c_v64>::value && is_same<c_arg1, c_v64>::value &&
        is_same<c_arg2, c_v64>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)v64_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned,
          (void *)c_v64_load_aligned, (void *)ref_simd, ref_d, s1, s2);
      // V128_V128V128
    } else if (is_same<c_ret, c_v128>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, c_v128>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)v128_store_aligned, (void *)v128_load_aligned,
          (void *)v128_load_aligned, (void *)simd, d,
          (void *)c_v128_store_aligned, (void *)c_v128_load_aligned,
          (void *)c_v128_load_aligned, (void *)ref_simd, ref_d, s1, s2);
      // U32_V64V64
    } else if (is_same<c_ret, uint32_t>::value &&
               is_same<c_arg1, c_v64>::value && is_same<c_arg2, c_v64>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)u32_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_u32_store_aligned, (void *)c_v64_load_aligned,
          (void *)c_v64_load_aligned, (void *)ref_simd, ref_d, s1, s2);
      // U32_V128V128
    } else if (is_same<c_ret, uint32_t>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, c_v128>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)u32_store_aligned, (void *)v128_load_aligned,
          (void *)v128_load_aligned, (void *)simd, d,
          (void *)c_u32_store_aligned, (void *)c_v128_load_aligned,
          (void *)c_v128_load_aligned, (void *)ref_simd, ref_d, s1, s2);
      // V128_V64V64
    } else if (is_same<c_ret, c_v128>::value && is_same<c_arg1, c_v64>::value &&
               is_same<c_arg2, c_v64>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)v128_store_aligned, (void *)v64_load_aligned,
          (void *)v64_load_aligned, (void *)simd, d,
          (void *)c_v128_store_aligned, (void *)c_v64_load_aligned,
          (void *)c_v64_load_aligned, (void *)ref_simd, ref_d, s1, s2);
      // V128_V128U32
    } else if (is_same<c_ret, c_v128>::value &&
               is_same<c_arg1, c_v128>::value &&
               is_same<c_arg2, uint32_t>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)v128_store_aligned, (void *)v128_load_aligned,
          (void *)u32_load_aligned, (void *)simd, d,
          (void *)c_v128_store_aligned, (void *)c_v128_load_aligned,
          (void *)c_u32_load_aligned, (void *)ref_simd, ref_d, s1, s2);
      // V64_V64U32
    } else if (is_same<c_ret, c_v64>::value && is_same<c_arg1, c_v64>::value &&
               is_same<c_arg2, uint32_t>::value) {
      error = compare_simd2<ret, arg1, arg2, c_ret, c_arg1, c_arg2>(
          (void *)v64_store_aligned, (void *)v64_load_aligned,
          (void *)u32_load_aligned, (void *)simd, d,
          (void *)c_v64_store_aligned, (void *)c_v64_load_aligned,
          (void *)c_u32_load_aligned, (void *)ref_simd, ref_d, s1, s2);
    } else {
      std::cout << "Internal errer: Unknown intrinsic function "
                << typeid(c_ret).name() << " f(" << typeid(c_arg1).name()
                << ", " << typeid(c_arg2).name() << ")" << std::endl;
      abort();
    }
  }

  EXPECT_EQ(0, error) << "Error: mismatch for " << name << "("
                      << print((uint8_t *)s1, sizeof(s1)) << ", "
                      << print((uint8_t *)s2, sizeof(s2)) << ") -> "
                      << print(d, sizeof(d)) << " (simd), "
                      << print(ref_d, sizeof(ref_d)) << " (ref)" << std::endl;
}

const int iterations = 65536;

TEST_P(V64_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(U64_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(S64_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(U64_V128, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V128, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V64_V64V64, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V128V128, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(U32_V64V64, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(U32_V128V128, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V64V64, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V128U32, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V64_V64U32, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

// Google Test allows up to 50 tests per case, so split the largest
TEST_P(V64_V64_Part2, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V64_V64V64_Part2, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V128V128_Part2, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V128_Part2, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

TEST_P(V128_V128_Part3, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, name, simd, ref_simd);
}

}  // namespace simd_test
