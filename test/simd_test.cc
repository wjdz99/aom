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

#include <cstdlib>
#include <string>
#include <typeinfo>

#include "third_party/googletest/src/include/gtest/gtest.h"

#include "./aom_config.h"
#include "./aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"
#include "aom_dsp/aom_simd.h"
#include "aom_dsp/simd/v128_intrinsics_c.h"

using libaom_test::ACMRandom;

namespace {

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
  }

  virtual void TearDown() { libaom_test::ClearSystemState(); }

 protected:
  simd_signature simd;
  ref_signature ref_simd;
  uint32_t mask, maskwidth;
};

typedef uint32_t c_uint32_t;
typedef uint64_t c_uint64_t;
typedef int64_t c_int64_t;

// Wrap templates around intrinsics using immediate values
template <int shift>
v64 imm_v64_shl_n_byte(v64 a) {
  return v64_shl_n_byte(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_byte(v64 a) {
  return v64_shr_n_byte(a, shift);
}
template <int shift>
v64 imm_v64_shl_n_8(v64 a) {
  return v64_shl_n_8(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_u8(v64 a) {
  return v64_shr_n_u8(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_s8(v64 a) {
  return v64_shr_n_s8(a, shift);
}
template <int shift>
v64 imm_v64_shl_n_16(v64 a) {
  return v64_shl_n_16(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_u16(v64 a) {
  return v64_shr_n_u16(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_s16(v64 a) {
  return v64_shr_n_s16(a, shift);
}
template <int shift>
v64 imm_v64_shl_n_32(v64 a) {
  return v64_shl_n_32(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_u32(v64 a) {
  return v64_shr_n_u32(a, shift);
}
template <int shift>
v64 imm_v64_shr_n_s32(v64 a) {
  return v64_shr_n_s32(a, shift);
}
template <int shift>
v64 imm_v64_align(v64 a, v64 b) {
  return v64_align(a, b, shift);
}

template <int shift>
v128 imm_v128_shl_n_byte(v128 a) {
  return v128_shl_n_byte(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_byte(v128 a) {
  return v128_shr_n_byte(a, shift);
}
template <int shift>
v128 imm_v128_shl_n_8(v128 a) {
  return v128_shl_n_8(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_u8(v128 a) {
  return v128_shr_n_u8(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_s8(v128 a) {
  return v128_shr_n_s8(a, shift);
}
template <int shift>
v128 imm_v128_shl_n_16(v128 a) {
  return v128_shl_n_16(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_u16(v128 a) {
  return v128_shr_n_u16(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_s16(v128 a) {
  return v128_shr_n_s16(a, shift);
}
template <int shift>
v128 imm_v128_shl_n_32(v128 a) {
  return v128_shl_n_32(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_u32(v128 a) {
  return v128_shr_n_u32(a, shift);
}
template <int shift>
v128 imm_v128_shr_n_s32(v128 a) {
  return v128_shr_n_s32(a, shift);
}
template <int shift>
v128 imm_v128_align(v128 a, v128 b) {
  return v128_align(a, b, shift);
}

template <int shift>
c_v64 c_imm_v64_shl_n_byte(c_v64 a) {
  return c_v64_shl_n_byte(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_byte(c_v64 a) {
  return c_v64_shr_n_byte(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shl_n_8(c_v64 a) {
  return c_v64_shl_n_8(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_u8(c_v64 a) {
  return c_v64_shr_n_u8(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_s8(c_v64 a) {
  return c_v64_shr_n_s8(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shl_n_16(c_v64 a) {
  return c_v64_shl_n_16(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_u16(c_v64 a) {
  return c_v64_shr_n_u16(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_s16(c_v64 a) {
  return c_v64_shr_n_s16(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shl_n_32(c_v64 a) {
  return c_v64_shl_n_32(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_u32(c_v64 a) {
  return c_v64_shr_n_u32(a, shift);
}
template <int shift>
c_v64 c_imm_v64_shr_n_s32(c_v64 a) {
  return c_v64_shr_n_s32(a, shift);
}
template <int shift>
c_v64 c_imm_v64_align(c_v64 a, c_v64 b) {
  return c_v64_align(a, b, shift);
}

template <int shift>
c_v128 c_imm_v128_shl_n_byte(c_v128 a) {
  return c_v128_shl_n_byte(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_byte(c_v128 a) {
  return c_v128_shr_n_byte(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shl_n_8(c_v128 a) {
  return c_v128_shl_n_8(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_u8(c_v128 a) {
  return c_v128_shr_n_u8(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_s8(c_v128 a) {
  return c_v128_shr_n_s8(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shl_n_16(c_v128 a) {
  return c_v128_shl_n_16(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_u16(c_v128 a) {
  return c_v128_shr_n_u16(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_s16(c_v128 a) {
  return c_v128_shr_n_s16(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shl_n_32(c_v128 a) {
  return c_v128_shl_n_32(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_u32(c_v128 a) {
  return c_v128_shr_n_u32(a, shift);
}
template <int shift>
c_v128 c_imm_v128_shr_n_s32(c_v128 a) {
  return c_v128_shr_n_s32(a, shift);
}
template <int shift>
c_v128 c_imm_v128_align(c_v128 a, c_v128 b) {
  return c_v128_align(a, b, shift);
}

// Wrappers around the the SAD and SSD functions
uint32_t v64_sad_u8(v64 a, v64 b) {
  return v64_sad_u8_sum(::v64_sad_u8(v64_sad_u8_init(), a, b));
}
uint32_t v64_ssd_u8(v64 a, v64 b) {
  return v64_ssd_u8_sum(::v64_ssd_u8(v64_ssd_u8_init(), a, b));
}
uint32_t v128_sad_u8(v128 a, v128 b) {
  return v128_sad_u8_sum(::v128_sad_u8(v128_sad_u8_init(), a, b));
}
uint32_t v128_ssd_u8(v128 a, v128 b) {
  return v128_ssd_u8_sum(::v128_ssd_u8(v128_ssd_u8_init(), a, b));
}

uint32_t c_v64_sad_u8(c_v64 a, c_v64 b) {
  return c_v64_sad_u8_sum(::c_v64_sad_u8(c_v64_sad_u8_init(), a, b));
}
uint32_t c_v64_ssd_u8(c_v64 a, c_v64 b) {
  return c_v64_ssd_u8_sum(::c_v64_ssd_u8(c_v64_ssd_u8_init(), a, b));
}
uint32_t c_v128_sad_u8(c_v128 a, c_v128 b) {
  return c_v128_sad_u8_sum(::c_v128_sad_u8(c_v128_sad_u8_init(), a, b));
}
uint32_t c_v128_ssd_u8(c_v128 a, c_v128 b) {
  return c_v128_ssd_u8_sum(::c_v128_ssd_u8(c_v128_ssd_u8_init(), a, b));
}

// Crate one typedef for each function signature
#define TYPEDEF_SIMD1(name, r, a)                                      \
  typedef TestIntrinsic<                                               \
      std::tr1::tuple<r (*)(a), c_##r (*)(c_##a), uint32_t, uint32_t>, \
      r (*)(a), c_##r (*)(c_##a)>                                      \
      name;

#define TYPEDEF_SIMD2(name, r, a, b)                                          \
  typedef TestIntrinsic<std::tr1::tuple<r (*)(a, b), c_##r (*)(c_##a, c_##b), \
                                        uint32_t, uint32_t>,                  \
                        r (*)(a, b), c_##r (*)(c_##a, c_##b)>                 \
      name;

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
                ret (*simd)(arg), c_ret (*ref_simd)(c_arg)) {
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

  EXPECT_EQ(0, error) << "Error: mismatch for f("
                      << print((uint8_t *)s, sizeof(s)) << ") -> "
                      << print(d, sizeof(d)) << " (simd), "
                      << print(ref_d, sizeof(ref_d)) << " (ref)" << std::endl;
}

template <typename ret, typename arg1, typename arg2, typename c_ret,
          typename c_arg1, typename c_arg2>
void test_simd2(uint32_t iterations, uint32_t mask, uint32_t maskwidth,
                ret (*simd)(arg1, arg2), c_ret (*ref_simd)(c_arg1, c_arg2)) {
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

  EXPECT_EQ(0, error) << "Error: mismatch for f("
                      << print((uint8_t *)s1, sizeof(s1)) << ", "
                      << print((uint8_t *)s2, sizeof(s2)) << ") -> "
                      << print(d, sizeof(d)) << " (simd), "
                      << print(ref_d, sizeof(ref_d)) << " (ref)" << std::endl;
}

const int iterations = 65536;

TEST_P(V64_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(U64_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(S64_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(U64_V128, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V128_V128, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V128_V64, TestIntrinsics) {
  test_simd1(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V64_V64V64, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V128_V128V128, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(U32_V64V64, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(U32_V128V128, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V128_V64V64, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V128_V128U32, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

TEST_P(V64_V64U32, TestIntrinsics) {
  test_simd2(iterations, mask, maskwidth, simd, ref_simd);
}

#define INSTANTIATE_SIMD_TEST(simd, test, mask, maskwidth)                    \
  INSTANTIATE_TEST_CASE_P(simd, test, ::testing::Values(std::tr1::make_tuple( \
                                          simd, c_##simd, mask, maskwidth)))

#define INSTANTIATE_SIMD_IMM_TEST(simd, test, imm)                \
  INSTANTIATE_TEST_CASE_P(simd##__##imm, test,                    \
                          ::testing::Values(std::tr1::make_tuple( \
                              imm_##simd<imm>, c_imm_##simd<imm>, 0U, 0U)))

INSTANTIATE_SIMD_TEST(v128_sad_u8, U32_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ssd_u8, U32_V128V128, 0U, 0U);

INSTANTIATE_SIMD_TEST(v64_sad_u8, U32_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ssd_u8, U32_V64V64, 0U, 0U);

INSTANTIATE_SIMD_TEST(v128_add_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_add_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_sadd_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_add_32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_sub_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ssub_u8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ssub_s8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_sub_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ssub_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_sub_32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziplo_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziphi_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziplo_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziphi_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziplo_32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziphi_32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziplo_64, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_ziphi_64, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unziphi_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unziplo_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unziphi_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unziplo_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unziphi_32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unziplo_32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_pack_s32_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_pack_s16_u8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_pack_s16_s8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_or, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_xor, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_and, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_andn, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_mullo_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_mulhi_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_mullo_s32, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_madd_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_madd_us8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_avg_u8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_rdavg_u8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_avg_u16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_min_u8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_max_u8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_min_s8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_max_s8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_min_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_max_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_cmpgt_s8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_cmplt_s8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_cmpeq_8, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_cmpgt_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_cmplt_s16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_cmpeq_16, V128_V128V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_shuffle_8, V128_V128V128, 15U, 8U);

INSTANTIATE_SIMD_TEST(v64_add_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_add_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_sadd_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_add_32, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_sub_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ssub_u8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ssub_s8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_sub_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ssub_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_sub_32, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ziplo_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ziphi_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ziplo_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ziphi_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ziplo_32, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_ziphi_32, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_pack_s32_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_pack_s16_u8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_pack_s16_s8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unziphi_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unziplo_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unziphi_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unziplo_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_or, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_xor, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_and, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_andn, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_mullo_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_mulhi_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_mullo_s32, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_madd_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_madd_us8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_avg_u8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_rdavg_u8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_avg_u16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_min_u8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_max_u8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_min_s8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_max_s8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_min_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_max_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_cmpgt_s8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_cmplt_s8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_cmpeq_8, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_cmpgt_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_cmplt_s16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_cmpeq_16, V64_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_shuffle_8, V64_V64V64, 7U, 8U);

INSTANTIATE_SIMD_TEST(v128_abs_s16, V128_V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_padd_s16, V128_V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unpacklo_u16_s32, V128_V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unpacklo_s16_s32, V128_V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unpackhi_u16_s32, V128_V128, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unpackhi_s16_s32, V128_V128, 0U, 0U);

INSTANTIATE_SIMD_TEST(v64_abs_s16, V64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unpacklo_u8_s16, V64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unpackhi_u8_s16, V64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unpacklo_u16_s32, V64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unpacklo_s16_s32, V64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unpackhi_u16_s32, V64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_unpackhi_s16_s32, V64_V64, 0U, 0U);

INSTANTIATE_SIMD_TEST(v128_from_v64, V128_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_zip_8, V128_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_zip_16, V128_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_zip_32, V128_V64V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_mul_s16, V128_V64V64, 0U, 0U);

INSTANTIATE_SIMD_TEST(v128_unpack_u8_s16, V128_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unpack_u16_s32, V128_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_unpack_s16_s32, V128_V64, 0U, 0U);

INSTANTIATE_SIMD_TEST(v128_shl_8, V128_V128U32, 7U, 32U);
INSTANTIATE_SIMD_TEST(v128_shr_u8, V128_V128U32, 7U, 32U);
INSTANTIATE_SIMD_TEST(v128_shr_s8, V128_V128U32, 7U, 32U);
INSTANTIATE_SIMD_TEST(v128_shl_16, V128_V128U32, 15U, 32U);
INSTANTIATE_SIMD_TEST(v128_shr_u16, V128_V128U32, 15U, 32U);
INSTANTIATE_SIMD_TEST(v128_shr_s16, V128_V128U32, 15U, 32U);
INSTANTIATE_SIMD_TEST(v128_shl_32, V128_V128U32, 31U, 32U);
INSTANTIATE_SIMD_TEST(v128_shr_u32, V128_V128U32, 31U, 32U);
INSTANTIATE_SIMD_TEST(v128_shr_s32, V128_V128U32, 31U, 32U);

INSTANTIATE_SIMD_TEST(v64_shl_8, V64_V64U32, 7U, 32U);
INSTANTIATE_SIMD_TEST(v64_shr_u8, V64_V64U32, 7U, 32U);
INSTANTIATE_SIMD_TEST(v64_shr_s8, V64_V64U32, 7U, 32U);
INSTANTIATE_SIMD_TEST(v64_shl_16, V64_V64U32, 15U, 32U);
INSTANTIATE_SIMD_TEST(v64_shr_u16, V64_V64U32, 15U, 32U);
INSTANTIATE_SIMD_TEST(v64_shr_s16, V64_V64U32, 15U, 32U);
INSTANTIATE_SIMD_TEST(v64_shl_32, V64_V64U32, 31U, 32U);
INSTANTIATE_SIMD_TEST(v64_shr_u32, V64_V64U32, 31U, 32U);
INSTANTIATE_SIMD_TEST(v64_shr_s32, V64_V64U32, 31U, 32U);

INSTANTIATE_SIMD_TEST(v64_hadd_u8, U64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v64_hadd_s16, S64_V64, 0U, 0U);
INSTANTIATE_SIMD_TEST(v128_hadd_u8, U64_V128, 0U, 0U);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 3);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 5);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 7);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 9);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 10);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 11);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 13);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 14);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_byte, V128_V128, 15);

INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 3);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 5);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 7);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 9);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 10);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 11);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 13);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 14);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_byte, V128_V128, 15);

INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 3);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 5);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_8, V128_V128, 7);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 3);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 5);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u8, V128_V128, 7);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 3);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 5);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s8, V128_V128, 7);

INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 10);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_16, V128_V128, 14);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 10);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u16, V128_V128, 14);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 10);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s16, V128_V128, 14);

INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 16);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 20);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 24);
INSTANTIATE_SIMD_IMM_TEST(v128_shl_n_32, V128_V128, 28);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 16);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 20);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 24);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_u32, V128_V128, 28);

INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 16);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 20);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 24);
INSTANTIATE_SIMD_IMM_TEST(v128_shr_n_s32, V128_V128, 28);

INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 1);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 2);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 3);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 4);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 5);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 6);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 7);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 8);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 9);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 10);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 11);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 12);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 13);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 14);
INSTANTIATE_SIMD_IMM_TEST(v128_align, V128_V128V128, 15);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 3);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 5);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_byte, V64_V64, 7);

INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 3);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 5);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_byte, V64_V64, 7);

INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 3);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 5);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_8, V64_V64, 7);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 3);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 5);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u8, V64_V64, 7);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 3);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 5);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s8, V64_V64, 7);

INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 8);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 10);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 12);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_16, V64_V64, 14);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 8);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 10);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 12);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u16, V64_V64, 14);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 8);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 10);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 12);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s16, V64_V64, 14);

INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 8);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 12);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 16);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 20);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 24);
INSTANTIATE_SIMD_IMM_TEST(v64_shl_n_32, V64_V64, 28);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 8);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 12);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 16);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 20);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 24);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_u32, V64_V64, 28);

INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 8);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 12);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 16);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 20);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 24);
INSTANTIATE_SIMD_IMM_TEST(v64_shr_n_s32, V64_V64, 28);

INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 1);
INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 2);
INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 3);
INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 4);
INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 5);
INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 6);
INSTANTIATE_SIMD_IMM_TEST(v64_align, V64_V64V64, 7);

}  // namespace
