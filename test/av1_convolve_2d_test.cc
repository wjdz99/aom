/*
 * Copyright (c) 2020, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include <memory>
#include <set>
#include <vector>
#include "config/av1_rtcd.h"
#include "config/aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

// We cannot use the BLOCK_SIZE enum because we need to simulate chroma block
// sizes, which can go down to 2xN and Nx2.
struct BlockSize {
  BlockSize() : width(0), height(0) {}
  BlockSize(int w, int h) : width(w), height(h) {}
  int width;
  int height;
  bool operator<(const BlockSize &other) const {
    if (width == other.width) {
      return height < other.height;
    }
    return width < other.width;
  }
};

// Generate the list of all block widths / heights that need to be tested.
std::set<BlockSize> GetBlockSizes() {
  std::set<BlockSize> sizes;
  for (int b = BLOCK_4X4; b < BLOCK_SIZES_ALL; ++b) {
    const int w = block_size_wide[b];
    const int h = block_size_high[b];
    sizes.insert(BlockSize(w, h));
    // Add in smaller chroma sizes as well.
    if (w == 4 || h == 4) {
      sizes.insert(BlockSize(w / 2, h / 2));
    }
  }
  return sizes;
}

::testing::internal::ParamGenerator<BlockSize> GenerateBlockSizes() {
  return ::testing::ValuesIn(GetBlockSizes());
}

// ConvolveTest is the base class that all convolve tests should derive from.
// It provides storage/methods for generating randomized buffers for both
// low bit-depth and high bit-depth.
template <typename T>
class ConvolveTest : public ::testing::TestWithParam<T> {
 public:
  virtual ~ConvolveTest() {}

  virtual void SetUp() override {
    rnd_.Reset(libaom_test::ACMRandom::DeterministicSeed());
  }

  virtual void TearDown() override {
    libaom_test::ClearSystemState();
    buffer8_.clear();
    buffer16_.clear();
  }

  // Randomizes the 8-bit input buffer and returns a pointer to it. Note that
  // the pointer is safe to use with an 8-tap filter. The stride can range
  // from width to (width + kPadding).
  static constexpr int kPadding = 8;

  const uint8_t *RandomInput8(int width, int height) {
    const int padded_width = width + kPadding;
    const int padded_height = height + kPadding;
    return RandomUint8(padded_width * padded_height) + 3 * padded_width + 3;
  }

  // Generate a random 16-bit input buffer, like above. Note that the
  // values are capped so they do not exceed the bit-depth.
  const uint16_t *RandomInput16(int width, int height, int bit_depth) {
    const int padded_width = width + kPadding;
    const int padded_height = height + kPadding;
    return RandomUint16(padded_width * padded_height, bit_depth) +
           3 * padded_width + 3;
  }

  // 8-bit output buffer of size width * height. It is aligned on a 32-byte
  // boundary.
  uint8_t *Output8(int width, int height) {
    size_t size = width * height + 32;
    buffer8_.emplace_back(new uint8_t[size]);
    uint8_t *p = buffer8_.back().get();
    void *vp = static_cast<void *>(p);
    return static_cast<uint8_t *>(std::align(32, sizeof(p[0]), vp, size));
  }

  // 16-bit output buffer of size width * height. It is aligned on a 32-byte
  // boundary.
  uint16_t *Output16(int width, int height) {
    size_t size = width * height + 16;
    buffer16_.emplace_back(new uint16_t[size]);
    uint16_t *p = buffer16_.back().get();
    void *vp = static_cast<void *>(p);
    return static_cast<uint16_t *>(std::align(32, sizeof(p[0]), vp, size));
  }

  // Check that two 8-bit buffers are identical.
  void AssertEq(const uint8_t *p1, const uint8_t *p2, int width, int height,
                int stride) {
    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        int idx = i + j * stride;
        ASSERT_EQ(p1[idx], p2[idx])
            << width << "x" << height << " Pixel mismatch at index " << idx
            << " = (" << i << ", " << j << ")";
      }
    }
  }

  // Check that two 16-bit buffers are identical.
  void AssertEq(const uint16_t *p1, const uint16_t *p2, int width, int height,
                int stride) {
    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        int idx = i + j * stride;
        ASSERT_EQ(p1[idx], p2[idx])
            << width << "x" << height << " Pixel mismatch at index " << idx
            << " = (" << i << ", " << j << ")";
      }
    }
  }

 private:
  const uint8_t *RandomUint8(int size) {
    buffer8_.emplace_back(new uint8_t[size]);
    uint8_t *p = buffer8_.back().get();
    for (int i = 0; i < size; ++i) {
      p[i] = rnd_.Rand8();
    }
    return p;
  }

  const uint16_t *RandomUint16(int size, int bit_depth) {
    buffer16_.emplace_back(new uint16_t[size]);
    uint16_t *p = buffer16_.back().get();
    for (int i = 0; i < size; ++i) {
      p[i] = rnd_.Rand16() & ((1 << bit_depth) - 1);
    }
    return p;
  }

  std::vector<std::unique_ptr<uint8_t[]>> buffer8_;
  std::vector<std::unique_ptr<uint16_t[]>> buffer16_;
  libaom_test::ACMRandom rnd_;
};

// Base structure for 1d-convolve parameters.
struct Convolve1DParam {
  BlockSize block_size;
  int sub_pixel;
  InterpFilter filter;
  int bit_depth;
};

::testing::internal::ParamGenerator<Convolve1DParam> GenerateConvolve1DParams(
    std::initializer_list<int> bit_depths) {
  std::vector<Convolve1DParam> result;
  for (const auto &block_size : GetBlockSizes()) {
    for (int i = 0; i < 16; ++i) {
      for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
        Convolve1DParam p;
        p.block_size = block_size;
        p.sub_pixel = i;
        p.filter = static_cast<InterpFilter>(f);
        for (const int bit_depth : bit_depths) {
          p.bit_depth = bit_depth;
          result.push_back(p);
        }
      }
    }
  }
  return ::testing::ValuesIn(result);
}

////////////////////////////////////////////////////////
// Single reference convolve-x functions (low bit-depth)
////////////////////////////////////////////////////////
typedef void (*convolve_x_func)(const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int w, int h,
                                const InterpFilterParams *filter_params_x,
                                const int subpel_x_qn,
                                ConvolveParams *conv_params);

class LowbdConvolveXTest : public ConvolveTest<Convolve1DParam> {
 public:
  void RunTest(convolve_x_func test_func) {
    BlockSize b = GetParam().block_size;
    int sub_x = GetParam().sub_pixel;
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(GetParam().filter,
                                                     b.width);
    ConvolveParams conv_params1 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    const uint8_t *input = RandomInput8(b.width, b.height);
    uint8_t *reference = Output8(b.width, b.height);
    av1_convolve_x_sr(input, b.width, reference, b.width, b.width, b.height,
                      filter_params_x, sub_x, &conv_params1);

    ConvolveParams conv_params2 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    uint8_t *test = Output8(b.width, b.height);
    test_func(input, b.width, test, b.width, b.width, b.height, filter_params_x,
              sub_x, &conv_params2);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(LowbdConvolveXTest, C) { RunTest(av1_convolve_x_sr_c); };

#if HAVE_SSE2
TEST_P(LowbdConvolveXTest, SSE2) { RunTest(av1_convolve_x_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(LowbdConvolveXTest, AVX2) { RunTest(av1_convolve_x_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(LowbdConvolveXTest, NEON) { RuNTest(av1_convolve_x_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdConvolveXTest,
                        GenerateConvolve1DParams({ 8 }));

/////////////////////////////////////////////////////////
// Single reference convolve-x functions (high bit-depth)
/////////////////////////////////////////////////////////
typedef void (*highbd_convolve_x_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn,
    ConvolveParams *conv_params, int bd);

class HighbdConvolveXTest : public ConvolveTest<Convolve1DParam> {
 public:
  void RunTest(highbd_convolve_x_func test_func) {
    BlockSize b = GetParam().block_size;
    int sub_x = GetParam().sub_pixel;
    int bit_depth = GetParam().bit_depth;
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(GetParam().filter,
                                                     b.width);
    ConvolveParams conv_params1 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    const uint16_t *input = RandomInput16(b.width, b.height, bit_depth);
    uint16_t *reference = Output16(b.width, b.height);
    av1_highbd_convolve_x_sr(input, b.width, reference, b.width, b.width,
                             b.height, filter_params_x, sub_x, &conv_params1,
                             bit_depth);

    ConvolveParams conv_params2 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    uint16_t *test = Output16(b.width, b.height);
    test_func(input, b.width, test, b.width, b.width, b.height, filter_params_x,
              sub_x, &conv_params2, bit_depth);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(HighbdConvolveXTest, C) { RunTest(av1_highbd_convolve_x_sr_c); };

#if HAVE_SSSE3
TEST_P(HighbdConvolveXTest, SSSE3) { RunTest(av1_highbd_convolve_x_sr_ssse3); }
#endif

#if HAVE_AVX2
TEST_P(HighbdConvolveXTest, AVX2) { RunTest(av1_highbd_convolve_x_sr_avx2); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdConvolveXTest,
                        GenerateConvolve1DParams({ 10, 12 }));

////////////////////////////////////////////////////////
// Single reference convolve-y functions (low bit-depth)
////////////////////////////////////////////////////////
typedef void (*convolve_y_func)(const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int w, int h,
                                const InterpFilterParams *filter_params_y,
                                const int subpel_y_qn);

class LowbdConvolveYTest : public ConvolveTest<Convolve1DParam> {
 public:
  void RunTest(convolve_y_func test_func) {
    BlockSize b = GetParam().block_size;
    int sub_y = GetParam().sub_pixel;
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(GetParam().filter,
                                                     b.height);
    const uint8_t *input = RandomInput8(b.width, b.height);
    uint8_t *reference = Output8(b.width, b.height);
    av1_convolve_y_sr(input, b.width, reference, b.width, b.width, b.height,
                      filter_params_y, sub_y);
    uint8_t *test = Output8(b.width, b.height);
    test_func(input, b.width, test, b.width, b.width, b.height, filter_params_y,
              sub_y);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(LowbdConvolveYTest, C) { RunTest(av1_convolve_y_sr_c); };

#if HAVE_SSE2
TEST_P(LowbdConvolveYTest, SSE2) { RunTest(av1_convolve_y_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(LowbdConvolveYTest, AVX2) { RunTest(av1_convolve_y_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(LowbdConvolveYTest, NEON) { RunTest(av1_convolve_y_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdConvolveYTest,
                        GenerateConvolve1DParams({ 8 }));

/////////////////////////////////////////////////////////
// Single reference convolve-y functions (high bit-depth)
/////////////////////////////////////////////////////////
typedef void (*highbd_convolve_y_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn,
    int bd);

class HighbdConvolveYTest : public ConvolveTest<Convolve1DParam> {
 public:
  void RunTest(highbd_convolve_y_func test_func) {
    BlockSize b = GetParam().block_size;
    const int bit_depth = GetParam().bit_depth;
    int sub_y = GetParam().sub_pixel;
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(GetParam().filter,
                                                     b.height);
    const uint16_t *input = RandomInput16(b.width, b.height, bit_depth);
    uint16_t *reference = Output16(b.width, b.height);
    av1_highbd_convolve_y_sr(input, b.width, reference, b.width, b.width,
                             b.height, filter_params_y, sub_y, bit_depth);
    uint16_t *test = Output16(b.width, b.height);
    test_func(input, b.width, test, b.width, b.width, b.height, filter_params_y,
              sub_y, bit_depth);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(HighbdConvolveYTest, C) { RunTest(av1_highbd_convolve_y_sr_c); };

#if HAVE_SSSE3
TEST_P(HighbdConvolveYTest, SSSE3) { RunTest(av1_highbd_convolve_y_sr_ssse3); }
#endif

#if HAVE_AVX2
TEST_P(HighbdConvolveYTest, AVX2) { RunTest(av1_highbd_convolve_y_sr_avx2); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdConvolveYTest,
                        GenerateConvolve1DParams({ 10, 12 }));

//////////////////////////////////////////////////////////////
// Single reference convolve-2d-copy functions (low bit-depth)
//////////////////////////////////////////////////////////////
typedef void (*convolve_2d_copy_func)(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride, int w,
                                      int h);

class LowbdConvolve2DCopyTest : public ConvolveTest<BlockSize> {
 public:
  void RunTest(convolve_2d_copy_func test_func) {
    BlockSize b = GetParam();
    const uint8_t *input = RandomInput8(b.width, b.height);
    uint8_t *reference = Output8(b.width, b.height);
    av1_convolve_2d_copy_sr(input, b.width, reference, b.width, b.width,
                            b.height);
    uint8_t *test = Output8(b.width, b.height);
    test_func(input, b.width, test, b.width, b.width, b.height);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(LowbdConvolve2DCopyTest, C) { RunTest(av1_convolve_2d_copy_sr_c); }

#if HAVE_SSE2
TEST_P(LowbdConvolve2DCopyTest, SSE2) { RunTest(av1_convolve_2d_copy_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(LowbdConvolve2DCopyTest, AVX2) { RunTest(av1_convolve_2d_copy_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(LowbdConvolve2DCopyTest, NEON) { RunTest(av1_convolve_2d_copy_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdConvolve2DCopyTest,
                        GenerateBlockSizes());

///////////////////////////////////////////////////////////////
// Single reference convolve-2d-copy functions (high bit-depth)
///////////////////////////////////////////////////////////////
typedef void (*highbd_convolve_2d_copy_func)(const uint16_t *src,
                                             int src_stride, uint16_t *dst,
                                             int dst_stride, int w, int h);

class HighbdConvolve2DCopyTest : public ConvolveTest<BlockSize> {
 public:
  void RunTest(highbd_convolve_2d_copy_func test_func) {
    RunTestAux(test_func, 10);
    RunTestAux(test_func, 12);
  }

 private:
  void RunTestAux(highbd_convolve_2d_copy_func test_func, int bit_depth) {
    BlockSize b = GetParam();
    const uint16_t *input = RandomInput16(b.width, b.height, bit_depth);
    uint16_t *reference = Output16(b.width, b.height);
    av1_highbd_convolve_2d_copy_sr(input, b.width, reference, b.width, b.width,
                                   b.height);
    uint16_t *test = Output16(b.width, b.height);
    test_func(input, b.width, test, b.width, b.width, b.height);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(HighbdConvolve2DCopyTest, C) {
  RunTest(av1_highbd_convolve_2d_copy_sr_c);
}

#if HAVE_SSE2
TEST_P(HighbdConvolve2DCopyTest, SSE2) {
  RunTest(av1_highbd_convolve_2d_copy_sr_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(HighbdConvolve2DCopyTest, AVX2) {
  RunTest(av1_highbd_convolve_2d_copy_sr_avx2);
}
#endif

#if HAVE_NEON
TEST_P(HighbdConvolve2DCopyTest, NEON) {
  RunTest(av1_highbd_convolve_2d_copy_sr_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdConvolve2DCopyTest,
                        GenerateBlockSizes());

// 2D convolve parameters.
struct Convolve2DParam {
  BlockSize block_size;
  InterpFilter h_filter;
  InterpFilter v_filter;
  int sub_x;
  int sub_y;
  int bit_depth;
};

::testing::internal::ParamGenerator<Convolve2DParam> GenerateConvolve2DParams(
    std::initializer_list<int> bit_depths) {
  std::vector<Convolve2DParam> result;
  Convolve2DParam p;
  for (const auto &block_size : GetBlockSizes()) {
    p.block_size = block_size;
    for (int x = 0; x < 16; ++x) {
      p.sub_x = x;
      for (int h_f = EIGHTTAP_REGULAR; h_f < INTERP_FILTERS_ALL; ++h_f) {
        p.h_filter = static_cast<InterpFilter>(h_f);
        for (int y = 0; y < 16; ++y) {
          p.sub_y = y;
          for (int v_f = EIGHTTAP_REGULAR; v_f < INTERP_FILTERS_ALL; ++v_f) {
            p.v_filter = static_cast<InterpFilter>(v_f);
            for (const int bit_depth : bit_depths) {
              p.bit_depth = bit_depth;
              result.push_back(p);
            }
          }
        }
      }
    }
  }
  return ::testing::ValuesIn(result);
}

/////////////////////////////////////////////////////////
// Single reference convolve-2D functions (low bit-depth)
/////////////////////////////////////////////////////////
typedef void (*convolve_2d_func)(const uint8_t *src, int src_stride,
                                 uint8_t *dst, int dst_stride, int w, int h,
                                 const InterpFilterParams *filter_params_x,
                                 const InterpFilterParams *filter_params_y,
                                 const int subpel_x_qn, const int subpel_y_qn,
                                 ConvolveParams *conv_params);

class LowbdConvolve2DTest : public ConvolveTest<Convolve2DParam> {
 public:
  void RunTest(convolve_2d_func test_func) {
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(
            GetParam().h_filter, GetParam().block_size.width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(
            GetParam().v_filter, GetParam().block_size.height);
    BlockSize b = GetParam().block_size;
    const int sub_x = GetParam().sub_x;
    const int sub_y = GetParam().sub_y;
    const uint8_t *input = RandomInput8(b.width, b.height);
    uint8_t *reference = Output8(b.width, b.height);
    ConvolveParams conv_params1 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    av1_convolve_2d_sr(input, b.width, reference, b.width, b.width, b.height,
                       filter_params_x, filter_params_y, sub_x, sub_y,
                       &conv_params1);
    uint8_t *test = Output8(b.width, b.height);
    ConvolveParams conv_params2 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    test_func(input, b.width, test, b.width, b.width, b.height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params2);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(LowbdConvolve2DTest, C) { RunTest(av1_convolve_2d_sr_c); }

#if HAVE_SSE2
TEST_P(LowbdConvolve2DTest, SSE2) { RunTest(av1_convolve_2d_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(LowbdConvolve2DTest, AVX2) { RunTest(av1_convolve_2d_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(LowbdConvolve2DTest, NEON) { RunTest(av1_convolve_2d_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdConvolve2DTest,
                        GenerateConvolve2DParams({ 8 }));

//////////////////////////////////////////////////////////
// Single reference convolve-2d functions (high bit-depth)
//////////////////////////////////////////////////////////

typedef void (*highbd_convolve_2d_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int subpel_x_qn,
    const int subpel_y_qn, ConvolveParams *conv_params, int bd);

class HighbdConvolve2DTest : public ConvolveTest<Convolve2DParam> {
 public:
  void RunTest(highbd_convolve_2d_func test_func) {
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(
            GetParam().h_filter, GetParam().block_size.width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(
            GetParam().v_filter, GetParam().block_size.height);
    BlockSize b = GetParam().block_size;
    const int sub_x = GetParam().sub_x;
    const int sub_y = GetParam().sub_y;
    const int bit_depth = GetParam().bit_depth;
    const uint16_t *input = RandomInput16(b.width, b.height, bit_depth);
    uint16_t *reference = Output16(b.width, b.height);
    ConvolveParams conv_params1 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    av1_highbd_convolve_2d_sr(input, b.width, reference, b.width, b.width,
                              b.height, filter_params_x, filter_params_y, sub_x,
                              sub_y, &conv_params1, bit_depth);
    uint16_t *test = Output16(b.width, b.height);
    ConvolveParams conv_params2 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    test_func(input, b.width, test, b.width, b.width, b.height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params2, bit_depth);
    AssertEq(reference, test, b.width, b.height, b.width);
  }
};

TEST_P(HighbdConvolve2DTest, C) { RunTest(av1_highbd_convolve_2d_sr_c); }

#if HAVE_SSSE3
TEST_P(HighbdConvolve2DTest, SSSE3) {
  RunTest(av1_highbd_convolve_2d_sr_ssse3);
}
#endif

#if HAVE_AVX2
TEST_P(HighbdConvolve2DTest, AVX2) { RunTest(av1_highbd_convolve_2d_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(HighbdConvolve2DTest, NEON) { RunTest(av1_highbd_convolve_2d_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdConvolve2DTest,
                        GenerateConvolve2DParams({ 10, 12 }));

//////////////////////////
// Compound Convolve Tests
//////////////////////////

// Possible variations in the distance weighting settings.
struct DistWtdSetting {
  bool use_dist_wtd_comp_avg;
  int fwd_offset;
  int bck_offset;
};

std::vector<DistWtdSetting> GetDistWtdSettings() {
  std::vector<DistWtdSetting> result;
  DistWtdSetting setting;
  for (int k = 0; k < 2; ++k) {
    for (int l = 0; l < 4; ++l) {
      setting.use_dist_wtd_comp_avg = true;
      setting.fwd_offset = quant_dist_lookup_table[k][l][0];
      setting.bck_offset = quant_dist_lookup_table[k][l][1];
      result.push_back(setting);
      setting.use_dist_wtd_comp_avg = false;
      result.push_back(setting);
    }
  }
  return result;
}

// Base structure for compound 1d-convolve parameters.
struct Compound1DParam {
  BlockSize block_size;
  DistWtdSetting dist_wtd_setting;
  int sub_pixel;
  InterpFilter filter;
  int bit_depth;
};

::testing::internal::ParamGenerator<Compound1DParam> GenerateCompound1DParams(
    std::initializer_list<int> bit_depths) {
  std::vector<Compound1DParam> result;
  Compound1DParam p;
  for (const auto &block_size : GetBlockSizes()) {
    p.block_size = block_size;
    for (const auto &setting : GetDistWtdSettings()) {
      p.dist_wtd_setting = setting;
      for (int i = 0; i < 16; ++i) {
        p.sub_pixel = i;
        for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
          p.filter = static_cast<InterpFilter>(f);
          for (const int bit_depth : bit_depths) {
            p.bit_depth = bit_depth;
            result.push_back(p);
          }
        }
      }
    }
  }
  return ::testing::ValuesIn(result);
}

////////////////////////////////////////////////
// Compound convolve-x functions (low bit-depth)
////////////////////////////////////////////////
typedef void (*compound_conv_1d_func)(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride, int w,
                                      int h,
                                      const InterpFilterParams *filter_params,
                                      const int subpel_qn,
                                      ConvolveParams *conv_params);

class LowbdCompoundConvolve1DTest : public ConvolveTest<Compound1DParam> {
 public:
  virtual ~LowbdCompoundConvolve1DTest() {}
  void RunTest(compound_conv_1d_func test_func) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;

    const uint8_t *input1 = RandomInput8(width, height);
    const uint8_t *input2 = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(ReferenceFunc(), input1, input2, reference, reference_conv_buf);

    uint8_t *test = Output8(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func, input1, input2, test, test_conv_buf);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 protected:
  virtual compound_conv_1d_func ReferenceFunc() const = 0;
  virtual const InterpFilterParams *FilterParams() const = 0;

 private:
  void Convolve(compound_conv_1d_func test_func, const uint8_t *src1,
                const uint8_t *src2, uint8_t *dst, uint16_t *conv_buf) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const int sub_pix = GetParam().sub_pixel;
    const InterpFilterParams *filter_params = FilterParams();

    const DistWtdSetting setting = GetParam().dist_wtd_setting;
    ConvolveParams conv_params =
        get_conv_params_no_round(0, 0, conv_buf, width, 1, 8);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src1, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params);

    conv_params = get_conv_params_no_round(1, 0, conv_buf, width, 1, 8);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src2, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params);
  }
};

class LowbdCompoundConvolveXTest : public LowbdCompoundConvolve1DTest {
  compound_conv_1d_func ReferenceFunc() const {
    return av1_dist_wtd_convolve_x;
  }
  const InterpFilterParams *FilterParams() const {
    return av1_get_interp_filter_params_with_block_size(
        GetParam().filter, GetParam().block_size.width);
  }
};

TEST_P(LowbdCompoundConvolveXTest, C) { RunTest(av1_dist_wtd_convolve_x_c); }

#if HAVE_SSE2
TEST_P(LowbdCompoundConvolveXTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_x_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(LowbdCompoundConvolveXTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_x_avx2);
}
#endif

#if HAVE_NEON
TEST_P(LowbdCompoundConvolveXTest, NEON) {
  RunTest(av1_dist_wtd_convolve_x_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdCompoundConvolveXTest,
                        GenerateCompound1DParams({ 8 }));

/////////////////////////////////////////////////
// Compound convolve-x functions (high bit-depth)
/////////////////////////////////////////////////
typedef void (*highbd_compound_conv_1d_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params, const int subpel_qn,
    ConvolveParams *conv_params, int bd);

class HighbdCompoundConvolve1DTest : public ConvolveTest<Compound1DParam> {
 public:
  virtual ~HighbdCompoundConvolve1DTest() {}

  void RunTest(highbd_compound_conv_1d_func test_func) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const int bit_depth = GetParam().bit_depth;
    const uint16_t *input1 = RandomInput16(width, height, bit_depth);
    const uint16_t *input2 = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(ReferenceFunc(), input1, input2, reference, reference_conv_buf);

    uint16_t *test = Output16(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func, input1, input2, test, test_conv_buf);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 protected:
  virtual highbd_compound_conv_1d_func ReferenceFunc() const = 0;
  virtual const InterpFilterParams *FilterParams() const = 0;

 private:
  void Convolve(highbd_compound_conv_1d_func test_func, const uint16_t *src1,
                const uint16_t *src2, uint16_t *dst, uint16_t *conv_buf) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const int sub_pix = GetParam().sub_pixel;
    const int bit_depth = GetParam().bit_depth;
    const InterpFilterParams *filter_params = FilterParams();

    const DistWtdSetting setting = GetParam().dist_wtd_setting;
    ConvolveParams conv_params =
        get_conv_params_no_round(0, 0, conv_buf, width, 1, bit_depth);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src1, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params, bit_depth);

    conv_params = get_conv_params_no_round(1, 0, conv_buf, width, 1, bit_depth);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src2, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params, bit_depth);
  }
};

class HighbdCompoundConvolveXTest : public HighbdCompoundConvolve1DTest {
  highbd_compound_conv_1d_func ReferenceFunc() const {
    return av1_highbd_dist_wtd_convolve_x;
  }
  const InterpFilterParams *FilterParams() const {
    return av1_get_interp_filter_params_with_block_size(
        GetParam().filter, GetParam().block_size.width);
  }
};

TEST_P(HighbdCompoundConvolveXTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_x_c);
}

#if HAVE_SSE4_1
TEST_P(HighbdCompoundConvolveXTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_x_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(HighbdCompoundConvolveXTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_x_avx2);
}
#endif

#if HAVE_NEON
TEST_P(HighbdCompoundConvolveXTest, NEON) {
  RunTest(av1_highbd_dist_wtd_convolve_x_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdCompoundConvolveXTest,
                        GenerateCompound1DParams({ 10, 12 }));

////////////////////////////////////////////////
// Compound convolve-y functions (low bit-depth)
////////////////////////////////////////////////

class LowbdCompoundConvolveYTest : public LowbdCompoundConvolve1DTest {
  compound_conv_1d_func ReferenceFunc() const {
    return av1_dist_wtd_convolve_y;
  }
  const InterpFilterParams *FilterParams() const {
    return av1_get_interp_filter_params_with_block_size(
        GetParam().filter, GetParam().block_size.height);
  }
};

TEST_P(LowbdCompoundConvolveYTest, C) { RunTest(av1_dist_wtd_convolve_y_c); }

#if HAVE_SSE2
TEST_P(LowbdCompoundConvolveYTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_y_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(LowbdCompoundConvolveYTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_y_avx2);
}
#endif

#if HAVE_NEON
TEST_P(LowbdCompoundConvolveYTest, NEON) {
  RunTest(av1_dist_wtd_convolve_y_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdCompoundConvolveYTest,
                        GenerateCompound1DParams({ 8 }));

/////////////////////////////////////////////////
// Compound convolve-y functions (high bit-depth)
/////////////////////////////////////////////////

class HighbdCompoundConvolveYTest : public HighbdCompoundConvolve1DTest {
  highbd_compound_conv_1d_func ReferenceFunc() const {
    return av1_highbd_dist_wtd_convolve_y;
  }
  const InterpFilterParams *FilterParams() const {
    return av1_get_interp_filter_params_with_block_size(
        GetParam().filter, GetParam().block_size.height);
  }
};

TEST_P(HighbdCompoundConvolveYTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_y_c);
}

#if HAVE_SSE4_1
TEST_P(HighbdCompoundConvolveYTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_y_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(HighbdCompoundConvolveYTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_y_avx2);
}
#endif

#if HAVE_NEON
TEST_P(HighbdCompoundConvolveYTest, NEON) {
  RunTest(av1_highbd_dist_wtd_convolve_y_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdCompoundConvolveYTest,
                        GenerateCompound1DParams({ 10, 12 }));

// Base structure for compound 2d-copy parameters.
struct Compound2DCopyParam {
  BlockSize block_size;
  DistWtdSetting dist_wtd_setting;
  int bit_depth;
};

::testing::internal::ParamGenerator<Compound2DCopyParam>
GenerateCompound2DCopyParams(std::initializer_list<int> bit_depths) {
  std::vector<Compound2DCopyParam> result;
  Compound2DCopyParam p;
  for (const int bit_depth : bit_depths) {
    p.bit_depth = bit_depth;
    for (const auto &block_size : GetBlockSizes()) {
      // Note that for 2d-copy, the width / height must be multiples of 4 due to
      // implementation details.
      if (block_size.width % 4 != 0 || block_size.height % 4 != 0) {
        continue;
      }
      p.block_size = block_size;
      for (const auto &setting : GetDistWtdSettings()) {
        p.dist_wtd_setting = setting;
        result.push_back(p);
      }
    }
  }
  return ::testing::ValuesIn(result);
}

//////////////////////////////////////////////////////
// Compound convolve-2d-copy functions (low bit-depth)
//////////////////////////////////////////////////////
typedef void (*compound_conv_2d_copy_func)(const uint8_t *src, int src_stride,
                                           uint8_t *dst, int dst_stride, int w,
                                           int h, ConvolveParams *conv_params);

class LowbdCompoundConvolve2DCopyTest
    : public ConvolveTest<Compound2DCopyParam> {
 public:
  void RunTest(compound_conv_2d_copy_func test_func) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;

    const uint8_t *input1 = RandomInput8(width, height);
    const uint8_t *input2 = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_dist_wtd_convolve_2d_copy, input1, input2, reference,
             reference_conv_buf);

    uint8_t *test = Output8(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func, input1, input2, test, test_conv_buf);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(compound_conv_2d_copy_func test_func, const uint8_t *src1,
                const uint8_t *src2, uint8_t *dst, uint16_t *conv_buf) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const DistWtdSetting setting = GetParam().dist_wtd_setting;
    ConvolveParams conv_params =
        get_conv_params_no_round(0, 0, conv_buf, width, 1, 8);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src1, width, dst, width, width, height, &conv_params);

    conv_params = get_conv_params_no_round(1, 0, conv_buf, width, 1, 8);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src2, width, dst, width, width, height, &conv_params);
  }
};

TEST_P(LowbdCompoundConvolve2DCopyTest, C) {
  RunTest(av1_dist_wtd_convolve_2d_copy_c);
}

#if HAVE_SSE2
TEST_P(LowbdCompoundConvolve2DCopyTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_2d_copy_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(LowbdCompoundConvolve2DCopyTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_2d_copy_avx2);
}
#endif

#if HAVE_NEON
TEST_P(LowbdCompoundConvolve2DCopyTest, NEON) {
  RunTest(av1_dist_wtd_convolve_2d_copy_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdCompoundConvolve2DCopyTest,
                        GenerateCompound2DCopyParams({ 8 }));

///////////////////////////////////////////////////////
// Compound convolve-2d-copy functions (high bit-depth)
///////////////////////////////////////////////////////
typedef void (*highbd_compound_conv_2d_copy_func)(const uint16_t *src,
                                                  int src_stride, uint16_t *dst,
                                                  int dst_stride, int w, int h,
                                                  ConvolveParams *conv_params,
                                                  int bd);

class HighbdCompoundConvolve2DCopyTest
    : public ConvolveTest<Compound2DCopyParam> {
 public:
  void RunTest(highbd_compound_conv_2d_copy_func test_func) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const int bit_depth = GetParam().bit_depth;

    const uint16_t *input1 = RandomInput16(width, height, bit_depth);
    const uint16_t *input2 = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_highbd_dist_wtd_convolve_2d_copy, input1, input2, reference,
             reference_conv_buf);

    uint16_t *test = Output16(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func, input1, input2, test, test_conv_buf);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(highbd_compound_conv_2d_copy_func test_func,
                const uint16_t *src1, const uint16_t *src2, uint16_t *dst,
                uint16_t *conv_buf) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const int bit_depth = GetParam().bit_depth;
    const DistWtdSetting setting = GetParam().dist_wtd_setting;
    ConvolveParams conv_params =
        get_conv_params_no_round(0, 0, conv_buf, width, 1, bit_depth);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src1, width, dst, width, width, height, &conv_params, bit_depth);

    conv_params = get_conv_params_no_round(1, 0, conv_buf, width, 1, bit_depth);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src2, width, dst, width, width, height, &conv_params, bit_depth);
  }
};

TEST_P(HighbdCompoundConvolve2DCopyTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_c);
}

#if HAVE_SSE4_1
TEST_P(HighbdCompoundConvolve2DCopyTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(HighbdCompoundConvolve2DCopyTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_avx2);
}
#endif

#if HAVE_NEON
TEST_P(HighbdCompoundConvolve2DCopyTest, NEON) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdCompoundConvolve2DCopyTest,
                        GenerateCompound2DCopyParams({ 10, 12 }));

// Base structure for compound 2d-convolve parameters.
struct Compound2DParam {
  BlockSize block_size;
  DistWtdSetting dist_wtd_setting;
  InterpFilter h_filter;
  InterpFilter v_filter;
  int sub_x;
  int sub_y;
  int bit_depth;
};

::testing::internal::ParamGenerator<Compound2DParam> GenerateCompound2DParams(
    std::initializer_list<int> bit_depths) {
  std::vector<Compound2DParam> result;
  Compound2DParam p;
  for (const auto &block_size : GetBlockSizes()) {
    p.block_size = block_size;
    for (const auto &setting : GetDistWtdSettings()) {
      p.dist_wtd_setting = setting;
      for (int h_f = EIGHTTAP_REGULAR; h_f < INTERP_FILTERS_ALL; ++h_f) {
        p.h_filter = static_cast<InterpFilter>(h_f);
        for (int v_f = EIGHTTAP_REGULAR; v_f < INTERP_FILTERS_ALL; ++v_f) {
          p.v_filter = static_cast<InterpFilter>(v_f);
          for (int x = 0; x < 16; ++x) {
            p.sub_x = x;
            for (int y = 0; y < 16; ++y) {
              p.sub_y = y;
              for (const int bit_depth : bit_depths) {
                p.bit_depth = bit_depth;
                result.push_back(p);
              }
            }
          }
        }
      }
    }
  }
  return ::testing::ValuesIn(result);
}

/////////////////////////////////////////////////
// Compound convolve-2d functions (low bit-depth)
/////////////////////////////////////////////////
typedef void (*compound_conv_2d_func)(
    const uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int subpel_x_qn,
    const int subpel_y_qn, ConvolveParams *conv_params);

class LowbdCompoundConvolve2DTest : public ConvolveTest<Compound2DParam> {
 public:
  void RunTest(compound_conv_2d_func test_func) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;

    const uint8_t *input1 = RandomInput8(width, height);
    const uint8_t *input2 = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_dist_wtd_convolve_2d, input1, input2, reference,
             reference_conv_buf);

    uint8_t *test = Output8(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func, input1, input2, test, test_conv_buf);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(compound_conv_2d_func test_func, const uint8_t *src1,
                const uint8_t *src2, uint8_t *dst, uint16_t *conv_buf) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(GetParam().h_filter,
                                                     width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(GetParam().v_filter,
                                                     height);
    const int sub_x = GetParam().sub_x;
    const int sub_y = GetParam().sub_y;
    const DistWtdSetting setting = GetParam().dist_wtd_setting;
    ConvolveParams conv_params =
        get_conv_params_no_round(0, 0, conv_buf, width, 1, 8);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src1, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params);

    conv_params = get_conv_params_no_round(1, 0, conv_buf, width, 1, 8);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src2, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params);
  }
};

TEST_P(LowbdCompoundConvolve2DTest, C) { RunTest(av1_dist_wtd_convolve_2d_c); }

#if HAVE_SSE2
TEST_P(LowbdCompoundConvolve2DTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_2d_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(LowbdCompoundConvolve2DTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_2d_avx2);
}
#endif

#if HAVE_NEON
TEST_P(LowbdCompoundConvolve2DTest, NEON) {
  RunTest(av1_dist_wtd_convolve_2d_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, LowbdCompoundConvolve2DTest,
                        GenerateCompound2DParams({ 8 }));

//////////////////////////////////////////////////
// Compound convolve-2d functions (high bit-depth)
//////////////////////////////////////////////////
typedef void (*highbd_compound_conv_2d_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int subpel_x_qn,
    const int subpel_y_qn, ConvolveParams *conv_params, int bd);

class HighbdCompoundConvolve2DTest : public ConvolveTest<Compound2DParam> {
 public:
  void RunTest(highbd_compound_conv_2d_func test_func) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const int bit_depth = GetParam().bit_depth;
    const uint16_t *input1 = RandomInput16(width, height, bit_depth);
    const uint16_t *input2 = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_highbd_dist_wtd_convolve_2d, input1, input2, reference,
             reference_conv_buf);

    uint16_t *test = Output16(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func, input1, input2, test, test_conv_buf);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(highbd_compound_conv_2d_func test_func, const uint16_t *src1,
                const uint16_t *src2, uint16_t *dst, uint16_t *conv_buf) {
    const int width = GetParam().block_size.width;
    const int height = GetParam().block_size.height;
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(GetParam().h_filter,
                                                     width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(GetParam().v_filter,
                                                     height);
    const int bit_depth = GetParam().bit_depth;
    const int sub_x = GetParam().sub_x;
    const int sub_y = GetParam().sub_y;
    const DistWtdSetting setting = GetParam().dist_wtd_setting;
    ConvolveParams conv_params =
        get_conv_params_no_round(0, 0, conv_buf, width, 1, bit_depth);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src1, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params, bit_depth);

    conv_params = get_conv_params_no_round(1, 0, conv_buf, width, 1, bit_depth);
    conv_params.use_dist_wtd_comp_avg = setting.use_dist_wtd_comp_avg;
    conv_params.fwd_offset = setting.fwd_offset;
    conv_params.bck_offset = setting.bck_offset;
    conv_params.dst = conv_buf;
    conv_params.dst_stride = width;

    test_func(src2, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params, bit_depth);
  }
};

TEST_P(HighbdCompoundConvolve2DTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_c);
}

#if HAVE_SSE4_1
TEST_P(HighbdCompoundConvolve2DTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(HighbdCompoundConvolve2DTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_avx2);
}
#endif

#if HAVE_NEON
TEST_P(HighbdCompoundConvolve2DTest, NEON) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(Convolve, HighbdCompoundConvolve2DTest,
                        GenerateCompound2DParams({ 10, 12 }));

}  // namespace
