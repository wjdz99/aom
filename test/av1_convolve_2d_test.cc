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

#include <set>
#include "config/av1_rtcd.h"
#include "config/aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

// Generate the list of all block widths / heights that need to be tested.
// We cannot use the BLOCK_SIZE enum because we need to simulate chroma block
// sizes, which can go down to 2xN and Nx2.
struct BlockSize {
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

// AV1ConvolveTest is the base class that all other tests should derive from.
// It is templated on the input type (e.g., uint8_t or uint16_t) and provides
// methods to randomize the input/output buffers before a run. Tests should
// call the RunCheckOutput, which runs a reference convolve function (e.g., the
// C implementation) against a test convolve function (e.g., an intrinsics
// implementation) and checks that the output is identical.
template <typename T, typename Param>
class AV1ConvolveTest : public ::testing::TestWithParam<Param> {
 public:
  virtual ~AV1ConvolveTest() {}

  virtual void SetUp() override {
    rnd_.Reset(libaom_test::ACMRandom::DeterministicSeed());
    RandomizeMemory(input_, kMaxSize * kMaxSize);
    RandomizeMemory(reference_output_, MAX_SB_SQUARE);
    RandomizeMemory(test_output_, MAX_SB_SQUARE);
  }

  virtual void TearDown() override { libaom_test::ClearSystemState(); }

  // Returns the pointer to a random area in the input. Guaranteed
  // to have enough space for an 8-tap filter applied to a block with
  // the given dimensions.
  const T *RandomInputArea(int w, int h) {
    const int offset_r = 3 + rnd_.PseudoUniform(kMaxSize - h - 7);
    const int offset_c = 3 + rnd_.PseudoUniform(kMaxSize - w - 7);
    return input_ + offset_r * kMaxSize + offset_c;
  }

  void RunCheckOutput() {
    // Check that the bit-depth / template type are reasonable.
    ASSERT_TRUE(BitDepth() == 8 || BitDepth() == 10 || BitDepth() == 12);
    ASSERT_EQ(sizeof(T) == sizeof(uint8_t), BitDepth() == 8);
    ASSERT_TRUE(sizeof(T) <= sizeof(uint16_t));
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    // Try 10 different random areas.
    for (int run = 0; run < 10; ++run) {
      const T *input = RandomInputArea(block_w, block_h);
      ReferenceImpl(input, kMaxSize, ReferenceOutput(), MAX_SB_SIZE);
      TestImpl(input, kMaxSize, TestOutput(), MAX_SB_SIZE);

      // Check that the test output and reference output match.
      for (int j = 0; j < block_h; ++j) {
        for (int i = 0; i < block_w; ++i) {
          int idx = j * MAX_SB_SIZE + i;
          ASSERT_EQ(ReferenceOutput()[idx], TestOutput()[idx])
              << block_w << "x" << block_h << " Pixel mismatch at index " << idx
              << " = (" << i << ", " << j << ")";
        }
      }
    }
  }

 protected:
  // When putting random values in the input buffer, they will be capped by
  // to fit within this many bits.
  virtual int BitDepth() const = 0;
  virtual int BlockWidth() const = 0;
  virtual int BlockHeight() const = 0;
  virtual void ReferenceImpl(const T *src, int src_stride, T *dst,
                             int dst_stride) = 0;
  virtual void TestImpl(const T *src, int src_stride, T *dst,
                        int dst_stride) = 0;

 private:
  static constexpr int kMaxSize = 128 + 32;  // padding

  // Get a reference to the test output buffer and reference output buffer,
  // when calling the implementations.
  T *TestOutput() { return test_output_; }
  T *ReferenceOutput() { return reference_output_; }

  // Put random values in the memory.
  void RandomizeMemory(T *p, int size) {
    for (int i = 0; i < size; ++i) {
      p[i] = rnd_.Rand31() & ((1 << BitDepth()) - 1);
    }
  }

  T input_[kMaxSize * kMaxSize];
  DECLARE_ALIGNED(32, T, reference_output_[MAX_SB_SQUARE]);
  DECLARE_ALIGNED(32, T, test_output_[MAX_SB_SQUARE]);
  libaom_test::ACMRandom rnd_;
};

////////////////////////////////////////////////////////
// Single reference convolve-x functions (low bit-depth)
////////////////////////////////////////////////////////
typedef void (*convolve_x_func)(const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int w, int h,
                                const InterpFilterParams *filter_params_x,
                                const int subpel_x_qn,
                                ConvolveParams *conv_params);

struct ConvolveXParam {
  ConvolveXParam(convolve_x_func func, BlockSize block,
                 InterpFilter horizontal_filter, int sub_x)
      : func(func), block(block), horizontal_filter(horizontal_filter),
        sub_x(sub_x) {}
  convolve_x_func func;
  BlockSize block;
  InterpFilter horizontal_filter;
  int sub_x;
};

class AV1ConvolveXTest : public AV1ConvolveTest<uint8_t, ConvolveXParam> {
 public:
  virtual ~AV1ConvolveXTest() {}

 protected:
  int BitDepth() const override { return 8; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }
  void ReferenceImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_convolve_x_sr_c);
  }

  void TestImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }

  void RunConvolve(const uint8_t *src, int src_stride, uint8_t *dst,
                   int dst_stride, convolve_x_func func) {
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    InterpFilter f = GetParam().horizontal_filter;
    int sub_x = GetParam().sub_x;
    ConvolveParams conv_params = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    func(src, src_stride, dst, dst_stride, block_w, block_h,
         av1_get_interp_filter_params_with_block_size(f, block_w), sub_x,
         &conv_params);
  }
};

::testing::internal::ParamGenerator<ConvolveXParam> GetConvolveXParams(
    convolve_x_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<ConvolveXParam> result;
  for (const BlockSize &b : sizes) {
    for (int sub_x = 0; sub_x < 16; ++sub_x) {
      for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
        result.push_back(
            ConvolveXParam(func, b, static_cast<InterpFilter>(f), sub_x));
      }
    }
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1ConvolveXTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(C_X, AV1ConvolveXTest,
                        GetConvolveXParams(av1_convolve_x_sr_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2_X, AV1ConvolveXTest,
                        GetConvolveXParams(av1_convolve_x_sr_sse2));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(AVX2_X, AV1ConvolveXTest,
                        GetConvolveXParams(av1_convolve_x_sr_avx2));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON_X, AV1ConvolveXTest,
                        GetConvolveXParams(av1_convolve_x_sr_neon));
#endif

/////////////////////////////////////////////////////////
// Single reference convolve-x functions (high bit-depth)
/////////////////////////////////////////////////////////
typedef void (*highbd_convolve_x_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn,
    ConvolveParams *conv_params, int bd);

struct HighbdConvolveXParam {
  HighbdConvolveXParam(highbd_convolve_x_func func, BlockSize block,
                       InterpFilter horizontal_filter, int sub_x, int bd)
      : func(func), block(block), horizontal_filter(horizontal_filter),
        sub_x(sub_x), bitdepth(bd) {}
  highbd_convolve_x_func func;
  BlockSize block;
  InterpFilter horizontal_filter;
  int sub_x;
  int bitdepth;
};

class AV1HighbdConvolveXTest
    : public AV1ConvolveTest<uint16_t, HighbdConvolveXParam> {
 public:
  virtual ~AV1HighbdConvolveXTest() {}

 protected:
  int BitDepth() const override { return GetParam().bitdepth; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }
  void ReferenceImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_highbd_convolve_x_sr_c);
  }

  void TestImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }

  void RunConvolve(const uint16_t *src, int src_stride, uint16_t *dst,
                   int dst_stride, highbd_convolve_x_func func) {
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    InterpFilter f = GetParam().horizontal_filter;
    int sub_x = GetParam().sub_x;
    ConvolveParams conv_params = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    func(src, src_stride, dst, dst_stride, block_w, block_h,
         av1_get_interp_filter_params_with_block_size(f, block_w), sub_x,
         &conv_params, BitDepth());
  }
};

::testing::internal::ParamGenerator<HighbdConvolveXParam>
GetHighbdConvolveXParams(highbd_convolve_x_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<HighbdConvolveXParam> result;
  for (const BlockSize &b : sizes) {
    for (int sub_x = 0; sub_x < 16; ++sub_x) {
      for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
        result.push_back(HighbdConvolveXParam(
            func, b, static_cast<InterpFilter>(f), sub_x, 10));
        result.push_back(HighbdConvolveXParam(
            func, b, static_cast<InterpFilter>(f), sub_x, 12));
      }
    }
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1HighbdConvolveXTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(HIGHBD_C_X, AV1HighbdConvolveXTest,
                        GetHighbdConvolveXParams(av1_highbd_convolve_x_sr_c));

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(
    HIGHBD_SSSE3_X, AV1HighbdConvolveXTest,
    GetHighbdConvolveXParams(av1_highbd_convolve_x_sr_ssse3));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(
    HIGHBD_AVX2_X, AV1HighbdConvolveXTest,
    GetHighbdConvolveXParams(av1_highbd_convolve_x_sr_avx2));
#endif

////////////////////////////////////////////////////////
// Single reference convolve-y functions (low bit-depth)
////////////////////////////////////////////////////////

typedef void (*convolve_y_func)(const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int w, int h,
                                const InterpFilterParams *filter_params_y,
                                const int subpel_y_qn);
struct ConvolveYParam {
  ConvolveYParam(convolve_y_func func, BlockSize block,
                 InterpFilter vertical_filter, int sub_y)
      : func(func), block(block), vertical_filter(vertical_filter),
        sub_y(sub_y) {}
  convolve_y_func func;
  BlockSize block;
  InterpFilter vertical_filter;
  int sub_y;
};

class AV1ConvolveYTest : public AV1ConvolveTest<uint8_t, ConvolveYParam> {
 public:
  virtual ~AV1ConvolveYTest() {}

 protected:
  int BitDepth() const override { return 8; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }

  void RunConvolve(const uint8_t *src, int src_stride, uint8_t *dst,
                   int dst_stride, convolve_y_func func) {
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    InterpFilter f = GetParam().vertical_filter;
    int sub_y = GetParam().sub_y;
    func(src, src_stride, dst, dst_stride, block_w, block_h,
         av1_get_interp_filter_params_with_block_size(f, block_w), sub_y);
  }

  void ReferenceImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_convolve_y_sr_c);
  }

  void TestImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }
};

::testing::internal::ParamGenerator<ConvolveYParam> GetConvolveYParams(
    convolve_y_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<ConvolveYParam> result;
  for (const BlockSize &b : sizes) {
    for (int sub_y = 0; sub_y < 16; ++sub_y) {
      for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
        result.push_back(
            ConvolveYParam(func, b, static_cast<InterpFilter>(f), sub_y));
      }
    }
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1ConvolveYTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(C_Y, AV1ConvolveYTest,
                        GetConvolveYParams(av1_convolve_y_sr_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2_Y, AV1ConvolveYTest,
                        GetConvolveYParams(av1_convolve_y_sr_sse2));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(AVX2_Y, AV1ConvolveYTest,
                        GetConvolveYParams(av1_convolve_y_sr_avx2));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON_Y, AV1ConvolveYTest,
                        GetConvolveYParams(av1_convolve_y_sr_neon));
#endif

/////////////////////////////////////////////////////////
// Single reference convolve-y functions (high bit-depth)
/////////////////////////////////////////////////////////
typedef void (*highbd_convolve_y_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn,
    int bd);

struct HighbdConvolveYParam {
  HighbdConvolveYParam(highbd_convolve_y_func func, BlockSize block,
                       InterpFilter vertical_filter, int sub_y, int bd)
      : func(func), block(block), vertical_filter(vertical_filter),
        sub_y(sub_y), bitdepth(bd) {}
  highbd_convolve_y_func func;
  BlockSize block;
  InterpFilter vertical_filter;
  int sub_y;
  int bitdepth;
};

class AV1HighbdConvolveYTest
    : public AV1ConvolveTest<uint16_t, HighbdConvolveYParam> {
 public:
  virtual ~AV1HighbdConvolveYTest() {}

 protected:
  int BitDepth() const override { return GetParam().bitdepth; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }
  void ReferenceImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_highbd_convolve_y_sr_c);
  }

  void TestImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }

  void RunConvolve(const uint16_t *src, int src_stride, uint16_t *dst,
                   int dst_stride, highbd_convolve_y_func func) {
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    InterpFilter f = GetParam().vertical_filter;
    int sub_y = GetParam().sub_y;
    func(src, src_stride, dst, dst_stride, block_w, block_h,
         av1_get_interp_filter_params_with_block_size(f, block_w), sub_y,
         BitDepth());
  }
};

::testing::internal::ParamGenerator<HighbdConvolveYParam>
GetHighbdConvolveYParams(highbd_convolve_y_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<HighbdConvolveYParam> result;
  for (const BlockSize &b : sizes) {
    for (int sub_x = 0; sub_x < 16; ++sub_x) {
      for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
        result.push_back(HighbdConvolveYParam(
            func, b, static_cast<InterpFilter>(f), sub_x, 10));
        result.push_back(HighbdConvolveYParam(
            func, b, static_cast<InterpFilter>(f), sub_x, 12));
      }
    }
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1HighbdConvolveYTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(HIGHBD_C_Y, AV1HighbdConvolveYTest,
                        GetHighbdConvolveYParams(av1_highbd_convolve_y_sr_c));

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(
    HIGHBD_SSSE3_Y, AV1HighbdConvolveYTest,
    GetHighbdConvolveYParams(av1_highbd_convolve_y_sr_ssse3));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(
    HIGHBD_AVX2_Y, AV1HighbdConvolveYTest,
    GetHighbdConvolveYParams(av1_highbd_convolve_y_sr_avx2));
#endif

//////////////////////////////////////////////////////////////
// Single reference convolve-2d-copy functions (low bit-depth)
//////////////////////////////////////////////////////////////

typedef void (*convolve_2d_copy_func)(const uint8_t *src, int src_stride,
                                      uint8_t *dst, int dst_stride, int w,
                                      int h);

struct Convolve2DCopyParam {
  Convolve2DCopyParam(convolve_2d_copy_func func, BlockSize block)
      : func(func), block(block) {}
  convolve_2d_copy_func func;
  BlockSize block;
};

class AV1Convolve2DCopyTest
    : public AV1ConvolveTest<uint8_t, Convolve2DCopyParam> {
 public:
  virtual ~AV1Convolve2DCopyTest() {}

 protected:
  int BitDepth() const override { return 8; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }

  void RunConvolve(const uint8_t *src, int src_stride, uint8_t *dst,
                   int dst_stride, convolve_2d_copy_func func) {
    func(src, src_stride, dst, dst_stride, BlockWidth(), BlockHeight());
  }

  void ReferenceImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_convolve_2d_copy_sr_c);
  }

  void TestImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }
};

::testing::internal::ParamGenerator<Convolve2DCopyParam>
GetConvolve2DCopyParams(convolve_2d_copy_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<Convolve2DCopyParam> result;
  for (const BlockSize &b : sizes) {
    result.push_back(Convolve2DCopyParam(func, b));
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1Convolve2DCopyTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(C_2D_COPY, AV1Convolve2DCopyTest,
                        GetConvolve2DCopyParams(av1_convolve_2d_copy_sr_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2_2D_COPY, AV1Convolve2DCopyTest,
                        GetConvolve2DCopyParams(av1_convolve_2d_copy_sr_sse2));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(AVX2_2D_COPY, AV1Convolve2DCopyTest,
                        GetConvolve2DCopyParams(av1_convolve_2d_copy_sr_avx2));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON_2D_COPY, AV1Convolve2DCopyTest,
                        GetConvolve2DCopyParams(av1_convolve_2d_copy_sr_neon));
#endif

///////////////////////////////////////////////////////////////
// Single reference convolve-2d-copy functions (high bit-depth)
///////////////////////////////////////////////////////////////

typedef void (*highbd_convolve_2d_copy_func)(const uint16_t *src,
                                             int src_stride, uint16_t *dst,
                                             int dst_stride, int w, int h);

struct HighbdConvolve2DCopyParam {
  HighbdConvolve2DCopyParam(highbd_convolve_2d_copy_func func, BlockSize block,
                            int bd)
      : func(func), block(block), bitdepth(bd) {}
  highbd_convolve_2d_copy_func func;
  BlockSize block;
  int bitdepth;
};

class AV1HighbdConvolve2DCopyTest
    : public AV1ConvolveTest<uint16_t, HighbdConvolve2DCopyParam> {
 public:
  virtual ~AV1HighbdConvolve2DCopyTest() {}

 protected:
  int BitDepth() const override { return GetParam().bitdepth; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }

  void RunConvolve(const uint16_t *src, int src_stride, uint16_t *dst,
                   int dst_stride, highbd_convolve_2d_copy_func func) {
    func(src, src_stride, dst, dst_stride, BlockWidth(), BlockHeight());
  }

  void ReferenceImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride,
                av1_highbd_convolve_2d_copy_sr_c);
  }

  void TestImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }
};

::testing::internal::ParamGenerator<HighbdConvolve2DCopyParam>
GetHighbdConvolve2DCopyParams(highbd_convolve_2d_copy_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<HighbdConvolve2DCopyParam> result;
  for (const BlockSize &b : sizes) {
    result.push_back(HighbdConvolve2DCopyParam(func, b, 10));
    result.push_back(HighbdConvolve2DCopyParam(func, b, 12));
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1HighbdConvolve2DCopyTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(
    HIGHBD_C_2D_COPY, AV1HighbdConvolve2DCopyTest,
    GetHighbdConvolve2DCopyParams(av1_highbd_convolve_2d_copy_sr_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(
    HIGHBD_SSE2_2D_COPY, AV1HighbdConvolve2DCopyTest,
    GetHighbdConvolve2DCopyParams(av1_highbd_convolve_2d_copy_sr_sse2));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(
    HIGHBD_AVX2_2D_COPY, AV1HighbdConvolve2DCopyTest,
    GetHighbdConvolve2DCopyParams(av1_highbd_convolve_2d_copy_sr_avx2));
#endif

/////////////////////////////////////////////////////////
// Single reference convolve-2d functions (low bit-depth)
/////////////////////////////////////////////////////////

typedef void (*convolve_2d_func)(const uint8_t *src, int src_stride,
                                 uint8_t *dst, int dst_stride, int w, int h,
                                 const InterpFilterParams *filter_params_x,
                                 const InterpFilterParams *filter_params_y,
                                 const int subpel_x_qn, const int subpel_y_qn,
                                 ConvolveParams *conv_params);

struct Convolve2DParam {
  Convolve2DParam(convolve_2d_func func, BlockSize block,
                  InterpFilter horizontal_filter, InterpFilter vertical_filter,
                  int sub_x, int sub_y)
      : func(func), block(block), horizontal_filter(horizontal_filter),
        vertical_filter(vertical_filter), sub_x(sub_x), sub_y(sub_y) {}
  convolve_2d_func func;
  BlockSize block;
  InterpFilter horizontal_filter;
  InterpFilter vertical_filter;
  int sub_x;
  int sub_y;
};

class AV1Convolve2DTest : public AV1ConvolveTest<uint8_t, Convolve2DParam> {
 public:
  virtual ~AV1Convolve2DTest() {}

 protected:
  int BitDepth() const override { return 8; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }

  void RunConvolve(const uint8_t *src, int src_stride, uint8_t *dst,
                   int dst_stride, convolve_2d_func func) {
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    InterpFilter horiz_f = GetParam().horizontal_filter;
    InterpFilter vert_f = GetParam().vertical_filter;
    int sub_x = GetParam().sub_x;
    int sub_y = GetParam().sub_y;
    ConvolveParams conv_params = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    func(src, src_stride, dst, dst_stride, block_w, block_h,
         av1_get_interp_filter_params_with_block_size(horiz_f, block_w),
         av1_get_interp_filter_params_with_block_size(vert_f, block_h), sub_x,
         sub_y, &conv_params);
  }

  void ReferenceImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_convolve_2d_sr_c);
  }

  void TestImpl(const uint8_t *src, int src_stride, uint8_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }
};

::testing::internal::ParamGenerator<Convolve2DParam> GetConvolve2DParams(
    convolve_2d_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<Convolve2DParam> result;
  for (const BlockSize &b : sizes) {
    for (int horiz_f = EIGHTTAP_REGULAR; horiz_f < INTERP_FILTERS_ALL;
         ++horiz_f) {
      for (int vert_f = EIGHTTAP_REGULAR; vert_f < INTERP_FILTERS_ALL;
           ++vert_f) {
        for (int sub_x = 0; sub_x < 16; ++sub_x) {
          for (int sub_y = 0; sub_y < 16; ++sub_y) {
            result.push_back(Convolve2DParam(
                func, b, static_cast<InterpFilter>(horiz_f),
                static_cast<InterpFilter>(vert_f), sub_x, sub_y));
          }
        }
      }
    }
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1Convolve2DTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(C_2D, AV1Convolve2DTest,
                        GetConvolve2DParams(av1_convolve_2d_sr_c));

#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2_2D, AV1Convolve2DTest,
                        GetConvolve2DParams(av1_convolve_2d_sr_sse2));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(AVX2_2D, AV1Convolve2DTest,
                        GetConvolve2DParams(av1_convolve_2d_sr_avx2));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_CASE_P(NEON_2D, AV1Convolve2DTest,
                        GetConvolve2DParams(av1_convolve_2d_sr_neon));
#endif

/////////////////////////////////////////////////////////
// Single reference convolve-2d functions (high bit-depth)
/////////////////////////////////////////////////////////

typedef void (*highbd_convolve_2d_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int subpel_x_qn,
    const int subpel_y_qn, ConvolveParams *conv_params, int bd);

struct HighbdConvolve2DParam {
  HighbdConvolve2DParam(highbd_convolve_2d_func func, BlockSize block,
                        InterpFilter horizontal_filter,
                        InterpFilter vertical_filter, int sub_x, int sub_y,
                        int bd)
      : func(func), block(block), horizontal_filter(horizontal_filter),
        vertical_filter(vertical_filter), sub_x(sub_x), sub_y(sub_y),
        bitdepth(bd) {}
  highbd_convolve_2d_func func;
  BlockSize block;
  InterpFilter horizontal_filter;
  InterpFilter vertical_filter;
  int sub_x;
  int sub_y;
  int bitdepth;
};

class AV1HighbdConvolve2DTest
    : public AV1ConvolveTest<uint16_t, HighbdConvolve2DParam> {
 public:
  virtual ~AV1HighbdConvolve2DTest() {}

 protected:
  int BitDepth() const override { return GetParam().bitdepth; }
  int BlockWidth() const override { return GetParam().block.width; }
  int BlockHeight() const override { return GetParam().block.height; }

  void RunConvolve(const uint16_t *src, int src_stride, uint16_t *dst,
                   int dst_stride, highbd_convolve_2d_func func) {
    int block_w = BlockWidth();
    int block_h = BlockHeight();
    InterpFilter horiz_f = GetParam().horizontal_filter;
    InterpFilter vert_f = GetParam().vertical_filter;
    int sub_x = GetParam().sub_x;
    int sub_y = GetParam().sub_y;
    ConvolveParams conv_params = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    func(src, src_stride, dst, dst_stride, block_w, block_h,
         av1_get_interp_filter_params_with_block_size(horiz_f, block_w),
         av1_get_interp_filter_params_with_block_size(vert_f, block_h), sub_x,
         sub_y, &conv_params, BitDepth());
  }

  void ReferenceImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                     int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, av1_highbd_convolve_2d_sr_c);
  }

  void TestImpl(const uint16_t *src, int src_stride, uint16_t *dst,
                int dst_stride) override {
    RunConvolve(src, src_stride, dst, dst_stride, GetParam().func);
  }
};

::testing::internal::ParamGenerator<HighbdConvolve2DParam>
GetHighbdConvolve2DParams(highbd_convolve_2d_func func) {
  std::set<BlockSize> sizes = GetBlockSizes();
  std::vector<HighbdConvolve2DParam> result;
  for (const BlockSize &b : sizes) {
    for (int horiz_f = EIGHTTAP_REGULAR; horiz_f < INTERP_FILTERS_ALL;
         ++horiz_f) {
      for (int vert_f = EIGHTTAP_REGULAR; vert_f < INTERP_FILTERS_ALL;
           ++vert_f) {
        for (int sub_x = 0; sub_x < 16; ++sub_x) {
          for (int sub_y = 0; sub_y < 16; ++sub_y) {
            result.push_back(HighbdConvolve2DParam(
                func, b, static_cast<InterpFilter>(horiz_f),
                static_cast<InterpFilter>(vert_f), sub_x, sub_y, 10));
            result.push_back(HighbdConvolve2DParam(
                func, b, static_cast<InterpFilter>(horiz_f),
                static_cast<InterpFilter>(vert_f), sub_x, sub_y, 12));
          }
        }
      }
    }
  }
  return ::testing::ValuesIn(result);
}

TEST_P(AV1HighbdConvolve2DTest, RunCheckOutput) { RunCheckOutput(); }

INSTANTIATE_TEST_CASE_P(HIGHBD_C_2D, AV1HighbdConvolve2DTest,
                        GetHighbdConvolve2DParams(av1_highbd_convolve_2d_sr_c));

#if HAVE_SSSE3
INSTANTIATE_TEST_CASE_P(
    HIGHBD_SSSE3_2D, AV1HighbdConvolve2DTest,
    GetHighbdConvolve2DParams(av1_highbd_convolve_2d_sr_ssse3));
#endif

#if HAVE_AVX2
INSTANTIATE_TEST_CASE_P(
    AVX2_2D, AV1HighbdConvolve2DTest,
    GetHighbdConvolve2DParams(av1_highbd_convolve_2d_sr_avx2));
#endif

}  // namespace
