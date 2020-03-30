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
#include <vector>
#include "config/av1_rtcd.h"
#include "config/aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

// All convolve tests are parameterized on block size. This can either be
// block sizes that include 2xN and Nx2 (to include chroma sizes), or just
// the block sizes listed in the BLOCK_SIZE enum.
//
// Note that parameterizing on just block size (and not other parameters, such
// as bit-depth, filter parameters, etc.) is a conscious decision - otherwise
// we would have several million parameters to pass in, and the gtest framework
// does not handle this well (increased overhead per test, huge amount of output
// to stdout, etc.). Alternatively, the tests could not be parameterized,
// but that reduces parallelization possibilities.
class BlockSize {
 public:
  BlockSize(int w, int h) : width_(w), height_(h) {}

  int Width() const { return width_; }
  int Height() const { return height_; }

  bool operator<(const BlockSize &other) const {
    if (Width() == other.Width()) {
      return Height() < other.Height();
    }
    return Width() < other.Width();
  }

  bool operator==(const BlockSize &other) const {
    return Width() == other.Width() && Height() == other.Height();
  }

 private:
  int width_;
  int height_;
};

// Generate the list of all block widths / heights that need to be tested,
// includes chroma and luma sizes.
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

// ConvolveTest is the base class that all convolve tests should derive from.
// It provides storage/methods for generating randomized buffers for both
// low bit-depth and high bit-depth. It is templated over the parameters
// that should be tested, which are iterated over in RunTest().
// Implementors should call that method and implement TestConvolve, which
// takes the current testing parameters. Note that BlockSize is available
// via GetParam(), but implementors should ignore this and get all information
// from the templatized parameter.

// Class to iterate over parameters. Implementors should create
// template-specialized versions for their specific parameters. Note that
// the class should be copy-able.
template <typename T>
class ParamIterator {
  // Implementors should specify a constructor:
  //
  // ParamIterator(const BlockSize &b) { ... }
  //
  // And implement two methods:
  //
  // bool HasNext() const { ... }
  // T Next() { ... }
};

// We often need a high bit-depth param iterator as well. This is a wrapper
// class that adds a bit_depth field to a parameter.
template <typename T>
class Highbd {
 public:
  Highbd(const T &param, int bit_depth)
      : param_(param), bit_depth_(bit_depth) {}

  int BitDepth() const { return bit_depth_; }
  const T &Param() const { return param_; }

 private:
  T param_;
  int bit_depth_;
};

template <typename T>
class ParamIterator<Highbd<T>> {
 public:
  explicit ParamIterator(const BlockSize &block) : iter_(block), values_() {}

  bool HasNext() const { return !values_.empty() || iter_.HasNext(); }
  Highbd<T> Next() {
    if (!values_.empty()) {
      auto r = values_.back();
      values_.pop_back();
      return r;
    }
    T param = iter_.Next();
    values_.push_back(Highbd<T>(param, 12));
    return Highbd<T>(param, 10);
  }

 private:
  ParamIterator<T> iter_;
  std::vector<Highbd<T>> values_;
};

template <typename T>
class AV1ConvolveTest : public ::testing::TestWithParam<BlockSize> {
 public:
  virtual ~AV1ConvolveTest() { ClearMemory(); }

  virtual void SetUp() override {
    rnd_.Reset(libaom_test::ACMRandom::DeterministicSeed());
  }

  virtual void TearDown() override {
    libaom_test::ClearSystemState();
    ClearMemory();
  }

  void RunTest() {
    BlockSize block = GetParam();
    ParamIterator<T> iter(block);
    while (iter.HasNext()) {
      TestConvolve(iter.Next());
    }
  }

  // Randomizes the 8-bit input buffer and returns a pointer to it. Note that
  // the pointer is safe to use with an 8-tap filter. The stride can range
  // from width to (width + kPadding).
  static constexpr int kInputPadding = 8;

  const uint8_t *RandomInput8(int width, int height) {
    const int padded_width = width + kInputPadding;
    const int padded_height = height + kInputPadding;
    return RandomUint8(padded_width * padded_height) + 3 * padded_width + 3;
  }

  // Generate a random 16-bit input buffer, like above. Note that the
  // values are capped so they do not exceed the bit-depth.
  const uint16_t *RandomInput16(int width, int height, int bit_depth) {
    const int padded_width = width + kInputPadding;
    const int padded_height = height + kInputPadding;
    return RandomUint16(padded_width * padded_height, bit_depth) +
           3 * padded_width + 3;
  }

  // Some of the intrinsics perform writes in 16 byte chunks. Make sure
  // this padding exists along the width.
  static constexpr int kOutputPadding = 16;

  // 8-bit output buffer of size width * height. It is aligned on a 16-byte
  // boundary and padded by kOutputPadding for the width.
  uint8_t *Output8(int width, int height) {
    const int padded_width = width + kOutputPadding;
    const size_t size = padded_width * height;
    buffer8_.push_back(reinterpret_cast<uint8_t *>(aom_memalign(16, size)));
    return buffer8_.back();
  }

  // 16-bit output buffer of size width * height. It is aligned on a 16-byte
  // boundary and padded by kOutputPadding for the width.
  uint16_t *Output16(int width, int height) {
    const int padded_width = width + kOutputPadding;
    const size_t size = padded_width * height;
    buffer16_.push_back(reinterpret_cast<uint16_t *>(
        aom_memalign(16, sizeof(uint16_t) * size)));
    return buffer16_.back();
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

 protected:
  virtual void TestConvolve(const T &param) = 0;

 private:
  const uint8_t *RandomUint8(int size) {
    unaligned_buffer8_.push_back(new uint8_t[size]);
    uint8_t *p = unaligned_buffer8_.back();
    for (int i = 0; i < size; ++i) {
      p[i] = rnd_.Rand8();
    }
    return p;
  }

  const uint16_t *RandomUint16(int size, int bit_depth) {
    unaligned_buffer16_.push_back(new uint16_t[size]);
    uint16_t *p = unaligned_buffer16_.back();
    for (int i = 0; i < size; ++i) {
      p[i] = rnd_.Rand16() & ((1 << bit_depth) - 1);
    }
    return p;
  }

  void ClearMemory() {
    for (uint8_t *ptr : buffer8_) {
      aom_free(ptr);
    }
    buffer8_.clear();
    for (uint16_t *ptr : buffer16_) {
      aom_free(ptr);
    }
    buffer16_.clear();
    for (uint8_t *ptr : unaligned_buffer8_) {
      delete[] ptr;
    }
    unaligned_buffer8_.clear();
    for (uint16_t *ptr : unaligned_buffer16_) {
      delete[] ptr;
    }
    unaligned_buffer16_.clear();
  }

  // We maintain separate pools for aligned and unaligned memory allocations -
  // aligned memory must be deallocated with aom_free. The two pools are used
  // to test that the implementation can deal with unaligned input, even if
  // output requires alignment.
  std::vector<uint8_t *> buffer8_;
  std::vector<uint16_t *> buffer16_;
  std::vector<uint8_t *> unaligned_buffer8_;
  std::vector<uint16_t *> unaligned_buffer16_;
  libaom_test::ACMRandom rnd_;
};

// Single reference 1D convolution parameters.
class Sr1DParam {
 public:
  Sr1DParam(const BlockSize &b, int s, InterpFilter f)
      : block_(b), sub_pixel_(s), filter_(f) {}

  const BlockSize &Block() const { return block_; }
  int SubPixel() const { return sub_pixel_; }
  InterpFilter Filter() const { return filter_; }

 private:
  BlockSize block_;
  int sub_pixel_;
  InterpFilter filter_;
};

template <>
class ParamIterator<Sr1DParam> {
 public:
  ParamIterator<Sr1DParam>(const BlockSize &b)
      : block_(b), sub_pix_(0), filter_(EIGHTTAP_REGULAR) {}

  bool HasNext() const { return filter_ < INTERP_FILTERS_ALL; }
  Sr1DParam Next() {
    int s = sub_pix_;
    InterpFilter f = static_cast<InterpFilter>(filter_);
    if (sub_pix_ < 15) {
      ++sub_pix_;
    } else {
      sub_pix_ = 0;
      ++filter_;
    }
    return Sr1DParam(block_, s, f);
  }

 private:
  const BlockSize block_;
  int sub_pix_;
  int filter_;
};

// "Iterator" over a block size. For a given block size, the set of blocks is
// just that block.
template <>
class ParamIterator<BlockSize> {
 public:
  ParamIterator<BlockSize>(const BlockSize &b) : block_(b), has_next_(true) {}
  bool HasNext() const { return has_next_; }
  BlockSize Next() {
    has_next_ = false;
    return block_;
  }

 private:
  const BlockSize block_;
  bool has_next_;
};

// Tests for the block sizes and iterators.
class AV1ConvolveIteratorTest : public ::testing::Test {};

TEST_F(AV1ConvolveIteratorTest, GetBlockSizes) {
  auto v = GetBlockSizes();
  ASSERT_EQ(27U, v.size());
}

TEST_F(AV1ConvolveIteratorTest, BlockSizeIterator) {
  BlockSize block(32, 64);
  ParamIterator<BlockSize> iterator(block);
  ASSERT_TRUE(iterator.HasNext());
  ASSERT_EQ(block, iterator.Next());
  ASSERT_FALSE(iterator.HasNext());
}

TEST_F(AV1ConvolveIteratorTest, Sr1DParam) {
  BlockSize block(32, 64);
  ParamIterator<Sr1DParam> iterator(block);
  for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
    for (int i = 0; i < 16; ++i) {
      ASSERT_TRUE(iterator.HasNext());
      Sr1DParam p = iterator.Next();
      ASSERT_EQ(block, p.Block());
      ASSERT_EQ(i, p.SubPixel());
      ASSERT_EQ(static_cast<InterpFilter>(f), p.Filter());
    }
  }
  ASSERT_FALSE(iterator.HasNext());
}

TEST_F(AV1ConvolveIteratorTest, HighbdSr1DParam) {
  BlockSize block(32, 64);
  ParamIterator<Highbd<Sr1DParam>> iterator(block);
  for (int f = EIGHTTAP_REGULAR; f < INTERP_FILTERS_ALL; ++f) {
    for (int i = 0; i < 16; ++i) {
      ASSERT_TRUE(iterator.HasNext());
      Highbd<Sr1DParam> p = iterator.Next();
      ASSERT_EQ(block, p.Param().Block());
      ASSERT_EQ(i, p.Param().SubPixel());
      ASSERT_EQ(static_cast<InterpFilter>(f), p.Param().Filter());
      ASSERT_EQ(10, p.BitDepth());

      ASSERT_TRUE(iterator.HasNext());
      p = iterator.Next();
      ASSERT_EQ(block, p.Param().Block());
      ASSERT_EQ(i, p.Param().SubPixel());
      ASSERT_EQ(static_cast<InterpFilter>(f), p.Param().Filter());
      ASSERT_EQ(12, p.BitDepth());
    }
  }
  ASSERT_FALSE(iterator.HasNext());
}

////////////////////////////////////////////////////////
// Single reference convolve-x functions (low bit-depth)
////////////////////////////////////////////////////////
typedef void (*convolve_x_func)(const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int w, int h,
                                const InterpFilterParams *filter_params_x,
                                const int subpel_x_qn,
                                ConvolveParams *conv_params);

class AV1ConvolveXTest : public AV1ConvolveTest<Sr1DParam> {
 public:
  void RunTest(convolve_x_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Sr1DParam &p) {
    const int width = p.Block().Width();
    const int height = p.Block().Height();
    const int sub_x = p.SubPixel();
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(p.Filter(), width);
    ConvolveParams conv_params1 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    const uint8_t *input = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    av1_convolve_x_sr(input, width, reference, width, width, height,
                      filter_params_x, sub_x, &conv_params1);

    ConvolveParams conv_params2 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    uint8_t *test = Output8(width, height);
    test_func_(input, width, test, width, width, height, filter_params_x, sub_x,
               &conv_params2);
    AssertEq(reference, test, width, height, width);
  }

  convolve_x_func test_func_;
};

TEST_P(AV1ConvolveXTest, C) { RunTest(av1_convolve_x_sr_c); };

#if HAVE_SSE2
TEST_P(AV1ConvolveXTest, SSE2) { RunTest(av1_convolve_x_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(AV1ConvolveXTest, AVX2) { RunTest(av1_convolve_x_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(AV1ConvolveXTest, NEON) { RunTest(av1_convolve_x_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1ConvolveXTest,
                        ::testing::ValuesIn(GetBlockSizes()));

/////////////////////////////////////////////////////////
// Single reference convolve-x functions (high bit-depth)
/////////////////////////////////////////////////////////
typedef void (*highbd_convolve_x_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x, const int subpel_x_qn,
    ConvolveParams *conv_params, int bd);

class AV1HighbdConvolveXTest : public AV1ConvolveTest<Highbd<Sr1DParam>> {
 public:
  void RunTest(highbd_convolve_x_func func) {
    test_func_ = func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Highbd<Sr1DParam> &highbd) override {
    const int width = highbd.Param().Block().Width();
    const int height = highbd.Param().Block().Height();
    const int sub_x = highbd.Param().SubPixel();
    const int bit_depth = highbd.BitDepth();
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(highbd.Param().Filter(),
                                                     width);
    ConvolveParams conv_params1 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    const uint16_t *input = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    av1_highbd_convolve_x_sr(input, width, reference, width, width, height,
                             filter_params_x, sub_x, &conv_params1, bit_depth);

    ConvolveParams conv_params2 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    uint16_t *test = Output16(width, height);
    test_func_(input, width, test, width, width, height, filter_params_x, sub_x,
               &conv_params2, bit_depth);
    AssertEq(reference, test, width, height, width);
  }

  highbd_convolve_x_func test_func_;
};

TEST_P(AV1HighbdConvolveXTest, C) { RunTest(av1_highbd_convolve_x_sr_c); };

#if HAVE_SSSE3
TEST_P(AV1HighbdConvolveXTest, SSSE3) {
  RunTest(av1_highbd_convolve_x_sr_ssse3);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdConvolveXTest, AVX2) { RunTest(av1_highbd_convolve_x_sr_avx2); }
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdConvolveXTest,
                        testing::ValuesIn(GetBlockSizes()));

////////////////////////////////////////////////////////
// Single reference convolve-y functions (low bit-depth)
////////////////////////////////////////////////////////
typedef void (*convolve_y_func)(const uint8_t *src, int src_stride,
                                uint8_t *dst, int dst_stride, int w, int h,
                                const InterpFilterParams *filter_params_y,
                                const int subpel_y_qn);

class AV1ConvolveYTest : public AV1ConvolveTest<Sr1DParam> {
 public:
  void RunTest(convolve_y_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Sr1DParam &param) override {
    const int width = param.Block().Width();
    const int height = param.Block().Height();
    int sub_y = param.SubPixel();
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(param.Filter(), height);
    const uint8_t *input = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    av1_convolve_y_sr(input, width, reference, width, width, height,
                      filter_params_y, sub_y);
    uint8_t *test = Output8(width, height);
    test_func_(input, width, test, width, width, height, filter_params_y,
               sub_y);
    AssertEq(reference, test, width, height, width);
  }

  convolve_y_func test_func_;
};

TEST_P(AV1ConvolveYTest, C) { RunTest(av1_convolve_y_sr_c); };

#if HAVE_SSE2
TEST_P(AV1ConvolveYTest, SSE2) { RunTest(av1_convolve_y_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(AV1ConvolveYTest, AVX2) { RunTest(av1_convolve_y_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(AV1ConvolveYTest, NEON) { RunTest(av1_convolve_y_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1ConvolveYTest,
                        ::testing::ValuesIn(GetBlockSizes()));

/////////////////////////////////////////////////////////
// Single reference convolve-y functions (high bit-depth)
/////////////////////////////////////////////////////////
typedef void (*highbd_convolve_y_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_y, const int subpel_y_qn,
    int bd);

class AV1HighbdConvolveYTest : public AV1ConvolveTest<Highbd<Sr1DParam>> {
 public:
  void RunTest(highbd_convolve_y_func func) {
    test_func_ = func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Highbd<Sr1DParam> &highbd) override {
    const int width = highbd.Param().Block().Width();
    const int height = highbd.Param().Block().Height();
    const int sub_y = highbd.Param().SubPixel();
    const int bit_depth = highbd.BitDepth();
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(highbd.Param().Filter(),
                                                     height);
    const uint16_t *input = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    av1_highbd_convolve_y_sr(input, width, reference, width, width, height,
                             filter_params_y, sub_y, bit_depth);
    uint16_t *test = Output16(width, height);
    test_func_(input, width, test, width, width, height, filter_params_y, sub_y,
               bit_depth);
    AssertEq(reference, test, width, height, width);
  }

  highbd_convolve_y_func test_func_;
};

TEST_P(AV1HighbdConvolveYTest, C) { RunTest(av1_highbd_convolve_y_sr_c); };

#if HAVE_SSSE3
TEST_P(AV1HighbdConvolveYTest, SSSE3) {
  RunTest(av1_highbd_convolve_y_sr_ssse3);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdConvolveYTest, AVX2) { RunTest(av1_highbd_convolve_y_sr_avx2); }
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdConvolveYTest,
                        ::testing::ValuesIn(GetBlockSizes()));

//////////////////////////////////////////////////////////////
// Single reference convolve-copy functions (low bit-depth)
//////////////////////////////////////////////////////////////
typedef void (*convolve_copy_func)(const uint8_t *src, ptrdiff_t src_stride,
                                   uint8_t *dst, ptrdiff_t dst_stride, int w,
                                   int h);

class AV1ConvolveCopyTest : public AV1ConvolveTest<BlockSize> {
 public:
  void RunTest(convolve_copy_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const BlockSize &block) {
    const int width = block.Width();
    const int height = block.Height();
    const uint8_t *input = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    aom_convolve_copy(input, width, reference, width, width, height);
    uint8_t *test = Output8(width, height);
    test_func_(input, width, test, width, width, height);
    AssertEq(reference, test, width, height, width);
  }

  convolve_copy_func test_func_;
};

// Note that even though these are AOM convolve functions, we are using the
// newer AV1 test framework.
TEST_P(AV1ConvolveCopyTest, C) { RunTest(aom_convolve_copy_c); }

#if HAVE_SSE2
TEST_P(AV1ConvolveCopyTest, SSE2) { RunTest(aom_convolve_copy_sse2); }
#endif

#if HAVE_AVX2
TEST_P(AV1ConvolveCopyTest, AVX2) { RunTest(aom_convolve_copy_avx2); }
#endif

#if HAVE_NEON
TEST_P(AV1ConvolveCopyTest, NEON) { RunTest(aom_convolve_copy_neon); }
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1ConvolveCopyTest,
                        ::testing::ValuesIn(GetBlockSizes()));

///////////////////////////////////////////////////////////////
// Single reference convolve-copy functions (high bit-depth)
///////////////////////////////////////////////////////////////
typedef void (*highbd_convolve_copy_func)(const uint16_t *src, int src_stride,
                                          uint16_t *dst, int dst_stride, int w,
                                          int h);

class AV1HighbdConvolveCopyTest : public AV1ConvolveTest<Highbd<BlockSize>> {
 public:
  void RunTest(highbd_convolve_copy_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Highbd<BlockSize> &highbd) {
    const BlockSize &block = highbd.Param();
    const int width = block.Width();
    const int height = block.Height();
    const int bit_depth = highbd.BitDepth();
    const uint16_t *input = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    av1_highbd_convolve_2d_copy_sr(input, width, reference, width, width,
                                   height);
    uint16_t *test = Output16(width, height);
    test_func_(input, width, test, width, width, height);
    AssertEq(reference, test, width, height, width);
  }

  highbd_convolve_copy_func test_func_;
};

TEST_P(AV1HighbdConvolveCopyTest, C) {
  RunTest(av1_highbd_convolve_2d_copy_sr_c);
}

#if HAVE_SSE2
TEST_P(AV1HighbdConvolveCopyTest, SSE2) {
  RunTest(av1_highbd_convolve_2d_copy_sr_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdConvolveCopyTest, AVX2) {
  RunTest(av1_highbd_convolve_2d_copy_sr_avx2);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdConvolveCopyTest,
                        ::testing::ValuesIn(GetBlockSizes()));

// Single reference 2D convolution parameters.
class Sr2DParam {
 public:
  Sr2DParam(const BlockSize &b, InterpFilter h_f, InterpFilter v_f, int sub_x,
            int sub_y)
      : block_(b), horizontal_filter_(h_f), vertical_filter_(v_f),
        sub_x_(sub_x), sub_y_(sub_y) {}

  const BlockSize &Block() const { return block_; }
  int SubX() const { return sub_x_; }
  int SubY() const { return sub_y_; }
  InterpFilter HorizontalFilter() const { return horizontal_filter_; }
  InterpFilter VerticalFilter() const { return vertical_filter_; }

 private:
  BlockSize block_;
  InterpFilter horizontal_filter_;
  InterpFilter vertical_filter_;
  int sub_x_;
  int sub_y_;
};

template <>
class ParamIterator<Sr2DParam> {
 public:
  ParamIterator<Sr2DParam>(const BlockSize &b)
      : block_(b), horizontal_filter_(EIGHTTAP_REGULAR),
        vertical_filter_(EIGHTTAP_REGULAR), sub_x_(0), sub_y_(0) {}

  bool HasNext() const { return horizontal_filter_ < INTERP_FILTERS_ALL; }
  Sr2DParam Next() {
    const InterpFilter h_f = static_cast<InterpFilter>(horizontal_filter_);
    const InterpFilter v_f = static_cast<InterpFilter>(vertical_filter_);
    const int s_x = sub_x_;
    const int s_y = sub_y_;

    if (sub_y_ < 15) {
      ++sub_y_;
    } else if (sub_x_ < 15) {
      sub_y_ = 0;
      ++sub_x_;
    } else if (vertical_filter_ < INTERP_FILTERS_ALL - 1) {
      sub_y_ = 0;
      sub_x_ = 0;
      ++vertical_filter_;
    } else {
      sub_y_ = 0;
      sub_x_ = 0;
      vertical_filter_ = 0;
      ++horizontal_filter_;
    }
    return Sr2DParam(block_, h_f, v_f, s_x, s_y);
  }

 private:
  const BlockSize block_;
  int horizontal_filter_;
  int vertical_filter_;
  int sub_x_;
  int sub_y_;
};

TEST_F(AV1ConvolveIteratorTest, Sr2DParam) {
  BlockSize block(16, 32);
  ParamIterator<Sr2DParam> iterator(block);
  for (int h_f = EIGHTTAP_REGULAR; h_f < INTERP_FILTERS_ALL; ++h_f) {
    for (int v_f = EIGHTTAP_REGULAR; v_f < INTERP_FILTERS_ALL; ++v_f) {
      for (int sub_x = 0; sub_x < 16; ++sub_x) {
        for (int sub_y = 0; sub_y < 16; ++sub_y) {
          ASSERT_TRUE(iterator.HasNext());
          Sr2DParam p = iterator.Next();
          ASSERT_EQ(block, p.Block());
          ASSERT_EQ(static_cast<InterpFilter>(h_f), p.HorizontalFilter());
          ASSERT_EQ(static_cast<InterpFilter>(v_f), p.VerticalFilter());
          ASSERT_EQ(sub_x, p.SubX());
          ASSERT_EQ(sub_y, p.SubY());
        }
      }
    }
  }
  ASSERT_FALSE(iterator.HasNext());
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

class AV1Convolve2DTest : public AV1ConvolveTest<Sr2DParam> {
 public:
  void RunTest(convolve_2d_func func) {
    test_func_ = func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Sr2DParam &param) override {
    const int width = param.Block().Width();
    const int height = param.Block().Height();
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(param.HorizontalFilter(),
                                                     width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(param.VerticalFilter(),
                                                     height);
    const int sub_x = param.SubX();
    const int sub_y = param.SubY();
    const uint8_t *input = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    ConvolveParams conv_params1 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    av1_convolve_2d_sr(input, width, reference, width, width, height,
                       filter_params_x, filter_params_y, sub_x, sub_y,
                       &conv_params1);
    uint8_t *test = Output8(width, height);
    ConvolveParams conv_params2 = get_conv_params_no_round(0, 0, NULL, 0, 0, 8);
    test_func_(input, width, test, width, width, height, filter_params_x,
               filter_params_y, sub_x, sub_y, &conv_params2);
    AssertEq(reference, test, width, height, width);
  }

  convolve_2d_func test_func_;
};

TEST_P(AV1Convolve2DTest, C) { RunTest(av1_convolve_2d_sr_c); }

#if HAVE_SSE2
TEST_P(AV1Convolve2DTest, SSE2) { RunTest(av1_convolve_2d_sr_sse2); }
#endif

#if HAVE_AVX2
TEST_P(AV1Convolve2DTest, AVX2) { RunTest(av1_convolve_2d_sr_avx2); }
#endif

#if HAVE_NEON
TEST_P(AV1Convolve2DTest, NEON) { RunTest(av1_convolve_2d_sr_neon); }
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1Convolve2DTest,
                        ::testing::ValuesIn(GetBlockSizes()));

//////////////////////////////////////////////////////////
// Single reference convolve-2d functions (high bit-depth)
//////////////////////////////////////////////////////////

typedef void (*highbd_convolve_2d_func)(
    const uint16_t *src, int src_stride, uint16_t *dst, int dst_stride, int w,
    int h, const InterpFilterParams *filter_params_x,
    const InterpFilterParams *filter_params_y, const int subpel_x_qn,
    const int subpel_y_qn, ConvolveParams *conv_params, int bd);

class AV1HighbdConvolve2DTest : public AV1ConvolveTest<Highbd<Sr2DParam>> {
 public:
  void RunTest(highbd_convolve_2d_func func) {
    test_func_ = func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Highbd<Sr2DParam> &highbd) override {
    const int width = highbd.Param().Block().Width();
    const int height = highbd.Param().Block().Height();
    const int bit_depth = highbd.BitDepth();
    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(
            highbd.Param().HorizontalFilter(), width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(
            highbd.Param().VerticalFilter(), height);
    const int sub_x = highbd.Param().SubX();
    const int sub_y = highbd.Param().SubY();
    const uint16_t *input = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    ConvolveParams conv_params1 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    av1_highbd_convolve_2d_sr(input, width, reference, width, width, height,
                              filter_params_x, filter_params_y, sub_x, sub_y,
                              &conv_params1, bit_depth);
    uint16_t *test = Output16(width, height);
    ConvolveParams conv_params2 =
        get_conv_params_no_round(0, 0, NULL, 0, 0, bit_depth);
    test_func_(input, width, test, width, width, height, filter_params_x,
               filter_params_y, sub_x, sub_y, &conv_params2, bit_depth);
    AssertEq(reference, test, width, height, width);
  }

  highbd_convolve_2d_func test_func_;
};

TEST_P(AV1HighbdConvolve2DTest, C) { RunTest(av1_highbd_convolve_2d_sr_c); }

#if HAVE_SSSE3
TEST_P(AV1HighbdConvolve2DTest, SSSE3) {
  RunTest(av1_highbd_convolve_2d_sr_ssse3);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdConvolve2DTest, AVX2) {
  RunTest(av1_highbd_convolve_2d_sr_avx2);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdConvolve2DTest,
                        ::testing::ValuesIn(GetBlockSizes()));

//////////////////////////
// Compound Convolve Tests
//////////////////////////

// The compound 2d-copy functions do not work for chroma block sizes. Only
// generate luma sizes.
std::set<BlockSize> GetLumaBlockSizes() {
  std::set<BlockSize> sizes;
  for (int b = BLOCK_4X4; b < BLOCK_SIZES_ALL; ++b) {
    const int w = block_size_wide[b];
    const int h = block_size_high[b];
    sizes.insert(BlockSize(w, h));
  }
  return sizes;
}

// Compound cases also need to test different frame offsets and weightings.
template <typename T>
class Compound {
 public:
  Compound(const T &param, bool use_dist_wtd_comp_avg, int fwd_offset,
           int bck_offset)
      : param_(param), use_dist_wtd_comp_avg_(use_dist_wtd_comp_avg),
        fwd_offset_(fwd_offset), bck_offset_(bck_offset) {}

  const T &Param() const { return param_; }
  bool UseDistWtdCompAvg() const { return use_dist_wtd_comp_avg_; }
  int FwdOffset() const { return fwd_offset_; }
  int BckOffset() const { return bck_offset_; }

 private:
  T param_;
  bool use_dist_wtd_comp_avg_;
  int fwd_offset_;
  int bck_offset_;
};

template <typename T>
class ParamIterator<Compound<T>> {
 public:
  explicit ParamIterator(const BlockSize &block) : iter_(block), values_() {}

  bool HasNext() const { return !values_.empty() || iter_.HasNext(); }
  Compound<T> Next() {
    if (!values_.empty()) {
      auto r = values_.back();
      values_.pop_back();
      return r;
    }
    T param = iter_.Next();
    for (int k = 1; k >= 0; --k) {
      for (int l = 3; l >= 0; --l) {
        int fwd_offset = quant_dist_lookup_table[k][l][0];
        int bck_offset = quant_dist_lookup_table[k][l][1];
        values_.push_back(Compound<T>(param, true, fwd_offset, bck_offset));
      }
    }
    values_.push_back(Compound<T>(param, false, 0, 0));
    return Next();
  }

 private:
  ParamIterator<T> iter_;
  std::vector<Compound<T>> values_;
};

TEST_F(AV1ConvolveIteratorTest, Compound) {
  BlockSize block(128, 128);
  ParamIterator<Compound<BlockSize>> iterator(block);
  ASSERT_TRUE(iterator.HasNext());
  Compound<BlockSize> c = iterator.Next();
  ASSERT_EQ(block, c.Param());
  ASSERT_FALSE(c.UseDistWtdCompAvg());
  for (int k = 0; k < 2; ++k) {
    for (int l = 0; l < 4; ++l) {
      ASSERT_TRUE(iterator.HasNext());
      c = iterator.Next();
      ASSERT_EQ(block, c.Param());
      ASSERT_TRUE(c.UseDistWtdCompAvg());
      ASSERT_EQ(quant_dist_lookup_table[k][l][0], c.FwdOffset());
      ASSERT_EQ(quant_dist_lookup_table[k][l][1], c.BckOffset());
    }
  }
  ASSERT_FALSE(iterator.HasNext());
}

////////////////////////////////////////////////
// Compound convolve-x functions (low bit-depth)
////////////////////////////////////////////////

template <typename T>
ConvolveParams GetConvolveParams(int do_average, CONV_BUF_TYPE *conv_buf,
                                 int width, int bit_depth,
                                 const Compound<T> &compound) {
  ConvolveParams conv_params =
      get_conv_params_no_round(do_average, 0 /* plane */, conv_buf, width,
                               1 /* is_compound */, bit_depth);
  conv_params.use_dist_wtd_comp_avg = compound.UseDistWtdCompAvg();
  conv_params.fwd_offset = compound.FwdOffset();
  conv_params.bck_offset = compound.BckOffset();
  return conv_params;
}

class AV1CompoundConvolveXTest : public AV1ConvolveTest<Compound<Sr1DParam>> {
 public:
  void RunTest(convolve_x_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Compound<Sr1DParam> &compound) override {
    const int width = compound.Param().Block().Width();
    const int height = compound.Param().Block().Height();

    const uint8_t *input1 = RandomInput8(width, height);
    const uint8_t *input2 = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(ReferenceFunc(), input1, input2, reference, reference_conv_buf,
             compound);

    uint8_t *test = Output8(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func_, input1, input2, test, test_conv_buf, compound);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

  virtual const InterpFilterParams *FilterParams(InterpFilter f,
                                                 const BlockSize &block) const {
    return av1_get_interp_filter_params_with_block_size(f, block.Width());
  }

  virtual convolve_x_func ReferenceFunc() const {
    return av1_dist_wtd_convolve_x;
  }

 private:
  void Convolve(convolve_x_func test_func, const uint8_t *src1,
                const uint8_t *src2, uint8_t *dst, uint16_t *conv_buf,
                const Compound<Sr1DParam> &compound) {
    const int width = compound.Param().Block().Width();
    const int height = compound.Param().Block().Height();
    const int sub_pix = compound.Param().SubPixel();
    const InterpFilterParams *filter_params =
        FilterParams(compound.Param().Filter(), compound.Param().Block());

    ConvolveParams conv_params =
        GetConvolveParams(0, conv_buf, width, 8, compound);
    test_func(src1, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params);

    conv_params = GetConvolveParams(1, conv_buf, width, 8, compound);
    test_func(src2, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params);
  }

  convolve_x_func test_func_;
};

TEST_P(AV1CompoundConvolveXTest, C) { RunTest(av1_dist_wtd_convolve_x_c); }

#if HAVE_SSE2
TEST_P(AV1CompoundConvolveXTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_x_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(AV1CompoundConvolveXTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_x_avx2);
}
#endif

#if HAVE_NEON
TEST_P(AV1CompoundConvolveXTest, NEON) {
  RunTest(av1_dist_wtd_convolve_x_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1CompoundConvolveXTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

/////////////////////////////////////////////////
// Compound convolve-x functions (high bit-depth)
/////////////////////////////////////////////////
class AV1HighbdCompoundConvolveXTest
    : public AV1ConvolveTest<Highbd<Compound<Sr1DParam>>> {
 public:
  void RunTest(highbd_convolve_x_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  virtual const InterpFilterParams *FilterParams(InterpFilter f,
                                                 const BlockSize &block) const {
    return av1_get_interp_filter_params_with_block_size(f, block.Width());
  }

  virtual highbd_convolve_x_func ReferenceFunc() const {
    return av1_highbd_dist_wtd_convolve_x;
  }

  void TestConvolve(const Highbd<Compound<Sr1DParam>> &highbd) override {
    const Compound<Sr1DParam> &compound = highbd.Param();
    const int width = compound.Param().Block().Width();
    const int height = compound.Param().Block().Height();
    const int bit_depth = highbd.BitDepth();

    const uint16_t *input1 = RandomInput16(width, height, bit_depth);
    const uint16_t *input2 = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);

    Convolve(ReferenceFunc(), input1, input2, reference, reference_conv_buf,
             highbd);

    uint16_t *test = Output16(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func_, input1, input2, test, test_conv_buf, highbd);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(highbd_convolve_x_func test_func, const uint16_t *src1,
                const uint16_t *src2, uint16_t *dst, uint16_t *conv_buf,
                const Highbd<Compound<Sr1DParam>> &highbd) {
    const Compound<Sr1DParam> &compound = highbd.Param();
    const int width = compound.Param().Block().Width();
    const int height = compound.Param().Block().Height();
    const int sub_pix = compound.Param().SubPixel();
    const int bit_depth = highbd.BitDepth();
    const InterpFilterParams *filter_params =
        FilterParams(compound.Param().Filter(), compound.Param().Block());
    ConvolveParams conv_params =
        GetConvolveParams(0, conv_buf, width, bit_depth, compound);
    test_func(src1, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params, bit_depth);

    conv_params = GetConvolveParams(1, conv_buf, width, bit_depth, compound);
    test_func(src2, width, dst, width, width, height, filter_params, sub_pix,
              &conv_params, bit_depth);
  }

  highbd_convolve_x_func test_func_;
};

TEST_P(AV1HighbdCompoundConvolveXTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_x_c);
}

#if HAVE_SSE4_1
TEST_P(AV1HighbdCompoundConvolveXTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_x_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdCompoundConvolveXTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_x_avx2);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdCompoundConvolveXTest,
                        testing::ValuesIn(GetLumaBlockSizes()));

////////////////////////////////////////////////
// Compound convolve-y functions (low bit-depth)
////////////////////////////////////////////////

// Note that the X and Y convolve functions have the same type signature and
// logic; they only differentiate the filter parameters and reference function.
class AV1CompoundConvolveYTest : public AV1CompoundConvolveXTest {
 protected:
  virtual const InterpFilterParams *FilterParams(
      InterpFilter f, const BlockSize &block) const override {
    return av1_get_interp_filter_params_with_block_size(f, block.Height());
  }

  virtual convolve_x_func ReferenceFunc() const override {
    return av1_dist_wtd_convolve_y;
  }
};

TEST_P(AV1CompoundConvolveYTest, C) { RunTest(av1_dist_wtd_convolve_y_c); }

#if HAVE_SSE2
TEST_P(AV1CompoundConvolveYTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_y_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(AV1CompoundConvolveYTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_y_avx2);
}
#endif

#if HAVE_NEON
TEST_P(AV1CompoundConvolveYTest, NEON) {
  RunTest(av1_dist_wtd_convolve_y_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1CompoundConvolveYTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

/////////////////////////////////////////////////
// Compound convolve-y functions (high bit-depth)
/////////////////////////////////////////////////

// Again, the X and Y convolve functions have the same type signature and logic.
class AV1HighbdCompoundConvolveYTest : public AV1HighbdCompoundConvolveXTest {
  virtual highbd_convolve_x_func ReferenceFunc() const override {
    return av1_highbd_dist_wtd_convolve_y;
  }
  virtual const InterpFilterParams *FilterParams(
      InterpFilter f, const BlockSize &block) const override {
    return av1_get_interp_filter_params_with_block_size(f, block.Height());
  }
};

TEST_P(AV1HighbdCompoundConvolveYTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_y_c);
}

#if HAVE_SSE4_1
TEST_P(AV1HighbdCompoundConvolveYTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_y_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdCompoundConvolveYTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_y_avx2);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdCompoundConvolveYTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

//////////////////////////////////////////////////////
// Compound convolve-2d-copy functions (low bit-depth)
//////////////////////////////////////////////////////
typedef void (*compound_conv_2d_copy_func)(const uint8_t *src, int src_stride,
                                           uint8_t *dst, int dst_stride, int w,
                                           int h, ConvolveParams *conv_params);

class AV1CompoundConvolve2DCopyTest
    : public AV1ConvolveTest<Compound<BlockSize>> {
 public:
  void RunTest(compound_conv_2d_copy_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Compound<BlockSize> &compound) override {
    const BlockSize &block = compound.Param();
    const int width = block.Width();
    const int height = block.Height();

    const uint8_t *input1 = RandomInput8(width, height);
    const uint8_t *input2 = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_dist_wtd_convolve_2d_copy, input1, input2, reference,
             reference_conv_buf, compound);

    uint8_t *test = Output8(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func_, input1, input2, test, test_conv_buf, compound);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(compound_conv_2d_copy_func test_func, const uint8_t *src1,
                const uint8_t *src2, uint8_t *dst, uint16_t *conv_buf,
                const Compound<BlockSize> &compound) {
    const BlockSize &block = compound.Param();
    const int width = block.Width();
    const int height = block.Height();
    ConvolveParams conv_params =
        GetConvolveParams(0, conv_buf, width, 8, compound);
    test_func(src1, width, dst, width, width, height, &conv_params);

    conv_params = GetConvolveParams(1, conv_buf, width, 8, compound);
    test_func(src2, width, dst, width, width, height, &conv_params);
  }

  compound_conv_2d_copy_func test_func_;
};

TEST_P(AV1CompoundConvolve2DCopyTest, C) {
  RunTest(av1_dist_wtd_convolve_2d_copy_c);
}

#if HAVE_SSE2
TEST_P(AV1CompoundConvolve2DCopyTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_2d_copy_sse2);
}
#endif

#if HAVE_AVX2
TEST_P(AV1CompoundConvolve2DCopyTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_2d_copy_avx2);
}
#endif

#if HAVE_NEON
TEST_P(AV1CompoundConvolve2DCopyTest, NEON) {
  RunTest(av1_dist_wtd_convolve_2d_copy_neon);
}
#endif

TEST_F(AV1ConvolveIteratorTest, GetLumaBlockSizes) {
  ASSERT_EQ(22U, GetLumaBlockSizes().size());
}

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1CompoundConvolve2DCopyTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

///////////////////////////////////////////////////////
// Compound convolve-2d-copy functions (high bit-depth)
///////////////////////////////////////////////////////
typedef void (*highbd_compound_conv_2d_copy_func)(const uint16_t *src,
                                                  int src_stride, uint16_t *dst,
                                                  int dst_stride, int w, int h,
                                                  ConvolveParams *conv_params,
                                                  int bd);

class AV1HighbdCompoundConvolve2DCopyTest
    : public AV1ConvolveTest<Highbd<Compound<BlockSize>>> {
 public:
  void RunTest(highbd_compound_conv_2d_copy_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Highbd<Compound<BlockSize>> &highbd) override {
    const Compound<BlockSize> &compound = highbd.Param();
    const BlockSize &block = compound.Param();
    const int width = block.Width();
    const int height = block.Height();
    const int bit_depth = highbd.BitDepth();

    const uint16_t *input1 = RandomInput16(width, height, bit_depth);
    const uint16_t *input2 = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_highbd_dist_wtd_convolve_2d_copy, input1, input2, reference,
             reference_conv_buf, highbd);

    uint16_t *test = Output16(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func_, input1, input2, test, test_conv_buf, highbd);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(highbd_compound_conv_2d_copy_func test_func,
                const uint16_t *src1, const uint16_t *src2, uint16_t *dst,
                uint16_t *conv_buf, const Highbd<Compound<BlockSize>> &highbd) {
    const Compound<BlockSize> &compound = highbd.Param();
    const BlockSize &block = compound.Param();
    const int width = block.Width();
    const int height = block.Height();
    const int bit_depth = highbd.BitDepth();

    ConvolveParams conv_params =
        GetConvolveParams(0, conv_buf, width, bit_depth, compound);
    test_func(src1, width, dst, width, width, height, &conv_params, bit_depth);

    conv_params = GetConvolveParams(1, conv_buf, width, bit_depth, compound);
    test_func(src2, width, dst, width, width, height, &conv_params, bit_depth);
  }

  highbd_compound_conv_2d_copy_func test_func_;
};

TEST_P(AV1HighbdCompoundConvolve2DCopyTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_c);
}

#if HAVE_SSE4_1
TEST_P(AV1HighbdCompoundConvolve2DCopyTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdCompoundConvolve2DCopyTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_copy_avx2);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdCompoundConvolve2DCopyTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

/////////////////////////////////////////////////
// Compound convolve-2d functions (low bit-depth)
/////////////////////////////////////////////////

class AV1CompoundConvolve2DTest : public AV1ConvolveTest<Compound<Sr2DParam>> {
 public:
  void RunTest(convolve_2d_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Compound<Sr2DParam> &compound) override {
    const Sr2DParam &param2d = compound.Param();
    const BlockSize &block = param2d.Block();
    const int width = block.Width();
    const int height = block.Height();

    const uint8_t *input1 = RandomInput8(width, height);
    const uint8_t *input2 = RandomInput8(width, height);
    uint8_t *reference = Output8(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_dist_wtd_convolve_2d, input1, input2, reference,
             reference_conv_buf, compound);

    uint8_t *test = Output8(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func_, input1, input2, test, test_conv_buf, compound);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(convolve_2d_func test_func, const uint8_t *src1,
                const uint8_t *src2, uint8_t *dst, uint16_t *conv_buf,
                const Compound<Sr2DParam> &compound) {
    const Sr2DParam &param2d = compound.Param();
    const BlockSize &block = param2d.Block();
    const int width = block.Width();
    const int height = block.Height();

    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(param2d.HorizontalFilter(),
                                                     width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(param2d.VerticalFilter(),
                                                     height);
    const int sub_x = param2d.SubX();
    const int sub_y = param2d.SubY();
    ConvolveParams conv_params =
        GetConvolveParams(0, conv_buf, width, 8, compound);

    test_func(src1, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params);

    conv_params = GetConvolveParams(1, conv_buf, width, 8, compound);
    test_func(src2, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params);
  }

  convolve_2d_func test_func_;
};

TEST_P(AV1CompoundConvolve2DTest, C) { RunTest(av1_dist_wtd_convolve_2d_c); }

#if HAVE_SSE2
TEST_P(AV1CompoundConvolve2DTest, SSE2) {
  RunTest(av1_dist_wtd_convolve_2d_sse2);
}
#endif

#if HAVE_SSSE3
TEST_P(AV1CompoundConvolve2DTest, SSSE3) {
  RunTest(av1_dist_wtd_convolve_2d_ssse3);
}
#endif

#if HAVE_AVX2
TEST_P(AV1CompoundConvolve2DTest, AVX2) {
  RunTest(av1_dist_wtd_convolve_2d_avx2);
}
#endif

#if HAVE_NEON
TEST_P(AV1CompoundConvolve2DTest, NEON) {
  RunTest(av1_dist_wtd_convolve_2d_neon);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1CompoundConvolve2DTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

//////////////////////////////////////////////////
// Compound convolve-2d functions (high bit-depth)
//////////////////////////////////////////////////

class AV1HighbdCompoundConvolve2DTest
    : public AV1ConvolveTest<Highbd<Compound<Sr2DParam>>> {
 public:
  void RunTest(highbd_convolve_2d_func test_func) {
    test_func_ = test_func;
    AV1ConvolveTest::RunTest();
  }

 protected:
  void TestConvolve(const Highbd<Compound<Sr2DParam>> &highbd) override {
    const Compound<Sr2DParam> &compound = highbd.Param();
    const Sr2DParam &param2d = compound.Param();
    const BlockSize &block = param2d.Block();
    const int width = block.Width();
    const int height = block.Height();
    const int bit_depth = highbd.BitDepth();
    const uint16_t *input1 = RandomInput16(width, height, bit_depth);
    const uint16_t *input2 = RandomInput16(width, height, bit_depth);
    uint16_t *reference = Output16(width, height);
    uint16_t *reference_conv_buf = Output16(width, height);
    Convolve(av1_highbd_dist_wtd_convolve_2d, input1, input2, reference,
             reference_conv_buf, highbd);

    uint16_t *test = Output16(width, height);
    uint16_t *test_conv_buf = Output16(width, height);
    Convolve(test_func_, input1, input2, test, test_conv_buf, highbd);

    AssertEq(reference_conv_buf, test_conv_buf, width, height, width);
    AssertEq(reference, test, width, height, width);
  }

 private:
  void Convolve(highbd_convolve_2d_func test_func, const uint16_t *src1,
                const uint16_t *src2, uint16_t *dst, uint16_t *conv_buf,
                const Highbd<Compound<Sr2DParam>> &highbd) {
    const Compound<Sr2DParam> &compound = highbd.Param();
    const Sr2DParam &param2d = compound.Param();
    const BlockSize &block = param2d.Block();
    const int width = block.Width();
    const int height = block.Height();

    const InterpFilterParams *filter_params_x =
        av1_get_interp_filter_params_with_block_size(param2d.HorizontalFilter(),
                                                     width);
    const InterpFilterParams *filter_params_y =
        av1_get_interp_filter_params_with_block_size(param2d.VerticalFilter(),
                                                     height);
    const int bit_depth = highbd.BitDepth();
    const int sub_x = param2d.SubX();
    const int sub_y = param2d.SubY();
    ConvolveParams conv_params =
        GetConvolveParams(0, conv_buf, width, bit_depth, compound);
    test_func(src1, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params, bit_depth);

    conv_params = GetConvolveParams(1, conv_buf, width, bit_depth, compound);
    test_func(src2, width, dst, width, width, height, filter_params_x,
              filter_params_y, sub_x, sub_y, &conv_params, bit_depth);
  }

  highbd_convolve_2d_func test_func_;
};

TEST_P(AV1HighbdCompoundConvolve2DTest, C) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_c);
}

#if HAVE_SSE4_1
TEST_P(AV1HighbdCompoundConvolve2DTest, SSE4_1) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_sse4_1);
}
#endif

#if HAVE_AVX2
TEST_P(AV1HighbdCompoundConvolve2DTest, AVX2) {
  RunTest(av1_highbd_dist_wtd_convolve_2d_avx2);
}
#endif

INSTANTIATE_TEST_CASE_P(AV1Convolve, AV1HighbdCompoundConvolve2DTest,
                        ::testing::ValuesIn(GetLumaBlockSizes()));

}  // namespace
