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

#include <vector>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "./av1_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

#define DO_PERF_TEST 0
#if DO_PERF_TEST
#define TEST_ITERS 1000
#else
#define TEST_ITERS 10
#endif

namespace {

const int vpad = 32;
const int hpad = 32;
const int x_step_qn = 16;
const int y_step_qn = 20;

using std::tr1::tuple;
using libaom_test::ACMRandom;

// A 16-bit filter with a configurable number of taps.
struct TestFilter {
  void set(int n, bool backwards);

  InterpFilterParams params_;

 private:
  std::vector<int16_t> coeffs_;
};

void TestFilter::set(int n, bool backwards) {
  // The filter has n * SUBPEL_SHIFTS proper elements and an extra 8 bogus
  // elements at the end so that convolutions can read off the end safely.
  coeffs_.resize(n * SUBPEL_SHIFTS + 8);

  // The coefficients are pretty much arbitrary, but convolutions shouldn't
  // over or underflow. For the first filter (subpels = 0), we use an
  // increasing or decreasing ramp (depending on the backwards parameter). We
  // don't want any zero coefficients, so we make it have an x-intercept at -1
  // or n.
  //
  // When increasing, the function has the form:
  //
  //   f(x) = A * (x + 1)
  //
  // Summing and rearranging for A gives A = 2 * I / (n * (n + 1)). If the
  // filter is reversed, we have the same A but with formula
  //
  //   g(x) = A * (n - x)
  const float I = 1 << FILTER_BITS;
  const float A = 2 * I / (n * (n + 1.f));
  for (int i = 0; i < n; ++i)
    coeffs_[i] = (int16_t)(A * (backwards ? (n - i) : (i + 1)));

  // For the other filters, make them slightly different by swapping two
  // columns. Filter k will have the columns (k % n) and (7 * k) % n swapped.
  size_t filter_size = sizeof(coeffs_[0] * n);
  int16_t *filter0 = &coeffs_[0];
  for (int k = 1; k < SUBPEL_SHIFTS; ++k) {
    int16_t *filterk = &coeffs_[k * n];
    memcpy(filterk, filter0, filter_size);

    int idx0 = k % n;
    int idx1 = (7 * k) % n;

    int16_t tmp = filterk[idx0];
    filterk[idx0] = filterk[idx1];
    filterk[idx1] = tmp;
  }

  // Finally, write some rubbish at the end to make sure we don't use it.
  for (int i = 0; i < 8; ++i) coeffs_[n * SUBPEL_SHIFTS + i] = 123 + i;

  // Fill in params
  params_.filter_ptr = &coeffs_[0];
  params_.taps = n;
  // These are ignored by the functions being tested. Set them to whatever.
  params_.subpel_shifts = SUBPEL_SHIFTS;
  params_.interp_filter = EIGHTTAP_REGULAR;
}

template <typename src_pixel_t>
class TestImage {
 public:
  TestImage(int w, int h, int bd) : w_(w), h_(h), bd_(bd) {
    assert(bd < 16);
    assert(bd <= 8 * (int)sizeof(src_pixel_t));

    // Pad width by 2*hpad and then round up to the next multiple of 16
    // to get src_stride_. Add another 16 for dst_stride_ (to make sure
    // something goes wrong if we use the wrong one)
    src_stride_ = (w_ + 2 * hpad + 15) & ~15;
    dst_stride_ = src_stride_ + 16;

    // Allocate image data
    src_data_.resize(2 * GetBlockSize(false));
    dst_data_.resize(2 * GetBlockSize(true));
  }

  void Initialize(ACMRandom *rnd);
  void Check(const char *test_name) const;

  int GetStride(bool dst) const { return dst ? dst_stride_ : src_stride_; }

  int GetBlockSize(bool dst) const { return (h_ + 2 * vpad) * GetStride(dst); }

  const src_pixel_t *GetSrcData(bool ref, bool borders) const {
    const src_pixel_t *block = &src_data_[ref ? 0 : GetBlockSize(false)];
    return borders ? block : block + hpad + src_stride_ * vpad;
  }

  int32_t *GetDstData(bool ref, bool borders) {
    int32_t *block = &dst_data_[ref ? 0 : GetBlockSize(true)];
    return borders ? block : block + hpad + dst_stride_ * vpad;
  }

 private:
  int w_, h_, bd_;
  int src_stride_, dst_stride_;

  std::vector<src_pixel_t> src_data_;
  std::vector<int32_t> dst_data_;
};

template <typename pixel_t>
void FillEdge(ACMRandom *rnd, int nelts, int bd, bool trash, pixel_t *data) {
  if (!trash) {
    memset(data, 0, sizeof(*data) * nelts);
    return;
  }
  pixel_t mask = (1 << bd) - 1;
  for (int i = 0; i < nelts; ++i) data[i] = rnd->Rand16() & mask;
}

template <typename pixel_t>
void PrepBuffers(ACMRandom *rnd, int w, int h, int stride, int bd,
                 bool trash_edges, pixel_t *data) {
  assert(rnd);
  pixel_t mask = (1 << bd) - 1;

  // Fill in the first buffer with random data
  // Top border
  FillEdge(rnd, stride * vpad, bd, trash_edges, data);
  for (int r = 0; r < h; ++r) {
    pixel_t *row_data = data + (vpad + r) * stride;
    // Left border, contents, right border
    FillEdge(rnd, hpad, bd, trash_edges, row_data);
    for (int c = 0; c < w; ++c) row_data[hpad + c] = rnd->Rand16() & mask;
    FillEdge(rnd, hpad, bd, trash_edges, row_data + hpad + w);
  }
  // Bottom border
  FillEdge(rnd, stride * vpad, bd, trash_edges, data + stride * (vpad + h));

  const int bpp = sizeof(*data);
  const int block_elts = stride * (h + 2 * vpad);
  const int block_size = bpp * block_elts;

  // Now copy that to the second buffer
  memcpy(data + block_elts, data, block_size);
}

template <typename src_pixel_t>
void TestImage<src_pixel_t>::Initialize(ACMRandom *rnd) {
  PrepBuffers(rnd, w_, h_, src_stride_, bd_, false, &src_data_[0]);
  PrepBuffers(rnd, w_, h_, dst_stride_, bd_, true, &dst_data_[0]);
}

template <typename src_pixel_t>
void TestImage<src_pixel_t>::Check(const char *test_name) const {
  // If memcmp returns 0, there's nothing to do.
  const int nelts = GetBlockSize(true);
  const int32_t *ref_dst = &dst_data_[0];
  const int32_t *tst_dst = &dst_data_[nelts];

  if (0 == memcmp(ref_dst, tst_dst, sizeof(*ref_dst) * nelts)) return;

  // Otherwise, iterate through the buffer looking for differences (including
  // the edges)
  const int stride = dst_stride_;
  for (int r = 0; r < h_ + 2 * vpad; ++r) {
    for (int c = 0; c < w_ + 2 * hpad; ++c) {
      int32_t ref_value = ref_dst[r * stride + c];
      int32_t tst_value = tst_dst[r * stride + c];

      EXPECT_EQ(tst_value, ref_value)
          << "Error running " << test_name << " at row: " << (r - vpad)
          << ", col: " << (c - hpad);
    }
  }
}

typedef tuple<int, int> BlockDimension;

// This structure is just used as a templated namespace and is never
// instantiatied (since this is C++03, there are no templated
// typedefs)
template <typename src_pixel_t>
struct TestTypes {
  typedef void (*convolve_fun_t)(const src_pixel_t *src, int src_stride,
                                 int32_t *dst, int dst_stride, int w, int h,
                                 InterpFilterParams *filter_params_x,
                                 InterpFilterParams *filter_params_y,
                                 const int subpel_x_qn, const int x_step_qn,
                                 const int subpel_y_qn, const int y_step_qn,
                                 ConvolveParams *conv_params);
  typedef tuple<convolve_fun_t, convolve_fun_t, BlockDimension, uint16_t,
                uint16_t, bool>
      Params;

  typedef ::testing::TestWithParam<Params> TestBase;
};

// Test parameter list:
//  <ref_fun, tst_fun, dims, ntaps_x, ntaps_y, avg>
template <typename src_pixel_t>
class ConvolveScaleTest : public TestTypes<src_pixel_t>::TestBase {
 public:
  ConvolveScaleTest() : image_(NULL) {}
  virtual ~ConvolveScaleTest() { delete image_; }

  // For some reason, this needs the explicit template type
  // argument. This just passes the call straight through to gtest.h
  const typename TestTypes<src_pixel_t>::Params &GetParam() const {
    return TestTypes<src_pixel_t>::TestBase::GetParam();
  }

  virtual void SetUp() {
    ref_fun_ = GET_PARAM(0);
    tst_fun_ = GET_PARAM(1);

    BlockDimension block = GET_PARAM(2);
    width_ = std::tr1::get<0>(block);
    height_ = std::tr1::get<1>(block);
    ntaps_x_ = GET_PARAM(3);
    ntaps_y_ = GET_PARAM(4);
    bd_ = 8;
    avg_ = GET_PARAM(5);

    filter_x_.set(ntaps_x_, false);
    filter_y_.set(ntaps_y_, true);
    convolve_params_ = get_conv_params(0, avg_ != false, 0);

    delete image_;
    image_ = new TestImage<src_pixel_t>(width_, height_, bd_);
  }

  virtual void TearDown() { libaom_test::ClearSystemState(); }

 protected:
  void Run() {
    ACMRandom rnd(ACMRandom::DeterministicSeed());

    const uint8_t *ref_src = image_->GetSrcData(true, false);
    const uint8_t *tst_src = image_->GetSrcData(false, false);
    CONV_BUF_TYPE *ref_dst = image_->GetDstData(true, false);
    CONV_BUF_TYPE *tst_dst = image_->GetDstData(false, false);
    int src_stride = image_->GetStride(false);
    int dst_stride = image_->GetStride(true);

    for (int i = 0; i < TEST_ITERS; ++i) {
      Prep(&rnd);

      ref_fun_(tst_src, src_stride, tst_dst, dst_stride, width_, height_,
               &filter_x_.params_, &filter_y_.params_, subpel_x_, x_step_qn,
               subpel_y_, y_step_qn, &convolve_params_);
      tst_fun_(ref_src, src_stride, ref_dst, dst_stride, width_, height_,
               &filter_x_.params_, &filter_y_.params_, subpel_x_, x_step_qn,
               subpel_y_, y_step_qn, &convolve_params_);

      image_->Check("av1_convolve_2d_scale");
    }
  }

  void Prep(ACMRandom *rnd) {
    assert(rnd);

    // Choose subpel_x_ and subpel_y_. They should be less than
    // SCALE_SUBPEL_SHIFTS; we also add extra weight to "interesting" values: 0
    // and SCALE_SUBPEL_SHIFTS - 1
    uint8_t subpel_mode = rnd->Rand8();
    if ((subpel_mode & 7) == 0)
      subpel_x_ = 0;
    else if ((subpel_mode & 7) == 1)
      subpel_x_ = SCALE_SUBPEL_SHIFTS - 1;
    else
      subpel_x_ = 1 + rnd->PseudoUniform(SCALE_SUBPEL_SHIFTS - 2);

    if (((subpel_mode >> 3) & 7) == 0)
      subpel_y_ = 0;
    else if (((subpel_mode >> 3) & 7) == 1)
      subpel_y_ = SCALE_SUBPEL_SHIFTS - 1;
    else
      subpel_y_ = 1 + rnd->PseudoUniform(SCALE_SUBPEL_SHIFTS - 2);

    image_->Initialize(rnd);
  }

  typename TestTypes<src_pixel_t>::convolve_fun_t ref_fun_, tst_fun_;
  int width_, height_, ntaps_x_, ntaps_y_, bd_;
  bool avg_;
  int subpel_x_, subpel_y_;
  TestFilter filter_x_, filter_y_;
  TestImage<src_pixel_t> *image_;
  ConvolveParams convolve_params_;
};

typedef ConvolveScaleTest<uint8_t> LowBDConvolveScaleTest;

using std::tr1::make_tuple;
const BlockDimension kBlockDim[] = {
  make_tuple(2, 2),    make_tuple(2, 4),    make_tuple(4, 4),
  make_tuple(4, 8),    make_tuple(8, 4),    make_tuple(8, 8),
  make_tuple(8, 16),   make_tuple(16, 8),   make_tuple(16, 16),
  make_tuple(16, 32),  make_tuple(32, 16),  make_tuple(32, 32),
  make_tuple(32, 64),  make_tuple(64, 32),  make_tuple(64, 64),
  make_tuple(64, 128), make_tuple(128, 64), make_tuple(128, 128),
};

const uint16_t kNTaps[] = { 8, 10, 12 };
const bool kAvg[] = { false, true };

TEST_P(LowBDConvolveScaleTest, Check) { Run(); }

INSTANTIATE_TEST_CASE_P(
    SSE4_1, LowBDConvolveScaleTest,
    ::testing::Combine(::testing::Values(av1_convolve_2d_scale_c),
                       ::testing::Values(av1_convolve_2d_scale_sse4_1),
                       ::testing::ValuesIn(kBlockDim),
                       ::testing::ValuesIn(kNTaps), ::testing::ValuesIn(kNTaps),
                       ::testing::ValuesIn(kAvg)));
}  // namespace
