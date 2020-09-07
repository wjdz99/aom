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

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "av1/common/mfqe.h"

#define MFQE_TEST_STRIDE 80
#define MFQE_TEST_HEIGHT 32
#define MFQE_TEST_WIDTH 32

namespace {

class MFQETest : public ::testing::Test {
 protected:
  Y_BUFFER_CONFIG *tmp_;
  RefCntBuffer **ref_frames_;

  void SetUp() {
    tmp_ = reinterpret_cast<Y_BUFFER_CONFIG *>(aom_memalign(32, sizeof(*tmp_)));
    tmp_->stride = MFQE_TEST_STRIDE;
    tmp_->height = MFQE_TEST_HEIGHT;
    tmp_->width = MFQE_TEST_WIDTH;

    int buf_bytes = (tmp_->height + 2 * MFQE_PADDING_SIZE) * tmp_->stride;
    tmp_->buffer_orig = reinterpret_cast<uint8_t *>(
        aom_memalign(32, sizeof(uint8_t) * buf_bytes));
    tmp_->buffer = tmp_->buffer_orig + MFQE_PADDING_SIZE * tmp_->stride;
    memset(tmp_->buffer_orig, 0, sizeof(uint8_t) * buf_bytes);

    ref_frames_ = reinterpret_cast<RefCntBuffer **>(
        aom_memalign(32, sizeof(*ref_frames_) * (MFQE_NUM_REFS)));

    for (int i = 0; i < MFQE_NUM_REFS; i++) {
      ref_frames_[i] = reinterpret_cast<RefCntBuffer *>(
          aom_memalign(32, sizeof(**ref_frames_)));

      ref_frames_[i]->buf.y_width = MFQE_TEST_WIDTH;
      ref_frames_[i]->buf.y_height = MFQE_TEST_HEIGHT;
      ref_frames_[i]->buf.y_stride = MFQE_TEST_STRIDE;

      buf_bytes = (MFQE_TEST_HEIGHT + 2 * MFQE_PADDING_SIZE) * MFQE_TEST_STRIDE;
      uint8_t *buffer = reinterpret_cast<uint8_t *>(
          aom_memalign(32, sizeof(uint8_t) * buf_bytes));
      ref_frames_[i]->buf.y_buffer =
          buffer + MFQE_PADDING_SIZE * MFQE_TEST_STRIDE;
      memset(buffer, 0, sizeof(uint8_t) * buf_bytes);
    }

    // Assert that pointers to data are properly set up.
    ASSERT_NE(tmp_, nullptr);
    ASSERT_NE(ref_frames_, nullptr);

    // Assert that stride, height, width are properly set up.
    ASSERT_EQ(tmp_->stride, MFQE_TEST_STRIDE);
    ASSERT_EQ(tmp_->height, MFQE_TEST_HEIGHT);
    ASSERT_EQ(tmp_->width, MFQE_TEST_WIDTH);
  }

  void TearDown() {
    aom_free(tmp_->buffer_orig);
    aom_free(tmp_);

    for (int i = 0; i < MFQE_NUM_REFS; i++) {
      uint8_t *buffer =
          ref_frames_[i]->buf.y_buffer - MFQE_PADDING_SIZE * MFQE_TEST_STRIDE;
      aom_free(buffer);
      aom_free(ref_frames_[i]);
    }
    aom_free(ref_frames_);
  }
};

}  // namespace

TEST_F(MFQETest, TestLowbitdepth1) {
  // Test when the current frame and reference frames are all zeros, then the
  // current frame should still be all zeros after apply MFQE.
  int buf_size = tmp_->stride * tmp_->height;
  for (int i = 0; i < buf_size; i++) ASSERT_EQ(tmp_->buffer[i], 0);

  for (int i = 0; i < MFQE_NUM_REFS; i++) {
    uint8_t *ref_y_buffer = ref_frames_[i]->buf.y_buffer;
    for (int j = 0; j < buf_size; j++) ASSERT_EQ(ref_y_buffer[i], 0);
  }

  int high_bd = 0;
  int bitdepth = 8;
  av1_apply_loop_mfqe(tmp_, ref_frames_, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_EQ(tmp_->buffer[i], 0);

  // Test when every 11th bytes is set to 1 in the current frame, then the
  // values must be either 0 or 1 in the current frame after applying the
  // MFQE, given that reference frames are all zeros.
  for (int i = 0; i < buf_size; i++) {
    if (i % 11 == 0) tmp_->buffer[i] = 1;
  }

  av1_apply_loop_mfqe(tmp_, ref_frames_, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_LE(tmp_->buffer[i], 1);
}

TEST_F(MFQETest, TestLowbitdepth2) {
  // Test when the bytes in the current frame are set to its index modulo 7,
  // and the reference frames are all zeros, then the bytes in the current
  // frame must still have a value between 0 and 7 after applying MFQE.
  const int high_bd = 0;
  const int bitdepth = 8;

  for (int i = 0; i < buf_size; i++) {
    tmp_->buffer[i] = (i % 7);
  }
  av1_apply_loop_mfqe(tmp_, ref_frames_, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_LE(tmp_->buffer[i], 7);

  // Test when the bytes in the current frame are set to its index modulo 9
  // plus 10 and the reference frames are all zeros, then the bytes in the
  // current frames must be greater or equal to 10 after applying MFQE.
  for (int i = 0; i < buf_size; i++) {
    tmp_->buffer[i] = (i % 9) + 10;
  }
  av1_apply_loop_mfqe(tmp_, ref_frames_, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_GE(tmp_->buffer[i], 10);

  // Test when the bytes in the current frame are set to 50 and the reference
  // frames are all zeros, then the bytes in the current frames must still
  // be 50 after applying MFQE.
  for (int i = 0; i < buf_size; i++) {
    tmp_->buffer[i] = 50;
  }
  av1_apply_loop_mfqe(tmp_, ref_frames_, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_EQ(tmp_->buffer[i], 50);

  // Test when the bytes in the current frame are set to 100 and the reference
  // frames are all zeros, then the bytes in the current frames must still
  // be 100 after applying MFQE.
  for (int i = 0; i < buf_size; i++) {
    tmp_->buffer[i] = 100;
  }
  av1_apply_loop_mfqe(tmp_, ref_frames_, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_EQ(tmp_->buffer[i], 100);
}
