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

#include <assert.h>
#include <stdio.h>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "av1/common/mfqe.h"

#define MFQE_TEST_STRIDE 80
#define MFQE_TEST_HEIGHT 32
#define MFQE_TEST_WIDTH 32

namespace {

class MFQETest : public ::testing::Test {
 protected:
  Y_BUFFER_CONFIG *tmp;
  RefCntBuffer **ref_frames;

  void SetUp() {
    tmp = (Y_BUFFER_CONFIG *)aom_memalign(32, sizeof(*tmp));
    tmp->stride = MFQE_TEST_STRIDE;
    tmp->height = MFQE_TEST_HEIGHT;
    tmp->width = MFQE_TEST_WIDTH;

    int buf_bytes = (tmp->height + 2 * MFQE_PADDING_SIZE) * tmp->stride;
    tmp->buffer_orig =
        (uint8_t *)aom_memalign(32, sizeof(uint8_t) * buf_bytes);
    tmp->buffer = tmp->buffer_orig + MFQE_PADDING_SIZE * tmp->stride;
    memset(tmp->buffer_orig, 0, sizeof(uint8_t) * buf_bytes);

    ref_frames = (RefCntBuffer **)aom_memalign(
        32, sizeof(*ref_frames) * (MFQE_NUM_REFS));

    for (int i = 0; i < MFQE_NUM_REFS; i++) {
      ref_frames[i] = (RefCntBuffer *)aom_memalign(32, sizeof(**ref_frames));

      ref_frames[i]->buf.y_width = MFQE_TEST_WIDTH;
      ref_frames[i]->buf.y_height = MFQE_TEST_HEIGHT;
      ref_frames[i]->buf.y_stride = MFQE_TEST_STRIDE;

      buf_bytes = (MFQE_TEST_HEIGHT + 2 * MFQE_PADDING_SIZE) * MFQE_TEST_STRIDE;
      uint8_t *buffer =
          (uint8_t *)aom_memalign(32, sizeof(uint8_t) * buf_bytes);
      ref_frames[i]->buf.y_buffer =
          buffer + MFQE_PADDING_SIZE * MFQE_TEST_STRIDE;
      memset(buffer, 0, sizeof(uint8_t) * buf_bytes);
    }
  }

  void TearDown() {
    aom_free(tmp->buffer_orig);
    aom_free(tmp);

    for (int i = 0; i < MFQE_NUM_REFS; i++) {
      uint8_t *buffer =
          ref_frames[i]->buf.y_buffer - MFQE_PADDING_SIZE * MFQE_TEST_STRIDE;
      aom_free(buffer);
      aom_free(ref_frames[i]);
    }
    aom_free(ref_frames);
  }
};

}  // namespace

TEST_F(MFQETest, TestLowbitdepth1) {
  ASSERT_NE(tmp, nullptr);
  ASSERT_NE(ref_frames, nullptr);

  ASSERT_EQ(tmp->stride, MFQE_TEST_STRIDE);
  ASSERT_EQ(tmp->height, MFQE_TEST_HEIGHT);
  ASSERT_EQ(tmp->width, MFQE_TEST_WIDTH);

  int buf_size = tmp->stride * tmp->height;
  for (int i = 0; i < buf_size; i++)
    ASSERT_EQ(tmp->buffer[i], 0);

  for (int i = 0; i < MFQE_NUM_REFS; i++) {
    uint8_t *ref_y_buffer = ref_frames[i]->buf.y_buffer;
    for (int j = 0; j < buf_size; j++)
      ASSERT_EQ(ref_y_buffer[i], 0);
  }

  int high_bd = 0;
  int bitdepth = 8;
  av1_apply_loop_mfqe(tmp, ref_frames, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_EQ(tmp->buffer[i], 0);

  for (int i = 0; i < buf_size; i++) {
    if (i % 11 == 0) tmp->buffer[i] = 1;
  }

  av1_apply_loop_mfqe(tmp, ref_frames, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_LE(tmp->buffer[i], 1);
}

TEST_F(MFQETest, TestLowbitdepth2) {
  int buf_size = tmp->stride * tmp->height;
  for (int i = 0; i < buf_size; i++)
    ASSERT_EQ(tmp->buffer[i], 0);

  for (int i = 0; i < MFQE_NUM_REFS; i++) {
    uint8_t *ref_y_buffer = ref_frames[i]->buf.y_buffer;
    for (int j = 0; j < buf_size; j++)
      ASSERT_EQ(ref_y_buffer[i], 0);
  }

  int high_bd = 0;
  int bitdepth = 8;

  for (int i = 0; i < buf_size; i++) {
    tmp->buffer[i] = (i % 7);
  }

  av1_apply_loop_mfqe(tmp, ref_frames, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_LE(tmp->buffer[i], 7);

  for (int i = 0; i < buf_size; i++) {
    tmp->buffer[i] = (i % 9) + 10;
  }

  av1_apply_loop_mfqe(tmp, ref_frames, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_GE(tmp->buffer[i], (i % 9) + 10);

  for (int i = 0; i < buf_size; i++) {
    tmp->buffer[i] = 100;
  }

  av1_apply_loop_mfqe(tmp, ref_frames, MFQE_BLOCK_SIZE, MFQE_SCALE_SIZE,
                      high_bd, bitdepth);
  for (int i = 0; i < buf_size; i++) ASSERT_EQ(tmp->buffer[i], 100);
}
