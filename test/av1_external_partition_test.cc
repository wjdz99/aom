/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "aom/aom_codec.h"
#include "aom/aom_external_partition.h"
#include "av1/common/blockd.h"
//#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/y4m_video_source.h"
#include "test/util.h"

#if CONFIG_AV1_ENCODER
namespace {

constexpr int kFrameNum = 8;

aom_ext_part_status_t ext_part_create_model(
    void *priv, const aom_ext_part_config_t *part_config,
    aom_ext_part_model_t *ext_part_model_ptr) {
  (void)priv;
  (void)part_config;
  (void)ext_part_model_ptr;
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_send_features(
    aom_ext_part_model_t *ext_part_model_ptr,
    const aom_partition_features_t *part_features) {
  (void)ext_part_model_ptr;
  (void)part_features;
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_get_partition_decision(
    aom_ext_part_model_t *ext_part_model_ptr,
    aom_partition_decision_t *ext_part_decision) {
  (void)ext_part_model_ptr;
  (void)ext_part_decision;
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_send_partition_stats(
    aom_ext_part_model_t *ext_part_model_ptr,
    const aom_partition_stats_t *ext_part_stats) {
  (void)ext_part_model_ptr;
  (void)ext_part_stats;
  return AOM_EXT_PART_OK;
}

aom_ext_part_status_t ext_part_delete_model(
    aom_ext_part_model_t *ext_part_model_ptr) {
  (void)ext_part_model_ptr;
  return AOM_EXT_PART_OK;
}

class ExternalPartitionTest : public ::testing::Test,
                              public ::libaom_test::EncoderTest {
 protected:
  ExternalPartitionTest() : EncoderTest(&::libaom_test::kAV1) {}
  virtual ~ExternalPartitionTest() {}

  virtual void SetUp() {
    InitializeConfig(::libaom_test::kTwoPassGood);
    const aom_rational timebase = { 1, 30 };
    cfg_.g_timebase = timebase;
    cfg_.rc_end_usage = AOM_VBR;
    cfg_.g_threads = 1;
    cfg_.g_lag_in_frames = 35;
    cfg_.rc_target_bitrate = 400;
  }

  virtual bool DoDecode() const { return 0; }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 0) {
      aom_ext_part_funcs_t ext_part_funcs;
      ext_part_funcs.create_model = ext_part_create_model;
      ext_part_funcs.send_features = ext_part_send_features;
      ext_part_funcs.get_partition_decision = ext_part_get_partition_decision;
      ext_part_funcs.send_partition_stats = ext_part_send_partition_stats;
      ext_part_funcs.delete_model = ext_part_delete_model;

      encoder->Control(AOME_SET_CPUUSED, 4);
      encoder->Control(AV1E_SET_EXTERNAL_PARTITION, &ext_part_funcs);
    }
  }
};

TEST_F(ExternalPartitionTest, EncodeTest) {
  ::libaom_test::Y4mVideoSource video("niklas_1280_720_30.y4m", 0, kFrameNum);
  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
}

}  // namespace
#endif
