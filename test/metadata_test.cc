#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "aom/aom_codec.h"
#include "aom/internal/aom_image_internal.h"
#include "aom_scale/yv12config.h"
#include "av1/encoder/bitstream.h"

#if CONFIG_AV1_ENCODER
#include "test/codec_factory.h"
#include "test/encode_test_driver.h"
#include "test/i420_video_source.h"
#include "test/util.h"
#include "test/video_source.h"

namespace {

class MetadataEncodeTest
    : public ::libaom_test::CodecTestWithParam<libaom_test::TestMode>,
      public ::libaom_test::EncoderTest {
 protected:
  MetadataEncodeTest() : EncoderTest(GET_PARAM(0)) {}

  virtual ~MetadataEncodeTest() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(GET_PARAM(1));
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video) {
    aom_image_t *current_frame = video->img();
    if (current_frame) {
      if (current_frame->metadata) aom_img_remove_metadata(current_frame);
      const size_t dataSize = 10;
      uint8_t metadata[dataSize] = { 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 };
      ASSERT_EQ(aom_img_add_metadata(current_frame, OBU_METADATA_TYPE_ITUT_T35,
                                     metadata, dataSize),
                0);
    }
  }

  virtual void FramePktHook(const aom_codec_cx_pkt_t *pkt) {
    if (pkt->kind == AOM_CODEC_CX_FRAME_PKT) {
      const size_t metadataObuSize = 14;
      const uint8_t metadataObu[metadataObuSize] = { 42, 12, 4, 7, 7, 7, 7,
                                                     7,  7,  7, 7, 7, 7, 128 };
      const size_t bitstreamSize = pkt->data.frame.sz;
      const uint8_t *bitstream = (const uint8_t *)pkt->data.frame.buf;
      size_t j = 0;
      // look for expected metadata in bitstream
      for (size_t i = 0; i < bitstreamSize && j < metadataObuSize; ++i, ++j) {
        if (bitstream[i] != metadataObu[j]) j = 0;
      }
      EXPECT_EQ(j, metadataObuSize);
    }
  }
};

TEST_P(MetadataEncodeTest, TestMetadataEncoding) {
  ::libaom_test::I420VideoSource video("hantro_collage_w352h288.yuv", 352, 288,
                                       30, 1, 0, 5);
  init_flags_ = AOM_CODEC_USE_PSNR;

  cfg_.g_w = 352;
  cfg_.g_h = 288;

  cfg_.rc_buf_initial_sz = 500;
  cfg_.rc_buf_optimal_sz = 600;
  cfg_.rc_buf_sz = 1000;
  cfg_.rc_min_quantizer = 2;
  cfg_.rc_max_quantizer = 56;
  cfg_.rc_undershoot_pct = 50;
  cfg_.rc_overshoot_pct = 50;
  cfg_.rc_end_usage = AOM_CBR;
  cfg_.kf_mode = AOM_KF_AUTO;
  cfg_.g_lag_in_frames = 1;
  cfg_.kf_min_dist = cfg_.kf_max_dist = 3000;
  // Enable dropped frames.
  cfg_.rc_dropframe_thresh = 1;
  // Disable error_resilience mode.
  cfg_.g_error_resilient = 0;
  // Run at low bitrate.
  cfg_.rc_target_bitrate = 40;

  ASSERT_NO_FATAL_FAILURE(RunLoop(&video));
}

AV1_INSTANTIATE_TEST_CASE(MetadataEncodeTest,
                          ::testing::Values(::libaom_test::kOnePassGood));

}  // namespace
#endif

TEST(MetadataTest, MetadataAllocation) {
  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  aom_metadata_t *metadata =
      aom_img_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);
  ASSERT_NE(metadata, nullptr);
  EXPECT_EQ(aom_img_metadata_free(metadata), 0);
}

TEST(MetadataTest, MetadataArrayAllocation) {
  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  aom_metadata_array_t *metadata_array = aom_img_metadata_array_alloc(2);
  ASSERT_NE(metadata_array, nullptr);

  metadata_array->metadata_array[0] =
      aom_img_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);
  metadata_array->metadata_array[1] =
      aom_img_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);

  EXPECT_EQ(aom_img_metadata_array_free(metadata_array), 2u);
}

TEST(MetadataTest, AddMetadataToImage) {
  aom_image_t image;
  image.metadata = NULL;

  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  ASSERT_EQ(aom_img_add_metadata(&image, OBU_METADATA_TYPE_ITUT_T35, data, 10),
            0);
  EXPECT_EQ(aom_img_metadata_array_free(image.metadata), 1u);
  EXPECT_EQ(aom_img_add_metadata(NULL, OBU_METADATA_TYPE_ITUT_T35, data, 10),
            -1);
}

TEST(MetadataTest, RemoveMetadataFromImage) {
  aom_image_t image;
  image.metadata = NULL;

  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  ASSERT_EQ(aom_img_add_metadata(&image, OBU_METADATA_TYPE_ITUT_T35, data, 10),
            0);
  EXPECT_EQ(aom_img_remove_metadata(&image), 1u);
  EXPECT_EQ(aom_img_remove_metadata(NULL), 0u);
}

TEST(MetadataTest, CopyMetadataToFrameBuffer) {
  YV12_BUFFER_CONFIG yvBuf;
  yvBuf.metadata = NULL;
  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  aom_metadata_array_t *metadata_array = aom_img_metadata_array_alloc(1);
  ASSERT_NE(metadata_array, nullptr);

  metadata_array->metadata_array[0] =
      aom_img_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);

  // Metadata_array
  int status = aom_copy_metadata_to_frame_buffer(&yvBuf, metadata_array);
  EXPECT_EQ(status, 0);
  status = aom_copy_metadata_to_frame_buffer(NULL, metadata_array);
  EXPECT_EQ(status, -1);
  EXPECT_EQ(aom_img_metadata_array_free(metadata_array), 1u);

  // Metadata_array_2
  aom_metadata_array_t *metadata_array_2 = aom_img_metadata_array_alloc(0);
  ASSERT_NE(metadata_array_2, nullptr);
  status = aom_copy_metadata_to_frame_buffer(&yvBuf, metadata_array_2);
  EXPECT_EQ(status, -1);
  EXPECT_EQ(aom_img_metadata_array_free(metadata_array_2), 0u);

  // YV12_BUFFER_CONFIG
  status = aom_copy_metadata_to_frame_buffer(&yvBuf, NULL);
  EXPECT_EQ(status, -1);
  EXPECT_EQ(aom_remove_metadata_from_frame_buffer(NULL), 0u);
  EXPECT_EQ(aom_remove_metadata_from_frame_buffer(&yvBuf), 1u);
}
