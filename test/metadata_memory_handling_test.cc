#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "aom/aom_codec.h"
#include "aom/internal/aom_image_internal.h"
#include "aom_scale/yv12config.h"

TEST(MetadataMemoryHandlingTest, MetadataAllocation) {
  aom_metadata_t *metadata;
  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  metadata = aom_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);
  ASSERT_TRUE(metadata != NULL);

  int status = aom_metadata_free(metadata);
  EXPECT_EQ(status, 0);
}

TEST(MetadataMemoryHandlingTest, MetadataArrayAllocation) {
  aom_metadata_array_t *metadata_array;
  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  metadata_array = aom_metadata_array_alloc(2);
  ASSERT_TRUE(metadata_array != NULL);

  metadata_array->buffer[0] =
      aom_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);
  metadata_array->buffer[1] =
      aom_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);

  int status = aom_metadata_array_free(metadata_array);
  EXPECT_EQ(status, 2);
}

TEST(MetadataMemoryHandlingTest, AddMetadataToImage) {
  aom_image_t image;
  image.metadata = NULL;

  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  int status =
      aom_add_metadata_to_img(&image, OBU_METADATA_TYPE_ITUT_T35, data, 10);
  ASSERT_EQ(status, 0);

  aom_metadata_array_t *array = (aom_metadata_array_t *)image.metadata;

  status = aom_metadata_array_free(array);
  EXPECT_GE(status, 0);
}

TEST(MetadataMemoryHandlingTest, RemoveMetadataFromImage) {
  aom_image_t image;
  image.metadata = NULL;

  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  int status =
      aom_add_metadata_to_img(&image, OBU_METADATA_TYPE_ITUT_T35, data, 10);
  ASSERT_EQ(status, 0);

  status = aom_remove_metadata_from_img(&image);
  EXPECT_EQ(status, 0);
}

TEST(MetadataMemoryHandlingTest, CopyMetadataToFrameBUffer) {
  YV12_BUFFER_CONFIG yvBuf;
  yvBuf.metadata = NULL;
  aom_metadata_array_t *metadata_array;
  uint8_t data[10] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

  metadata_array = aom_metadata_array_alloc(1);
  ASSERT_TRUE(metadata_array != NULL);

  metadata_array->buffer[0] =
      aom_metadata_alloc(OBU_METADATA_TYPE_ITUT_T35, data, 10);

  int status = aom_copy_metadata_to_frame_buffer(&yvBuf, metadata_array);
  EXPECT_EQ(status, 0);

  status = aom_remove_metadata_from_frame_buffer(&yvBuf);
  EXPECT_EQ(status, 0);
}
