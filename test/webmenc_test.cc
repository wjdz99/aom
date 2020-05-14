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

#include <string>
#include "common/webmenc.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

#if CONFIG_WEBM_IO

class WebmencTest : public ::testing::Test {};

// All of these variations on output should be identical.
TEST(WebmencTest, ExtractEncoderSettingsOutput1) {
  const char *argv[] = { "aomenc", "-o", "out.webm", "input.y4m",
                         "--target-bitrate=300" };
  int argc = 5;
  const std::string expected("version:1.2.3 --target-bitrate=300");
  char *result = extract_encoder_settings("1.2.3", argv, argc, "input.y4m");
  ASSERT_EQ(expected, std::string(result));
  free(result);
}

TEST(WebmencTest, ExtractEncoderSettingsOutput2) {
  const char *argv[] = { "aomenc", "--output", "out.webm", "input.y4m",
                         "--cpu-used=3" };
  int argc = 5;
  const std::string expected("version:abc --cpu-used=3");
  char *result = extract_encoder_settings("abc", argv, argc, "input.y4m");
  ASSERT_EQ(expected, std::string(result));
  free(result);
}

TEST(WebmencTest, ExtractEncoderSettingsOutput3) {
  const char *argv[] = { "aomenc", "--cq-level=63", "--end-usage=q",
                         "--output=out.webm", "input.y4m" };
  int argc = 5;
  const std::string expected("version:23 --cq-level=63 --end-usage=q");
  char *result = extract_encoder_settings("23", argv, argc, "input.y4m");
  ASSERT_EQ(expected, std::string(result));
  free(result);
}

TEST(WebmencTest, ExtractEncoderSettingsInput) {
  // Check that input filename is filtered regardless of position.
  const char *argv[] = { "aomenc", "-o", "out.webm", "input.y4m", "-p", "2" };
  int argc = 6;
  const char version[] = "1.0.0";
  const std::string expected("version:1.0.0 -p 2");
  char *result = extract_encoder_settings(version, argv, argc, "input.y4m");
  ASSERT_EQ(expected, std::string(result));
  free(result);

  const char *argv2[] = { "aomenc", "input.y4m", "-o", "out.webm", "-p", "2" };
  result = extract_encoder_settings(version, argv2, argc, "input.y4m");
  ASSERT_EQ(expected, std::string(result));
  free(result);
}

#endif  // CONFIG_WEBM_IO
}  // namespace
