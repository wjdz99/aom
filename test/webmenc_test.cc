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
#include <string>
#include "common/webmenc.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

#if CONFIG_WEBM_IO

class WebmencTest : public ::testing::Test {};

TEST(WebmencTest, ExtractEncoderSettingsOutput) {
  // All of these variations on output should be identical.
  const char *argv1[] = { "aomenc",        "-o",
                          "out.webm",      "input.y4m",
                          "--cq-level=20", "--end-usage=q" };
  int argc1 = 6;
  const char *argv2[] = { "aomenc",    "--output",      "out.webm",
                          "input.y4m", "--cq-level=20", "--end-usage=q" };
  int argc2 = 6;
  const char *argv3[] = { "aomenc", "--output=out.webm", "input.y4m",
                          "--cq-level=20", "--end-usage=q" };
  int argc3 = 5;
  const char version[] = "1.0.0";
  const std::string expected("version:1.0.0 --cq-level=20 --end-usage=q");
  std::unique_ptr<char> result(
      extract_encoder_settings(version, argv1, argc1, "input.y4m"));
  ASSERT_EQ(expected, std::string(result.get()));
  result.reset(extract_encoder_settings(version, argv2, argc2, "input.y4m"));
  ASSERT_EQ(expected, std::string(result.get()));
  result.reset(extract_encoder_settings(version, argv3, argc3, "input.y4m"));
  ASSERT_EQ(expected, std::string(result.get()));
}

TEST(WebmencTest, ExtractEncoderSettingsInput) {
  // Check that input filename is filtered.
  const char *argv[] = { "aomenc", "-o", "out.webm", "input.y4m", "-p", "2" };
  int argc = 6;
  const char version[] = "1.0.0";
  const std::string expected("version:1.0.0 -p 2");
  std::unique_ptr<char> result(
      extract_encoder_settings(version, argv, argc, "input.y4m"));
  ASSERT_EQ(expected, std::string(result.get()));
}

#endif  // CONFIG_WEBM_IO
}  // namespace
