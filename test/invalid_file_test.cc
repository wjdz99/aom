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

#include <string>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"
#include "test/codec_factory.h"
#include "test/ivf_video_source.h"
#include "test/util.h"
#include "test/video_source.h"

namespace {

struct DecodeParam {
  int threads;
  const char *filename;
};

std::ostream &operator<<(std::ostream &os, const DecodeParam &dp) {
  return os << "threads: " << dp.threads << " file: " << dp.filename;
}

class InvalidFileTest : public ::libaom_test::DecoderTest,
                        public ::libaom_test::CodecTestWithParam<DecodeParam> {
 protected:
  InvalidFileTest() : DecoderTest(GET_PARAM(0)) {}

  virtual ~InvalidFileTest() {}

  bool HandleDecodeResult(const aom_codec_err_t res_dec,
                          libaom_test::Decoder *decoder) override {
    EXPECT_EQ(expected_res_, res_dec) << decoder->DecodeError();
    return expected_res_ == res_dec;
  }

  void HandlePeekResult(libaom_test::Decoder *const /*decoder*/,
                        libaom_test::CompressedVideoSource * /*video*/,
                        const aom_codec_err_t /*res_peek*/) override {}

  void RunTest() {
    const DecodeParam input = GET_PARAM(1);
    aom_codec_dec_cfg_t cfg = { 0, 0, 0, CONFIG_LOWBITDEPTH, { 1 } };
    cfg.threads = input.threads;
    const std::string filename = input.filename;
    libaom_test::IVFVideoSource decode_video(filename);
    decode_video.Init();

    expected_res_ = AOM_CODEC_CORRUPT_FRAME;
    ASSERT_NO_FATAL_FAILURE(RunLoop(&decode_video, cfg));
  }

 private:
  int expected_res_;
};

TEST_P(InvalidFileTest, ReturnCode) { RunTest(); }

const DecodeParam kAV1InvalidFileTests[] = {
  { 1, "invalid-bug-1814.ivf" },
};

AV1_INSTANTIATE_TEST_CASE(InvalidFileTest,
                          ::testing::ValuesIn(kAV1InvalidFileTests));

}  // namespace
