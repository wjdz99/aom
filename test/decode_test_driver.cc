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

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "test/codec_factory.h"
#include "test/decode_test_driver.h"
#include "test/register_state_check.h"
#include "test/video_source.h"

namespace libaom_test {

const char kVP8Name[] = "WebM Project VP8";
const char kAV1Name[] = "AOMedia Project AV1 Decoder";

AomCodecErrT Decoder::PeekStream(const uint8_t *cxdata, size_t size,
                                 AomCodecStreamInfoT *stream_info) {
  return aom_codec_peek_stream_info(
      CodecInterface(), cxdata, static_cast<unsigned int>(size), stream_info);
}

AomCodecErrT Decoder::DecodeFrame(const uint8_t *cxdata, size_t size) {
  return DecodeFrame(cxdata, size, NULL);
}

AomCodecErrT Decoder::DecodeFrame(const uint8_t *cxdata, size_t size,
                                  void *user_priv) {
  AomCodecErrT res_dec;
  InitOnce();
  API_REGISTER_STATE_CHECK(
      res_dec = aom_codec_decode(
          &decoder_, cxdata, static_cast<unsigned int>(size), user_priv, 0));
  return res_dec;
}

bool Decoder::IsVP8() const {
  const char *codec_name = GetDecoderName();
  return strncmp(kVP8Name, codec_name, sizeof(kVP8Name) - 1) == 0;
}

bool Decoder::IsAV1() const {
  const char *codec_name = GetDecoderName();
  return strncmp(kAV1Name, codec_name, sizeof(kAV1Name) - 1) == 0;
}

void DecoderTest::HandlePeekResult(Decoder *const decoder,
                                   CompressedVideoSource *video,
                                   const AomCodecErrT res_peek) {
  const bool is_vp8 = decoder->IsVP8();
  if (is_vp8) {
    /* Vp8's implementation of PeekStream returns an error if the frame you
     * pass it is not a keyframe, so we only expect AOM_CODEC_OK on the first
     * frame, which must be a keyframe. */
    if (video->frame_number() == 0)
      ASSERT_EQ(AOM_CODEC_OK, res_peek)
          << "Peek return failed: " << aom_codec_err_to_string(res_peek);
  } else {
    /* The Av1 implementation of PeekStream returns an error only if the
     * data passed to it isn't a valid Av1 chunk. */
    ASSERT_EQ(AOM_CODEC_OK, res_peek)
        << "Peek return failed: " << aom_codec_err_to_string(res_peek);
  }
}

void DecoderTest::RunLoop(CompressedVideoSource *video,
                          const AomCodecDecCfgT &dec_cfg) {
  Decoder *const decoder = codec_->CreateDecoder(dec_cfg, flags_);
  ASSERT_TRUE(decoder != NULL);
  bool end_of_file = false;

  // Decode frames.
  for (video->Begin(); !::testing::Test::HasFailure() && !end_of_file;
       video->Next()) {
    PreDecodeFrameHook(*video, decoder);

    AomCodecStreamInfoT stream_info;
    if (video->cxdata() != NULL) {
      const AomCodecErrT res_peek = decoder->PeekStream(
          video->cxdata(), video->frame_size(), &stream_info);
      HandlePeekResult(decoder, video, res_peek);
      ASSERT_FALSE(::testing::Test::HasFailure());

      AomCodecErrT res_dec =
          decoder->DecodeFrame(video->cxdata(), video->frame_size());
      if (!HandleDecodeResult(res_dec, decoder)) break;
    } else {
      // Signal end of the file to the decoder.
      const AomCodecErrT res_dec = decoder->DecodeFrame(NULL, 0);
      ASSERT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
      end_of_file = true;
    }

    DxDataIterator dec_iter = decoder->GetDxData();
    const AomImageT *img = NULL;

    // Get decompressed data
    while ((img = dec_iter.Next()))
      DecompressedFrameHook(*img, video->frame_number());
  }
  delete decoder;
}

void DecoderTest::RunLoop(CompressedVideoSource *video) {
  AomCodecDecCfgT dec_cfg = AomCodecDecCfgT();
  RunLoop(video, dec_cfg);
}

void DecoderTest::set_cfg(const AomCodecDecCfgT &dec_cfg) {
  memcpy(&cfg_, &dec_cfg, sizeof(cfg_));
}

void DecoderTest::set_flags(const AomCodecFlagsT flags) { flags_ = flags; }

}  // namespace libaom_test
