/*
 * Copyright (c) 2018, Alliance for Open Media. All rights reserved
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
#include "test/encode_test_driver.h"
#include "test/util.h"
#include "test/y4m_video_source.h"

#include "aom/internal/aom_codec_internal.h"
#include "aom/aom_decoder.h"

#include "av1/common/onyxc_int.h"
#include "av1/decoder/decoder.h"

namespace {

typedef struct aom_codec_alg_priv {
  aom_codec_priv_t base;
  aom_codec_dec_cfg_t cfg;
  aom_codec_stream_info_t si;
  int postproc_cfg_set;
  aom_postproc_cfg_t postproc_cfg;
  aom_image_t img;
  int img_avail;
  int flushed;
  int invert_tile_order;
  int last_show_frame;  // Index of last output frame.
  int byte_alignment;
  int skip_loop_filter;
  int decode_tile_row;
  int decode_tile_col;
  unsigned int tile_mode;
  unsigned int ext_tile_debug;
  unsigned int row_mt;
  EXTERNAL_REFERENCES ext_refs;
  unsigned int is_annexb;
  int operating_point;
  int output_all_layers;

  AVxWorker *frame_workers;
  int num_frame_workers;
  int next_submit_worker_id;
  int last_submit_worker_id;
  int next_output_worker_id;
  int available_threads;
  aom_image_t *image_with_grain[MAX_NUM_SPATIAL_LAYERS];
  int need_resync;  // wait for key/intra-only frame
  // BufferPool that holds all reference frames. Shared by all the FrameWorkers.
  BufferPool *buffer_pool;

  // External frame buffer info to save for AV1 common.
  void *ext_priv;  // Private data associated with the external frame buffers.
  aom_get_frame_buffer_cb_fn_t get_ext_fb_cb;
  aom_release_frame_buffer_cb_fn_t release_ext_fb_cb;

#if CONFIG_INSPECTION
  aom_inspect_cb inspect_cb;
  void *inspect_ctx;
#endif
} aom_codec_alg_priv_t;

const int kCpuUsed = 2;

struct EncodePerfTestVideo {
  const char *name;
  uint32_t width;
  uint32_t height;
  uint32_t bitrate;
  int frames;
};

const EncodePerfTestVideo kAV1EncodePerfTestVectors[] = {
  { "niklas_1280_720_30.y4m", 1280, 720, 600, 3 },
};

struct EncodeParameters {
  int32_t tile_rows;
  int32_t tile_cols;
  int32_t lossless;
  int32_t error_resilient;
  int32_t frame_parallel;
  /*
  aom_color_range_t color_range;
  aom_color_space_t cs;
  */
  int render_size[2];
};

const EncodeParameters kAV1EncodeParameterSet[] = {
  { 0, 0, 0, 1, 0, /*AOM_CR_STUDIO_RANGE, AOM_CS_BT_601,*/ { 0, 0 } },
  { 0, 0, 0, 0, 0, /*AOM_CR_FULL_RANGE, AOM_CS_BT_709,*/ { 0, 0 } },
  { 0, 0, 1, 0, 0, /*AOM_CR_FULL_RANGE, AOM_CS_BT_2020,*/ { 0, 0 } },
  { 0, 2, 0, 0, 1, /*AOM_CR_STUDIO_RANGE, AOM_CS_UNKNOWN,*/ { 640, 480 } },
};

class AVxEncoderParmsGetToDecoder
    : public ::libaom_test::EncoderTest,
      public ::libaom_test::CodecTestWith2Params<EncodeParameters,
                                                 EncodePerfTestVideo> {
 protected:
  AVxEncoderParmsGetToDecoder()
      : EncoderTest(GET_PARAM(0)), encode_parms(GET_PARAM(1)) {}

  virtual ~AVxEncoderParmsGetToDecoder() {}

  virtual void SetUp() {
    InitializeConfig();
    SetMode(::libaom_test::kTwoPassGood);
    cfg_.g_lag_in_frames = 25;
    cfg_.g_error_resilient = encode_parms.error_resilient;
    test_video_ = GET_PARAM(2);
    cfg_.rc_target_bitrate = test_video_.bitrate;
  }

  virtual void PreEncodeFrameHook(::libaom_test::VideoSource *video,
                                  ::libaom_test::Encoder *encoder) {
    if (video->frame() == 1) {
      /*
      encoder->Control(AV1E_SET_COLOR_SPACE, encode_parms.cs);
      encoder->Control(AV1E_SET_COLOR_RANGE, encode_parms.color_range);
      */
      encoder->Control(AV1E_SET_LOSSLESS, encode_parms.lossless);
      encoder->Control(AV1E_SET_FRAME_PARALLEL_DECODING,
                       encode_parms.frame_parallel);
      encoder->Control(AV1E_SET_TILE_ROWS, encode_parms.tile_rows);
      encoder->Control(AV1E_SET_TILE_COLUMNS, encode_parms.tile_cols);
      encoder->Control(AOME_SET_CPUUSED, kCpuUsed);
      encoder->Control(AOME_SET_ENABLEAUTOALTREF, 1);
      encoder->Control(AOME_SET_ARNR_MAXFRAMES, 7);
      encoder->Control(AOME_SET_ARNR_STRENGTH, 5);
      if (encode_parms.render_size[0] > 0 && encode_parms.render_size[1] > 0) {
        encoder->Control(AV1E_SET_RENDER_SIZE, encode_parms.render_size);
      }
    }
  }

  virtual bool HandleDecodeResult(const aom_codec_err_t res_dec,
                                  libaom_test::Decoder *decoder) {
    aom_codec_ctx_t *const av1_decoder = decoder->GetDecoder();
    aom_codec_alg_priv_t *const priv =
        reinterpret_cast<aom_codec_alg_priv_t *>(av1_decoder->priv);
    AVxWorker *const worker = priv->frame_workers;
    FrameWorkerData *const frame_worker_data = (FrameWorkerData *)worker->data1;
    const AV1Decoder *pbi = frame_worker_data->pbi;
    const AV1_COMMON *common = &pbi->common;

    if (encode_parms.lossless) {
      EXPECT_EQ(0, common->base_qindex);
      EXPECT_EQ(0, common->y_dc_delta_q);
      EXPECT_EQ(0, common->u_dc_delta_q);
      EXPECT_EQ(0, common->u_ac_delta_q);
      EXPECT_EQ(0, common->v_dc_delta_q);
      EXPECT_EQ(0, common->v_ac_delta_q);
      EXPECT_EQ(ONLY_4X4, common->tx_mode);
    }
    EXPECT_EQ(encode_parms.error_resilient, common->error_resilient_mode);
    if (encode_parms.error_resilient) {
      // EXPECT_EQ(1, common->frame_parallel_decoding_mode);
      // EXPECT_EQ(0, common->use_prev_frame_mvs);
    } else {
      // EXPECT_EQ(encode_parms.frame_parallel,
      // common->frame_parallel_decoding_mode);
    }
    // EXPECT_EQ(encode_parms.color_range, common->color_range);
    // EXPECT_EQ(encode_parms.cs, common->color_space);
    if (encode_parms.render_size[0] > 0 && encode_parms.render_size[1] > 0) {
      EXPECT_EQ(encode_parms.render_size[0], common->render_width);
      EXPECT_EQ(encode_parms.render_size[1], common->render_height);
    }
    EXPECT_EQ(encode_parms.tile_cols, common->log2_tile_cols);
    EXPECT_EQ(encode_parms.tile_rows, common->log2_tile_rows);

    EXPECT_EQ(AOM_CODEC_OK, res_dec) << decoder->DecodeError();
    return AOM_CODEC_OK == res_dec;
  }

  EncodePerfTestVideo test_video_;

 private:
  EncodeParameters encode_parms;
};

TEST_P(AVxEncoderParmsGetToDecoder, BitstreamParms) {
  init_flags_ = AOM_CODEC_USE_PSNR;

  testing::internal::scoped_ptr<libaom_test::VideoSource> video(
      new libaom_test::Y4mVideoSource(test_video_.name, 0, test_video_.frames));
  ASSERT_TRUE(video.get() != NULL);

  ASSERT_NO_FATAL_FAILURE(RunLoop(video.get()));
}

AV1_INSTANTIATE_TEST_CASE(AVxEncoderParmsGetToDecoder,
                          ::testing::ValuesIn(kAV1EncodeParameterSet),
                          ::testing::ValuesIn(kAV1EncodePerfTestVectors));
}  // namespace
