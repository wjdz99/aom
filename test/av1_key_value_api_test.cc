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

#include <cstdlib>
#include <tuple>

#include "aom/aom_codec.h"
#include "aom/aom_decoder.h"
#include "aom/aomcx.h"
#include "aom/aomdx.h"
#include "common/args_helper.h"
#include "config/aom_config.h"
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {
typedef std::tuple<const char *, const char *> key_val_param;

class BaseKeyValAPI : public testing::Test {
 protected:
#if CONFIG_AV1_ENCODER
  aom_codec_iface_t *iface_cx;
  aom_codec_ctx_t enc;
  aom_codec_enc_cfg_t enc_cfg;
#endif
#if CONFIG_AV1_DECODER
  aom_codec_iface_t *iface_dx;
  aom_codec_ctx_t dec;
  aom_codec_dec_cfg_t dec_cfg;
#endif
  char err_msg[ARG_ERR_MSG_MAX_LEN];

 public:
  virtual void SetUp() override {
#if CONFIG_AV1_ENCODER
    iface_cx = aom_codec_av1_cx();
    EXPECT_EQ(AOM_CODEC_OK,
              aom_codec_enc_config_default(iface_cx, &enc_cfg, 0));
    EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface_cx, &enc_cfg, 0));
#endif
#if CONFIG_AV1_DECODER
    iface_dx = aom_codec_av1_dx();
    dec_cfg = { 0, 0, 0, !FORCE_HIGHBITDEPTH_DECODING };
    EXPECT_EQ(AOM_CODEC_OK, aom_codec_dec_init(&dec, iface_dx, &dec_cfg, 0));
#endif
  }

  virtual void TearDown() override {
#if CONFIG_AV1_ENCODER
    EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));
#endif
#if CONFIG_AV1_DECODER
    EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&dec));
#endif
  }
};

// Tests on encoder options.
// Need to add ones for the decoder in the future if it is also supported in the
// key & value API.
#if CONFIG_AV1_ENCODER
class EncValidTest : public BaseKeyValAPI,
                     public testing::WithParamInterface<key_val_param> {};
class EncInvalidTest : public BaseKeyValAPI,
                       public testing::WithParamInterface<key_val_param> {};

TEST_P(EncValidTest, WithMsg) {
  const char *key = std::get<0>(GetParam());
  const char *val = std::get<1>(GetParam());
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_set_option(&enc, key, val, err_msg));
}
TEST_P(EncValidTest, NullMsg) {
  const char *key = std::get<0>(GetParam());
  const char *val = std::get<1>(GetParam());
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_set_option(&enc, key, val, NULL));
}

TEST_P(EncInvalidTest, NullCtx) {
  const char *key = std::get<0>(GetParam());
  const char *val = std::get<1>(GetParam());
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_set_option(NULL, key, val, NULL));
}
TEST_P(EncInvalidTest, WithMsg) {
  const char *key = std::get<0>(GetParam());
  const char *val = std::get<1>(GetParam());
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_set_option(&enc, key, val, err_msg));
}
TEST_P(EncInvalidTest, NullMsg) {
  const char *key = std::get<0>(GetParam());
  const char *val = std::get<1>(GetParam());
  EXPECT_EQ(AOM_CODEC_INVALID_PARAM,
            aom_codec_set_option(&enc, key, val, NULL));
}

// No test for ratio / list for now since the API does not support any of the
// parameters of these kinds.
const key_val_param enc_valid_params[] = {
  std::make_tuple("min-gf-interval", "10"),          // uint
  std::make_tuple("min-partition-size", "4"),        // int
  std::make_tuple("tune", "psnr"),                   // enum
  std::make_tuple("film-grain-table", "table.txt"),  // string
};

const key_val_param enc_invalid_params[] = {
  // NULL
  std::make_tuple("min-gf-interval", (const char *)NULL),
  std::make_tuple((const char *)NULL, "0"),
  std::make_tuple((const char *)NULL, (const char *)NULL),
  // no match
  std::make_tuple("a-b-c", "10"),
  // uint
  std::make_tuple("min-gf-interval", "-1"),
  std::make_tuple("min-gf-interval", "1.1"),
  std::make_tuple("min-gf-interval", "abc"),
  // int
  std::make_tuple("min-partition-size", "1.1"),
  std::make_tuple("min-partition-size", "abc"),
  // enum
  std::make_tuple("tune", "PsnR1"),
  // out of range
  std::make_tuple("cq-level", "1000"),
};

INSTANTIATE_TEST_SUITE_P(KeyValAPI, EncValidTest,
                         testing::ValuesIn(enc_valid_params));
INSTANTIATE_TEST_SUITE_P(KeyValAPI, EncInvalidTest,
                         testing::ValuesIn(enc_invalid_params));

GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(KeyValInvalidTest);
GTEST_ALLOW_UNINSTANTIATED_PARAMETERIZED_TEST(KeyValValidTest);
#endif  // CONFIG_AV1_ENCODER

}  // namespace
