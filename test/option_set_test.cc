/*
 * Copyright (c) 2024, Alliance for Open Media. All rights reserved.
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "gtest/gtest.h"

#include "aom/aomcx.h"
#include "aom/aom_encoder.h"
#include "config/aom_config.h"

namespace {

TEST(OptionSetTest, Ssimulacra2) {
  aom_codec_iface_t *iface = aom_codec_av1_cx();
  aom_codec_enc_cfg_t cfg;
  EXPECT_EQ(AOM_CODEC_OK,
            aom_codec_enc_config_default(iface, &cfg, AOM_USAGE_ALL_INTRA));
  aom_codec_ctx_t enc;
  EXPECT_EQ(AOM_CODEC_OK, aom_codec_enc_init(&enc, iface, &cfg, 0));
  EXPECT_EQ(AOM_CODEC_OK,
            aom_codec_set_option(&enc, "option-set", "ssimulacra2"));

  EXPECT_EQ(AOM_CODEC_OK, aom_codec_destroy(&enc));
}

}  // namespace
