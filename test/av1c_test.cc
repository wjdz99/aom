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
#include "common/av1c.h"

#include <string.h>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

//
// Input buffers containing exactly one Sequence Header OBU.
//
// Each buffer is named according to the OBU storage format (Annex-B vs Low
// Overhead Bitstream Format) and the type of Sequence Header OBU ("Full"
// Sequence Header OBUs vs Sequence Header OBUs with the
// reduced_still_image_flag set).
//
const uint8_t kAnnexBSequenceHeaderObu[] = {
  0x0c, 0x08, 0x00, 0x00, 0x00, 0x04, 0x45, 0x7e, 0x3e, 0xff, 0xfc, 0xc0, 0x20
};
const uint8_t kAnnexBReducedStillImageSequenceHeaderObu[] = {
  0x08, 0x08, 0x18, 0x22, 0x2b, 0xf1, 0xfe, 0xc0, 0x20
};

const uint8_t kLobfReducedStillImageSequenceHeaderObu[] = {
  0x0a, 0x07, 0x18, 0x22, 0x2b, 0xf1, 0xfe, 0xc0, 0x20
};

const uint8_t kLobfSequenceHeaderObu[] = {
  0x0a, 0x0b, 0x00, 0x00, 0x00, 0x04, 0x45, 0x7e, 0x3e, 0xff, 0xfc, 0xc0, 0x20
};

const uint8_t kAv1cAllZero[] = { 0, 0, 0, 0 };

// The size of AV1 config when no configOBUs are present at the end of the
// configuration structure.
const size_t kAv1cNoConfigObusSize = 4;

TEST(Av1c, ObuInvalidInputs) {
  av1c av1_config;
  memset(&av1_config, 0, sizeof(av1_config));
  ASSERT_EQ(-1, get_av1c_from_obu(NULL, 0, 0, NULL));
  ASSERT_EQ(-1, get_av1c_from_obu(&kLobfSequenceHeaderObu[0], 0, 0, NULL));
  ASSERT_EQ(-1, get_av1c_from_obu(&kLobfSequenceHeaderObu[0],
                                  sizeof(kLobfSequenceHeaderObu), 0, NULL));
  ASSERT_EQ(-1,
            get_av1c_from_obu(NULL, sizeof(kLobfSequenceHeaderObu), 0, NULL));
  ASSERT_EQ(-1,
            get_av1c_from_obu(&kLobfSequenceHeaderObu[0], 0, 0, &av1_config));
}

TEST(Av1c, ReadInvalidInputs) {
  av1c av1_config;
  memset(&av1_config, 0, sizeof(av1_config));
  size_t bytes_read = 0;
  ASSERT_EQ(-1, read_av1c(NULL, 0, NULL, NULL));
  ASSERT_EQ(-1, read_av1c(NULL, 4, NULL, NULL));
  ASSERT_EQ(-1, read_av1c(&kAv1cAllZero[0], 0, NULL, NULL));
  ASSERT_EQ(-1, read_av1c(&kAv1cAllZero[0], 4, &bytes_read, NULL));
  ASSERT_EQ(-1, read_av1c(NULL, 4, &bytes_read, &av1_config));
}

TEST(Av1c, WriteInvalidInputs) {
  av1c av1_config;
  memset(&av1_config, 0, sizeof(av1_config));
  size_t bytes_written = 0;
  uint8_t av1c_buffer[4] = { 0 };
  ASSERT_EQ(-1, write_av1c(NULL, 0, NULL, NULL));
  ASSERT_EQ(-1, write_av1c(&av1_config, 0, NULL, NULL));
  ASSERT_EQ(-1, write_av1c(&av1_config, 0, &bytes_written, NULL));

  ASSERT_EQ(-1, write_av1c(&av1_config, 0, &bytes_written, &av1c_buffer[0]));
  ASSERT_EQ(-1, write_av1c(&av1_config, 4, &bytes_written, NULL));
}

TEST(Av1c, GetAv1cFromLobfObu) {
  av1c av1_config;
  memset(&av1_config, 0, sizeof(av1_config));

  // Test parsing of a Sequence Header OBU with the reduced_still_picture_header
  // unset-- aka a full Sequence Header OBU.
  ASSERT_EQ(0,
            get_av1c_from_obu(&kLobfSequenceHeaderObu[0],
                              sizeof(kLobfSequenceHeaderObu), 0, &av1_config));
  ASSERT_EQ(1, av1_config.marker);
  ASSERT_EQ(1, av1_config.version);
  ASSERT_EQ(0, av1_config.seq_profile);
  ASSERT_EQ(0, av1_config.seq_level_idx_0);
  ASSERT_EQ(0, av1_config.seq_tier_0);
  ASSERT_EQ(0, av1_config.high_bitdepth);
  ASSERT_EQ(0, av1_config.twelve_bit);
  ASSERT_EQ(0, av1_config.monochrome);
  ASSERT_EQ(1, av1_config.chroma_subsampling_x);
  ASSERT_EQ(1, av1_config.chroma_subsampling_y);
  ASSERT_EQ(0, av1_config.chroma_sample_position);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_present);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_minus_one);

  // Test parsing of a reduced still image Sequence Header OBU.
  memset(&av1_config, 0, sizeof(av1_config));
  ASSERT_EQ(
      0, get_av1c_from_obu(&kLobfReducedStillImageSequenceHeaderObu[0],
                           sizeof(kLobfReducedStillImageSequenceHeaderObu), 0,
                           &av1_config));
  ASSERT_EQ(1, av1_config.marker);
  ASSERT_EQ(1, av1_config.version);
  ASSERT_EQ(0, av1_config.seq_profile);
  ASSERT_EQ(0, av1_config.seq_level_idx_0);
  ASSERT_EQ(0, av1_config.seq_tier_0);
  ASSERT_EQ(0, av1_config.high_bitdepth);
  ASSERT_EQ(0, av1_config.twelve_bit);
  ASSERT_EQ(0, av1_config.monochrome);
  ASSERT_EQ(1, av1_config.chroma_subsampling_x);
  ASSERT_EQ(1, av1_config.chroma_subsampling_y);
  ASSERT_EQ(0, av1_config.chroma_sample_position);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_present);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_minus_one);
}

TEST(Av1c, GetAv1cFromAnnexBObu) {
  av1c av1_config;
  memset(&av1_config, 0, sizeof(av1_config));

  // Test parsing of a Sequence Header OBU with the reduced_still_picture_header
  // unset-- aka a full Sequence Header OBU.
  ASSERT_EQ(
      0, get_av1c_from_obu(&kAnnexBSequenceHeaderObu[0],
                           sizeof(kAnnexBSequenceHeaderObu), 1, &av1_config));
  ASSERT_EQ(1, av1_config.marker);
  ASSERT_EQ(1, av1_config.version);
  ASSERT_EQ(0, av1_config.seq_profile);
  ASSERT_EQ(0, av1_config.seq_level_idx_0);
  ASSERT_EQ(0, av1_config.seq_tier_0);
  ASSERT_EQ(0, av1_config.high_bitdepth);
  ASSERT_EQ(0, av1_config.twelve_bit);
  ASSERT_EQ(0, av1_config.monochrome);
  ASSERT_EQ(1, av1_config.chroma_subsampling_x);
  ASSERT_EQ(1, av1_config.chroma_subsampling_y);
  ASSERT_EQ(0, av1_config.chroma_sample_position);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_present);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_minus_one);

  // Test parsing of a reduced still image Sequence Header OBU.
  memset(&av1_config, 0, sizeof(av1_config));
  ASSERT_EQ(
      0, get_av1c_from_obu(&kAnnexBReducedStillImageSequenceHeaderObu[0],
                           sizeof(kAnnexBReducedStillImageSequenceHeaderObu), 1,
                           &av1_config));
  ASSERT_EQ(1, av1_config.marker);
  ASSERT_EQ(1, av1_config.version);
  ASSERT_EQ(0, av1_config.seq_profile);
  ASSERT_EQ(0, av1_config.seq_level_idx_0);
  ASSERT_EQ(0, av1_config.seq_tier_0);
  ASSERT_EQ(0, av1_config.high_bitdepth);
  ASSERT_EQ(0, av1_config.twelve_bit);
  ASSERT_EQ(0, av1_config.monochrome);
  ASSERT_EQ(1, av1_config.chroma_subsampling_x);
  ASSERT_EQ(1, av1_config.chroma_subsampling_y);
  ASSERT_EQ(0, av1_config.chroma_sample_position);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_present);
  ASSERT_EQ(0, av1_config.initial_presentation_delay_minus_one);
}

TEST(Av1c, ReadWriteConfig) {
  av1c av1_config;
  memset(&av1_config, 0, sizeof(av1_config));

  // Test writing out the AV1 config.
  size_t bytes_written = 0;
  uint8_t av1c_buffer[4] = { 0 };
  ASSERT_EQ(0, write_av1c(&av1_config, sizeof(av1c_buffer), &bytes_written,
                          &av1c_buffer[0]));
  ASSERT_EQ(0, memcmp(&av1c_buffer[0], &kAv1cAllZero[0], sizeof(kAv1cAllZero)));
  ASSERT_EQ(kAv1cNoConfigObusSize, bytes_written);

  // Test reading the AV1 config.
  size_t bytes_read = 0;
  ASSERT_EQ(0, read_av1c(&kAv1cAllZero[0], sizeof(kAv1cAllZero), &bytes_read,
                         &av1_config));
  ASSERT_EQ(kAv1cNoConfigObusSize, bytes_read);
  ASSERT_EQ(0, write_av1c(&av1_config, sizeof(av1c_buffer), &bytes_written,
                          &av1c_buffer[0]));
  ASSERT_EQ(0, memcmp(&av1c_buffer[0], &kAv1cAllZero[0], sizeof(kAv1cAllZero)));
}

}  // namespace
