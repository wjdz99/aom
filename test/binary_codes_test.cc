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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "config/aom_config.h"

#include "test/acm_random.h"
#include "aom/aom_integer.h"
#include "aom_dsp/bitreader.h"
#include "aom_dsp/bitwriter.h"
#include "aom_dsp/binary_codes_reader.h"
#include "aom_dsp/binary_codes_writer.h"

#define ACCT_STR __func__

using libaom_test::ACMRandom;

namespace {

// Test for Finite subexponential code with reference
TEST(AV1, TestPrimitiveRefsubexpfin) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  const unsigned int kBufferSize = 65536;
  aom_writer bw;
  uint8_t bw_buffer[kBufferSize];
  const uint16_t kRanges = 8;
  const uint16_t kSubexpParams = 6;
  const uint16_t kReferences = 8;
  const uint16_t kValues = 16;
  uint16_t enc_values[kRanges][kSubexpParams][kReferences][kValues][4];
  const uint16_t range_vals[kRanges] = { 1, 13, 64, 120, 230, 420, 1100, 8000 };
  aom_start_encode(&bw, bw_buffer);
  for (int n = 0; n < kRanges; ++n) {
    const uint16_t range = range_vals[n];
    for (int k = 0; k < kSubexpParams; ++k) {
      for (int r = 0; r < kReferences; ++r) {
        const uint16_t ref = rnd(range);
        for (int v = 0; v < kValues; ++v) {
          const uint16_t value = rnd(range);
          enc_values[n][k][r][v][0] = range;
          enc_values[n][k][r][v][1] = k;
          enc_values[n][k][r][v][2] = ref;
          enc_values[n][k][r][v][3] = value;
          aom_write_primitive_refsubexpfin(&bw, range, k, ref, value);
        }
      }
    }
  }
  aom_stop_encode(&bw);
  GTEST_ASSERT_LE(bw.pos, kBufferSize);
  aom_reader br;
  aom_reader_init(&br, bw_buffer, bw.pos);
  GTEST_ASSERT_GE(aom_reader_tell(&br), 0u);
  GTEST_ASSERT_LE(aom_reader_tell(&br), 1u);
  for (int n = 0; n < kRanges; ++n) {
    for (int k = 0; k < kSubexpParams; ++k) {
      for (int r = 0; r < kReferences; ++r) {
        for (int v = 0; v < kValues; ++v) {
          const uint16_t range = enc_values[n][k][r][v][0];
          assert(k == enc_values[n][k][r][v][1]);
          const uint16_t ref = enc_values[n][k][r][v][2];
          const uint16_t value =
              aom_read_primitive_refsubexpfin(&br, range, k, ref, ACCT_STR);
          GTEST_ASSERT_EQ(value, enc_values[n][k][r][v][3]);
        }
      }
    }
  }
}

// Test for aom_write_symbol with nsymbs=2.
TEST(AV1, TestSymbolBoolean) {
  aom_cdf_prob cdf[4][3] = {
    { 16384, 0, 0 },
    { 32768 - 8386, 0, 0 },
    { 32768 - 24312, 0, 0 },
    { 16384, 0, 0 },
  };
  constexpr int kSymbols[4][4] = { { 0, 0, 1, 1 },  //
                                   { 0, 1, 1, 0 },  //
                                   { 1, 0, 1, 0 },  //
                                   { 1, 0, 0, 1 } };
  const unsigned int kBufferSize = 65536;
  uint8_t bw_buffer[kBufferSize];
  aom_writer bw;
  bw.allow_update_cdf = 1;
  aom_start_encode(&bw, bw_buffer);
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        aom_write_symbol(&bw, kSymbols[j][k], cdf[k], 2);
      }
    }
  }
  aom_stop_encode(&bw);
  GTEST_ASSERT_LE(bw.pos, kBufferSize);
  printf("  constexpr size_t kNumBytes = %u;\n", bw.pos);
  printf("  constexpr uint8_t kBytes[] = {");
  int count = 0;
  for (unsigned int i = 0; i < bw.pos; ++i) {
    if (count++ % 12 == 0) {
      printf("\n      ");
    } else {
      printf(" ");
    }
    printf("0x%02x,", bw_buffer[i]);
  }
  printf("\n  };\n");
  aom_cdf_prob cdf0[4][3] = {
    { 16384, 0, 0 },
    { 32768 - 8386, 0, 0 },
    { 32768 - 24312, 0, 0 },
    { 16384, 0, 0 },
  };
  aom_reader br;
  br.allow_update_cdf = 1;
  aom_reader_init(&br, bw_buffer, bw.pos);
  int symbol;
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        symbol = aom_read_symbol(&br, cdf0[k], 2, ACCT_STR);
        EXPECT_EQ(symbol, kSymbols[j][k]);
      }
    }
  }
}

// Test for aom_write_symbol with nsymbs=3.
TEST(AV1, TestSymbol3) {
  aom_cdf_prob cdf[4][4] = {
    // pdf: 1/3, 1/3, 1/3
    { 32768 - 10923, 32768 - 21845, 0, 0 },
    // pdf: 1/6, 2/6, 3/6
    { 32768 - 5461, 32768 - 16384, 0, 0 },
    // pdf: 2/6, 3/6, 1/6
    { 32768 - 10923, 32768 - 27307, 0, 0 },
    // pdf: 3/6, 1/6, 2/6
    { 32768 - 16384, 32768 - 21845, 0, 0 },
  };
  constexpr int kSymbols[6][4] = { { 0, 2, 1, 2 },  //
                                   { 1, 1, 2, 1 },  //
                                   { 2, 0, 0, 0 },  //
                                   { 0, 2, 0, 2 },  //
                                   { 1, 2, 1, 0 },  //
                                   { 2, 1, 1, 0 } };
  const unsigned int kBufferSize = 65536;
  uint8_t bw_buffer[kBufferSize];
  aom_writer bw;
  bw.allow_update_cdf = 1;
  aom_start_encode(&bw, bw_buffer);
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 4; ++k) {
        aom_write_symbol(&bw, kSymbols[j][k], cdf[k], 3);
      }
    }
  }
  aom_stop_encode(&bw);
  GTEST_ASSERT_LE(bw.pos, kBufferSize);
  printf("  constexpr size_t kNumBytes = %u;\n", bw.pos);
  printf("  constexpr uint8_t kBytes[] = {");
  int count = 0;
  for (unsigned int i = 0; i < bw.pos; ++i) {
    if (count++ % 12 == 0) {
      printf("\n      ");
    } else {
      printf(" ");
    }
    printf("0x%02x,", bw_buffer[i]);
  }
  printf("\n  };\n");
  aom_cdf_prob cdf0[4][4] = {
    // pdf: 1/3, 1/3, 1/3
    { 32768 - 10923, 32768 - 21845, 0, 0 },
    // pdf: 1/6, 2/6, 3/6
    { 32768 - 5461, 32768 - 16384, 0, 0 },
    // pdf: 2/6, 3/6, 1/6
    { 32768 - 10923, 32768 - 27307, 0, 0 },
    // pdf: 3/6, 1/6, 2/6
    { 32768 - 16384, 32768 - 21845, 0, 0 },
  };
  aom_reader br;
  br.allow_update_cdf = 1;
  aom_reader_init(&br, bw_buffer, bw.pos);
  int symbol;
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 4; ++k) {
        symbol = aom_read_symbol(&br, cdf0[k], 3, ACCT_STR);
        EXPECT_EQ(symbol, kSymbols[j][k]);
      }
    }
  }
}

// Test for aom_write_symbol with nsymbs=4.
TEST(AV1, TestSymbol4) {
  aom_cdf_prob cdf[4][5] = {
    // pdf: 1/4, 1/4, 1/4, 1/4
    { 32768 - 8192, 32768 - 16384, 32768 - 24576, 0, 0 },
    // pdf: 2/8, 1/8, 2/8, 3/8
    { 32768 - 8192, 32768 - 12288, 32768 - 20480, 0, 0 },
    // pdf: 1/4, 1/4, 1/4, 1/4
    { 32768 - 8192, 32768 - 16384, 32768 - 24576, 0, 0 },
    // pdf: 2/8, 3/8, 2/8, 1/8
    { 32768 - 8192, 32768 - 20480, 32768 - 28672, 0, 0 },
  };
  constexpr int kSymbols[8][4] = { { 0, 0, 3, 3 },  //
                                   { 0, 0, 2, 2 },  //
                                   { 1, 1, 0, 0 },  //
                                   { 1, 2, 1, 1 },  //
                                   { 2, 2, 3, 2 },  //
                                   { 2, 3, 2, 1 },  //
                                   { 3, 3, 0, 0 },  //
                                   { 3, 3, 1, 1 } };
  const unsigned int kBufferSize = 65536;
  uint8_t bw_buffer[kBufferSize];
  aom_writer bw;
  bw.allow_update_cdf = 1;
  aom_start_encode(&bw, bw_buffer);
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int k = 0; k < 4; ++k) {
        aom_write_symbol(&bw, kSymbols[j][k], cdf[k], 4);
      }
    }
  }
  aom_stop_encode(&bw);
  GTEST_ASSERT_LE(bw.pos, kBufferSize);
  printf("  constexpr size_t kNumBytes = %u;\n", bw.pos);
  printf("  constexpr uint8_t kBytes[] = {");
  int count = 0;
  for (unsigned int i = 0; i < bw.pos; ++i) {
    if (count++ % 12 == 0) {
      printf("\n      ");
    } else {
      printf(" ");
    }
    printf("0x%02x,", bw_buffer[i]);
  }
  printf("\n  };\n");
  aom_cdf_prob cdf0[4][5] = {
    // pdf: 1/4, 1/4, 1/4, 1/4
    { 32768 - 8192, 32768 - 16384, 32768 - 24576, 0, 0 },
    // pdf: 2/8, 1/8, 2/8, 3/8
    { 32768 - 8192, 32768 - 12288, 32768 - 20480, 0, 0 },
    // pdf: 1/4, 1/4, 1/4, 1/4
    { 32768 - 8192, 32768 - 16384, 32768 - 24576, 0, 0 },
    // pdf: 2/8, 3/8, 2/8, 1/8
    { 32768 - 8192, 32768 - 20480, 32768 - 28672, 0, 0 },
  };
  aom_reader br;
  br.allow_update_cdf = 1;
  aom_reader_init(&br, bw_buffer, bw.pos);
  int symbol;
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 8; ++j) {
      for (int k = 0; k < 4; ++k) {
        symbol = aom_read_symbol(&br, cdf0[k], 4, ACCT_STR);
        EXPECT_EQ(symbol, kSymbols[j][k]);
      }
    }
  }
}

// Test for aom_write_symbol with nsymbs=7.
TEST(AV1, TestSymbol7) {
  aom_cdf_prob cdf[4][8] = {
    // pdf: 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7
    { 32768 - 4681, 32768 - 9362, 32768 - 14043, 32768 - 18725, 32768 - 23406,
      32768 - 28087, 0, 0 },
    // pdf: 3/14, 2/14, 2/14, 2/14, 2/14, 2/14, 1/14
    { 32768 - 7022, 32768 - 11703, 32768 - 16384, 32768 - 21065, 32768 - 25746,
      32768 - 30427, 0, 0 },
    // pdf: 1/14, 1/14, 2/14, 2/14, 2/14, 3/14, 3/14
    { 32768 - 2341, 32768 - 4681, 32768 - 9362, 32768 - 14043, 32768 - 18725,
      32768 - 25746, 0, 0 },
    // pdf: 1/14, 2/14, 3/14, 3/14, 2/14, 2/14, 1/14
    { 32768 - 2341, 32768 - 7022, 32768 - 14043, 32768 - 21065, 32768 - 25746,
      32768 - 30427, 0, 0 },
  };
  constexpr int kSymbols[14][4] = { { 0, 4, 6, 3 },  //
                                    { 1, 5, 5, 2 },  //
                                    { 2, 6, 4, 1 },  //
                                    { 3, 0, 3, 0 },  //
                                    { 4, 1, 2, 6 },  //
                                    { 5, 2, 1, 5 },  //
                                    { 6, 3, 0, 4 },  //
                                    { 0, 0, 6, 5 },  //
                                    { 2, 1, 4, 3 },  //
                                    { 4, 3, 6, 1 },  //
                                    { 6, 5, 2, 4 },  //
                                    { 1, 0, 5, 2 },  //
                                    { 3, 2, 3, 2 },  //
                                    { 5, 4, 5, 3 } };
  const unsigned int kBufferSize = 65536;
  uint8_t bw_buffer[kBufferSize];
  aom_writer bw;
  bw.allow_update_cdf = 1;
  aom_start_encode(&bw, bw_buffer);
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 14; ++j) {
      for (int k = 0; k < 4; ++k) {
        aom_write_symbol(&bw, kSymbols[j][k], cdf[k], 7);
      }
    }
  }
  aom_stop_encode(&bw);
  GTEST_ASSERT_LE(bw.pos, kBufferSize);
  printf("  constexpr size_t kNumBytes = %u;\n", bw.pos);
  printf("  constexpr uint8_t kBytes[] = {");
  int count = 0;
  for (unsigned int i = 0; i < bw.pos; ++i) {
    if (count++ % 12 == 0) {
      printf("\n      ");
    } else {
      printf(" ");
    }
    printf("0x%02x,", bw_buffer[i]);
  }
  printf("\n  };\n");
  aom_cdf_prob cdf0[4][8] = {
    // pdf: 1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7
    { 32768 - 4681, 32768 - 9362, 32768 - 14043, 32768 - 18725, 32768 - 23406,
      32768 - 28087, 0, 0 },
    // pdf: 3/14, 2/14, 2/14, 2/14, 2/14, 2/14, 1/14
    { 32768 - 7022, 32768 - 11703, 32768 - 16384, 32768 - 21065, 32768 - 25746,
      32768 - 30427, 0, 0 },
    // pdf: 1/14, 1/14, 2/14, 2/14, 2/14, 3/14, 3/14
    { 32768 - 2341, 32768 - 4681, 32768 - 9362, 32768 - 14043, 32768 - 18725,
      32768 - 25746, 0, 0 },
    // pdf: 1/14, 2/14, 3/14, 3/14, 2/14, 2/14, 1/14
    { 32768 - 2341, 32768 - 7022, 32768 - 14043, 32768 - 21065, 32768 - 25746,
      32768 - 30427, 0, 0 },
  };
  aom_reader br;
  br.allow_update_cdf = 1;
  aom_reader_init(&br, bw_buffer, bw.pos);
  int symbol;
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 14; ++j) {
      for (int k = 0; k < 4; ++k) {
        symbol = aom_read_symbol(&br, cdf0[k], 7, ACCT_STR);
        EXPECT_EQ(symbol, kSymbols[j][k]);
      }
    }
  }
}

// Test for aom_write_symbol with nsymbs=8.
TEST(AV1, TestSymbol8) {
  aom_cdf_prob cdf[4][9] = {
    // pdf: 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8
    { 32768 - 4096, 32768 - 8192, 32768 - 12288, 32768 - 16384, 32768 - 20480,
      32768 - 24576, 32768 - 28672, 0, 0 },
    // pdf: 3/16, 2/16, 2/16, 2/16, 2/16, 2/16, 2/16, 1/16
    { 32768 - 6144, 32768 - 10240, 32768 - 14336, 32768 - 18432, 32768 - 22528,
      32768 - 26624, 32768 - 30720, 0, 0 },
    // pdf: 1/16, 1/16, 2/16, 2/16, 2/16, 2/16, 3/16, 3/16
    { 32768 - 2048, 32768 - 4096, 32768 - 8192, 32768 - 12288, 32768 - 16384,
      32768 - 20480, 32768 - 26624, 0, 0 },
    // pdf: 1/16, 1/16, 3/16, 3/16, 3/16, 3/16, 1/16, 1/16
    { 32768 - 2048, 32768 - 4096, 32768 - 10240, 32768 - 16384, 32768 - 22528,
      32768 - 28672, 32768 - 30720, 0, 0 },
  };
  constexpr int kSymbols[16][4] = { { 0, 4, 7, 3 },  //
                                    { 1, 5, 6, 2 },  //
                                    { 2, 6, 5, 1 },  //
                                    { 3, 7, 4, 0 },  //
                                    { 4, 0, 3, 7 },  //
                                    { 5, 1, 2, 6 },  //
                                    { 6, 2, 1, 5 },  //
                                    { 7, 3, 0, 4 },  //
                                    { 0, 0, 6, 5 },  //
                                    { 2, 1, 4, 3 },  //
                                    { 4, 3, 6, 4 },  //
                                    { 6, 5, 2, 2 },  //
                                    { 1, 0, 7, 3 },  //
                                    { 3, 2, 5, 5 },  //
                                    { 5, 4, 7, 2 },  //
                                    { 7, 6, 3, 4 } };
  const unsigned int kBufferSize = 65536;
  uint8_t bw_buffer[kBufferSize];
  aom_writer bw;
  bw.allow_update_cdf = 1;
  aom_start_encode(&bw, bw_buffer);
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 16; ++j) {
      for (int k = 0; k < 4; ++k) {
        aom_write_symbol(&bw, kSymbols[j][k], cdf[k], 8);
      }
    }
  }
  aom_stop_encode(&bw);
  GTEST_ASSERT_LE(bw.pos, kBufferSize);
  printf("  constexpr size_t kNumBytes = %u;\n", bw.pos);
  printf("  constexpr uint8_t kBytes[] = {");
  int count = 0;
  for (unsigned int i = 0; i < bw.pos; ++i) {
    if (count++ % 12 == 0) {
      printf("\n      ");
    } else {
      printf(" ");
    }
    printf("0x%02x,", bw_buffer[i]);
  }
  printf("\n  };\n");
  aom_cdf_prob cdf0[4][9] = {
    // pdf: 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8
    { 32768 - 4096, 32768 - 8192, 32768 - 12288, 32768 - 16384, 32768 - 20480,
      32768 - 24576, 32768 - 28672, 0, 0 },
    // pdf: 3/16, 2/16, 2/16, 2/16, 2/16, 2/16, 2/16, 1/16
    { 32768 - 6144, 32768 - 10240, 32768 - 14336, 32768 - 18432, 32768 - 22528,
      32768 - 26624, 32768 - 30720, 0, 0 },
    // pdf: 1/16, 1/16, 2/16, 2/16, 2/16, 2/16, 3/16, 3/16
    { 32768 - 2048, 32768 - 4096, 32768 - 8192, 32768 - 12288, 32768 - 16384,
      32768 - 20480, 32768 - 26624, 0, 0 },
    // pdf: 1/16, 1/16, 3/16, 3/16, 3/16, 3/16, 1/16, 1/16
    { 32768 - 2048, 32768 - 4096, 32768 - 10240, 32768 - 16384, 32768 - 22528,
      32768 - 28672, 32768 - 30720, 0, 0 },
  };
  aom_reader br;
  br.allow_update_cdf = 1;
  aom_reader_init(&br, bw_buffer, bw.pos);
  int symbol;
  for (int i = 0; i < 1024; ++i) {
    for (int j = 0; j < 16; ++j) {
      for (int k = 0; k < 4; ++k) {
        symbol = aom_read_symbol(&br, cdf0[k], 8, ACCT_STR);
        EXPECT_EQ(symbol, kSymbols[j][k]);
      }
    }
  }
}

// TODO(debargha): Adds tests for other primitives
}  // namespace
