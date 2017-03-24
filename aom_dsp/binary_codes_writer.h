/*
 * Copyright (c) 2017, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_DSP_BINARY_CODES_WRITER_H_
#define AOM_DSP_BINARY_CODES_WRITER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include "./aom_config.h"
#include "aom/aom_integer.h"
#include "aom_dsp/bitwriter.h"

// mag_bits is number of bits for magnitude. The alphabet is of size
// 2 * 2^mag_bits + 1, symmetric around 0, where one bit is used to
// indicate 0 or non-zero, mag_bits bits are used to indicate magnitide
// and 1 more bit for the sign if non-zero.
void aom_write_primitive_symmetric(aom_writer *w, int word,
                                   unsigned int mag_bits);

// Encodes a value v in [0, n-1] based on a reference ref also in [0, n-1]
// The closest p values of v from ref are coded using a p-ary quasi-unoform
// short code while the remaining n-p values are coded with a longer code.
void aom_write_primitive_refbilevel(aom_writer *w, uint16_t n, uint16_t p,
                                    uint16_t ref, uint16_t v);
#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // AOM_DSP_BINARY_CODES_WRITER_H_
