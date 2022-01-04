/*
 * Copyright (c) 2022, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#ifndef AOM_AOM_PORTS_MEM_BYTES_H_
#define AOM_AOM_PORTS_MEM_BYTES_H_

#include "aom/aom_integer.h"

/*!\brief force enum to be unsigned 1 byte*/
#define UENUM1BYTE(enumvar) \
  ;                         \
  typedef uint8_t enumvar

/*!\brief force enum to be signed 1 byte*/
#define SENUM1BYTE(enumvar) \
  ;                         \
  typedef int8_t enumvar

/*!\brief force enum to be unsigned 2 byte*/
#define UENUM2BYTE(enumvar) \
  ;                         \
  typedef uint16_t enumvar

/*!\brief force enum to be signed 2 byte*/
#define SENUM2BYTE(enumvar) \
  ;                         \
  typedef int16_t enumvar

/*!\brief force enum to be unsigned 4 byte*/
#define UENUM4BYTE(enumvar) \
  ;                         \
  typedef uint32_t enumvar

/*!\brief force enum to be unsigned 4 byte*/
#define SENUM4BYTE(enumvar) \
  ;                         \
  typedef int32_t enumvar

#endif  // AOM_AOM_PORTS_MEM_BYTES_H_
