/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 *
 * Based on code from the OggTheora software codec source code,
 * Copyright (C) 2002-2010 The Xiph.Org Foundation and contributors.
 */
#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "aom/aom_integer.h"
#include "aom_ports/msvc.h"
#include "y4minput.h"

// Reads 'size' bytes from 'file' into 'buf' with some fault tolerance.
// Returns true on success.
static int file_read(void *buf, size_t size, FILE *file) {
  const int kMaxRetries = 5;
  int retry_count = 0;
  int file_error;
  size_t len = 0;
  do {
    const size_t n = fread((uint8_t *)buf + len, 1, size - len, file);
    len += n;
    file_error = ferror(file);
    if (file_error) {
      if (errno == EINTR || errno == EAGAIN) {
        clearerr(file);
        continue;
      } else {
        fprintf(stderr, "Error reading file: %u of %u bytes read, %d: %s\n",
                (uint32_t)len, (uint32_t)size, errno, strerror(errno));
        return 0;
      }
    }
  } while (!feof(file) && len < size && ++retry_count < kMaxRetries);

  if (!feof(file) && len != size) {
    fprintf(stderr,
            "Error reading file: %u of %u bytes read,"
            " error: %d, retries: %d, %d: %s\n",
            (uint32_t)len, (uint32_t)size, file_error, retry_count, errno,
            strerror(errno));
  }
  return len == size;
}

// Parses an int from the string, storing it in the destination.
static bool parse_int(int *dst, const char *buf) {
  if (sscanf(buf, "%d", dst) != 1) {
    fprintf(stderr, "Error parsing integer: %s\n", buf);
    return false;
  }
  return true;
}

// Parses a rational number (represented as X:Y) from the string,
// storing the numerator and denominator.
static bool parse_rational(int *n, int *d, const char *buf) {
  if (sscanf(buf, "%d:%d", n, d) != 2) {
    fprintf(stderr, "Error parsing rational: %s\n", buf);
    return false;
  }
  return true;
}

// Checks that the buffer is 1-character long, and assigns it into the
// destination buffer.
static bool parse_char(char *dst, const char *buf) {
  if (strlen(buf) != 1) {
    fprintf(stderr, "Too many characters (expected 1): %s\n", buf);
    return false;
  }
  *dst = buf[0];
  return true;
}

// Checks that the string + null-char do not exceed the length, and stores them
// in the destination buffer.
static bool parse_string(char *dst, size_t len, const char *buf) {
  if (len < strlen(buf) + 1) {
    fprintf(stderr, "Buffer not large enough (%lu bytes) for value: %s\n",
            (unsigned long)len, buf);
    return false;
  }
  strncpy(dst, buf, strlen(buf) + 1);
  return true;
}

static bool parse_color_range(y4m_input *_y4m, const char *buf) {
  // Note that default is studio range.
  if (strcmp(buf, "LIMITED") == 0) {
    return true;
  }
  if (strcmp(buf, "FULL") == 0) {
    _y4m->color_range = AOM_CR_FULL_RANGE;
    return true;
  }
  fprintf(stderr, "Unknown color range value: %s\n", buf);
  return false;
}

static bool parse_metadata(y4m_input *_y4m, const char *buf) {
  if (strncmp(buf, "COLORRANGE=", 11) == 0) {
    return parse_color_range(_y4m, buf + 11);
  }
  return true;  // No support for other metadata, just ignore them.
}

// Returns true if there is no tag, or if the tag is parsed successfully;
// false otherwise. Note that the tag ends at the first space or newline
// character.
static bool parse_single_tag(const char *buf, y4m_input *_y4m) {
  // Empty tag.
  if (buf[0] == '\0') {
    return true;
  }
  // Since leading space characters are consumed, the first character
  // should never be a space.
  assert(buf[0] != ' ' && buf[0] != '\n');
  switch (buf[0]) {
    case 'W': return parse_int(&_y4m->pic_w, buf + 1);
    case 'H': return parse_int(&_y4m->pic_h, buf + 1);
    case 'F': return parse_rational(&_y4m->fps_n, &_y4m->fps_d, buf + 1);
    case 'I': return parse_char(&_y4m->interlace, buf + 1);
    case 'A': return parse_rational(&_y4m->par_n, &_y4m->par_d, buf + 1);
    case 'C':
      return parse_string(_y4m->chroma_type, sizeof(_y4m->chroma_type),
                          buf + 1);
    case 'X': return parse_metadata(_y4m, buf + 1);
    default: return true;  // Skip the tag.
  }
}

// Copy a single tag into the buffer, along with a null character.
// Returns false if any file IO errors occur.
static bool copy_tag(char *buf, size_t buf_len, char *end_tag, FILE *_fin) {
  assert(buf_len >= 1);
  // Skip leading space characters.
  do {
    if (!file_read(buf, 1, _fin)) {
      return false;
    }
  } while (buf[0] == ' ');

  // If we hit the newline, treat this as the "empty" tag.
  if (buf[0] == '\n') {
    buf[0] = '\0';
    *end_tag = '\n';
    return true;
  }

  // Copy over characters until a space is hit, or the buffer is exhausted.
  size_t i;
  for (i = 1; i < buf_len; ++i) {
    if (!file_read(buf + i, 1, _fin)) {
      return false;
    }
    if (buf[i] == ' ' || buf[i] == '\n') {
      break;
    }
  }
  if (i == buf_len) {
    fprintf(stderr, "Error: Y4M header tags must be less than %lu characters\n",
            (unsigned long)i);
    return false;
  }
  *end_tag = buf[i];
  buf[i] = '\0';
  return true;
}

// Returns true if tags were parsed successfully and a newline was encountered,
// false otherwise.
static bool parse_tags(y4m_input *_y4m, FILE *_fin) {
  // Set Y4M tags to defaults, updating them as processing occurs. Mandatory
  // fields are marked with -1 and will be checked after the tags are parsed.
  _y4m->pic_w = -1;
  _y4m->pic_h = -1;
  _y4m->fps_n = -1;  // Also serves as marker for fps_d
  _y4m->par_n = 0;
  _y4m->par_d = 0;
  _y4m->interlace = '?';
  snprintf(_y4m->chroma_type, sizeof(_y4m->chroma_type), "420");
  _y4m->color_range = AOM_CR_STUDIO_RANGE;

  // Find one tag at a time. Buffer contains the tag.
  char buf[256];
  char end_tag;  // Character denoting the end of the tag, ' ' or '\n'.
  do {
    if (!copy_tag(buf, sizeof(buf), &end_tag, _fin)) {
      return false;
    }
    if (!parse_single_tag(buf, _y4m)) {
      return false;
    }
  } while (end_tag != '\n');

  // Check the mandatory fields.
  if (_y4m->pic_w == -1) {
    fprintf(stderr, "Width field missing\n");
    return false;
  }
  if (_y4m->pic_h == -1) {
    fprintf(stderr, "Height field missing\n");
    return false;
  }
  if (_y4m->fps_n == -1) {
    fprintf(stderr, "FPS field missing\n");
    return false;
  }
  return true;
}

/*All anti-aliasing filters in the following conversion functions are based on
   one of two window functions:
  The 6-tap Lanczos window (for down-sampling and shifts):
   sinc(\pi*t)*sinc(\pi*t/3), |t|<3  (sinc(t)==sin(t)/t)
   0,                         |t|>=3
  The 4-tap Mitchell window (for up-sampling):
   7|t|^3-12|t|^2+16/3,             |t|<1
   -(7/3)|x|^3+12|x|^2-20|x|+32/3,  |t|<2
   0,                               |t|>=2
  The number of taps is intentionally kept small to reduce computational
   overhead and limit ringing.

  The taps from these filters are scaled so that their sum is 1, and the
  result is scaled by 128 and rounded to integers to create a filter whose
   intermediate values fit inside 16 bits.
  Coefficients are rounded in such a way as to ensure their sum is still 128,
   which is usually equivalent to normal rounding.

  Conversions which require both horizontal and vertical filtering could
   have these steps pipelined, for less memory consumption and better cache
   performance, but we do them separately for simplicity.*/
#define OC_MINI(_a, _b) ((_a) > (_b) ? (_b) : (_a))
#define OC_MAXI(_a, _b) ((_a) < (_b) ? (_b) : (_a))
#define OC_CLAMPI(_a, _b, _c) (OC_MAXI(_a, OC_MINI(_b, _c)))

/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  420mpeg2 chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  BR      |       BR      |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  BR      |       BR      |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  We use a resampling filter to shift the site locations one quarter pixel (at
   the chroma plane's resolution) to the right.
  The 4:2:2 modes look exactly the same, except there are twice as many chroma
   lines, and they are vertically co-sited with the luma samples in both the
   mpeg2 and jpeg cases (thus requiring no vertical resampling).*/
static void y4m_42xmpeg2_42xjpeg_helper(unsigned char *_dst,
                                        const unsigned char *_src, int _c_w,
                                        int _c_h) {
  int y;
  int x;
  for (y = 0; y < _c_h; y++) {
    /*Filter: [4 -17 114 35 -9 1]/128, derived from a 6-tap Lanczos
       window.*/
    for (x = 0; x < OC_MINI(_c_w, 2); x++) {
      _dst[x] = (unsigned char)OC_CLAMPI(
          0,
          (4 * _src[0] - 17 * _src[OC_MAXI(x - 1, 0)] + 114 * _src[x] +
           35 * _src[OC_MINI(x + 1, _c_w - 1)] -
           9 * _src[OC_MINI(x + 2, _c_w - 1)] + _src[OC_MINI(x + 3, _c_w - 1)] +
           64) >>
              7,
          255);
    }
    for (; x < _c_w - 3; x++) {
      _dst[x] = (unsigned char)OC_CLAMPI(
          0,
          (4 * _src[x - 2] - 17 * _src[x - 1] + 114 * _src[x] +
           35 * _src[x + 1] - 9 * _src[x + 2] + _src[x + 3] + 64) >>
              7,
          255);
    }
    for (; x < _c_w; x++) {
      _dst[x] = (unsigned char)OC_CLAMPI(
          0,
          (4 * _src[x - 2] - 17 * _src[x - 1] + 114 * _src[x] +
           35 * _src[OC_MINI(x + 1, _c_w - 1)] -
           9 * _src[OC_MINI(x + 2, _c_w - 1)] + _src[_c_w - 1] + 64) >>
              7,
          255);
    }
    _dst += _c_w;
    _src += _c_w;
  }
}

/*Handles both 422 and 420mpeg2 to 422jpeg and 420jpeg, respectively.*/
static void y4m_convert_42xmpeg2_42xjpeg(y4m_input *_y4m, unsigned char *_dst,
                                         unsigned char *_aux) {
  int c_w;
  int c_h;
  int c_sz;
  int pli;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  for (pli = 1; pli < 3; pli++) {
    y4m_42xmpeg2_42xjpeg_helper(_dst, _aux, c_w, c_h);
    _dst += c_sz;
    _aux += c_sz;
  }
}

/*This format is only used for interlaced content, but is included for
   completeness.

  420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  420paldv chroma samples are sited like:
  YR------Y-------YR------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YB------Y-------YB------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YR------Y-------YR------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YB------Y-------YB------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  We use a resampling filter to shift the site locations one quarter pixel (at
   the chroma plane's resolution) to the right.
  Then we use another filter to move the C_r location down one quarter pixel,
   and the C_b location up one quarter pixel.*/
static void y4m_convert_42xpaldv_42xjpeg(y4m_input *_y4m, unsigned char *_dst,
                                         unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int pli;
  int y;
  int x;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + 1) / 2;
  c_h = (_y4m->pic_h + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  c_sz = c_w * c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*First do the horizontal re-sampling.
      This is the same as the mpeg2 case, except that after the horizontal
       case, we need to apply a second vertical filter.*/
    y4m_42xmpeg2_42xjpeg_helper(tmp, _aux, c_w, c_h);
    _aux += c_sz;
    switch (pli) {
      case 1: {
        /*Slide C_b up a quarter-pel.
          This is the same filter used above, but in the other order.*/
        for (x = 0; x < c_w; x++) {
          for (y = 0; y < OC_MINI(c_h, 3); y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (tmp[0] - 9 * tmp[OC_MAXI(y - 2, 0) * c_w] +
                 35 * tmp[OC_MAXI(y - 1, 0) * c_w] + 114 * tmp[y * c_w] -
                 17 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] +
                 4 * tmp[OC_MINI(y + 2, c_h - 1) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h - 2; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (tmp[(y - 3) * c_w] - 9 * tmp[(y - 2) * c_w] +
                 35 * tmp[(y - 1) * c_w] + 114 * tmp[y * c_w] -
                 17 * tmp[(y + 1) * c_w] + 4 * tmp[(y + 2) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (tmp[(y - 3) * c_w] - 9 * tmp[(y - 2) * c_w] +
                 35 * tmp[(y - 1) * c_w] + 114 * tmp[y * c_w] -
                 17 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] +
                 4 * tmp[(c_h - 1) * c_w] + 64) >>
                    7,
                255);
          }
          _dst++;
          tmp++;
        }
        _dst += c_sz - c_w;
        tmp -= c_w;
      } break;
      case 2: {
        /*Slide C_r down a quarter-pel.
          This is the same as the horizontal filter.*/
        for (x = 0; x < c_w; x++) {
          for (y = 0; y < OC_MINI(c_h, 2); y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (4 * tmp[0] - 17 * tmp[OC_MAXI(y - 1, 0) * c_w] +
                 114 * tmp[y * c_w] + 35 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] -
                 9 * tmp[OC_MINI(y + 2, c_h - 1) * c_w] +
                 tmp[OC_MINI(y + 3, c_h - 1) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h - 3; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (4 * tmp[(y - 2) * c_w] - 17 * tmp[(y - 1) * c_w] +
                 114 * tmp[y * c_w] + 35 * tmp[(y + 1) * c_w] -
                 9 * tmp[(y + 2) * c_w] + tmp[(y + 3) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (4 * tmp[(y - 2) * c_w] - 17 * tmp[(y - 1) * c_w] +
                 114 * tmp[y * c_w] + 35 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] -
                 9 * tmp[OC_MINI(y + 2, c_h - 1) * c_w] + tmp[(c_h - 1) * c_w] +
                 64) >>
                    7,
                255);
          }
          _dst++;
          tmp++;
        }
      } break;
    }
    /*For actual interlaced material, this would have to be done separately on
       each field, and the shift amounts would be different.
      C_r moves down 1/8, C_b up 3/8 in the top field, and C_r moves down 3/8,
       C_b up 1/8 in the bottom field.
      The corresponding filters would be:
       Down 1/8 (reverse order for up): [3 -11 125 15 -4 0]/128
       Down 3/8 (reverse order for up): [4 -19 98 56 -13 2]/128*/
  }
}

/*Perform vertical filtering to reduce a single plane from 4:2:2 to 4:2:0.
  This is used as a helper by several conversion routines.*/
static void y4m_422jpeg_420jpeg_helper(unsigned char *_dst,
                                       const unsigned char *_src, int _c_w,
                                       int _c_h) {
  int y;
  int x;
  /*Filter: [3 -17 78 78 -17 3]/128, derived from a 6-tap Lanczos window.*/
  for (x = 0; x < _c_w; x++) {
    for (y = 0; y < OC_MINI(_c_h, 2); y += 2) {
      _dst[(y >> 1) * _c_w] =
          OC_CLAMPI(0,
                    (64 * _src[0] + 78 * _src[OC_MINI(1, _c_h - 1) * _c_w] -
                     17 * _src[OC_MINI(2, _c_h - 1) * _c_w] +
                     3 * _src[OC_MINI(3, _c_h - 1) * _c_w] + 64) >>
                        7,
                    255);
    }
    for (; y < _c_h - 3; y += 2) {
      _dst[(y >> 1) * _c_w] =
          OC_CLAMPI(0,
                    (3 * (_src[(y - 2) * _c_w] + _src[(y + 3) * _c_w]) -
                     17 * (_src[(y - 1) * _c_w] + _src[(y + 2) * _c_w]) +
                     78 * (_src[y * _c_w] + _src[(y + 1) * _c_w]) + 64) >>
                        7,
                    255);
    }
    for (; y < _c_h; y += 2) {
      _dst[(y >> 1) * _c_w] = OC_CLAMPI(
          0,
          (3 * (_src[(y - 2) * _c_w] + _src[(_c_h - 1) * _c_w]) -
           17 * (_src[(y - 1) * _c_w] + _src[OC_MINI(y + 2, _c_h - 1) * _c_w]) +
           78 * (_src[y * _c_w] + _src[OC_MINI(y + 1, _c_h - 1) * _c_w]) +
           64) >>
              7,
          255);
    }
    _src++;
    _dst++;
  }
}

/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  422jpeg chroma samples are sited like:
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  We use a resampling filter to decimate the chroma planes by two in the
   vertical direction.*/
static void y4m_convert_422jpeg_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                        unsigned char *_aux) {
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_w;
  int dst_c_h;
  int dst_c_sz;
  int pli;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = dst_c_w * dst_c_h;
  for (pli = 1; pli < 3; pli++) {
    y4m_422jpeg_420jpeg_helper(_dst, _aux, c_w, c_h);
    _aux += c_sz;
    _dst += dst_c_sz;
  }
}

/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  422 chroma samples are sited like:
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  We use a resampling filter to shift the original site locations one quarter
   pixel (at the original chroma resolution) to the right.
  Then we use a second resampling filter to decimate the chroma planes by two
   in the vertical direction.*/
static void y4m_convert_422_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                    unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_h;
  int dst_c_sz;
  int pli;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = c_w * dst_c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*In reality, the horizontal and vertical steps could be pipelined, for
       less memory consumption and better cache performance, but we do them
       separately for simplicity.*/
    /*First do horizontal filtering (convert to 422jpeg)*/
    y4m_42xmpeg2_42xjpeg_helper(tmp, _aux, c_w, c_h);
    /*Now do the vertical filtering.*/
    y4m_422jpeg_420jpeg_helper(_dst, tmp, c_w, c_h);
    _aux += c_sz;
    _dst += dst_c_sz;
  }
}

/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  411 chroma samples are sited like:
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |

  We use a filter to resample at site locations one eighth pixel (at the source
   chroma plane's horizontal resolution) and five eighths of a pixel to the
   right.
  Then we use another filter to decimate the planes by 2 in the vertical
   direction.*/
static void y4m_convert_411_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                    unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_w;
  int dst_c_h;
  int dst_c_sz;
  int tmp_sz;
  int pli;
  int y;
  int x;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = dst_c_w * dst_c_h;
  tmp_sz = dst_c_w * c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*In reality, the horizontal and vertical steps could be pipelined, for
       less memory consumption and better cache performance, but we do them
       separately for simplicity.*/
    /*First do horizontal filtering (convert to 422jpeg)*/
    for (y = 0; y < c_h; y++) {
      /*Filters: [1 110 18 -1]/128 and [-3 50 86 -5]/128, both derived from a
         4-tap Mitchell window.*/
      for (x = 0; x < OC_MINI(c_w, 1); x++) {
        tmp[x << 1] = (unsigned char)OC_CLAMPI(
            0,
            (111 * _aux[0] + 18 * _aux[OC_MINI(1, c_w - 1)] -
             _aux[OC_MINI(2, c_w - 1)] + 64) >>
                7,
            255);
        tmp[x << 1 | 1] = (unsigned char)OC_CLAMPI(
            0,
            (47 * _aux[0] + 86 * _aux[OC_MINI(1, c_w - 1)] -
             5 * _aux[OC_MINI(2, c_w - 1)] + 64) >>
                7,
            255);
      }
      for (; x < c_w - 2; x++) {
        tmp[x << 1] =
            (unsigned char)OC_CLAMPI(0,
                                     (_aux[x - 1] + 110 * _aux[x] +
                                      18 * _aux[x + 1] - _aux[x + 2] + 64) >>
                                         7,
                                     255);
        tmp[x << 1 | 1] = (unsigned char)OC_CLAMPI(
            0,
            (-3 * _aux[x - 1] + 50 * _aux[x] + 86 * _aux[x + 1] -
             5 * _aux[x + 2] + 64) >>
                7,
            255);
      }
      for (; x < c_w; x++) {
        tmp[x << 1] = (unsigned char)OC_CLAMPI(
            0,
            (_aux[x - 1] + 110 * _aux[x] + 18 * _aux[OC_MINI(x + 1, c_w - 1)] -
             _aux[c_w - 1] + 64) >>
                7,
            255);
        if ((x << 1 | 1) < dst_c_w) {
          tmp[x << 1 | 1] = (unsigned char)OC_CLAMPI(
              0,
              (-3 * _aux[x - 1] + 50 * _aux[x] +
               86 * _aux[OC_MINI(x + 1, c_w - 1)] - 5 * _aux[c_w - 1] + 64) >>
                  7,
              255);
        }
      }
      tmp += dst_c_w;
      _aux += c_w;
    }
    tmp -= tmp_sz;
    /*Now do the vertical filtering.*/
    y4m_422jpeg_420jpeg_helper(_dst, tmp, dst_c_w, c_h);
    _dst += dst_c_sz;
  }
}

/*Convert 444 to 420jpeg.*/
static void y4m_convert_444_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                    unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_w;
  int dst_c_h;
  int dst_c_sz;
  int tmp_sz;
  int pli;
  int y;
  int x;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = dst_c_w * dst_c_h;
  tmp_sz = dst_c_w * c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*Filter: [3 -17 78 78 -17 3]/128, derived from a 6-tap Lanczos window.*/
    for (y = 0; y < c_h; y++) {
      for (x = 0; x < OC_MINI(c_w, 2); x += 2) {
        tmp[x >> 1] = OC_CLAMPI(0,
                                (64 * _aux[0] + 78 * _aux[OC_MINI(1, c_w - 1)] -
                                 17 * _aux[OC_MINI(2, c_w - 1)] +
                                 3 * _aux[OC_MINI(3, c_w - 1)] + 64) >>
                                    7,
                                255);
      }
      for (; x < c_w - 3; x += 2) {
        tmp[x >> 1] = OC_CLAMPI(0,
                                (3 * (_aux[x - 2] + _aux[x + 3]) -
                                 17 * (_aux[x - 1] + _aux[x + 2]) +
                                 78 * (_aux[x] + _aux[x + 1]) + 64) >>
                                    7,
                                255);
      }
      for (; x < c_w; x += 2) {
        tmp[x >> 1] =
            OC_CLAMPI(0,
                      (3 * (_aux[x - 2] + _aux[c_w - 1]) -
                       17 * (_aux[x - 1] + _aux[OC_MINI(x + 2, c_w - 1)]) +
                       78 * (_aux[x] + _aux[OC_MINI(x + 1, c_w - 1)]) + 64) >>
                          7,
                      255);
      }
      tmp += dst_c_w;
      _aux += c_w;
    }
    tmp -= tmp_sz;
    /*Now do the vertical filtering.*/
    y4m_422jpeg_420jpeg_helper(_dst, tmp, dst_c_w, c_h);
    _dst += dst_c_sz;
  }
}

/*The image is padded with empty chroma components at 4:2:0.*/
static void y4m_convert_mono_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                     unsigned char *_aux) {
  int c_sz;
  (void)_aux;
  _dst += _y4m->pic_w * _y4m->pic_h;
  c_sz = ((_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h) *
         ((_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v);
  memset(_dst, 128, c_sz * 2);
}

/*No conversion function needed.*/
static void y4m_convert_null(y4m_input *_y4m, unsigned char *_dst,
                             unsigned char *_aux) {
  (void)_y4m;
  (void)_dst;
  (void)_aux;
}

static const char TAG[] = "YUV4MPEG2";

int y4m_input_open(y4m_input *_y4m, FILE *_fin, char *_skip, int _nskip,
                   aom_chroma_sample_position_t csp, int only_420) {
  // File must start with TAG.
  char buffer[9];  // 9 == strlen(TAG)
  // Read as much as possible from the skip-buffer, which were characters
  // that were previously read from the file to do input-type detection.
  assert(_nskip >= 0 && _nskip <= 8);
  if (_nskip > 0) {
    memcpy(buffer, _skip, _nskip);
  }
  // Start reading from the file now that the skip buffer is depleted.
  if (!file_read(buffer + _nskip, 9 - _nskip, _fin)) {
    return -1;
  }
  if (memcmp("YUV4MPEG2", buffer, 9) != 0) {
    fprintf(stderr, "Error parsing header: must start with %s\n", TAG);
    return -1;
  }

  // Next character must be a space.
  if (!file_read(buffer, 1, _fin) || buffer[0] != ' ') {
    fprintf(stderr, "Error parsing header: space must follow %s\n", TAG);
    return -1;
  }

  if (!parse_tags(_y4m, _fin)) {
    fprintf(stderr, "Error parsing %s header.\n", TAG);
    return -1;
  }

  if (_y4m->interlace == '?') {
    fprintf(stderr,
            "Warning: Input video interlacing format unknown; "
            "assuming progressive scan.\n");
  } else if (_y4m->interlace != 'p') {
    fprintf(stderr,
            "Input video is interlaced; "
            "Only progressive scan handled.\n");
    return -1;
  }
  /* Only support vertical chroma sample position if the input format is
   * already 420mpeg2. Colocated is not supported in Y4M.
   */
  if (csp == AOM_CSP_VERTICAL && strcmp(_y4m->chroma_type, "420mpeg2") != 0) {
    fprintf(stderr,
            "Vertical chroma sample position only supported "
            "for 420mpeg2 input\n");
    return -1;
  }
  if (csp == AOM_CSP_COLOCATED) {
    fprintf(stderr, "Colocated chroma sample position not supported in Y4M\n");
    return -1;
  }
  _y4m->aux_buf = NULL;
  _y4m->dst_buf = NULL;
  _y4m->aom_fmt = AOM_IMG_FMT_I420;
  _y4m->bps = 12;
  _y4m->bit_depth = 8;
  if (strcmp(_y4m->chroma_type, "420") == 0 ||
      strcmp(_y4m->chroma_type, "420jpeg") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = _y4m->src_c_dec_v =
        _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz =
        _y4m->pic_w * _y4m->pic_h +
        2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    /* Natively supported: no conversion required. */
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
  } else if (strcmp(_y4m->chroma_type, "420p10") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 2;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz =
        2 * (_y4m->pic_w * _y4m->pic_h +
             2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2));
    /* Natively supported: no conversion required. */
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    _y4m->bit_depth = 10;
    _y4m->bps = 15;
    _y4m->aom_fmt = AOM_IMG_FMT_I42016;
    if (only_420) {
      fprintf(stderr, "Unsupported conversion from 420p10 to 420jpeg\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "420p12") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 2;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz =
        2 * (_y4m->pic_w * _y4m->pic_h +
             2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2));
    /* Natively supported: no conversion required. */
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    _y4m->bit_depth = 12;
    _y4m->bps = 18;
    _y4m->aom_fmt = AOM_IMG_FMT_I42016;
    if (only_420) {
      fprintf(stderr, "Unsupported conversion from 420p12 to 420jpeg\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "420mpeg2") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = _y4m->src_c_dec_v =
        _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.*/
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz =
        2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    _y4m->convert = y4m_convert_null;
    if (csp != AOM_CSP_VERTICAL) {
      _y4m->convert = y4m_convert_42xmpeg2_42xjpeg;
      snprintf(_y4m->chroma_type, sizeof(_y4m->chroma_type), "420");
    }
  } else if (strcmp(_y4m->chroma_type, "420paldv") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = _y4m->src_c_dec_v =
        _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.
      We need to make two filter passes, so we need some extra space in the
       aux buffer.*/
    _y4m->aux_buf_sz = 3 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    _y4m->aux_buf_read_sz =
        2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    _y4m->convert = y4m_convert_42xpaldv_42xjpeg;
  } else if (strcmp(_y4m->chroma_type, "422jpeg") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.*/
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz =
        2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
    _y4m->convert = y4m_convert_422jpeg_420jpeg;
  } else if (strcmp(_y4m->chroma_type, "422") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    if (only_420) {
      _y4m->dst_c_dec_h = 2;
      _y4m->dst_c_dec_v = 2;
      _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
      /*Chroma filter required: read into the aux buf first.
        We need to make two filter passes, so we need some extra space in the
         aux buffer.*/
      _y4m->aux_buf_read_sz = 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      _y4m->aux_buf_sz =
          _y4m->aux_buf_read_sz + ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      _y4m->convert = y4m_convert_422_420jpeg;
    } else {
      _y4m->aom_fmt = AOM_IMG_FMT_I422;
      _y4m->bps = 16;
      _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
      _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
      _y4m->dst_buf_read_sz =
          _y4m->pic_w * _y4m->pic_h + 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      /*Natively supported: no conversion required.*/
      _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
      _y4m->convert = y4m_convert_null;
    }
  } else if (strcmp(_y4m->chroma_type, "422p10") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->aom_fmt = AOM_IMG_FMT_I42216;
    _y4m->bps = 20;
    _y4m->bit_depth = 10;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * (_y4m->pic_w * _y4m->pic_h +
                                 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h);
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      fprintf(stderr, "Unsupported conversion from 422p10 to 420jpeg\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "422p12") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->aom_fmt = AOM_IMG_FMT_I42216;
    _y4m->bps = 24;
    _y4m->bit_depth = 12;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * (_y4m->pic_w * _y4m->pic_h +
                                 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h);
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      fprintf(stderr, "Unsupported conversion from 422p12 to 420jpeg\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "411") == 0) {
    _y4m->src_c_dec_h = 4;
    _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.
      We need to make two filter passes, so we need some extra space in the
       aux buffer.*/
    _y4m->aux_buf_read_sz = 2 * ((_y4m->pic_w + 3) / 4) * _y4m->pic_h;
    _y4m->aux_buf_sz =
        _y4m->aux_buf_read_sz + ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
    _y4m->convert = y4m_convert_411_420jpeg;
  } else if (strcmp(_y4m->chroma_type, "444") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    if (only_420) {
      _y4m->dst_c_dec_h = 2;
      _y4m->dst_c_dec_v = 2;
      _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
      /*Chroma filter required: read into the aux buf first.
        We need to make two filter passes, so we need some extra space in the
         aux buffer.*/
      _y4m->aux_buf_read_sz = 2 * _y4m->pic_w * _y4m->pic_h;
      _y4m->aux_buf_sz =
          _y4m->aux_buf_read_sz + ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      _y4m->convert = y4m_convert_444_420jpeg;
    } else {
      _y4m->aom_fmt = AOM_IMG_FMT_I444;
      _y4m->bps = 24;
      _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
      _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
      _y4m->dst_buf_read_sz = 3 * _y4m->pic_w * _y4m->pic_h;
      /*Natively supported: no conversion required.*/
      _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
      _y4m->convert = y4m_convert_null;
    }
  } else if (strcmp(_y4m->chroma_type, "444p10") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    _y4m->aom_fmt = AOM_IMG_FMT_I44416;
    _y4m->bps = 30;
    _y4m->bit_depth = 10;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * 3 * _y4m->pic_w * _y4m->pic_h;
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      fprintf(stderr, "Unsupported conversion from 444p10 to 420jpeg\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "444p12") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    _y4m->aom_fmt = AOM_IMG_FMT_I44416;
    _y4m->bps = 36;
    _y4m->bit_depth = 12;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * 3 * _y4m->pic_w * _y4m->pic_h;
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      fprintf(stderr, "Unsupported conversion from 444p12 to 420jpeg\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "444alpha") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    if (only_420) {
      _y4m->dst_c_dec_h = 2;
      _y4m->dst_c_dec_v = 2;
      _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
      /*Chroma filter required: read into the aux buf first.
        We need to make two filter passes, so we need some extra space in the
         aux buffer.
        The extra plane also gets read into the aux buf.
        It will be discarded.*/
      _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 3 * _y4m->pic_w * _y4m->pic_h;
      _y4m->convert = y4m_convert_444_420jpeg;
    } else {
      fprintf(stderr, "Unsupported format: 444A\n");
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "mono") == 0) {
    _y4m->src_c_dec_h = _y4m->src_c_dec_v = 0;
    _y4m->dst_c_dec_h = _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*No extra space required, but we need to clear the chroma planes.*/
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_mono_420jpeg;
  } else {
    fprintf(stderr, "Unknown chroma sampling type: %s\n", _y4m->chroma_type);
    return -1;
  }
  /*The size of the final frame buffers is always computed from the
     destination chroma decimation type.*/
  _y4m->dst_buf_sz =
      _y4m->pic_w * _y4m->pic_h +
      2 * ((_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h) *
          ((_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v);
  if (_y4m->bit_depth == 8)
    _y4m->dst_buf = (unsigned char *)malloc(_y4m->dst_buf_sz);
  else
    _y4m->dst_buf = (unsigned char *)malloc(2 * _y4m->dst_buf_sz);

  if (_y4m->aux_buf_sz > 0)
    _y4m->aux_buf = (unsigned char *)malloc(_y4m->aux_buf_sz);
  return 0;
}

void y4m_input_close(y4m_input *_y4m) {
  free(_y4m->dst_buf);
  free(_y4m->aux_buf);
}

int y4m_input_fetch_frame(y4m_input *_y4m, FILE *_fin, aom_image_t *_img) {
  char frame[6];
  int pic_sz;
  int c_w;
  int c_h;
  int c_sz;
  int bytes_per_sample = _y4m->bit_depth > 8 ? 2 : 1;
  /*Read and skip the frame header.*/
  if (!file_read(frame, 6, _fin)) return 0;
  if (memcmp(frame, "FRAME", 5)) {
    fprintf(stderr, "Loss of framing in Y4M input data\n");
    return -1;
  }
  if (frame[5] != '\n') {
    char c;
    int j;
    for (j = 0; j < 79 && file_read(&c, 1, _fin) && c != '\n'; j++) {
    }
    if (j == 79) {
      fprintf(stderr, "Error parsing Y4M frame header\n");
      return -1;
    }
  }
  /*Read the frame data that needs no conversion.*/
  if (!file_read(_y4m->dst_buf, _y4m->dst_buf_read_sz, _fin)) {
    fprintf(stderr, "Error reading Y4M frame data.\n");
    return -1;
  }
  /*Read the frame data that does need conversion.*/
  if (!file_read(_y4m->aux_buf, _y4m->aux_buf_read_sz, _fin)) {
    fprintf(stderr, "Error reading Y4M frame data.\n");
    return -1;
  }
  /*Now convert the just read frame.*/
  (*_y4m->convert)(_y4m, _y4m->dst_buf, _y4m->aux_buf);
  /*Fill in the frame buffer pointers.
    We don't use aom_img_wrap() because it forces padding for odd picture
     sizes, which would require a separate fread call for every row.*/
  memset(_img, 0, sizeof(*_img));
  /*Y4M has the planes in Y'CbCr order, which libaom calls Y, U, and V.*/
  _img->fmt = _y4m->aom_fmt;
  _img->w = _img->d_w = _y4m->pic_w;
  _img->h = _img->d_h = _y4m->pic_h;
  _img->x_chroma_shift = _y4m->dst_c_dec_h >> 1;
  _img->y_chroma_shift = _y4m->dst_c_dec_v >> 1;
  _img->bps = _y4m->bps;

  /*Set up the buffer pointers.*/
  pic_sz = _y4m->pic_w * _y4m->pic_h * bytes_per_sample;
  c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  c_w *= bytes_per_sample;
  c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  _img->stride[AOM_PLANE_Y] = _y4m->pic_w * bytes_per_sample;
  _img->stride[AOM_PLANE_U] = _img->stride[AOM_PLANE_V] = c_w;
  _img->planes[AOM_PLANE_Y] = _y4m->dst_buf;
  _img->planes[AOM_PLANE_U] = _y4m->dst_buf + pic_sz;
  _img->planes[AOM_PLANE_V] = _y4m->dst_buf + pic_sz + c_sz;
  return 1;
}
