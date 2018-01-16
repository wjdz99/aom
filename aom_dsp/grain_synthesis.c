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

/*!\file
 * \brief Describes film grain parameters and film grain synthesis
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "grain_synthesis.h"
//#include "aom_dsp_common.h"

// static const int elem_in_gauss_sequence = 2048;
static const int gaussian_sequence[2048] = {
  -53,  -122, -52,  -50,  -83,  -70,  -84,  -121, 16,   -128, -182, -69,  -332,
  112,  54,   -190, -237, 138,  80,   -209, -266, 118,  132,  -88,  70,   137,
  33,   -392, 18,   -5,   -4,   -283, -84,  58,   88,   43,   209,  -218, -34,
  144,  -431, 18,   51,   -3,   114,  -103, -120, 290,  -141, -157, 116,  -9,
  9,    111,  65,   89,   53,   -9,   -191, -72,  -88,  111,  29,   -186, 237,
  70,   175,  -68,  -96,  112,  -43,  7,    -59,  -183, 113,  118,  -50,  -72,
  88,   -14,  -79,  -255, -177, -163, 16,   -349, -175, -9,   193,  68,   51,
  14,   -13,  63,   -163, -155, -100, 74,   -32,  71,   -11,  -61,  14,   -6,
  -124, -7,   -73,  -176, 152,  -205, 35,   -30,  79,   63,   27,   -25,  32,
  -82,  263,  -110, 119,  -26,  9,    108,  -81,  93,   -31,  -201, 49,   -125,
  30,   65,   141,  -366, 64,   -17,  216,  -119, 138,  -160, -67,  -151, -27,
  93,   -237, 17,   99,   167,  208,  -62,  118,  -66,  25,   -56,  -25,  0,
  -25,  86,   43,   -24,  -131, -24,  -180, -221, -98,  -148, -102, -12,  -220,
  173,  -77,  -28,  114,  60,   -155, 4,    151,  -59,  -58,  85,   -57,  -298,
  91,   310,  12,   30,   86,   265,  108,  -138, -115, -113, -131, 109,  -169,
  55,   148,  -37,  62,   41,   -3,   56,   -1,   155,  232,  216,  -4,   116,
  -179, -140, -13,  232,  103,  68,   -235, 310,  5,    38,   98,   90,   -46,
  -114, 145,  14,   -41,  171,  -115, 142,  -114, -214, -15,  169,  214,  304,
  33,   -163, -247, 218,  -172, -138, -29,  -83,  -104, -94,  99,   44,   -17,
  266,  -106, 102,  -9,   -14,  -84,  54,   -31,  -132, -126, -194, -200, -90,
  47,   -35,  -110, 159,  -188, 10,   99,   34,   141,  -8,   102,  7,    -76,
  182,  -7,   3,    193,  -66,  -138, -38,  -89,  58,   60,   -115, -113, 61,
  0,    -140, 8,    -20,  -126, 41,   -146, -25,  182,  -95,  -93,  211,  -8,
  -173, 153,  -76,  -10,  54,   5,    -60,  190,  -60,  -242, -78,  77,   189,
  39,   -74,  51,   -42,  -3,   -206, 26,   37,   187,  -86,  20,   -60,  227,
  -153, -235, 1,    -26,  -73,  -51,  24,   1,    -44,  10,   118,  -55,  94,
  144,  -101, 129,  20,   -158, 364,  95,   52,   -164, 69,   -90,  -109, -82,
  -142, -62,  4,    141,  -77,  53,   10,   -219, 96,   -29,  -107, -59,  112,
  146,  157,  344,  139,  -17,  -59,  -93,  -1,   164,  -64,  16,   -82,  -193,
  -182, 37,   153,  194,  57,   0,    33,   2,    53,   -136, -100, 24,   -44,
  -241, -22,  -91,  -3,   -90,  -68,  135,  -46,  -105, -113, -22,  -45,  88,
  -47,  98,   220,  -144, 43,   147,  -342, -171, 74,   128,  84,   -53,  -102,
  -286, -143, -30,  366,  7,    -47,  264,  84,   44,   -121, -2,   -96,  -141,
  155,  -2,   -145, -117, 158,  -235, 56,   -216, -19,  -32,  -159, -191, 8,
  -68,  53,   -49,  -28,  80,   5,    -90,  -97,  153,  -89,  -229, 73,   18,
  -110, -214, 274,  170,  177,  -141, 15,   -17,  14,   63,   -22,  -47,  33,
  -25,  12,   -29,  -7,   312,  -107, 104,  -104, 68,   -287, 234,  -8,   -192,
  -1,   -34,  106,  29,   -46,  100,  -151, -61,  -26,  86,   -28,  8,    -15,
  183,  79,   -176, -113, -141, 233,  -57,  -142, 126,  -104, 7,    -31,  116,
  -98,  -1,   74,   153,  -102, 14,   -189, 206,  195,  -102, 133,  19,   16,
  90,   31,   -4,   -117, -128, -7,   206,  -15,  274,  -18,  -37,  113,  32,
  -116, -112, -5,   135,  264,  -140, -66,  -84,  65,   13,   -36,  54,   -226,
  -217, -239, -98,  104,  -84,  -136, 210,  27,   -94,  -25,  73,   -116, -106,
  -183, 106,  -129, -172, -9,   -182, 51,   -88,  359,  -60,  -69,  -13,  99,
  -215, 17,   163,  -152, -30,  25,   -34,  -45,  39,   -12,  107,  13,   107,
  180,  -67,  185,  -136, -206, -34,  -166, 137,  -2,   -44,  79,   150,  67,
  65,   54,   -11,  -35,  -140, -214, -19,  65,   -177, 133,  -186, -4,   -206,
  81,   -126, -69,  9,    -211, -7,   -176, 87,   167,  -73,  147,  102,  143,
  -15,  81,   21,   54,   215,  36,   40,   356,  -67,  148,  34,   -225, -168,
  188,  -301, -204, -141, -68,  -99,  -90,  21,   258,  216,  -51,  -87,  171,
  185,  156,  39,   178,  71,   -54,  46,   -115, -196, -8,   -16,  -9,   -14,
  64,   -86,  -5,   136,  -15,  16,   -249, 37,   -78,  -9,   100,  -327, 220,
  107,  -40,  -119, -100, 167,  37,   -24,  -69,  -323, 219,  2,    77,   135,
  118,  -105, -62,  -66,  98,   -17,  -179, 47,   -88,  -56,  -208, 119,  29,
  -22,  -163, -90,  83,   -143, 15,   0,    -113, -19,  -259, 14,   321,  14,
  0,    55,   178,  -80,  121,  -2,   -104, 60,   95,   30,   -147, 90,   -104,
  194,  14,   -136, -194, 69,   -153, 30,   -339, -112, 178,  -105, 86,   237,
  -17,  225,  -208, -95,  214,  62,   -25,  -109, 170,  6,    153,  -38,  -58,
  -98,  107,  -38,  35,   28,   69,   18,   92,   -28,  -147, 107,  -192, -169,
  46,   23,   165,  -61,  187,  -136, 3,    -41,  -55,  23,   116,  41,   -18,
  198,  132,  -233, -55,  -32,  90,   167,  40,   -33,  114,  66,   109,  164,
  68,   -8,   25,   -178, -26,  20,   47,   6,    19,   14,   -121, -172, 49,
  44,   137,  -65,  -27,  -54,  278,  113,  197,  44,   117,  -80,  -33,  -222,
  90,   -27,  98,   -47,  -257, 28,   -50,  -51,  28,   248,  35,   144,  254,
  -163, 56,   111,  125,  -70,  -41,  16,   22,   163,  -254, -241, -51,  -52,
  -48,  -187, -165, -29,  171,  85,   94,   211,  229,  41,   -134, 7,    -49,
  141,  -62,  -364, -97,  -109, 12,   28,   18,   52,   110,  -289, 68,   -96,
  -2,   -60,  -54,  2,    -183, -151, -27,  -56,  103,  53,   96,   94,   102,
  224,  128,  -135, 162,  -59,  -122, 10,   -3,   -50,  -65,  -26,  -155, 33,
  -166, -92,  134,  49,   171,  -67,  -177, -250, -46,  -77,  171,  18,   -85,
  -108, 52,   122,  -137, 24,   204,  -92,  -106, -103, -21,  56,   77,   -50,
  109,  262,  -98,  244,  -13,  -26,  -50,  149,  40,   -27,  -28,  25,   93,
  66,   -166, 1,    212,  -17,  237,  93,   86,   -25,  89,   -124, -132, -189,
  37,   -101, 29,   -42,  129,  282,  127,  -100, 22,   -26,  3,    -114, 144,
  -23,  -65,  161,  -247, -131, 81,   23,   139,  -85,  -67,  124,  -63,  101,
  -135, 65,   47,   9,    -103, -79,  131,  -15,  11,   -5,   -508, 297,  -33,
  -116, -48,  -83,  151,  88,   22,   114,  31,   -3,   67,   42,   190,  -69,
  30,   -75,  -99,  205,  -61,  -29,  63,   -25,  205,  -182, -220, -133, -41,
  -184, -34,  -16,  47,   -185, 69,   -8,   -33,  -36,  -134, -44,  112,  53,
  52,   83,   -54,  43,   -119, -2,   9,    -43,  -12,  -243, -24,  -196, 183,
  -26,  -32,  23,   -141, -80,  -75,  57,   -61,  22,   56,   163,  -169, -18,
  248,  66,   4,    -26,  359,  -46,  -99,  -33,  -251, -106, -35,  -155, -60,
  112,  -65,  -52,  -127, 78,   -40,  56,   -37,  194,  301,  -99,  -21,  198,
  -51,  -115, 83,   212,  210,  47,   135,  -34,  74,   -79,  -50,  110,  -73,
  57,   310,  38,   -111, 57,   -215, -237, -212, 21,   148,  -102, -47,  -3,
  -88,  -125, -3,   -135, 175,  3,    -32,  58,   230,  -204, 157,  -149, 38,
  86,   -160, 74,   73,   -7,   -146, 205,  3,    -48,  -87,  135,  -233, -192,
  -31,  107,  -80,  12,   37,   243,  -65,  18,   -181, -70,  48,   112,  76,
  154,  -173, -370, -115, -212, -93,  7,    221,  -171, 178,  44,   232,  32,
  233,  -56,  15,   -55,  141,  -18,  3,    -7,   140,  24,   -48,  70,   -147,
  -40,  107,  -61,  -14,  -10,  6,    154,  87,   93,   129,  -206, 71,   -224,
  -79,  195,  -101, 14,   -27,  149,  -59,  -202, -56,  111,  -32,  20,   -113,
  114,  -169, -216, 135,  -98,  290,  92,   -61,  35,   35,   274,  -79,  101,
  -248, 23,   279,  248,  52,   79,   63,   77,   86,   149,  -49,  195,  104,
  98,   51,   81,   -34,  53,   -22,  -21,  -44,  76,   -142, -110, -34,  -205,
  -124, -44,  -5,   72,   -76,  -49,  16,   -252, -93,  -203, 58,   -312, 11,
  -1,   -136, -67,  -151, 237,  189,  36,   20,   72,   -230, -51,  -93,  -103,
  182,  -165, -40,  33,   92,   -206, 156,  233,  -151, -106, -29,  -292, 133,
  -138, 43,   128,  63,   97,   153,  55,   28,   -59,  -20,  -121, -104, -104,
  23,   1,    -59,  -266, -100, -183, 101,  -175, -131, -35,  77,   -15,  -54,
  45,   -84,  -163, 23,   -112, 73,   41,   123,  130,  281,  55,   -14,  119,
  283,  208,  -39,  134,  -28,  -249, -29,  39,   -161, 126,  -29,  -3,   99,
  158,  -198, 225,  -30,  15,   18,   -38,  -47,  43,   82,   -88,  -141, -104,
  -58,  104,  -139, -62,  103,  105,  -207, -189, 64,   -104, -184, -176, 35,
  -94,  -139, -53,  137,  219,  59,   -14,  -109, -237, -301, 91,   14,   15,
  321,  -40,  -116, 207,  -79,  -131, 156,  73,   26,   -35,  73,   124,  -173,
  11,   -141, -97,  76,   -105, 276,  0,    227,  78,   29,   -114, -35,  -79,
  -24,  -221, -85,  -63,  -5,   87,   155,  -108, 36,   102,  72,   -20,  238,
  4,    137,  -94,  109,  176,  -32,  166,  189,  -153, 43,   183,  -39,  91,
  153,  22,   -13,  100,  -73,  -29,  -26,  6,    26,   155,  -254, -45,  -104,
  46,   -119, -76,  -24,  54,   70,   -62,  102,  7,    207,  76,   85,   -45,
  375,  -52,  217,  162,  -104, 106,  196,  -141, 65,   -28,  -65,  -201, 7,
  52,   -142, -55,  39,   -96,  -60,  -121, -21,  -104, 77,   216,  70,   -15,
  139,  -129, -148, 23,   -270, 106,  8,    0,    -6,   273,  -202, 26,   -48,
  87,   -12,  -83,  50,   31,   187,  -24,  235,  -56,  35,   -257, 75,   185,
  -93,  16,   391,  -88,  -91,  -282, 150,  -300, 162,  83,   61,   167,  -203,
  -144, -177, 20,   -78,  -270, 94,   -77,  67,   -115, 20,   222,  -21,  -87,
  -60,  -11,  48,   5,    -112, -46,  206,  169,  -118, -16,  99,   52,   50,
  -227, 160,  -92,  -86,  92,   78,   57,   209,  -144, -170, 83,   82,   -29,
  59,   -114, 58,   -18,  176,  -49,  -89,  115,  40,   -92,  120,  -99,  28,
  -28,  -127, 145,  140,  -92,  -86,  15,   58,   -127, -108, -43,  -10,  -58,
  23,   124,  177,  -110, -77,  154,  39,   -64,  -28,  226,  -158, -193, 57,
  8,    199,  -2,   44,   -27,  55,   -104, -28,  -350, 39,   -36,  -20,  -84,
  105,  -153, 104,  31,   75,   76,   -167, -172, -11,  -70,  -131, 131,  -104,
  37,   -250, -12,  169,  -53,  228,  -7,   -5,   176,  119,  104,  12,   144,
  41,   105,  14,   182,  -49,  -7,   -106, 3,    155,  30,   -26,  241,  -80,
  -78,  11,   188,  -84,  -70,  -105, 91,   120,  186,  54,   -14,  -216, 242,
  216,  193,  152,  -55,  202,  -112, 119,  26,   70,   -51,  -17,  199,  -264,
  -70,  -126, 10,   -23,  -80,  -221, -24,  -63,  -52,  -273, 301,  -35,  103,
  -310, 40,   23,   -57,  -247, 2,    243,  42,   -162, -241, -78,  210,  -100,
  -39,  -16,  84,   -59,  170,  -24,  27,   -200, -22,  34,   3,    123,  117,
  266,  5,    -1,   -20,  96,   -32,  -60,  -2,   98,   -6,   3,    115,  5,
  140,  -177, 249,  200,  58,   205,  55,   188,  -112, 133,  -78,  139,  120,
  -170, -124, 211,  78,   -135, -94,  101,  31,   -29,  -150, 15,   23,   -34,
  -185, 40,   113,  -199, 19,   145,  -24,  -143, -52,  -8,   -183, -100, 121,
  144,  17,   -349, 2,    -89,  4,    76,   -104, 188,  52,   57,   251,  -55,
  -42,  -138, -83,  -341, -82,  32,   11,   228,  36,   -53,  -180, 65,   -183,
  48,   -250, 91,   -91,  205,  151,  -121, 42,   -50,  -78,  -18,  -71,  26,
  -10,  121,  -121, 77,   193,  -30,  -107, 96,   93,   -194, 130,  140,  178,
  -43,  -110, -241, 2,    34,   -74,  151,  -65,  -230, 221,  96,   -154, 177,
  -120, 48,   6,    -160, 185,  -123, -26,  -120, -223, 187,  -32,  61,   -170,
  -158, -205, 125,  -214, 31,   81,   40,   -59,  222,  192,  202,  -27,  134,
  0,    -119, -124, 221,  -27,  -115, 46,   -51,  -5,   14,   -135, 81,   81,
  177,  138,  -31,  -59,  75,   80,   31,   89,   162,  -160, 129,  -140, 40,
  -96,  105,  -180, 227,  -42,  52,   -11,  123,  154,  -207, 64,   -14,  -89,
  114,  -161, 263,  20,   79,   -72,  -46,  61,   112,  -71,  104,  81,   52,
  86,   -42,  -156, 183,  -165, -243, 222,  210,  -23,  34,   -59,  44,   52,
  -299, -294, 31,   -165, -28,  80,   111,  -15,  -115, -16,  -113, -200, -62,
  95,   -70,  -30,  -11,  -150, 279,  -76,  -97,  81,   32,   -92,  240,  -96,
  -43,  141,  -34,  -15,  69,   -64,  77,   93,   143,  270,  123,  -69,  218,
  -217, 103,  -21,  50,   -129, -46,  65,   388,  -196, 103,  46,   -49,  -111,
  177,  178,  124,  69,   64,   -177, -66,  -83,  -161, 157,  -179, 47,   252,
  90,   -308, 186,  -100, -8,   174,  42,   -46,  -38,  159,  -113, -36,  198,
  54,   -84,  130,  75,   20,   -231, -267, 3,    -237, 75,   63,   123,  11,
  174,  168,  15,   -20,  -62,  280,  -1,   -7,   149,  75,   237,  -159, -242,
  76,   -126, -62,  188,  -89,  107,  -84,  -181, -167, -38,  -73,  144,  3,
  332,  -93,  157,  -79,  142,  154,  -34,  25,   -41,  -34,  -89,  -46,  296,
  235,  135,  44,   192,  -55,  -221, 99
};

static const int gauss_bits = 11;

static const int luma_subblock_size = 32;
static const int chroma_subblock_size = 16;

static int min_luma_legal_range = 16;
static int max_luma_legal_range = 235;

static int min_chroma_legal_range = 16;
static int max_chroma_legal_range = 240;

static int scaling_lut_y[256];
static int scaling_lut_cb[256];
static int scaling_lut_cr[256];

static int grain_center;
static int grain_min;
static int grain_max;

static uint16_t random_register = 0;  // random number generator register

void assign_default(aom_film_grain_t *params) {
  params->apply_grain = 1;
  params->update_parameters = 1;
  
  int scaling_points_y[][2] = { { 16, 0 },   { 25, 17 },  { 33, 18 },
                                { 41, 20 },  { 48, 21 },  { 56, 17 },
                                { 67, 16 },  { 82, 18 },  { 97, 19 },
                                { 113, 18 }, { 128, 22 }, { 143, 21 },
                                { 158, 22 }, { 178, 23 } };  // 36

  params->num_y_points = 14;

  int scaling_points_cb[][2] = { { 16, 0 },   { 20, 8 },  { 28, 11 },
                                 { 60, 13 },  { 90, 17 }, { 105, 20 },
                                 { 134, 21 }, { 168, 26 } };

  params->num_cb_points = 8;

  int scaling_points_cr[][2] = { { 16, 0 },   { 28, 12 },  { 56, 10 },
                                 { 66, 12 },  { 80, 13 },  { 108, 12 },
                                 { 122, 14 }, { 137, 14 }, { 169, 22 } };

  params->num_cr_points = 9;

  memcpy(params->scaling_points_y, scaling_points_y, sizeof(int) * 14 * 2);
  memcpy(params->scaling_points_cb, scaling_points_cb, sizeof(int) * 8 * 2);
  memcpy(params->scaling_points_cr, scaling_points_cr, sizeof(int) * 9 * 2);

  static const int ar_coeff_lag = 2;

  static const int num_pos_luma = 2 * ar_coeff_lag * (ar_coeff_lag + 1);
  static const int num_pos_chroma = num_pos_luma + 1;

  int arc_y[12] = {
    0, 0, -58, 0, 0, 0, -76, 100, -43, 0, -51, 82
  };  //  ARcoeffs * 256;

  int arc_cb[13] = {
    0, 0, -49, 0, 0, 0, -36, 22, -30, 0, -38, 7, 39
  };  //  ARcoeffs * 256;

  int arc_cr[13] = {
    0, 0, -47, 0, 0, 0, -31, 31, -25, 0, -32, 13, -100
  };  //  ARcoeffs * 256;

  memcpy(params->ar_coeffs_y, arc_y, sizeof(int) * num_pos_luma);
  memcpy(params->ar_coeffs_cb, arc_cb, sizeof(int) * num_pos_chroma);
  memcpy(params->ar_coeffs_cr, arc_cr, sizeof(int) * num_pos_chroma);

  params->ar_coeff_lag = ar_coeff_lag;

  params->cb_mult = 247;       // 8 bit
  params->cb_luma_mult = 192;  // 8bit
  params->cb_offset = 18;      // 9bit

  params->cr_mult = 229;       // 8 bit
  params->cr_luma_mult = 192;  // 8bit
  params->cr_offset = 54;      // 9bit

  params->ar_coeff_scale = 8;

  params->overlap_flag = 1;

  params->full_range = 0;

  params->random_seed = 45231;

  params->bit_depth = 10;
}

void init_arrays(aom_film_grain_t *params, int luma_stride, int chroma_stride,
                 int ***pred_pos_luma, int ***pred_pos_chroma,
                 int **luma_grain_block, int **cb_grain_block,
                 int **cr_grain_block, int **y_line_buf, int **cb_line_buf,
                 int **cr_line_buf, int **y_col_buf, int **cb_col_buf,
                 int **cr_col_buf, int luma_grain_samples,
                 int chroma_grain_samples) {
  int num_pos_luma = 2 * params->ar_coeff_lag * (params->ar_coeff_lag + 1);
  int num_pos_chroma = num_pos_luma + 1;

  (*pred_pos_luma) = (int **)malloc(sizeof(int *) * num_pos_luma);

  for (int row = 0; row < num_pos_luma; row++) {
    (*pred_pos_luma)[row] = (int *)malloc(sizeof(int) * 3);
  }

  (*pred_pos_chroma) = (int **)malloc(sizeof(int *) * num_pos_chroma);
  for (int row = 0; row < num_pos_chroma; row++) {
    (*pred_pos_chroma)[row] = (int *)malloc(sizeof(int) * 3);
  }

  int pos_ar_index = 0;

  for (int row = -params->ar_coeff_lag; row < 0; row++) {
    for (int col = -params->ar_coeff_lag; col < params->ar_coeff_lag + 1;
         col++) {
      (*pred_pos_luma)[pos_ar_index][0] = row;
      (*pred_pos_luma)[pos_ar_index][1] = col;
      (*pred_pos_luma)[pos_ar_index][2] = 0;

      (*pred_pos_chroma)[pos_ar_index][0] = row;
      (*pred_pos_chroma)[pos_ar_index][1] = col;
      (*pred_pos_chroma)[pos_ar_index][2] = 0;
      ++pos_ar_index;
    }
  }

  for (int col = -params->ar_coeff_lag; col < 0; col++) {
    (*pred_pos_luma)[pos_ar_index][0] = 0;
    (*pred_pos_luma)[pos_ar_index][1] = col;
    (*pred_pos_luma)[pos_ar_index][2] = 0;

    (*pred_pos_chroma)[pos_ar_index][0] = 0;
    (*pred_pos_chroma)[pos_ar_index][1] = col;
    (*pred_pos_chroma)[pos_ar_index][2] = 0;

    ++pos_ar_index;
  }

  (*pred_pos_chroma)[pos_ar_index][0] = 0;
  (*pred_pos_chroma)[pos_ar_index][1] = 0;
  (*pred_pos_chroma)[pos_ar_index][2] = 1;

  *y_line_buf = (int *)malloc(sizeof(int) * luma_stride * 2);
  *cb_line_buf = (int *)malloc(sizeof(int) * chroma_stride);
  *cr_line_buf = (int *)malloc(sizeof(int) * chroma_stride);

  *y_col_buf = (int *)malloc(sizeof(int) * (luma_subblock_size + 2) * 2);
  *cb_col_buf = (int *)malloc(sizeof(int) * (chroma_subblock_size + 1));
  *cr_col_buf = (int *)malloc(sizeof(int) * (chroma_subblock_size + 1));

  *luma_grain_block = (int *)malloc(sizeof(int) * luma_grain_samples);
  *cb_grain_block = (int *)malloc(sizeof(int) * chroma_grain_samples);
  *cr_grain_block = (int *)malloc(sizeof(int) * chroma_grain_samples);
}

void dealloc_arrays(aom_film_grain_t *params, int ***pred_pos_luma,
                    int ***pred_pos_chroma, int **luma_grain_block,
                    int **cb_grain_block, int **cr_grain_block,
                    int **y_line_buf, int **cb_line_buf, int **cr_line_buf,
                    int **y_col_buf, int **cb_col_buf, int **cr_col_buf) {
  int num_pos_luma = 2 * params->ar_coeff_lag * (params->ar_coeff_lag + 1);
  int num_pos_chroma = num_pos_luma + 1;

  for (int row = 0; row < num_pos_luma; row++) {
    free((*pred_pos_luma)[row]);
  }
  free(*pred_pos_luma);

  for (int row = 0; row < num_pos_chroma; row++) {
    free((*pred_pos_chroma)[row]);
  }
  free((*pred_pos_chroma));

  free(*y_line_buf);

  free(*cb_line_buf);

  free(*cr_line_buf);

  free(*y_col_buf);

  free(*cb_col_buf);

  free(*cr_col_buf);

  free(*luma_grain_block);

  free(*cb_grain_block);

  free(*cr_grain_block);
}

int get_random_number(int bits)  // between 0 and 2^bits - 1
{
  uint16_t bit;
  bit = ((random_register >> 0) ^ (random_register >> 1) ^
         (random_register >> 3) ^ (random_register >> 12)) &
        1;
  random_register = (random_register >> 1) | (bit << 15);
  return (random_register >> (16 - bits)) & ((1 << bits) - 1);
}

void init_random_generator(int lumaLine, uint16_t seed) {
  // same for the picture

  uint16_t msb = (seed >> 8) & 255;
  uint16_t lsb = seed & 255;

  random_register = (msb << 8) + lsb;

  //  changes for each row
  int luma_num = lumaLine >> 5;

  random_register ^= ((luma_num * 37 + 178) & 255) << 8;
  random_register ^= ((luma_num * 173 + 105) & 255);
}

void generate_luma_grain_block(aom_film_grain_t *params, int **pred_pos_luma,
                               int *luma_grain_block, int luma_block_size_y,
                               int luma_block_size_x, int luma_grain_stride,
                               int left_pad, int top_pad, int right_pad,
                               int bottom_pad) {
  int bit_depth = params->bit_depth;

  int num_pos_luma = 2 * params->ar_coeff_lag * (params->ar_coeff_lag + 1);

  if (bit_depth > 10)  // can be removed if support for 12 bits in not needed
                       // or gaussian table precision increased
  {
    for (int i = 0; i < luma_block_size_y; i++)
      for (int j = 0; j < luma_block_size_x; j++)
        luma_grain_block[i * luma_grain_stride + j] =
            (gaussian_sequence[get_random_number(gauss_bits)]
             << (bit_depth - 10));
  } else {
    for (int i = 0; i < luma_block_size_y; i++)
      for (int j = 0; j < luma_block_size_x; j++)
        luma_grain_block[i * luma_grain_stride + j] =
            (gaussian_sequence[get_random_number(gauss_bits)] +
             ((1 << (10 - bit_depth)) >> 1)) >>
            (10 - bit_depth);
  }

  for (int i = top_pad; i < luma_block_size_y - bottom_pad; i++)
    for (int j = left_pad; j < luma_block_size_x - right_pad; j++) {
      int wsum = 0;
      for (int pos = 0; pos < num_pos_luma; pos++) {
        wsum = wsum + params->ar_coeffs_y[pos] *
                          luma_grain_block[(i + pred_pos_luma[pos][0]) *
                                               luma_grain_stride +
                                           j + pred_pos_luma[pos][1]];
      }
      luma_grain_block[i * luma_grain_stride + j] =
          clamp(luma_grain_block[i * luma_grain_stride + j] +
                    (wsum >> params->ar_coeff_scale),
                grain_min, grain_max);
    }
}

void generate_chroma_grain_blocks(
    aom_film_grain_t *params,
    //                                  int** pred_pos_luma,
    int **pred_pos_chroma, int *luma_grain_block, int *cb_grain_block,
    int *cr_grain_block, int luma_grain_stride, int chroma_block_size_y,
    int chroma_block_size_x, int chroma_grain_stride, int left_pad, int top_pad,
    int right_pad, int bottom_pad) {
  int bit_depth = params->bit_depth;

  int num_pos_chroma =
      2 * params->ar_coeff_lag * (params->ar_coeff_lag + 1) + 1;

  if (bit_depth > 10)  // can be removed if support for 12 bits in not needed
                       // or gaussian table precision increased
  {
    for (int i = 0; i < chroma_block_size_y; i++)
      for (int j = 0; j < chroma_block_size_x; j++) {
        cb_grain_block[i * chroma_grain_stride + j] =
            gaussian_sequence[get_random_number(gauss_bits)]
            << (bit_depth - 10);
        cr_grain_block[i * chroma_grain_stride + j] =
            gaussian_sequence[get_random_number(gauss_bits)]
            << (bit_depth - 10);
      }
  } else {
    for (int i = 0; i < chroma_block_size_y; i++)
      for (int j = 0; j < chroma_block_size_x; j++) {
        cb_grain_block[i * chroma_grain_stride + j] =
            (gaussian_sequence[get_random_number(gauss_bits)] +
             ((1 << (10 - bit_depth)) >> 1)) >>
            (10 - bit_depth);
        cr_grain_block[i * chroma_grain_stride + j] =
            (gaussian_sequence[get_random_number(gauss_bits)] +
             ((1 << (10 - bit_depth)) >> 1)) >>
            (10 - bit_depth);
      }
  }

  for (int i = top_pad; i < chroma_block_size_y - bottom_pad; i++)
    for (int j = left_pad; j < chroma_block_size_x - right_pad; j++) {
      int wsum_cb = 0;
      int wsum_cr = 0;
      for (int pos = 0; pos < num_pos_chroma; pos++) {
        if (pred_pos_chroma[pos][2] == 0) {
          wsum_cb = wsum_cb + params->ar_coeffs_cb[pos] *
                                  cb_grain_block[(i + pred_pos_chroma[pos][0]) *
                                                     chroma_grain_stride +
                                                 j + pred_pos_chroma[pos][1]];
          wsum_cr = wsum_cr + params->ar_coeffs_cr[pos] *
                                  cr_grain_block[(i + pred_pos_chroma[pos][0]) *
                                                     chroma_grain_stride +
                                                 j + pred_pos_chroma[pos][1]];
        } else if (pred_pos_chroma[pos][2] == 1) {
          int av_luma =
              (luma_grain_block[(((i - top_pad) << 1) + top_pad) *
                                    luma_grain_stride +
                                ((j - left_pad) << 1) + left_pad] +
               luma_grain_block[(((i - top_pad) << 1) + top_pad + 1) *
                                    luma_grain_stride +
                                ((j - left_pad) << 1) + left_pad] +
               luma_grain_block[(((i - top_pad) << 1) + top_pad) *
                                    luma_grain_stride +
                                ((j - left_pad) << 1) + left_pad + 1] +
               luma_grain_block[(((i - top_pad) << 1) + top_pad + 1) *
                                    luma_grain_stride +
                                ((j - left_pad) << 1) + left_pad + 1] +
               2) >>
              2;

          wsum_cb = wsum_cb + params->ar_coeffs_cb[pos] * av_luma;
          wsum_cr = wsum_cr + params->ar_coeffs_cr[pos] * av_luma;
        } else {
          printf(
              "Grain synthesis: prediction between two chroma components is "
              "not supported!");
          exit(1);
        }
      }
      cb_grain_block[i * chroma_grain_stride + j] =
          clamp(cb_grain_block[i * chroma_grain_stride + j] +
                    (wsum_cb >> params->ar_coeff_scale),
                grain_min, grain_max);
      cr_grain_block[i * chroma_grain_stride + j] =
          clamp(cr_grain_block[i * chroma_grain_stride + j] +
                    (wsum_cr >> params->ar_coeff_scale),
                grain_min, grain_max);
    }
}

void init_scaling_function(int scaling_points[][2], int num_points,
                           int scaling_lut[]) {
  for (int i = 0; i < scaling_points[0][0]; i++)
    scaling_lut[i] = scaling_points[0][1];

  for (int point = 0; point < num_points - 1; point++) {
    int delta_y = scaling_points[point + 1][1] - scaling_points[point][1];
    int delta_x = scaling_points[point + 1][0] - scaling_points[point][0];

    int delta = delta_y * ((16384 + (delta_x >> 1)) / delta_x);

    for (int x = 0; x < delta_x; x++) {
      scaling_lut[scaling_points[point][0] + x] =
          scaling_points[point][1] + ((x * delta + 8192) >> 14);
    }
  }

  for (int i = scaling_points[num_points - 1][0]; i < 256; i++)
    scaling_lut[i] = scaling_points[num_points - 1][1];
}

// function that extracts samples from a LUT (and interpolates intemediate
// frames for 10 bit video)
int scale_LUT(int *scaling_lut, int index, int bit_depth) {
  int x = index >> (bit_depth - 8);

  if (!(bit_depth - 8) || x >= 255)
    return scaling_lut[x];
  else
    return scaling_lut[x] +
           (((scaling_lut[x + 1] - scaling_lut[x]) * (index & 3) + 2) >> 2);
}

void add_noise_to_block(aom_film_grain_t *params, uint8_t *luma, uint8_t *cb,
                        uint8_t *cr, int luma_stride, int chroma_stride,
                        int *luma_grain, int *cb_grain, int *cr_grain,
                        int luma_grain_stride, int chroma_grain_stride,
                        int chroma_height, int chroma_width, int bit_depth) {
  int cb_mult = params->cb_mult - 128;            // fixed scale
  int cb_luma_mult = params->cb_luma_mult - 128;  // fixed scale
  int cb_offset = params->cb_offset - 256;

  int cr_mult = params->cr_mult - 128;            // fixed scale
  int cr_luma_mult = params->cr_luma_mult - 128;  // fixed scale
  int cr_offset = params->cr_offset - 256;

  int min_luma, max_luma, min_chroma, max_chroma;

  if (params->full_range) {
    min_luma = min_chroma = 0;
    max_luma = max_chroma = 255;
  } else {
    min_luma = min_luma_legal_range;
    max_luma = max_luma_legal_range;

    min_chroma = min_chroma_legal_range;
    max_chroma = max_chroma_legal_range;
  }

  for (int i = 0; i < chroma_height; i++) {
    for (int j = 0; j < chroma_width; j++) {
      int average_luma = (luma[(i << 1) * luma_stride + (j << 1)] +
                          luma[((i << 1)) * luma_stride + (j << 1) + 1] + 1) >>
                         1;

      luma[((i) << 1) * luma_stride + ((j) << 1)] = clamp(
          luma[((i) << 1) * luma_stride + ((j) << 1)] +
              ((scale_LUT(scaling_lut_y,
                          luma[((i) << 1) * luma_stride + ((j) << 1)], 8) *
                luma_grain[(i << 1) * luma_grain_stride + (j << 1)]) >>
               8),
          min_luma, max_luma);
      luma[(((i) << 1) + 1) * luma_stride + ((j) << 1)] = clamp(
          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1)] +
              ((scale_LUT(scaling_lut_y,
                          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1)],
                          8) *
                luma_grain[((i << 1) + 1) * luma_grain_stride + (j << 1)]) >>
               8),
          min_luma, max_luma);
      luma[(((i) << 1)) * luma_stride + ((j) << 1) + 1] = clamp(
          luma[(((i) << 1)) * luma_stride + ((j) << 1) + 1] +
              ((scale_LUT(scaling_lut_y,
                          luma[(((i) << 1)) * luma_stride + ((j) << 1) + 1],
                          8) *
                luma_grain[(i << 1) * luma_grain_stride + (j << 1) + 1]) >>
               8),
          min_luma, max_luma);
      luma[(((i) << 1) + 1) * luma_stride + ((j) << 1) + 1] = clamp(
          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1) + 1] +
              ((scale_LUT(scaling_lut_y,
                          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1) + 1],
                          8) *
                luma_grain[((i << 1) + 1) * luma_grain_stride + (j << 1) +
                           1]) >>
               8),
          min_luma, max_luma);

      cb[i * chroma_stride + j] =
          clamp(cb[i * chroma_stride + j] +
                    ((scale_LUT(scaling_lut_cb,
                                clamp(((average_luma * cb_luma_mult +
                                        cb_mult * cb[i * chroma_stride + j]) >>
                                       6) +
                                          cb_offset,
                                      0, (256 << (bit_depth - 8)) - 1),
                                8) *
                      cb_grain[i * chroma_grain_stride + j]) >>
                     8),
                min_chroma, max_chroma);

      cr[i * chroma_stride + j] =
          clamp(cr[i * chroma_stride + j] +
                    ((scale_LUT(scaling_lut_cr,
                                clamp(((average_luma * cr_luma_mult +
                                        cr_mult * cr[i * chroma_stride + j]) >>
                                       6) +
                                          cr_offset,
                                      0, (256 << (bit_depth - 8)) - 1),
                                8) *
                      cr_grain[i * chroma_grain_stride + j]) >>
                     8),
                min_chroma, max_chroma);
    }
  }
}

void add_noise_to_block_hbd(aom_film_grain_t *params, uint16_t *luma,
                            uint16_t *cb, uint16_t *cr, int luma_stride,
                            int chroma_stride, int *luma_grain, int *cb_grain,
                            int *cr_grain, int luma_grain_stride,
                            int chroma_grain_stride, int chroma_height,
                            int chroma_width, int bit_depth) {
  int cb_mult = params->cb_mult - 128;            // fixed scale
  int cb_luma_mult = params->cb_luma_mult - 128;  // fixed scale
  // offset value depends on the bit depth
  int cb_offset = (params->cb_offset << (bit_depth - 8)) - (1 << bit_depth);

  int cr_mult = params->cr_mult - 128;            // fixed scale
  int cr_luma_mult = params->cr_luma_mult - 128;  // fixed scale
  // offset value depends on the bit depth
  int cr_offset = (params->cr_offset << (bit_depth - 8)) - (1 << bit_depth);

  int min_luma, max_luma, min_chroma, max_chroma;

  if (params->full_range) {
    min_luma = min_chroma = 0;
    max_luma = max_chroma = (256 << (bit_depth - 8)) - 1;
  } else {
    min_luma = min_luma_legal_range << (bit_depth - 8);
    max_luma = max_luma_legal_range << (bit_depth - 8);

    min_chroma = min_chroma_legal_range << (bit_depth - 8);
    max_chroma = max_chroma_legal_range << (bit_depth - 8);
  }

  for (int i = 0; i < chroma_height; i++) {
    for (int j = 0; j < chroma_width; j++) {
      int average_luma = (luma[(i << 1) * luma_stride + (j << 1)] +
                          luma[((i << 1)) * luma_stride + (j << 1) + 1] + 1) >>
                         1;

      luma[((i) << 1) * luma_stride + ((j) << 1)] =
          clamp(luma[((i) << 1) * luma_stride + ((j) << 1)] +
                    ((scale_LUT(scaling_lut_y,
                                luma[((i) << 1) * luma_stride + ((j) << 1)],
                                bit_depth) *
                      luma_grain[(i << 1) * luma_grain_stride + (j << 1)]) >>
                     8),
                min_luma, max_luma);
      luma[(((i) << 1) + 1) * luma_stride + ((j) << 1)] = clamp(
          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1)] +
              ((scale_LUT(scaling_lut_y,
                          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1)],
                          bit_depth) *
                luma_grain[((i << 1) + 1) * luma_grain_stride + (j << 1)]) >>
               8),
          min_luma, max_luma);
      luma[(((i) << 1)) * luma_stride + ((j) << 1) + 1] = clamp(
          luma[(((i) << 1)) * luma_stride + ((j) << 1) + 1] +
              ((scale_LUT(scaling_lut_y,
                          luma[(((i) << 1)) * luma_stride + ((j) << 1) + 1],
                          bit_depth) *
                luma_grain[(i << 1) * luma_grain_stride + (j << 1) + 1]) >>
               8),
          min_luma, max_luma);
      luma[(((i) << 1) + 1) * luma_stride + ((j) << 1) + 1] = clamp(
          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1) + 1] +
              ((scale_LUT(scaling_lut_y,
                          luma[(((i) << 1) + 1) * luma_stride + ((j) << 1) + 1],
                          bit_depth) *
                luma_grain[((i << 1) + 1) * luma_grain_stride + (j << 1) +
                           1]) >>
               8),
          min_luma, max_luma);

      cb[i * chroma_stride + j] =
          clamp(cb[i * chroma_stride + j] +
                    ((scale_LUT(scaling_lut_cb,
                                clamp(((average_luma * cb_luma_mult +
                                        cb_mult * cb[i * chroma_stride + j]) >>
                                       6) +
                                          cb_offset,
                                      0, (256 << (bit_depth - 8)) - 1),
                                bit_depth) *
                      cb_grain[i * chroma_grain_stride + j]) >>
                     8),
                min_chroma, max_chroma);

      cr[i * chroma_stride + j] =
          clamp(cr[i * chroma_stride + j] +
                    ((scale_LUT(scaling_lut_cr,
                                clamp(((average_luma * cr_luma_mult +
                                        cr_mult * cr[i * chroma_stride + j]) >>
                                       6) +
                                          cr_offset,
                                      0, (256 << (bit_depth - 8)) - 1),
                                bit_depth) *
                      cr_grain[i * chroma_grain_stride + j]) >>
                     8),
                min_chroma, max_chroma);
    }
  }
}

void copy_rect(uint8_t *src, int src_stride, uint8_t *dst, int dst_stride, int width, int height, int high_bit_depth)
{
  int hbd_coeff = high_bit_depth ? 2 : 0;
  while(height)
  {
    memcpy(dst, src, width * sizeof(uint8_t) * hbd_coeff);
    src += src_stride;
    dst += dst_stride;
    --height;
  }
  return;
}


void add_film_grain(aom_film_grain_t *params, aom_image_t* src, aom_image_t* dst)
{
  uint8_t *luma, *cb,  *cr;
  int height,  width, luma_stride, chroma_stride;
  
  if ( !(src->fmt == AOM_IMG_FMT_I42016) && !(src->fmt == AOM_IMG_FMT_I420))
  {
    printf("Film grain error: only 4:2:0 is currently supported!");
    exit(1);
  }
  
  dst->r_w = src->r_w;
  dst->r_h = src->r_h;
    
  copy_rect(src->planes[AOM_PLANE_Y], src->stride[AOM_PLANE_Y],
                 dst->planes[AOM_PLANE_Y], dst->stride[AOM_PLANE_Y],
                 dst->d_w, dst->d_h, src->bit_depth);
  
  copy_rect(src->planes[AOM_PLANE_U], src->stride[AOM_PLANE_U],
            dst->planes[AOM_PLANE_U], dst->stride[AOM_PLANE_U],
            dst->d_w/2, dst->d_h/2, src->bit_depth);

  copy_rect(src->planes[AOM_PLANE_V], src->stride[AOM_PLANE_V],
            dst->planes[AOM_PLANE_V], dst->stride[AOM_PLANE_V],
            dst->d_w/2, dst->d_h/2, src->bit_depth);
  
  luma = dst->planes[AOM_PLANE_Y];
  cb = dst->planes[AOM_PLANE_U];
  cr = dst->planes[AOM_PLANE_V];
  luma_stride = (dst->bit_depth > 8) ? dst->stride[AOM_PLANE_Y]/2 : dst->stride[AOM_PLANE_Y];
  chroma_stride = (dst->bit_depth > 8) ? dst->stride[AOM_PLANE_U]/2 : dst->stride[AOM_PLANE_U];
  width = dst->d_w;
  height = dst->d_h;
  params->bit_depth = dst->bit_depth;
  
  add_film_grain_run(params, luma, cb,
                 cr, height, width, luma_stride,
                 chroma_stride);
  return;
}

  
void add_film_grain_run(aom_film_grain_t *params, uint8_t *luma, uint8_t *cb,
                    uint8_t *cr, int height, int width, int luma_stride,
                    int chroma_stride) {
  uint16_t *luma_hbd;
  uint16_t *cb_hbd;
  uint16_t *cr_hbd;
  int use_high_bit_depth = 0;

  if (params->bit_depth > 8) {
    luma_hbd = (uint16_t *)luma;
    cb_hbd = (uint16_t *)cb;
    cr_hbd = (uint16_t *)cr;
    use_high_bit_depth = 1;
  }

  int **pred_pos_luma;
  int **pred_pos_chroma;
  int *luma_grain_block;
  int *cb_grain_block;
  int *cr_grain_block;

  int *y_line_buf;
  int *cb_line_buf;
  int *cr_line_buf;

  int *y_col_buf;
  int *cb_col_buf;
  int *cr_col_buf;

  random_register = params->random_seed;

  int left_pad = 3;
  int right_pad = 3;  // padding to offset for AR coefficients
  int top_pad = 3;
  int bottom_pad = 0;

  int ar_padding = 3;  // maximum lag used for stabilization of AR coefficients

  // Initial padding is only needed for generation of
  // film grain templates (to stabilize the AR process)
  // Only a 64x64 luma and 32x32 chroma part of a template
  // is used later for adding grain, padding can be discarded

  int luma_block_size_y =
      top_pad + 2 * ar_padding + luma_subblock_size * 2 + bottom_pad;
  int luma_block_size_x = left_pad + 2 * ar_padding + luma_subblock_size * 2 +
                          2 * ar_padding + right_pad;

  int chroma_block_size_y =
      top_pad + ar_padding + chroma_subblock_size * 2 + bottom_pad;
  int chroma_block_size_x =
      left_pad + ar_padding + chroma_subblock_size * 2 + ar_padding + right_pad;

  int luma_grain_stride = luma_block_size_x;
  int chroma_grain_stride = chroma_block_size_x;

  int overlap = params->overlap_flag;
  int bit_depth = params->bit_depth;

  grain_center = 128 << (bit_depth - 8);
  grain_min = 0 - grain_center;
  grain_max = (256 << (bit_depth - 8)) - 1 - grain_center;

  init_arrays(params, luma_stride, chroma_stride, &pred_pos_luma,
              &pred_pos_chroma, &luma_grain_block, &cb_grain_block,
              &cr_grain_block, &y_line_buf, &cb_line_buf, &cr_line_buf,
              &y_col_buf, &cb_col_buf, &cr_col_buf,
              luma_block_size_y * luma_block_size_x,
              chroma_block_size_y * chroma_block_size_x);

  generate_luma_grain_block(params, pred_pos_luma, luma_grain_block,
                            luma_block_size_y, luma_block_size_x,
                            luma_grain_stride, left_pad, top_pad, right_pad,
                            bottom_pad);

  generate_chroma_grain_blocks(
      params,
      //                               pred_pos_luma,
      pred_pos_chroma, luma_grain_block, cb_grain_block, cr_grain_block,
      luma_grain_stride, chroma_block_size_y, chroma_block_size_x,
      chroma_grain_stride, left_pad, top_pad, right_pad, bottom_pad);

  init_scaling_function(params->scaling_points_y, params->num_y_points,
                        scaling_lut_y);
  init_scaling_function(params->scaling_points_cb, params->num_cb_points,
                        scaling_lut_cb);
  init_scaling_function(params->scaling_points_cr, params->num_cr_points,
                        scaling_lut_cr);

  for (int y = 0; y < height / 2; y += chroma_subblock_size) {
    init_random_generator(y * 2, params->random_seed);

    for (int x = 0; x < width / 2; x += chroma_subblock_size) {
      int offset_y = get_random_number(8);
      int offset_x = (offset_y >> 4) & 15;
      offset_y &= 15;

      int luma_offset_y = left_pad + 2 * ar_padding + (offset_y << 1);
      int luma_offset_x = top_pad + 2 * ar_padding + (offset_x << 1);

      int chroma_offset_y = top_pad + ar_padding + offset_y;
      int chroma_offset_x = left_pad + ar_padding + offset_x;

      if (overlap && x) {
        for (int i = 0; i < AOMMIN(chroma_subblock_size, height / 2 - y); i++) {
          y_col_buf[(i << 1) * 2] =
              clamp((27 * y_col_buf[(i << 1) * 2] +
                     17 * luma_grain_block[(luma_offset_y + (i << 1)) *
                                               luma_grain_stride +
                                           luma_offset_x] +
                     16) >>
                        5,
                    grain_min, grain_max);
          ;
          y_col_buf[(i << 1) * 2 + 1] =
              clamp((17 * y_col_buf[(i << 1) * 2 + 1] +
                     27 * luma_grain_block[(luma_offset_y + (i << 1)) *
                                               luma_grain_stride +
                                           luma_offset_x + 1] +
                     16) >>
                        5,
                    grain_min, grain_max);

          y_col_buf[(i << 1) * 2 + 2] =
              clamp((27 * y_col_buf[(i << 1) * 2 + 2] +
                     17 * luma_grain_block[(luma_offset_y + (i << 1) + 1) *
                                               luma_grain_stride +
                                           luma_offset_x] +
                     16) >>
                        5,
                    grain_min, grain_max);
          y_col_buf[(i << 1) * 2 + 3] =
              clamp((17 * y_col_buf[(i << 1) * 2 + 3] +
                     27 * luma_grain_block[(luma_offset_y + (i << 1) + 1) *
                                               luma_grain_stride +
                                           luma_offset_x + 1] +
                     16) >>
                        5,
                    grain_min, grain_max);

          cb_col_buf[i] = clamp(
              (23 * cb_col_buf[i] +
               22 * cb_grain_block[(chroma_offset_y + i) * chroma_grain_stride +
                                   chroma_offset_x] +
               16) >>
                  5,
              grain_min, grain_max);
          cr_col_buf[i] = clamp(
              (23 * cr_col_buf[i] +
               22 * cr_grain_block[(chroma_offset_y + i) * chroma_grain_stride +
                                   chroma_offset_x] +
               16) >>
                  5,
              grain_min, grain_max);
        }

        int i = y ? 1 : 0;

        if (use_high_bit_depth) {
          add_noise_to_block_hbd(
              params,
              (uint16_t *)luma + ((y + i) << 1) * luma_stride + (x << 1),
              (uint16_t *)cb + (y + i) * chroma_stride + x,
              (uint16_t *)cr + (y + i) * chroma_stride + x, luma_stride,
              chroma_stride, y_col_buf + i * 4, cb_col_buf + i, cr_col_buf + i,
              2, 1, AOMMIN(chroma_subblock_size, height / 2 - y) - i, 1,
              bit_depth);
        } else {
          add_noise_to_block(
              params, luma + ((y + i) << 1) * luma_stride + (x << 1),
              cb + (y + i) * chroma_stride + x,
              cr + (y + i) * chroma_stride + x, luma_stride, chroma_stride,
              y_col_buf + i * 4, cb_col_buf + i, cr_col_buf + i, 2, 1,
              AOMMIN(chroma_subblock_size, height / 2 - y) - i, 1, bit_depth);
        }
      }

      if (overlap && y) {
        if (x) {
          y_line_buf[x << 1] =
              clamp((17 * y_col_buf[0] + 27 * y_line_buf[x << 1] + 16) >> 5,
                    grain_min, grain_max);
          y_line_buf[(x << 1) + 1] = clamp(
              (17 * y_col_buf[1] + 27 * y_line_buf[(x << 1) + 1] + 16) >> 5,
              grain_min, grain_max);
          y_line_buf[luma_stride + (x << 1)] =
              clamp((27 * y_col_buf[2] +
                     17 * y_line_buf[luma_stride + (x << 1)] + 16) >>
                        5,
                    grain_min, grain_max);
          y_line_buf[luma_stride + (x << 1) + 1] =
              clamp((27 * y_col_buf[3] +
                     17 * y_line_buf[luma_stride + (x << 1) + 1] + 16) >>
                        5,
                    grain_min, grain_max);

          cb_line_buf[x] =
              clamp((22 * cb_col_buf[0] + 23 * cb_line_buf[x] + 16) >> 5,
                    grain_min, grain_max);
          cr_line_buf[x] =
              clamp((22 * cr_col_buf[0] + 23 * cr_line_buf[x] + 16) >> 5,
                    grain_min, grain_max);
        }

        for (int j = x ? 1 : 0; j < AOMMIN(chroma_subblock_size, width / 2 - x);
             j++) {
          y_line_buf[(x + j) << 1] =
              clamp((27 * y_line_buf[(x + j) << 1] +
                     17 * luma_grain_block[luma_offset_y * luma_grain_stride +
                                           luma_offset_x + (j << 1)] +
                     16) >>
                        5,
                    grain_min, grain_max);
          y_line_buf[((x + j) << 1) + 1] =
              clamp((27 * y_line_buf[((x + j) << 1) + 1] +
                     17 * luma_grain_block[luma_offset_y * luma_grain_stride +
                                           luma_offset_x + (j << 1) + 1] +
                     16) >>
                        5,
                    grain_min, grain_max);

          y_line_buf[luma_stride + ((x + j) << 1)] = clamp(
              (17 * y_line_buf[luma_stride + ((x + j) << 1)] +
               27 * luma_grain_block[(luma_offset_y + 1) * luma_grain_stride +
                                     luma_offset_x + (j << 1)] +
               16) >>
                  5,
              grain_min, grain_max);
          y_line_buf[luma_stride + ((x + j) << 1) + 1] = clamp(
              (17 * y_line_buf[luma_stride + ((x + j) << 1) + 1] +
               27 * luma_grain_block[(luma_offset_y + 1) * luma_grain_stride +
                                     luma_offset_x + (j << 1) + 1] +
               16) >>
                  5,
              grain_min, grain_max);

          cb_line_buf[x + j] =
              clamp((23 * cb_line_buf[x + j] +
                     22 * cb_grain_block[chroma_offset_y * chroma_grain_stride +
                                         chroma_offset_x + j] +
                     16) >>
                        5,
                    grain_min, grain_max);
          cr_line_buf[x + j] =
              clamp((23 * cr_line_buf[x + j] +
                     22 * cr_grain_block[chroma_offset_y * chroma_grain_stride +
                                         chroma_offset_x + j] +
                     16) >>
                        5,
                    grain_min, grain_max);
        }

        if (use_high_bit_depth) {
          add_noise_to_block_hbd(
              params, (uint16_t *)luma + (y << 1) * luma_stride + (x << 1),
              (uint16_t *)cb + y * chroma_stride + x,
              (uint16_t *)cr + y * chroma_stride + x, luma_stride,
              chroma_stride, y_line_buf + (x << 1), cb_line_buf + x,
              cr_line_buf + x, luma_stride, chroma_stride, 1,
              AOMMIN(chroma_subblock_size, width / 2 - x), bit_depth);
        } else {
          add_noise_to_block(
              params, luma + (y << 1) * luma_stride + (x << 1),
              cb + y * chroma_stride + x, cr + y * chroma_stride + x,
              luma_stride, chroma_stride, y_line_buf + (x << 1),
              cb_line_buf + x, cr_line_buf + x, luma_stride, chroma_stride, 1,
              AOMMIN(chroma_subblock_size, width / 2 - x), bit_depth);
        }
      }

      int i = overlap && y ? 1 : 0;
      int j = overlap && x ? 1 : 0;

      if (use_high_bit_depth) {
        add_noise_to_block_hbd(
            params,
            (uint16_t *)luma + ((y + i) << 1) * luma_stride + ((x + j) << 1),
            (uint16_t *)cb + (y + i) * chroma_stride + x + j,
            (uint16_t *)cr + (y + i) * chroma_stride + x + j, luma_stride,
            chroma_stride,
            luma_grain_block + (luma_offset_y + (i << 1)) * luma_grain_stride +
                luma_offset_x + (j << 1),
            cb_grain_block + (chroma_offset_y + i) * chroma_grain_stride +
                chroma_offset_x + j,
            cr_grain_block + (chroma_offset_y + i) * chroma_grain_stride +
                chroma_offset_x + j,
            luma_grain_stride, chroma_grain_stride,
            AOMMIN(chroma_subblock_size, height / 2 - y) - i,
            AOMMIN(chroma_subblock_size, width / 2 - x) - j, bit_depth);
      } else {
        add_noise_to_block(
            params, luma + ((y + i) << 1) * luma_stride + ((x + j) << 1),
            cb + (y + i) * chroma_stride + x + j,
            cr + (y + i) * chroma_stride + x + j, luma_stride, chroma_stride,
            luma_grain_block + (luma_offset_y + (i << 1)) * luma_grain_stride +
                luma_offset_x + (j << 1),
            cb_grain_block + (chroma_offset_y + i) * chroma_grain_stride +
                chroma_offset_x + j,
            cr_grain_block + (chroma_offset_y + i) * chroma_grain_stride +
                chroma_offset_x + j,
            luma_grain_stride, chroma_grain_stride,
            AOMMIN(chroma_subblock_size, height / 2 - y) - i,
            AOMMIN(chroma_subblock_size, width / 2 - x) - j, bit_depth);
      }

      if (overlap) {
        if (x) {
          y_line_buf[x << 1] = clamp(
              (27 * y_line_buf[x << 1] +
               17 * luma_grain_block[(luma_offset_y + luma_subblock_size) *
                                         luma_grain_stride +
                                     luma_offset_x] +
               16) >>
                  5,
              grain_min, grain_max);
          y_line_buf[(x << 1) + 1] = clamp(
              (17 * y_line_buf[(x << 1) + 1] +
               27 * luma_grain_block[(luma_offset_y + luma_subblock_size) *
                                         luma_grain_stride +
                                     luma_offset_x + 1] +
               16) >>
                  5,
              grain_min, grain_max);
          y_line_buf[luma_stride + (x << 1)] = clamp(
              (27 * y_line_buf[luma_stride + (x << 1)] +
               17 * luma_grain_block[(luma_offset_y + luma_subblock_size + 1) *
                                         luma_grain_stride +
                                     luma_offset_x] +
               16) >>
                  5,
              grain_min, grain_max);
          y_line_buf[luma_stride + (x << 1) + 1] = clamp(
              (17 * y_line_buf[luma_stride + (x << 1) + 1] +
               27 * luma_grain_block[(luma_offset_y + luma_subblock_size + 1) *
                                         luma_grain_stride +
                                     luma_offset_x + 1] +
               16) >>
                  5,
              grain_min, grain_max);

          cb_line_buf[x] = clamp(
              (23 * cb_col_buf[x] +
               22 * cb_grain_block[(chroma_offset_y + chroma_subblock_size) *
                                       chroma_grain_stride +
                                   chroma_offset_x] +
               16) >>
                  5,
              grain_min, grain_max);
          cr_line_buf[x] = clamp(
              (23 * cr_col_buf[x] +
               22 * cr_grain_block[(chroma_offset_y + chroma_subblock_size) *
                                       chroma_grain_stride +
                                   chroma_offset_x] +
               16) >>
                  5,
              grain_min, grain_max);
        }

        for (int m = overlap && x ? 1 : 0;
             m < AOMMIN(chroma_subblock_size + 1, width / 2 - x); m++) {
          y_line_buf[(x + m) << 1] =
              luma_grain_block[(luma_offset_y + luma_subblock_size) *
                                   luma_grain_stride +
                               luma_offset_x + (m << 1)];
          y_line_buf[((x + m) << 1) + 1] =
              luma_grain_block[(luma_offset_y + luma_subblock_size) *
                                   luma_grain_stride +
                               luma_offset_x + (m << 1) + 1];
          y_line_buf[luma_stride + ((x + m) << 1)] =
              luma_grain_block[(luma_offset_y + luma_subblock_size + 1) *
                                   luma_grain_stride +
                               luma_offset_x + (m << 1)];
          y_line_buf[luma_stride + ((x + m) << 1) + 1] =
              luma_grain_block[(luma_offset_y + luma_subblock_size + 1) *
                                   luma_grain_stride +
                               luma_offset_x + (m << 1) + 1];

          cb_line_buf[x + m] =
              cb_grain_block[(chroma_offset_y + chroma_subblock_size) *
                                 chroma_grain_stride +
                             chroma_offset_x + m];
          cr_line_buf[x + m] =
              cr_grain_block[(chroma_offset_y + chroma_subblock_size) *
                                 chroma_grain_stride +
                             chroma_offset_x + m];
        }

        for (int n = 0; n < AOMMIN(chroma_subblock_size + 1, height / 2 - y);
             n++) {
          y_col_buf[(n << 1) * 2] =
              luma_grain_block[(luma_offset_y + (n << 1)) * luma_grain_stride +
                               luma_offset_x + luma_subblock_size];
          y_col_buf[(n << 1) * 2 + 1] =
              luma_grain_block[(luma_offset_y + (n << 1)) * luma_grain_stride +
                               luma_offset_x + luma_subblock_size + 1];

          y_col_buf[((n << 1) + 1) * 2] =
              luma_grain_block[(luma_offset_y + (n << 1) + 1) *
                                   luma_grain_stride +
                               luma_offset_x + luma_subblock_size];
          y_col_buf[((n << 1) + 1) * 2 + 1] =
              luma_grain_block[(luma_offset_y + (n << 1) + 1) *
                                   luma_grain_stride +
                               luma_offset_x + luma_subblock_size + 1];

          cb_col_buf[n] =
              cb_grain_block[(chroma_offset_y + n) * chroma_grain_stride +
                             chroma_offset_x + chroma_subblock_size];
          cr_col_buf[n] =
              cr_grain_block[(chroma_offset_y + n) * chroma_grain_stride +
                             chroma_offset_x + chroma_subblock_size];
        }
      }
    }
  }

  dealloc_arrays(params, &pred_pos_luma, &pred_pos_chroma, &luma_grain_block,
                 &cb_grain_block, &cr_grain_block, &y_line_buf, &cb_line_buf,
                 &cr_line_buf, &y_col_buf, &cb_col_buf, &cr_col_buf);
}
