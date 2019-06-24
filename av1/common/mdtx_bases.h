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

#ifndef MDTX_BASES_H_
#define MDTX_BASES_H_

#ifdef __cplusplus
extern "C" {
#endif

static const int32_t klt4_inter[16] = {
  -436, -213, 1439, 2466, -710, 1313,  2138, -1260,
  2147, 1852, -79,  586,  1756, -1786, 1319, -613,
};

static const int32_t klt8_inter[64] = {
  -2,    24,    154,   253,  111,   811,   2450, 3165, 71,    211,   -106,
  -381,  1392,  3127,  1272, -1801, -181,  75,   478,  -1131, -3117, -263,
  1957,  -1272, 116,   428,  -1369, -3425, -185, 674,  -1160, 1069,  116,
  -1857, -3056, 81,    647,  -1162, 1350,  -614, 3126, 2235,  -366,  -45,
  520,   -968,  744,   -339, 2317,  -1713, -451, 1013, -1685, 1676,  -1166,
  488,   1255,  -2273, 2227, -1589, 1245,  -850, 497,  -174,
};

static const int32_t mdt4_mode0[16] = {
  -520, -337, 1393, 2462, -617, 1398,  2167, -1165,
  2160, 1769, -116, 763,  1753, -1787, 1318, -620,
};

static const int32_t mdt4_mode1[16] = {
  -237, -113, 1502, 2462, -671, 885,   2273, -1410,
  1896, 2163, -78,  330,  2071, -1707, 980,  -477,
};

static const int32_t mdt4_mode2[16] = {
  -266, -166, 1359, 2538, -675, 800,   2373, -1289,
  1728, 2293, -79,  373,  2208, -1570, 951,  -381,
};

static const int32_t mdt4_mode3[16] = {
  -338, -62,  1463, 2476, -541, 1471,  2081, -1266,
  2199, 1770, -346, 548,  1774, -1757, 1342, -594,
};

static const int32_t mdt4_mode4[16] = {
  -291, -278, 1232, 2590, -567, 1146,  2369, -1068,
  1551, 2315, -460, 641,  2362, -1279, 1023, -359,
};

static const int32_t mdt4_mode5[16] = {
  -304, -275, 1230, 2590, -671, 1015,  2385, -1104,
  1723, 2252, -221, 546,  2208, -1487, 1067, -405,
};

static const int32_t mdt4_mode6[16] = {
  -331, -425, 1099, 2625, -830, 929,   2428, -971,
  1887, 2125, 55,   558,  2008, -1682, 1133, -494,
};

static const int32_t mdt4_mode7[16] = {
  -287, -148, 1249, 2593, -912, 871,   2328, -1173,
  2133, 1931, 231,  235,  1710, -1969, 1163, -483,
};

static const int32_t mdt4_mode8[16] = {
  -388, -250, 1302, 2546, -785, 1267,  2213, -1127,
  2132, 1885, -50,  536,  1754, -1780, 1339, -592,
};

static const int32_t mdt4_mode9[16] = {
  1486, 2237, 491, -968, 1371, -708,  -2343, -718,
  1527, 19,   134, 2458, 1403, -1698, 1625,  -947,
};

static const int32_t mdt4_mode10[16] = {
  -532, -410, 1445, 2419, -474, 1586,  2109, -1096,
  2143, 1683, -299, 935,  1814, -1695, 1328, -682,
};

static const int32_t mdt4_mode11[16] = {
  -621, -562, 1297, 2450, -561, 1548,  2187, -945,
  2040, 1746, -257, 1054, 1878, -1620, 1363, -617,
};

static const int32_t mdt4_mode12[16] = {
  -289, -240, 1284, 2569, -650, 623,   2456, -1243,
  2182, 1864, 253,  294,  1767, -2113, 801,  -400,
};

static const int32_t mdt8_mode0[64] = {
  23,   -34,   101,   786,   1582,  2075,  2322,  1986, -68,   -364,  50,
  1773, 2593,  874,   -1490, -1947, -312,  -361,  1220, 2472,  101,   -2287,
  -508, 1858,  -518,  525,   2682,  1066,  -1503, 758,  1395,  -1768, 1839,
  3108, 1238,  -504,  658,   171,   -1024, 663,   2558, 269,   -1562, 1300,
  -202, -1178, 1798,  -1170, 2167,  -1745, 383,   681,  -1554, 1796,  -1487,
  756,  1337,  -1860, 1991,  -1833, 1548,  -1123, 686,  -285,
};

static const int32_t mdt8_mode1[64] = {
  70,   73,    -8,    164,  841,   1613,  2531, 2651, 134,   22,    -294,
  396,  2323,  2526,  -136, -2173, -253,  137,  345,  -1842, -2400, 1263,
  1800, -1607, -226,  807,  -661,  -3051, 523,  913,  -1903, 1266,  1047,
  -851, -3241, -847,  625,  -1266, 1205,  -539, 3567, 1908,  408,   215,
  -391, 159,   -135,  -2,   1193,  -2707, -623, 687,  -1473, 1669,  -1311,
  702,  1180,  -2100, 2253, -1629, 1388,  -979, 552,  -239,
};

static const int32_t mdt8_mode2[64] = {
  15,   96,    -46,   4,    624,   1507,  2583,  2726, 157,   107,   -297,
  -44,  1655,  2706,  636,  -2487, -101,  142,   433,  -946,  -2930, -22,
  2240, -1436, 201,   427,  160,   -1882, -1481, 2309, -2094, 1037,  932,
  461,  -2149, -2880, 1041, -1201, 498,   -100,  2746, 2818,  788,   803,
  -64,  -163,  30,    -26,  1719,  -1370, -2711, 1624, -1249, 597,   -180,
  117,  2311,  -2554, 1972, -874,  469,   -162,  117,  -17,
};

static const int32_t mdt8_mode3[64] = {
  65,   23,   61,    459,   1187,  1965,  2448,  2301, 59,    -300,  -188,
  1346, 2727, 1587,  -983,  -1979, -221,  -538,  697,  2524,  937,   -2084,
  -989, 1838, 691,   -569,  -2925, -1506, 1219,  -353, -1270, 1389,  1875,
  3108, 1142, -521,  657,   138,   -1065, 666,   2925, -160,  -1253, 1361,
  -458, -925, 1629,  -1027, 1835,  -2018, 1087,  148,  -1297, 1771,  -1634,
  804,  896,  -1522, 1908,  -1992, 1801,  -1385, 876,  -343,
};

static const int32_t mdt8_mode4[64] = {
  37,   -38,   -164,  50,   795,   1832,  2561,  2489, 69,    -309,  -603,
  708,  2595,  2161,  -407, -2060, 377,   737,   -676, -2737, -1607, 1424,
  1027, -1574, -679,  657,  2938,  1519,  -974,  587,  1363,  -1340, 1620,
  3138, 1082,  -500,  946,  261,   -1143, 787,   2854, 200,   -1223, 1500,
  -305, -1103, 1637,  -926, 2051,  -1816, 933,   337,  -1375, 1848,  -1579,
  701,  1094,  -1588, 1940, -1966, 1731,  -1297, 812,  -308,
};

static const int32_t mdt8_mode5[64] = {
  0,    -87,   208,   964,   1541,  1908,  2266,  2156, -93,   -284,  302,
  1850, 2429,  997,   -1250, -2172, -223,  -309,  1304, 2434,  11,    -2371,
  -769, 1672,  -444,  381,   2558,  887,   -1814, 611,  1671,  -1629, 302,
  2510, 2016,  -1110, 560,   892,   -1717, 1018,  2897, 2006,  -769,  718,
  341,  -1100, 1201,  -697,  2502,  -1581, 204,   733,  -1528, 1698,  -1362,
  608,  1334,  -1905, 1924,  -1852, 1566,  -1122, 675,  -269,
};

static const int32_t mdt8_mode6[64] = {
  42,   51,    -18,   -275, -329,  917,   2613, 2986, 132,   -34,   -662,
  -635, 1529,  3059,  1058, -1761, 77,    436,  143,  -2083, -2998, -102,
  1075, -1439, 708,   -287, -2766, -2292, 434,  -229, -1369, 1083,  853,
  2791, 1782,  -1012, 242,  1050,  -1474, 853,  2834, 1560,  -1070, 1012,
  433,  -1243, 1363,  -744, 2224,  -1597, 101,  1221, -1745, 1730,  -1236,
  467,  1595,  -1930, 2083, -1780, 1366,  -936, 489,  -131,
};

static const int32_t mdt8_mode7[64] = {
  97,   -2,    -43,   126,  570,   1278,  2411, 2996, -17,   -177,  -220,
  854,  2483,  2503,  337,  -1850, 48,    353,  -240, -2133, -1976, 1201,
  1991, -1653, 273,   182,  -1705, -2562, 709,  1129, -2048, 1106,  661,
  -847, -2780, -152,  1132, -2013, 1476,  -600, 2925, 2716,  529,   -61,
  581,  -423,  160,   -141, 2494,  -1919, -796, 1306, -1639, 1195,  -631,
  161,  1217,  -2193, 2263, -1785, 1256,  -707, 300,  -112,
};

static const int32_t mdt8_mode8[64] = {
  39,   -19,  89,    602,   1342,  2040,  2413,  2150, -23,   -312,  -69,
  1493, 2663, 1347,  -1175, -2040, -265,  -489,  829,  2518,  709,   -2204,
  -830, 1841, -666,  352,   2724,  1449,  -1462, 386,  1418,  -1547, 1890,
  3093, 1377, -397,  576,   174,   -942,  579,   2888, -157,  -1439, 1295,
  -357, -953, 1603,  -1030, 1890,  -2076, 998,   204,  -1273, 1723,  -1612,
  818,  882,  -1553, 1923,  -1967, 1780,  -1374, 907,  -384,
};

static const int32_t mdt8_mode9[64] = {
  -97,   108,   949,   1875,  2257,  2002,  1498,  997,  -450,  6,     1953,
  2483,  79,    -2009, -1551, -384,  18,    1498,  2223, -449,  -2244, -57,
  1765,  1112,  993,   2392,  1053,  -735,  734,   1270, -1192, -2202, 2291,
  1746,  -1020, 124,   922,   -1412, -389,  2105,  2309, -554,  -784,  1675,
  -1047, -217,  1721,  -1901, 1999,  -1831, 1140,  -117, -966,  1652,  -1800,
  1102,  1003,  -1446, 1749,  -1876, 1795,  -1531, 1113, -543,
};

static const int32_t mdt8_mode10[64] = {
  -56,  189,  1247,  2233,  2265,  1687,  1229,  843,  -376,  342,   2181,
  2041, -742, -2303, -1298, -237,  90,    1474,  1906, -950,  -2141, 597,
  2070, 913,  752,   2259,  1012,  -834,  811,   1241, -1494, -2227, 2260,
  1939, -936, -54,   836,   -1317, -331,  2118,  2470, -436,  -751,  1527,
  -940, -206, 1582,  -2039, 1976,  -1820, 1124,  -90,  -966,  1621,  -1842,
  1154, 970,  -1444, 1762,  -1920, 1827,  -1524, 1029, -492,
};

static const int32_t mdt8_mode11[64] = {
  -66,  42,   990,   1975,  2223,  1926,  1569,  880,  -267,  277,   2076,
  2366, -13,  -2024, -1564, -424,  -149,  885,   1943, -370,  -2409, -188,
  2100, 1347, 324,   2034,  1590,  -1009, 91,    1737, -728,  -2331, 1862,
  2471, -395, -559,  1110,  -814,  -765,  2063,  2771, 92,    -933,  1344,
  -613, -620, 1617,  -1743, 2091,  -1743, 871,   190,  -1156, 1736,  -1750,
  1077, 1027, -1621, 1844,  -1993, 1754,  -1363, 815,  -352,
};

static const int32_t mdt8_mode12[64] = {
  26,    44,    -27,   353,  1184,  1594,  2373,  2659,  52,    1,     -187,
  847,   2598,  2061,  -325, -2217, 96,    -80,   -230,  1491,  2008,  -1579,
  -2176, 1795,  -146,  345,  -557,  -2685, 78,    1825,  -2129, 1118,  -510,
  -632,  -2233, -1931, 1592, -1756, 1237,  -510,  -1631, -2978, -1432, 1057,
  -994,  904,   -470,  231,  2775,  530,   -2641, 887,   -887,  475,   -133,
  48,    2474,  -2665, 1539, -934,  527,   -189,  34,    8,
};

static INLINE const int32_t *mdt4_arr(int mode) {
  if (mode >= INTER_MODE_START && mode < INTER_MODE_END) {
    return klt4_inter;
  } else {
    switch (mode) {
      case 0: return mdt4_mode0;
      case 1: return mdt4_mode1;
      case 2: return mdt4_mode2;
      case 3: return mdt4_mode3;
      case 4: return mdt4_mode4;
      case 5: return mdt4_mode5;
      case 6: return mdt4_mode6;
      case 7: return mdt4_mode7;
      case 8: return mdt4_mode8;
      case 9: return mdt4_mode9;
      case 10: return mdt4_mode10;
      case 11: return mdt4_mode11;
      case 12: return mdt4_mode12;
      default: assert(0); return 0;
    }
  }
}

static INLINE const int32_t *mdt8_arr(int mode) {
  if (mode >= INTER_MODE_START && mode < INTER_MODE_END) {
    return klt8_inter;
  } else {
    switch (mode) {
      case 0: return mdt8_mode0;
      case 1: return mdt8_mode1;
      case 2: return mdt8_mode2;
      case 3: return mdt8_mode3;
      case 4: return mdt8_mode4;
      case 5: return mdt8_mode5;
      case 6: return mdt8_mode6;
      case 7: return mdt8_mode7;
      case 8: return mdt8_mode8;
      case 9: return mdt8_mode9;
      case 10: return mdt8_mode10;
      case 11: return mdt8_mode11;
      case 12: return mdt8_mode12;
      default: assert(0); return 0;
    }
  }
}

#ifdef __cplusplus
}
#endif

#endif  // MDTX_BASES_H_
