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

//  Test and time AOM intra-predictor functions

#include <stdio.h>
#include <string>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "./aom_dsp_rtcd.h"
#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/md5_helper.h"
#include "aom/aom_integer.h"
#include "aom_ports/mem.h"
#include "aom_ports/aom_timer.h"
#include "av1/common/common_data.h"

// -----------------------------------------------------------------------------

namespace {

// Note:
// APPLY_UNIT_TESTS
// 1: Do unit tests
// 0: Generate MD5 array as required
#define APPLY_UNIT_TESTS 1

typedef void (*AvxPredFunc)(uint8_t *dst, ptrdiff_t y_stride,
                            const uint8_t *above, const uint8_t *left);

const int kBPS = 64;
const int kTotalPixels = kBPS * kBPS;

const int kNumAv1IntraFuncs = INTRA_MODES + 3;  // 4 DC predictor variants.

#if APPLY_UNIT_TESTS
const char *kAv1IntraPredNames[kNumAv1IntraFuncs] = {
  "DC_PRED",    "DC_LEFT_PRED", "DC_TOP_PRED",   "DC_128_PRED",
  "V_PRED",     "H_PRED",       "D45_PRED",      "D135_PRED",
  "D113_PRED",  "D157_PRED",    "D203_PRED",     "D67_PRED",
  "PAETH_PRED", "SMOOTH_PRED",  "SMOOTH_V_PRED", "SMOOTH_H_PRED",
};
#endif  // APPLY_UNIT_TESTS

template <typename Pixel>
struct IntraPredTestMem {
  void Init(int block_width, int block_height, int bd) {
    ASSERT_LE(block_width, kBPS);
    ASSERT_LE(block_height, kBPS);
    // Note: for blocks having width <= 32 and height <= 32, we generate 32x32
    // random pixels as before to avoid having to recalculate all hashes again.
    const int block_size_upto_32 = (block_width <= 32) && (block_height <= 32);
    stride = block_size_upto_32 ? 32 : kBPS;
    num_pixels = stride * stride;
    libaom_test::ACMRandom rnd(libaom_test::ACMRandom::DeterministicSeed());
    above = above_mem + 16;
    const int mask = (1 << bd) - 1;
    for (int i = 0; i < num_pixels; ++i) ref_src[i] = rnd.Rand16() & mask;
    for (int i = 0; i < stride; ++i) left[i] = rnd.Rand16() & mask;
    for (int i = -1; i < stride; ++i) above[i] = rnd.Rand16() & mask;

    for (int i = stride; i < 2 * stride; ++i) {
      left[i] = rnd.Rand16() & mask;
      above[i] = rnd.Rand16() & mask;
    }
  }

  DECLARE_ALIGNED(16, Pixel, src[kTotalPixels]);
  DECLARE_ALIGNED(16, Pixel, ref_src[kTotalPixels]);
  DECLARE_ALIGNED(16, Pixel, left[2 * kBPS]);
  Pixel *above;
  int stride;
  int num_pixels;

 private:
  DECLARE_ALIGNED(16, Pixel, above_mem[2 * kBPS + 16]);
};

// -----------------------------------------------------------------------------
// Low Bittdepth

typedef IntraPredTestMem<uint8_t> Av1IntraPredTestMem;

static const char *const kTxSizeStrings[TX_SIZES_ALL] = {
  "4X4",  "8X8",  "16X16", "32X32", "64X64", "4X8",   "8X4",
  "8X16", "16X8", "16X32", "32X16", "32X64", "64X32", "4X16",
  "16X4", "8X32", "32X8",  "16X64", "64X16",
};

void CheckMd5Signature(TX_SIZE tx_size, bool is_hbd,
                       const char *const signatures[], const void *data,
                       size_t data_size, int elapsed_time, int idx) {
  const std::string hbd_str = is_hbd ? "Hbd " : "";
  const std::string name_str = hbd_str + "Intra" + kTxSizeStrings[tx_size];
  libaom_test::MD5 md5;
  md5.Add(reinterpret_cast<const uint8_t *>(data), data_size);
#if APPLY_UNIT_TESTS
  printf("Mode %s[%13s]: %5d ms     MD5: %s\n", name_str.c_str(),
         kAv1IntraPredNames[idx], elapsed_time, md5.Get());
  EXPECT_STREQ(signatures[idx], md5.Get());
#else
  (void)signatures;
  (void)elapsed_time;
  (void)idx;
  printf("\"%s\",\n", md5.Get());
#endif
}

void TestIntraPred(TX_SIZE tx_size, AvxPredFunc const *pred_funcs,
                   const char *const signatures[]) {
  const int block_width = tx_size_wide[tx_size];
  const int block_height = tx_size_high[tx_size];
  const int num_pixels_per_test =
      block_width * block_height * kNumAv1IntraFuncs;
  const int kNumTests = static_cast<int>(2.e10 / num_pixels_per_test);
  Av1IntraPredTestMem intra_pred_test_mem;
  intra_pred_test_mem.Init(block_width, block_height, 8);

  for (int k = 0; k < kNumAv1IntraFuncs; ++k) {
    if (pred_funcs[k] == NULL) continue;
    memcpy(intra_pred_test_mem.src, intra_pred_test_mem.ref_src,
           sizeof(intra_pred_test_mem.src));
    aom_usec_timer timer;
    aom_usec_timer_start(&timer);
    for (int num_tests = 0; num_tests < kNumTests; ++num_tests) {
      pred_funcs[k](intra_pred_test_mem.src, intra_pred_test_mem.stride,
                    intra_pred_test_mem.above, intra_pred_test_mem.left);
    }
    libaom_test::ClearSystemState();
    aom_usec_timer_mark(&timer);
    const int elapsed_time =
        static_cast<int>(aom_usec_timer_elapsed(&timer) / 1000);
    CheckMd5Signature(
        tx_size, false, signatures, intra_pred_test_mem.src,
        intra_pred_test_mem.num_pixels * sizeof(*intra_pred_test_mem.src),
        elapsed_time, k);
  }
}

static const char *const kSignatures[TX_SIZES_ALL][kNumAv1IntraFuncs] = {
  {
      // 4X4
      "e7ed7353c3383fff942e500e9bfe82fe",
      "2a4a26fcc6ce005eadc08354d196c8a9",
      "269d92eff86f315d9c38fe7640d85b15",
      "ae2960eea9f71ee3dabe08b282ec1773",
      "6c1abcc44e90148998b51acd11144e9c",
      "f7bb3186e1ef8a2b326037ff898cad8e",
      "87e72798518d62e84bcc77dcb17d0f3b",
      "141624072a4a56773f68fadbdd07c4a7",
      "7be49b08687a5f24df3a2c612fca3876",
      "459bb5d9fd5b238348179c9a22108cd6",
      "3d98810f418a9de92acfe2c68909c61c",
      "6310eecda3cc9496987ca10186255558",
      "59fc0e923a08cfac0a493fb38988e2bb",
      "9ff8bb37d9c830e6ab8ecb0c435d3c91",
      "de6937fca02354f2874dbc5dbec5d5b3",
      "723cf948137f7d8c7860d814e55ae67d",
  },
  {
      // 8X8
      "d8bbae5d6547cfc17e4f5f44c8730e88",
      "373bab6d931868d41a601d9d88ce9ac3",
      "6fdd5ff4ff79656c14747598ca9e3706",
      "d9661c2811d6a73674f40ffb2b841847",
      "7c722d10b19ccff0b8c171868e747385",
      "f81dd986eb2b50f750d3a7da716b7e27",
      "e0b1292448f3350bf1c92ca283ca872a",
      "0e3523f9cab2142dd37fd07ec0760bce",
      "79ac4efe907f0a0f1885d43066cfedee",
      "19ecf2432ac305057de3b6578474eec6",
      "7ae38292cbe47b4aa0807c3bd5a543df",
      "d0ecffec1bb01f4b61ab5738164695c4",
      "064404361748dd111a890a1470d7f0ea",
      "dc29b7e1f78cc8e7525d5ea4c0ab9b78",
      "97111eb1bc26bade6272015df829f1ae",
      "d19a8a73cc46b807f2c5e817576cc1e1",
  },
  {
      // 16X16
      "50971c07ce26977d30298538fffec619",
      "527a6b9e0dc5b21b98cf276305432bef",
      "7eff2868f80ebc2c43a4f367281d80f7",
      "67cd60512b54964ef6aff1bd4816d922",
      "48371c87dc95c08a33b2048f89cf6468",
      "b0acf2872ee411d7530af6d2625a7084",
      "31d901ab2289d1e61e704e40240382a7",
      "dae208f3dca583529cff49b73f7c4183",
      "7af66a2f4c8e0b4908e40f047e60c47c",
      "125e3ab6ab9bc961f183ec366a7afa88",
      "ff230677e800977757d14b85a9eba404",
      "eb42dc39140515dd4f3ab1afe6c3e71b",
      "93d6b5352b571805ab16a55e1bbed86a",
      "03764e4c0aebbc180e4e2c68fb06df2b",
      "bb6c74c9076c9f266ab11fb57060d8e6",
      "0c5162bc28489756ddb847b5678e6f07",
  },
  {
      // 32X32
      "a0a618c900e65ae521ccc8af789729f2",
      "985aaa7c72b4a6c2fb431d32100cf13a",
      "10662d09febc3ca13ee4e700120daeb5",
      "b3b01379ba08916ef6b1b35f7d9ad51c",
      "9f4261755795af97e34679c333ec7004",
      "bc2c9da91ad97ef0d1610fb0a9041657",
      "f524b1a7e31c7bb9bfb2487fac3e16d8",
      "4039bb7da0f6860090d3c57b5c85468f",
      "b29fff7b61804e68383e3a609b33da58",
      "e1aa5e49067fd8dba66c2eb8d07b7a89",
      "db217e7891581cf93895ef5974bebb21",
      "beb6cdc52b52c8976b4d2407ec8d2313",
      "ef1653982b69e1f64bee3759f3e1ec45",
      "1a51a675deba2c83282142eb48d3dc3d",
      "866c224746dc260cda861a7b1b383fb3",
      "cea23799fc3526e1b6a6ff02b42b82af",
  },
  {
      // 64X64
      "6e1094fa7b50bc813aa2ba29f5df8755",
      "afe020786b83b793c2bbd9468097ff6e",
      "be91585259bc37bf4dc1651936e90b3e",
      "a1650dbcd56e10288c3e269eca37967d",
      "9e5c34f3797e0cdd3cd9d4c05b0d8950",
      "bc87be7ac899cc6a28f399d7516c49fe",
      "4212fc4c10b35d8007f3107c54baf288",
      "cb6fa17f4c4075e4604f2d0ce210cb0a",
      "0febfec10abe189787c9eca333ec3499",
      "fde42ff7e6802b7e539843fe9a1ef4f6",
      "3c12c9ef5eab1d27302dbf4c18400650",
      "cfca1c2f4fb65c6e5db9e09b7ee83f0f",
      "9811fd0d2dd515f06122f5d1bd18b784",
      "3c140e466f2c2c0d9cb7d2157ab8dc27",
      "9543de76c925a8f6adc884cc7f98dc91",
      "df1df0376cc944afe7e74e94f53e575a",
  },
  {
      // 4X8
      "d9fbebdc85f71ab1e18461b2db4a2adc",
      "5ccb2a68284bc9714d94b8a06ccadbb2",
      "735d059abc2744f3ff3f9590f7191b37",
      "d9fbebdc85f71ab1e18461b2db4a2adc",
      "6819497c44cd0ace120add83672996ee",
      "7e3244f5a2d3edf81c7e962a842b97f9",
      "3fa52ee9acf5a25594cac684be263f32",
      "c18dd23d57def4df4c6147c572dfc827",
      "d007fbf7e43cb8f49702daa20f0c9153",
      "5c0226c44c5df285728296b80cc6de4b",
      "b55d7b558bebc8c2042dfac58b3c4688",
      "6549362baa389b8faa2d954926b64e2f",
      "809350f164cd4d1650850bb0f59c3260",
      "1b60a394331eeab6927a6f8aaff57040",
      "5307de1bd7329ba6b281d2c1b0b457f9",
      "24c58a8138339846d95568efb91751db",
  },
  {
      // 8X4
      "23f9fc11344426c9bee2e06d57dfd628",
      "2d71a26d1bae1fb34734de7b42fc5eb7",
      "5af9c1b2fd9d5721fad67b67b3f7c816",
      "00d71b17be662753813d515f197d145e",
      "bef10ec984427e28f4390f43809d10af",
      "77773cdfb7ed6bc882ab202a64b0a470",
      "cba356970f6b9a1b6024e1dbe4a66f9b",
      "c58c21efc804242848e6f29a93a7984d",
      "dc92cc45a51c7a397506cab19f74e66d",
      "391f6a12224f81a3719ea09a2cf7a5ad",
      "b74b8b11f7eb2bbf723b25f381104ca9",
      "2234aaa06ca245624211cf53a0261017",
      "2cc48bd66d6b0121b5221d52ccd732af",
      "b302155e1c9eeeafe2ba2bf68e807a46",
      "561bc8d0e76d5041ebd5168fc6a115e1",
      "81d0113fb1d0a9a24ffd6f1987b77948",
  },
  {
      // 8X16
      "c849de88b24f773dfcdd1d48d1209796",
      "6cb807c1897b94866a0f3d3c56ed8695",
      "d56db05a8ac7981762f5b877f486c4ef",
      "b4bc01eb6e59a40922ad17715cafb04b",
      "09d178439534f4062ae687c351f66d64",
      "644501399cf73080ac606e5cef7ca09b",
      "0e8e968fa177204d7e73d7e04ce69ebb",
      "1d25f9287fdf7ba48a5105f1529b7e75",
      "02cacccf3752451763a6a6e2e784494f",
      "6044a1416d53e324ddc012d2e7763339",
      "57ac6e8f3ab5e943c9280043eeb174b8",
      "d51b9d65471194d9caebc7d67e75ef10",
      "278076495180e17c065a95ab7278539a",
      "9dd7f324816f242be408ffeb0c673732",
      "f520c4a20acfa0bea1d253c6f0f040fd",
      "85f38df809df2c2d7c8b4a157a65cd44",
  },
  {
      // 16X8
      "b4cbdbdf10ce13300b4063a3daf99e04",
      "3731e1e6202064a9d0604d7c293ecee4",
      "6c856188c4256a06452f0d5d70cac436",
      "1f2192b4c8c497589484ea7bf9c944e8",
      "84011bd4b7f565119d06787840e333a0",
      "0e48949f7a6aa36f0d76b5d01f91124a",
      "58114c06f6b9d8285e5020c7afd834ab",
      "e37afe84a8b3c5e0f048d4652ecbe09e",
      "c216348473fb029b45f8fb4f2862a7bd",
      "0b7385155dcef742cc456d5741ae93a3",
      "d55fadb221f0ea20266e57cd413e7b94",
      "9bd6eb226c7e169b8d53cf70aea98b3a",
      "60eff8064634b6c73b10681356baeee9",
      "1559aeb081a9c0c71111d6093c2ff9fd",
      "c15479b739713773e5cabb748451987b",
      "72e33ec12c9b67aea26d8d005fb82de2",
  },
  {
      // 16X32
      "abe5233d189cdbf79424721571bbaa7b",
      "282759f81e3cfb2e2d396fe406b72a8b",
      "e2224926c264f6f174cbc3167a233168",
      "6814e85c2b33f8c9415d62e80394b47b",
      "99cbbb60459c08a3061d72c4e4f6276a",
      "1d1567d40b8e816f8c1f71e576fe0f87",
      "5e989f9c748a0d2cd8c4ebf9d3fe1278",
      "7135a2f419452a3a192a35156f68b019",
      "06e10af5a726d2c81b8f8c708204f9fb",
      "c0882f0e7ba1ffa0aeef6d5c751df6de",
      "8477429e17d39a423f30e2082f651549",
      "ba35068a30c2d1d10901e4bfabd02a11",
      "36fdd371b624a075814d497c4832ec85",
      "8ab8da61b727442b6ff692b40d0df018",
      "e35a10ad7fdf2327e821504a90f6a6eb",
      "1f7211e727dc1de7d6a55d082fbdd821",
  },
  {
      // 32X16
      "d1aeb8d5fdcfd3307922af01a798a4dc",
      "b0bcb514ebfbee065faea9d34c12ae75",
      "d6a18c63b4e909871c0137ca652fad23",
      "fd047f2fc1b8ffb95d0eeef3e8796a45",
      "645ab60779ea348fd93c81561c31bab9",
      "4409633c9db8dff41ade4292a3a56e7f",
      "b9b2935b2287a9a461ac5c11251ac706",
      "43b05f808c0ac4fe8accd84d293b0488",
      "1d2cb43872d20c205ffb185102bcd22a",
      "2c1551b5e99592fd21053b5d14e397d9",
      "cd499ef0dd41e2e38d5dac3319dfdd97",
      "cd2610426637003f3b5d3984cb3320d5",
      "5e36a11e069b31c2a739f3a9c7b37c24",
      "e83b9483d702cfae496991c3c7fa92c0",
      "12f6ddf98c7f30a277307f1ea935b030",
      "354321d6c32bbdb0739e4fa2acbf41e1",
  },
  {
      // 32X64
      "0ce332b343934b34cd4417725faa85cb",
      "4e2a2cfd8f56f15939bdfc753145b303",
      "0f46d124ba9f48cdd5d5290acf786d6d",
      "e1e8ed803236367821981500a3d9eebe",
      "1d2f8e48e3adb7c448be05d9f66f4954",
      "9fb2e176636a5689b26f73ca73fcc512",
      "cfbb498ee3bfb6847a128d0fdfcdbefc",
      "27513c8e559c497e8ccc3c634ec256d1",
      "00e8b7f593174c18061b2a626b403192",
      "378e0ca2459ff4f97f05e0c75c6bacec",
      "1bab785e7bee41fff4cc5b85303d9b21",
      "f6a06ea228676890b98ac430fb8162cb",
      "e720ebccae7e25e36f23da53ae5b5d6a",
      "86fe4364734169aaa4520d799890d530",
      "b1870290764bb1b100d1974e2bd70f1d",
      "ce5b238e19d85ef69d85badfab4e63ae",
  },
  {
      // 64X32
      "a6c5aeb722615089efbca80b02951ceb",
      "538424b24bd0830f21788e7238ca762f",
      "80c15b303235f9bc2259027bb92dfdc4",
      "e48e1ac15e97191a8fda08d62fff343e",
      "12604b37875533665078405ef4582e35",
      "0048afa17bd3e1632d68b96048836530",
      "ef0ae3280e9510dc13cc7c3b773d919a",
      "6af8e4812cb6c44f60ae3c12d50eefad",
      "265875770abece235f58c6642d87d050",
      "d4f2822a5a6458c98bba826ad9ed0255",
      "5b9d388dfa74e4eaf556b205d6aeb810",
      "8e35512af78ba9454d8e603ef855c043",
      "07a0cfcb56a5eed50c4bd6c26814336b",
      "529d8a070de5bc6531fa3ee8f450c233",
      "33c50a11c7d78f72434064f634305e95",
      "e0ef7f0559c1a50ec5a8c12011b962f7",
  },
  {
      // 4X16
      "750491056568eb8fe15387b86bdf06b8",
      "3a52dae9f599f08cfb3bd1b910dc0e11",
      "af79f71e3e03dbeca44e2e13561f70c7",
      "ca7dfd7624afc0c06fb5552f44398535",
      "b591af115444bf43140c29c269f68fb2",
      "483d942ae36e69e62f31eb215331416f",
      "c07af04bd17e41b27f065f437930e025",
      "1a53c3af65f4cae813564eb77bc3959d",
      "ced23810f1f39af2f83c3e69619e0244",
      "9b4f4eef2a18ee8473406b5967a1ea4c",
      "7821f3be5787b7e921981ce2b4254d66",
      "8dd781f95acfa22a9b1b9e876980076f",
      "f14b58525e81870bc5d95c7ac71a347f",
      "371208bb4027d9badb04095d1590bbc4",
      "c7049c21b2924d70c7c12784d6b6b796",
      "7d87233f4b5b0f12086045e5d7b2d4c2",
  },
  {
      // 16X4
      "7c6e325a65e77e732b3adbe237e045e4",
      "24478f93ffcec47852e004d0fe948464",
      "258d042c67d4ba3ecfa667f0adc9aebf",
      "b2cd21d06959f159a1f3c4d9768ee7fb",
      "b4e1f38157bf8410e7c3da02f687a343",
      "869e703729eb0fc0711c254944ff5d5a",
      "76209c66c080e9832e6a6c25d67354ab",
      "0dab0b6b2d8a3e9567147d21700cb79e",
      "aadbe4b6a322e966287d93a4348ae2d5",
      "51afd7106766bd5dc6af8d0170e6527c",
      "1bfa41b2e50842850bab5a19460525e9",
      "9c3ba328e3c207983f0bc44ceed73926",
      "9638dd77105a640b146a8201ea7a0801",
      "919d932c6af8a1cc7486e8ce996dd487",
      "e1c9be493b6714c7ae48f30044c43140",
      "bf0fe3889d654b2f6eb98c8fc751f9e4",
  },
  {
      // 8X32
      "8dfac4319fe0bd40013ffb3102da8c72",
      "feb46b6dc4e2ca0a09533bfc51d4dcb0",
      "850837ec714c37262216527aaf4cbbe9",
      "4603c7800fb08361f163daca876e8bda",
      "1ff95e7d2debc27b05806fb25abfd624",
      "d81b9a51a062b23ca7823804cb7bec22",
      "08ed89d15a510e01e8f4bc7f2a72c936",
      "79fe2296d52cb614da518f20ab29b0bb",
      "fb799a6c6953818b9a0d24b947c998db",
      "00bfe3b6481cc125c95050124f41f312",
      "0ef7b0fbf647e5257040c8cc8e56e021",
      "cf130a14a52ddf5c02f264edd36d8f40",
      "f1d8978158766f46335203608cb807e7",
      "f3527096256258c0878d644a9d7d53ca",
      "cbde98ac8b009953eb112807ad2ea29e",
      "654fb1153415747feae599f538122af5",
  },
  {
      // 32X8
      "3d4ee16fab374357474f60b845327bc7",
      "bc17c5059473a476df4e85f56395ad55",
      "3d4ee16fab374357474f60b845327bc7",
      "c14b8db34dc2355b84e3735c9ba16c7f",
      "a71d25b5d47a92a8b9223c98f18458ee",
      "6c1cfe2b1893f4576a80675687cb6426",
      "fb4046f3fa67ddc7dcf14d71ad478251",
      "650c8e037ec0d93f97fa9bdb47f94086",
      "f609ec7f3f9607d7a587559f2c0a0fc5",
      "553f8f47b0fe69c402569b86ff72528a",
      "92f1ccb5356e3764f8038b95c1eb63a9",
      "897f466630b32dc9f4695a9d5ff1abe2",
      "92d11bbef8b85bb48d799bb055de3514",
      "bcf81d1db8ae5cc03360467f44f498ec",
      "79f8c564163555592e808e145eaf5c60",
      "46fff139cef2ef773938bcc8b0e5abb8",
  },
  {
      // 16X64
      "3b2a053ee8b05a8ac35ad23b0422a151",
      "12b0c69595328c465e0b25e0c9e3e9fc",
      "f77c544ac8035e01920deae40cee7b07",
      "727797ef15ccd8d325476fe8f12006a3",
      "f3be77c0fe67eb5d9d515e92bec21eb7",
      "f1ece6409e01e9dd98b800d49628247d",
      "cd3d9177d365a036227ff3059ada4f40",
      "3bee921a7745eac39e967d0b6a718e60",
      "a270f027dc4eac72d36298ecdf558a11",
      "f2793274301a1ba435e332dee0379813",
      "5a2682f1495e80707fe072a8e84e500a",
      "5d2fecd1a51bf12799bc17dba7c67978",
      "efd2ec9bfbbd4fd1f6604ea369df1894",
      "ec703de918422b9e03197ba0ed60a199",
      "739418efb89c07f700895deaa5d0b3e3",
      "9943ae1bbeeebfe1d3a92dc39e049d63",
  },
  {
      // 64X16
      "821b76b1494d4f84d20817840f719a1a",
      "69e462c3338a9aaf993c3f7cfbc15649",
      "516d8f6eb054d74d150e7b444185b6b9",
      "de1b736e9d99129609d6ef3a491507a0",
      "fd9b4276e7affe1e0e4ce4f428058994",
      "cd82fd361a4767ac29a9f406b480b8f3",
      "fe0efad08c4a63db47c36bbd3c70c2ae",
      "8151af4ec36b24b139aaad5a07cad024",
      "54b13f16b5dcd5a994b1ebdf08dc3e78",
      "49b54062c995de2cc3dee9028c653cdb",
      "13f6f9e45995c7d994909969bcc50bf9",
      "6f6f36a928b478ea2dc7c897e8983602",
      "2792c2f810157a4a6cb13c28529ff779",
      "1220442d90c4255ba0969d28b91e93a6",
      "c7253e10b45f7f67dfee3256c9b94825",
      "879792198071c7e0b50b9b5010d8c18f",
  },
};

}  // namespace

// Defines a test case for |arch| (e.g., C, SSE2, ...) passing the predictors
// to TestIntraPred. The test name is 'arch.TestIntraPred.tx_size', e.g.,
// C.TestIntraPred.0
#define INTRA_PRED_TEST(arch, tx_size, dc, dc_left, dc_top, dc_128, v, h,   \
                        d45e, d135, d117, d153, d207e, d63e, paeth, smooth, \
                        smooth_v, smooth_h)                                 \
  TEST(arch, DISABLED_##TestIntraPred_##tx_size) {                          \
    static const AvxPredFunc aom_intra_pred[] = {                           \
      dc,   dc_left, dc_top, dc_128, v,     h,      d45e,     d135,         \
      d117, d153,    d207e,  d63e,   paeth, smooth, smooth_v, smooth_h      \
    };                                                                      \
    TestIntraPred(tx_size, aom_intra_pred, kSignatures[tx_size]);           \
  }

// -----------------------------------------------------------------------------
// 4x4, 4x8, 4x16

INTRA_PRED_TEST(C_1, TX_4X4, aom_dc_predictor_4x4_c,
                aom_dc_left_predictor_4x4_c, aom_dc_top_predictor_4x4_c,
                aom_dc_128_predictor_4x4_c, aom_v_predictor_4x4_c,
                aom_h_predictor_4x4_c, aom_d45e_predictor_4x4_c,
                aom_d135_predictor_4x4_c, aom_d117_predictor_4x4_c,
                aom_d153_predictor_4x4_c, aom_d207e_predictor_4x4_c,
                aom_d63e_predictor_4x4_c, aom_paeth_predictor_4x4_c,
                aom_smooth_predictor_4x4_c, aom_smooth_v_predictor_4x4_c,
                aom_smooth_h_predictor_4x4_c)

INTRA_PRED_TEST(C_2, TX_4X8, aom_dc_predictor_4x8_c,
                aom_dc_left_predictor_4x8_c, aom_dc_top_predictor_4x8_c,
                aom_dc_128_predictor_4x8_c, aom_v_predictor_4x8_c,
                aom_h_predictor_4x8_c, aom_d45e_predictor_4x8_c,
                aom_d135_predictor_4x8_c, aom_d117_predictor_4x8_c,
                aom_d153_predictor_4x8_c, aom_d207e_predictor_4x8_c,
                aom_d63e_predictor_4x8_c, aom_paeth_predictor_4x8_c,
                aom_smooth_predictor_4x8_c, aom_smooth_v_predictor_4x8_c,
                aom_smooth_h_predictor_4x8_c)

INTRA_PRED_TEST(C_3, TX_4X16, aom_dc_predictor_4x16_c,
                aom_dc_left_predictor_4x16_c, aom_dc_top_predictor_4x16_c,
                aom_dc_128_predictor_4x16_c, aom_v_predictor_4x16_c,
                aom_h_predictor_4x16_c, aom_d45e_predictor_4x16_c,
                aom_d135_predictor_4x16_c, aom_d117_predictor_4x16_c,
                aom_d153_predictor_4x16_c, aom_d207e_predictor_4x16_c,
                aom_d63e_predictor_4x16_c, aom_paeth_predictor_4x16_c,
                aom_smooth_predictor_4x16_c, aom_smooth_v_predictor_4x16_c,
                aom_smooth_h_predictor_4x16_c)

#if HAVE_SSE2
INTRA_PRED_TEST(SSE2_1, TX_4X4, aom_dc_predictor_4x4_sse2,
                aom_dc_left_predictor_4x4_sse2, aom_dc_top_predictor_4x4_sse2,
                aom_dc_128_predictor_4x4_sse2, aom_v_predictor_4x4_sse2,
                aom_h_predictor_4x4_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_4x4_sse2, aom_d63e_predictor_4x4_sse2, NULL,
                NULL, NULL, NULL)
INTRA_PRED_TEST(SSE2_2, TX_4X8, aom_dc_predictor_4x8_sse2,
                aom_dc_left_predictor_4x8_sse2, aom_dc_top_predictor_4x8_sse2,
                aom_dc_128_predictor_4x8_sse2, aom_v_predictor_4x8_sse2,
                aom_h_predictor_4x8_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_4x8_sse2, aom_d63e_predictor_4x8_sse2, NULL,
                NULL, NULL, NULL)
#endif  // HAVE_SSE2

#if HAVE_SSSE3
INTRA_PRED_TEST(SSSE3_1, TX_4X4, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_4x4_ssse3, NULL, NULL,
                aom_d153_predictor_4x4_ssse3, NULL,
                aom_d63e_predictor_4x4_ssse3, aom_paeth_predictor_4x4_ssse3,
                aom_smooth_predictor_4x4_ssse3, NULL, NULL)
INTRA_PRED_TEST(SSSE3_2, TX_4X8, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_4x8_ssse3, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_4x8_ssse3, aom_smooth_predictor_4x8_ssse3,
                NULL, NULL)
#endif  // HAVE_SSSE3

#if HAVE_DSPR2
INTRA_PRED_TEST(DSPR2, TX_4X4, aom_dc_predictor_4x4_dspr2, NULL, NULL, NULL,
                NULL, aom_h_predictor_4x4_dspr2, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL)
#endif  // HAVE_DSPR2

#if HAVE_NEON
INTRA_PRED_TEST(NEON, TX_4X4, aom_dc_predictor_4x4_neon,
                aom_dc_left_predictor_4x4_neon, aom_dc_top_predictor_4x4_neon,
                aom_dc_128_predictor_4x4_neon, aom_v_predictor_4x4_neon,
                aom_h_predictor_4x4_neon, NULL, aom_d135_predictor_4x4_neon,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
#endif  // HAVE_NEON

#if HAVE_MSA
INTRA_PRED_TEST(MSA, TX_4X4, aom_dc_predictor_4x4_msa,
                aom_dc_left_predictor_4x4_msa, aom_dc_top_predictor_4x4_msa,
                aom_dc_128_predictor_4x4_msa, aom_v_predictor_4x4_msa,
                aom_h_predictor_4x4_msa, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_MSA

// -----------------------------------------------------------------------------
// 8x8, 8x4, 8x16, 8x32

INTRA_PRED_TEST(C_1, TX_8X8, aom_dc_predictor_8x8_c,
                aom_dc_left_predictor_8x8_c, aom_dc_top_predictor_8x8_c,
                aom_dc_128_predictor_8x8_c, aom_v_predictor_8x8_c,
                aom_h_predictor_8x8_c, aom_d45e_predictor_8x8_c,
                aom_d135_predictor_8x8_c, aom_d117_predictor_8x8_c,
                aom_d153_predictor_8x8_c, aom_d207e_predictor_8x8_c,
                aom_d63e_predictor_8x8_c, aom_paeth_predictor_8x8_c,
                aom_smooth_predictor_8x8_c, aom_smooth_v_predictor_8x8_c,
                aom_smooth_h_predictor_8x8_c)

INTRA_PRED_TEST(C_2, TX_8X4, aom_dc_predictor_8x4_c,
                aom_dc_left_predictor_8x4_c, aom_dc_top_predictor_8x4_c,
                aom_dc_128_predictor_8x4_c, aom_v_predictor_8x4_c,
                aom_h_predictor_8x4_c, aom_d45e_predictor_8x4_c,
                aom_d135_predictor_8x4_c, aom_d117_predictor_8x4_c,
                aom_d153_predictor_8x4_c, aom_d207e_predictor_8x4_c,
                aom_d63e_predictor_8x4_c, aom_paeth_predictor_8x4_c,
                aom_smooth_predictor_8x4_c, aom_smooth_v_predictor_8x4_c,
                aom_smooth_h_predictor_8x4_c)

INTRA_PRED_TEST(C_3, TX_8X16, aom_dc_predictor_8x16_c,
                aom_dc_left_predictor_8x16_c, aom_dc_top_predictor_8x16_c,
                aom_dc_128_predictor_8x16_c, aom_v_predictor_8x16_c,
                aom_h_predictor_8x16_c, aom_d45e_predictor_8x16_c,
                aom_d135_predictor_8x16_c, aom_d117_predictor_8x16_c,
                aom_d153_predictor_8x16_c, aom_d207e_predictor_8x16_c,
                aom_d63e_predictor_8x16_c, aom_paeth_predictor_8x16_c,
                aom_smooth_predictor_8x16_c, aom_smooth_v_predictor_8x16_c,
                aom_smooth_h_predictor_8x16_c)

INTRA_PRED_TEST(C_4, TX_8X32, aom_dc_predictor_8x32_c,
                aom_dc_left_predictor_8x32_c, aom_dc_top_predictor_8x32_c,
                aom_dc_128_predictor_8x32_c, aom_v_predictor_8x32_c,
                aom_h_predictor_8x32_c, aom_d45e_predictor_8x32_c,
                aom_d135_predictor_8x32_c, aom_d117_predictor_8x32_c,
                aom_d153_predictor_8x32_c, aom_d207e_predictor_8x32_c,
                aom_d63e_predictor_8x32_c, aom_paeth_predictor_8x32_c,
                aom_smooth_predictor_8x32_c, aom_smooth_v_predictor_8x32_c,
                aom_smooth_h_predictor_8x32_c)

#if HAVE_SSE2
INTRA_PRED_TEST(SSE2_1, TX_8X8, aom_dc_predictor_8x8_sse2,
                aom_dc_left_predictor_8x8_sse2, aom_dc_top_predictor_8x8_sse2,
                aom_dc_128_predictor_8x8_sse2, aom_v_predictor_8x8_sse2,
                aom_h_predictor_8x8_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_8x8_sse2, aom_d63e_predictor_8x8_sse2, NULL,
                NULL, NULL, NULL)
INTRA_PRED_TEST(SSE2_2, TX_8X4, aom_dc_predictor_8x4_sse2,
                aom_dc_left_predictor_8x4_sse2, aom_dc_top_predictor_8x4_sse2,
                aom_dc_128_predictor_8x4_sse2, aom_v_predictor_8x4_sse2,
                aom_h_predictor_8x4_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_8x4_sse2, aom_d63e_predictor_8x4_sse2, NULL,
                NULL, NULL, NULL)
INTRA_PRED_TEST(SSE2_3, TX_8X16, aom_dc_predictor_8x16_sse2,
                aom_dc_left_predictor_8x16_sse2, aom_dc_top_predictor_8x16_sse2,
                aom_dc_128_predictor_8x16_sse2, aom_v_predictor_8x16_sse2,
                aom_h_predictor_8x16_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_8x16_sse2, aom_d63e_predictor_8x16_sse2,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_SSE2

#if HAVE_SSSE3
INTRA_PRED_TEST(SSSE3_1, TX_8X8, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_8x8_ssse3, NULL, NULL,
                aom_d153_predictor_8x8_ssse3, NULL, NULL,
                aom_paeth_predictor_8x8_ssse3, aom_smooth_predictor_8x8_ssse3,
                NULL, NULL)
INTRA_PRED_TEST(SSSE3_2, TX_8X4, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_8x4_ssse3, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_8x4_ssse3, aom_smooth_predictor_8x4_ssse3,
                NULL, NULL)
INTRA_PRED_TEST(SSSE3_3, TX_8X16, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_8x16_ssse3, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_8x16_ssse3, aom_smooth_predictor_8x16_ssse3,
                NULL, NULL)
#endif  // HAVE_SSSE3

#if HAVE_DSPR2
INTRA_PRED_TEST(DSPR2, TX_8X8, aom_dc_predictor_8x8_dspr2, NULL, NULL, NULL,
                NULL, aom_h_predictor_8x8_dspr2, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL)
#endif  // HAVE_DSPR2

#if HAVE_NEON
INTRA_PRED_TEST(NEON, TX_8X8, aom_dc_predictor_8x8_neon,
                aom_dc_left_predictor_8x8_neon, aom_dc_top_predictor_8x8_neon,
                aom_dc_128_predictor_8x8_neon, aom_v_predictor_8x8_neon,
                aom_h_predictor_8x8_neon, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_NEON

#if HAVE_MSA
INTRA_PRED_TEST(MSA, TX_8X8, aom_dc_predictor_8x8_msa,
                aom_dc_left_predictor_8x8_msa, aom_dc_top_predictor_8x8_msa,
                aom_dc_128_predictor_8x8_msa, aom_v_predictor_8x8_msa,
                aom_h_predictor_8x8_msa, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_MSA

// -----------------------------------------------------------------------------
// 16x16, 16x8, 16x32, 16x4, 16x64

INTRA_PRED_TEST(C_1, TX_16X16, aom_dc_predictor_16x16_c,
                aom_dc_left_predictor_16x16_c, aom_dc_top_predictor_16x16_c,
                aom_dc_128_predictor_16x16_c, aom_v_predictor_16x16_c,
                aom_h_predictor_16x16_c, aom_d45e_predictor_16x16_c,
                aom_d135_predictor_16x16_c, aom_d117_predictor_16x16_c,
                aom_d153_predictor_16x16_c, aom_d207e_predictor_16x16_c,
                aom_d63e_predictor_16x16_c, aom_paeth_predictor_16x16_c,
                aom_smooth_predictor_16x16_c, aom_smooth_v_predictor_16x16_c,
                aom_smooth_h_predictor_16x16_c)

INTRA_PRED_TEST(C_2, TX_16X8, aom_dc_predictor_16x8_c,
                aom_dc_left_predictor_16x8_c, aom_dc_top_predictor_16x8_c,
                aom_dc_128_predictor_16x8_c, aom_v_predictor_16x8_c,
                aom_h_predictor_16x8_c, aom_d45e_predictor_16x8_c,
                aom_d135_predictor_16x8_c, aom_d117_predictor_16x8_c,
                aom_d153_predictor_16x8_c, aom_d207e_predictor_16x8_c,
                aom_d63e_predictor_16x8_c, aom_paeth_predictor_16x8_c,
                aom_smooth_predictor_16x8_c, aom_smooth_v_predictor_16x8_c,
                aom_smooth_h_predictor_16x8_c)

INTRA_PRED_TEST(C_3, TX_16X32, aom_dc_predictor_16x32_c,
                aom_dc_left_predictor_16x32_c, aom_dc_top_predictor_16x32_c,
                aom_dc_128_predictor_16x32_c, aom_v_predictor_16x32_c,
                aom_h_predictor_16x32_c, aom_d45e_predictor_16x32_c,
                aom_d135_predictor_16x32_c, aom_d117_predictor_16x32_c,
                aom_d153_predictor_16x32_c, aom_d207e_predictor_16x32_c,
                aom_d63e_predictor_16x32_c, aom_paeth_predictor_16x32_c,
                aom_smooth_predictor_16x32_c, aom_smooth_v_predictor_16x32_c,
                aom_smooth_h_predictor_16x32_c)

INTRA_PRED_TEST(C_4, TX_16X4, aom_dc_predictor_16x4_c,
                aom_dc_left_predictor_16x4_c, aom_dc_top_predictor_16x4_c,
                aom_dc_128_predictor_16x4_c, aom_v_predictor_16x4_c,
                aom_h_predictor_16x4_c, aom_d45e_predictor_16x4_c,
                aom_d135_predictor_16x4_c, aom_d117_predictor_16x4_c,
                aom_d153_predictor_16x4_c, aom_d207e_predictor_16x4_c,
                aom_d63e_predictor_16x4_c, aom_paeth_predictor_16x4_c,
                aom_smooth_predictor_16x4_c, aom_smooth_v_predictor_16x4_c,
                aom_smooth_h_predictor_16x4_c)

INTRA_PRED_TEST(C_5, TX_16X64, aom_dc_predictor_16x64_c,
                aom_dc_left_predictor_16x64_c, aom_dc_top_predictor_16x64_c,
                aom_dc_128_predictor_16x64_c, aom_v_predictor_16x64_c,
                aom_h_predictor_16x64_c, aom_d45e_predictor_16x64_c,
                aom_d135_predictor_16x64_c, aom_d117_predictor_16x64_c,
                aom_d153_predictor_16x64_c, aom_d207e_predictor_16x64_c,
                aom_d63e_predictor_16x64_c, aom_paeth_predictor_16x64_c,
                aom_smooth_predictor_16x64_c, aom_smooth_v_predictor_16x64_c,
                aom_smooth_h_predictor_16x64_c)

#if HAVE_SSE2
INTRA_PRED_TEST(SSE2_1, TX_16X16, aom_dc_predictor_16x16_sse2,
                aom_dc_left_predictor_16x16_sse2,
                aom_dc_top_predictor_16x16_sse2,
                aom_dc_128_predictor_16x16_sse2, aom_v_predictor_16x16_sse2,
                aom_h_predictor_16x16_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_16x16_sse2, aom_d63e_predictor_16x16_sse2,
                NULL, NULL, NULL, NULL)
INTRA_PRED_TEST(SSE2_2, TX_16X8, aom_dc_predictor_16x8_sse2,
                aom_dc_left_predictor_16x8_sse2, aom_dc_top_predictor_16x8_sse2,
                aom_dc_128_predictor_16x8_sse2, aom_v_predictor_16x8_sse2,
                aom_h_predictor_16x8_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_16x8_sse2, aom_d63e_predictor_16x8_sse2,
                NULL, NULL, NULL, NULL)
INTRA_PRED_TEST(SSE2_3, TX_16X32, aom_dc_predictor_16x32_sse2,
                aom_dc_left_predictor_16x32_sse2,
                aom_dc_top_predictor_16x32_sse2,
                aom_dc_128_predictor_16x32_sse2, aom_v_predictor_16x32_sse2,
                aom_h_predictor_16x32_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_16x32_sse2, aom_d63e_predictor_16x32_sse2,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_SSE2

#if HAVE_SSSE3
INTRA_PRED_TEST(SSSE3_1, TX_16X16, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_16x16_ssse3, NULL, NULL,
                aom_d153_predictor_16x16_ssse3, NULL, NULL,
                aom_paeth_predictor_16x16_ssse3,
                aom_smooth_predictor_16x16_ssse3, NULL, NULL)
INTRA_PRED_TEST(SSSE3_2, TX_16X8, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_16x8_ssse3, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_16x8_ssse3, aom_smooth_predictor_16x8_ssse3,
                NULL, NULL)
INTRA_PRED_TEST(SSSE3_3, TX_16X32, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_16x32_ssse3, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_16x32_ssse3,
                aom_smooth_predictor_16x32_ssse3, NULL, NULL)
#endif  // HAVE_SSSE3

#if HAVE_AVX2
INTRA_PRED_TEST(AVX2_1, TX_16X16, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, aom_paeth_predictor_16x16_avx2,
                NULL, NULL, NULL)
INTRA_PRED_TEST(AVX2_2, TX_16X8, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, aom_paeth_predictor_16x8_avx2, NULL,
                NULL, NULL)
INTRA_PRED_TEST(AVX2_3, TX_16X32, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL, aom_paeth_predictor_16x32_avx2,
                NULL, NULL, NULL)
#endif  // HAVE_AVX2

#if HAVE_DSPR2
INTRA_PRED_TEST(DSPR2, TX_16X16, aom_dc_predictor_16x16_dspr2, NULL, NULL, NULL,
                NULL, aom_h_predictor_16x16_dspr2, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL, NULL)
#endif  // HAVE_DSPR2

#if HAVE_NEON
INTRA_PRED_TEST(NEON, TX_16X16, aom_dc_predictor_16x16_neon,
                aom_dc_left_predictor_16x16_neon,
                aom_dc_top_predictor_16x16_neon,
                aom_dc_128_predictor_16x16_neon, aom_v_predictor_16x16_neon,
                aom_h_predictor_16x16_neon, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_NEON

#if HAVE_MSA
INTRA_PRED_TEST(MSA, TX_16X16, aom_dc_predictor_16x16_msa,
                aom_dc_left_predictor_16x16_msa, aom_dc_top_predictor_16x16_msa,
                aom_dc_128_predictor_16x16_msa, aom_v_predictor_16x16_msa,
                aom_h_predictor_16x16_msa, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_MSA

// -----------------------------------------------------------------------------
// 32x32, 32x16, 32x64, 32x8

INTRA_PRED_TEST(C_1, TX_32X32, aom_dc_predictor_32x32_c,
                aom_dc_left_predictor_32x32_c, aom_dc_top_predictor_32x32_c,
                aom_dc_128_predictor_32x32_c, aom_v_predictor_32x32_c,
                aom_h_predictor_32x32_c, aom_d45e_predictor_32x32_c,
                aom_d135_predictor_32x32_c, aom_d117_predictor_32x32_c,
                aom_d153_predictor_32x32_c, aom_d207e_predictor_32x32_c,
                aom_d63e_predictor_32x32_c, aom_paeth_predictor_32x32_c,
                aom_smooth_predictor_32x32_c, aom_smooth_v_predictor_32x32_c,
                aom_smooth_h_predictor_32x32_c)

INTRA_PRED_TEST(C_2, TX_32X16, aom_dc_predictor_32x16_c,
                aom_dc_left_predictor_32x16_c, aom_dc_top_predictor_32x16_c,
                aom_dc_128_predictor_32x16_c, aom_v_predictor_32x16_c,
                aom_h_predictor_32x16_c, aom_d45e_predictor_32x16_c,
                aom_d135_predictor_32x16_c, aom_d117_predictor_32x16_c,
                aom_d153_predictor_32x16_c, aom_d207e_predictor_32x16_c,
                aom_d63e_predictor_32x16_c, aom_paeth_predictor_32x16_c,
                aom_smooth_predictor_32x16_c, aom_smooth_v_predictor_32x16_c,
                aom_smooth_h_predictor_32x16_c)

INTRA_PRED_TEST(C_3, TX_32X64, aom_dc_predictor_32x64_c,
                aom_dc_left_predictor_32x64_c, aom_dc_top_predictor_32x64_c,
                aom_dc_128_predictor_32x64_c, aom_v_predictor_32x64_c,
                aom_h_predictor_32x64_c, aom_d45e_predictor_32x64_c,
                aom_d135_predictor_32x64_c, aom_d117_predictor_32x64_c,
                aom_d153_predictor_32x64_c, aom_d207e_predictor_32x64_c,
                aom_d63e_predictor_32x64_c, aom_paeth_predictor_32x64_c,
                aom_smooth_predictor_32x64_c, aom_smooth_v_predictor_32x64_c,
                aom_smooth_h_predictor_32x64_c)

INTRA_PRED_TEST(C_4, TX_32X8, aom_dc_predictor_32x8_c,
                aom_dc_left_predictor_32x8_c, aom_dc_top_predictor_32x8_c,
                aom_dc_128_predictor_32x8_c, aom_v_predictor_32x8_c,
                aom_h_predictor_32x8_c, aom_d45e_predictor_32x8_c,
                aom_d135_predictor_32x8_c, aom_d117_predictor_32x8_c,
                aom_d153_predictor_32x8_c, aom_d207e_predictor_32x8_c,
                aom_d63e_predictor_32x8_c, aom_paeth_predictor_32x8_c,
                aom_smooth_predictor_32x8_c, aom_smooth_v_predictor_32x8_c,
                aom_smooth_h_predictor_32x8_c)

#if HAVE_SSE2
INTRA_PRED_TEST(SSE2_1, TX_32X32, aom_dc_predictor_32x32_sse2,
                aom_dc_left_predictor_32x32_sse2,
                aom_dc_top_predictor_32x32_sse2,
                aom_dc_128_predictor_32x32_sse2, aom_v_predictor_32x32_sse2,
                aom_h_predictor_32x32_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_32x32_sse2, aom_d63e_predictor_32x32_sse2,
                NULL, NULL, NULL, NULL)
INTRA_PRED_TEST(SSE2_2, TX_32X16, aom_dc_predictor_32x16_sse2,
                aom_dc_left_predictor_32x16_sse2,
                aom_dc_top_predictor_32x16_sse2,
                aom_dc_128_predictor_32x16_sse2, aom_v_predictor_32x16_sse2,
                aom_h_predictor_32x16_sse2, NULL, NULL, NULL, NULL,
                aom_d207e_predictor_32x16_sse2, aom_d63e_predictor_32x16_sse2,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_SSE2

#if HAVE_SSSE3
INTRA_PRED_TEST(SSSE3_1, TX_32X32, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_32x32_ssse3, NULL, NULL,
                aom_d153_predictor_32x32_ssse3, NULL, NULL,
                aom_paeth_predictor_32x32_ssse3,
                aom_smooth_predictor_32x32_ssse3, NULL, NULL)
INTRA_PRED_TEST(SSSE3_2, TX_32X16, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_d45e_predictor_32x16_ssse3, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_32x16_ssse3,
                aom_smooth_predictor_32x16_ssse3, NULL, NULL)
#endif  // HAVE_SSSE3

#if HAVE_AVX2
INTRA_PRED_TEST(AVX2_1, TX_32X32, aom_dc_predictor_32x32_avx2,
                aom_dc_left_predictor_32x32_avx2,
                aom_dc_top_predictor_32x32_avx2,
                aom_dc_128_predictor_32x32_avx2, aom_v_predictor_32x32_avx2,
                aom_h_predictor_32x32_avx2, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_32x32_avx2, NULL, NULL, NULL)
INTRA_PRED_TEST(AVX2_2, TX_32X16, aom_dc_predictor_32x16_avx2,
                aom_dc_left_predictor_32x16_avx2,
                aom_dc_top_predictor_32x16_avx2,
                aom_dc_128_predictor_32x16_avx2, aom_v_predictor_32x16_avx2,
                NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                aom_paeth_predictor_32x16_avx2, NULL, NULL, NULL)
#endif  // HAVE_AVX2

#if HAVE_NEON
INTRA_PRED_TEST(NEON, TX_32X32, aom_dc_predictor_32x32_neon,
                aom_dc_left_predictor_32x32_neon,
                aom_dc_top_predictor_32x32_neon,
                aom_dc_128_predictor_32x32_neon, aom_v_predictor_32x32_neon,
                aom_h_predictor_32x32_neon, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_NEON

#if HAVE_MSA
INTRA_PRED_TEST(MSA, TX_32X32, aom_dc_predictor_32x32_msa,
                aom_dc_left_predictor_32x32_msa, aom_dc_top_predictor_32x32_msa,
                aom_dc_128_predictor_32x32_msa, aom_v_predictor_32x32_msa,
                aom_h_predictor_32x32_msa, NULL, NULL, NULL, NULL, NULL, NULL,
                NULL, NULL, NULL, NULL)
#endif  // HAVE_MSA

// -----------------------------------------------------------------------------
// 64x64, 64x32, 64x16

INTRA_PRED_TEST(C_1, TX_64X64, aom_dc_predictor_64x64_c,
                aom_dc_left_predictor_64x64_c, aom_dc_top_predictor_64x64_c,
                aom_dc_128_predictor_64x64_c, aom_v_predictor_64x64_c,
                aom_h_predictor_64x64_c, aom_d45e_predictor_64x64_c,
                aom_d135_predictor_64x64_c, aom_d117_predictor_64x64_c,
                aom_d153_predictor_64x64_c, aom_d207e_predictor_64x64_c,
                aom_d63e_predictor_64x64_c, aom_paeth_predictor_64x64_c,
                aom_smooth_predictor_64x64_c, aom_smooth_v_predictor_64x64_c,
                aom_smooth_h_predictor_64x64_c)

INTRA_PRED_TEST(C_2, TX_64X32, aom_dc_predictor_64x32_c,
                aom_dc_left_predictor_64x32_c, aom_dc_top_predictor_64x32_c,
                aom_dc_128_predictor_64x32_c, aom_v_predictor_64x32_c,
                aom_h_predictor_64x32_c, aom_d45e_predictor_64x32_c,
                aom_d135_predictor_64x32_c, aom_d117_predictor_64x32_c,
                aom_d153_predictor_64x32_c, aom_d207e_predictor_64x32_c,
                aom_d63e_predictor_64x32_c, aom_paeth_predictor_64x32_c,
                aom_smooth_predictor_64x32_c, aom_smooth_v_predictor_64x32_c,
                aom_smooth_h_predictor_64x32_c)

INTRA_PRED_TEST(C_3, TX_64X16, aom_dc_predictor_64x16_c,
                aom_dc_left_predictor_64x16_c, aom_dc_top_predictor_64x16_c,
                aom_dc_128_predictor_64x16_c, aom_v_predictor_64x16_c,
                aom_h_predictor_64x16_c, aom_d45e_predictor_64x16_c,
                aom_d135_predictor_64x16_c, aom_d117_predictor_64x16_c,
                aom_d153_predictor_64x16_c, aom_d207e_predictor_64x16_c,
                aom_d63e_predictor_64x16_c, aom_paeth_predictor_64x16_c,
                aom_smooth_predictor_64x16_c, aom_smooth_v_predictor_64x16_c,
                aom_smooth_h_predictor_64x16_c)

// -----------------------------------------------------------------------------
// High Bitdepth
namespace {

typedef void (*AvxHighbdPredFunc)(uint16_t *dst, ptrdiff_t y_stride,
                                  const uint16_t *above, const uint16_t *left,
                                  int bd);

typedef IntraPredTestMem<uint16_t> Av1HighbdIntraPredTestMem;

void TestHighbdIntraPred(TX_SIZE tx_size, AvxHighbdPredFunc const *pred_funcs,
                         const char *const signatures[]) {
  const int block_width = tx_size_wide[tx_size];
  const int block_height = tx_size_high[tx_size];
  const int num_pixels_per_test =
      block_width * block_height * kNumAv1IntraFuncs;
  const int kNumTests = static_cast<int>(2.e10 / num_pixels_per_test);
  Av1HighbdIntraPredTestMem intra_pred_test_mem;
  const int bd = 12;
  intra_pred_test_mem.Init(block_width, block_height, bd);

  for (int k = 0; k < kNumAv1IntraFuncs; ++k) {
    if (pred_funcs[k] == NULL) continue;
    memcpy(intra_pred_test_mem.src, intra_pred_test_mem.ref_src,
           sizeof(intra_pred_test_mem.src));
    aom_usec_timer timer;
    aom_usec_timer_start(&timer);
    for (int num_tests = 0; num_tests < kNumTests; ++num_tests) {
      pred_funcs[k](intra_pred_test_mem.src, intra_pred_test_mem.stride,
                    intra_pred_test_mem.above, intra_pred_test_mem.left, bd);
    }
    libaom_test::ClearSystemState();
    aom_usec_timer_mark(&timer);
    const int elapsed_time =
        static_cast<int>(aom_usec_timer_elapsed(&timer) / 1000);
    CheckMd5Signature(
        tx_size, true, signatures, intra_pred_test_mem.src,
        intra_pred_test_mem.num_pixels * sizeof(*intra_pred_test_mem.src),
        elapsed_time, k);
  }
}

static const char *const kHighbdSignatures[TX_SIZES_ALL][kNumAv1IntraFuncs] = {
  {
      // 4X4
      "11f74af6c5737df472f3275cbde062fa",
      "51bea056b6447c93f6eb8f6b7e8f6f71",
      "27e97f946766331795886f4de04c5594",
      "53ab15974b049111fb596c5168ec7e3f",
      "f0b640bb176fbe4584cf3d32a9b0320a",
      "729783ca909e03afd4b47111c80d967b",
      "d631a8544ccc87702db3e98fac494657",
      "293fc903254a33754133314c6cdba81f",
      "f8074d704233e73dfd35b458c6092374",
      "aa6363d08544a1ec4da33d7a0be5640d",
      "0bdc21a3acdebc393bc2c22e71bbeada",
      "a48f7a484ba4ad3916055c7160665b56",
      "6e30009c45474a22032678b1bd579c8f",
      "e57cba016d808aa8a35619df2a65f049",
      "55a6c37f39afcbbf5abca4a985b96459",
      "a623d45b37dafec1f8a75c4c5218913d",
  },
  {
      // 8X8
      "03da8829fe94663047fd108c5fcaa71d",
      "ecdb37b8120a2d3a4c706b016bd1bfd7",
      "1d4543ed8d2b9368cb96898095fe8a75",
      "f791c9a67b913cbd82d9da8ecede30e2",
      "065c70646f4dbaff913282f55a45a441",
      "51f87123616662ef7c35691497dfd0ba",
      "4f53cf8e5f43894dc0759f43c7081f60",
      "9ffe186a6bc7db95275f1bbddd6f7aba",
      "a3258a2eae2e2bd55cb8f71351b22998",
      "8d909f0a2066e39b3216092c6289ece4",
      "6751f60655aba44aff78aaaf4e967377",
      "d31a449872fab968a8d41de578338780",
      "85c01ba03df68f9ece7bd3fa0f8980e6",
      "ad19b7dac092f56df6d054e1f67f21e7",
      "0edc415b5dd7299f7a34fb9f71d31d78",
      "2bc8ec19e9f4b77a64b8a0a1f6aec7e7",
  },
  {
      // 16X16
      "e33cb3f56a878e2fddb1b2fc51cdd275",
      "c7bff6f04b6052c8ab335d726dbbd52d",
      "d0b0b47b654a9bcc5c6008110a44589b",
      "78f5da7b10b2b9ab39f114a33b6254e9",
      "c78e31d23831abb40d6271a318fdd6f3",
      "90d1347f4ec9198a0320daecb6ff90b8",
      "e38e12830e2ee5a01a064ec5998d5948",
      "cf28bd387b81ad3e5f1a1c779a4b70a0",
      "24c304330431ddeaf630f6ce94af2eac",
      "91a329798036bf64e8e00a87b131b8b1",
      "e536338d1a8ee192b9e591855db1a222",
      "54ecd47737f71c62d24e3779585113f2",
      "e63ded54ab3d0e8728b6f24d4f01e53f",
      "35ce21fbe0ea114c089fc3489a78155d",
      "f277f6ef8e4d717f1f0dfe2706ac197d",
      "e8014d3f41256976c02e0f1e622ba2b9",
  },
  {
      // 32X32
      "a3e8056ba7e36628cce4917cd956fedd",
      "cc7d3024fe8748b512407edee045377e",
      "2aab0a0f330a1d3e19b8ecb8f06387a3",
      "a547bc3fb7b06910bf3973122a426661",
      "26f712514da95042f93d6e8dc8e431dc",
      "bb08c6e16177081daa3d936538dbc2e3",
      "4e10f10b082a5b4265080c102d34eb47",
      "42867c8553285e94ee8e4df7abafbda8",
      "6496bdee96100667833f546e1be3d640",
      "2ebfa25bf981377e682e580208504300",
      "1788695b10a6f82ae1a56686dcbcd0a9",
      "c3b9c506604a7132bbb5f4e97bdb03f0",
      "84bf83f94a51b33654ca940c6f8bc057",
      "7168b03fc31bf29596a344d6a35d007c",
      "b073a70d3672f1282236994f5d12e94b",
      "c51607aebad5dcb3c1e3b58ef9e5b84e",
  },
  {
      // 64X64
      "a6baa0d4bfb2269a94c7a38f86a4bccf",
      "3f1ef5f473a49eba743f17a3324adf9d",
      "12ac11889ae5f55b7781454efd706a6a",
      "d9a906c0e692b22e1b4414e71a704b7e",
      "47d4cadd56f70c11ff8f3e5d8df81161",
      "de997744cf24c16c5ac2a36b02b351cc",
      "6e9ba318431d50763be1be37a426d016",
      "88b8bfe302def5c82e4bff52c72dfd6e",
      "ad1690f58d04d7e88d3a2c463f9ee040",
      "e1a50991fef4d7d22fbe019b00e1d7aa",
      "52e0429b079c9722f550edfa8f014c2a",
      "0c4a2ef0e786716c2055352a070d881c",
      "23781211ae178ddeb6c4bb97a6bd7d83",
      "a79d2e28340ca34b9e37daabbf030f63",
      "0372bd3ddfc258750a6ac106b70587f4",
      "228ef625d9460cbf6fa253a16a730976",
  },
  {
      // 4X8
      "22d519b796d59644043466320e4ccd14",
      "09513a738c49b3f9542d27f34abbe1d5",
      "807ae5e8813443ff01e71be6efacfb69",
      "cbfa18d0293430b6e9708b0be1fd2394",
      "346c354c34ec7fa780b576db355dab88",
      "f97dae85c35359632380b09ca98d611e",
      "aed1beef71de33856c814ff7d63dd9db",
      "49c47c04dd3d23d6fc5cc32bf9d40ae4",
      "a24aade6e22b323ee28c8bf08aa2d234",
      "aefef502f9e144e71cd27dc7383b3c28",
      "b284ae5277b85ebdd16b5952149f7458",
      "8dc5791167271f6f347582e07379f580",
      "698ae351d8896d89ed9e4e67b6e53eda",
      "dcc197034a9c45a3d8238bf085835f4e",
      "7a35e2c42ffdc2efc2d6d1d75a100fc7",
      "41ab6cebd4516c87a91b2a593e2c2506",
  },
  {
      // 8X4
      "d58cd4c4bf3b7bbaa5db5e1a5622ec78",
      "6e572c35aa782d00cafcb99e9ea047ea",
      "e8c22a3702b416dc9ab974505afbed09",
      "aaa4e4762a795aad7ad74de0c662c4e4",
      "a19f9101967383c3dcbd516dc317a291",
      "9ab8cb91f1a595b9ebe3fe8de58031aa",
      "c6c7d65264397d4d31e378e1f1cfd921",
      "5804158e463ff794b6b8a623f5d2c10d",
      "c342cdeb39aae4c4f7be10e057029298",
      "c1bbbcfe4b25f6b8eca6ad2f7ee793d3",
      "98d1dab8b949859b9c65298ee9f105f8",
      "396e803aaf6d7a03a231edc48b396051",
      "2cf9021d5f1169268699807ee118b65f",
      "ee9605fcbd6fb871f1c5cd81a6989327",
      "b4871af8316089e3e23522175df7e93f",
      "d33301e1c2cb173be46792a22d19881a",
  },
  {
      // 8X16
      "4562de1d0336610880fdd5685498a9ec",
      "16310fa7076394f16fc85c4b149d89c9",
      "0e94af88e1dc573b6f0f499cddd1f530",
      "dfd245ee20d091c67809160340365aa9",
      "d3562504327f70c096c5be23fd8a3747",
      "601b853558502acbb5135eadd2da117a",
      "e83f9a8bc16b507d2ed0b6b31a25d6f5",
      "fc8427d942246e8cba81247bb294afb5",
      "89cde712e4c1ef675ea156ad679c62c7",
      "0a68c2b28c3b171ad797cf76a7058f10",
      "e70724010e12d8f374cedd3910ceb0d5",
      "ad7987e91267503ba6fd3e8be42eb48c",
      "3c624345a723a1b2b1bea05a6a08bc99",
      "2a9c781de609e0184cc7ab442050f4e5",
      "0ddc5035c22252747126b61fc238c74d",
      "e43f5d83bab759af69c7b6773fc8f9b2",
  },
  {
      // 16X8
      "a57d6b5a9bfd30c29591d8717ace9c51",
      "f5907ba97ee6c53e339e953fc8d845ee",
      "ea3aa727913ce45af06f89dd1808db5f",
      "408af4f23e48d14b48ee35ae094fcd18",
      "85c41cbcb5d744f7961e8950026fbffe",
      "8a4e588a837638887ba671f8d4910485",
      "caae3cc3d419bbd28aa389dbe4febee1",
      "ea67fb80d71b6471467c79662af1186c",
      "c83f7252412dd1ad2fc6af848e7f6be8",
      "f45af3d697f42f1b9b8def4e46bac78c",
      "dca4a2aaf5f63db387e264ba5963943a",
      "d01b1bcc50b4b66c1231142eae628cd3",
      "b792d8826b67a21757ea7097cff9e05b",
      "f94ce7101bb87fd3bb9312112527dbf4",
      "688c6660a6dc6fa61fa1aa38e708c209",
      "0cdf641b4f81d69509c92ae0b93ef5ff",
  },
  {
      // 16X32
      "aee4b3b0e3cc02d48e2c40d77f807927",
      "8baef2b2e789f79c8df9d90ad10f34a4",
      "038c38ee3c4f090bb8d736eab136aafc",
      "1a3de2aaeaffd68a9fd6c7f6557b83f3",
      "385c6e0ea29421dd81011a2934641e26",
      "6cf96c285d1a2d4787f955dad715b08c",
      "21f82421fda1c3afca8baca0dc048a52",
      "eac3734852c99a051f6d15a921d9e7b9",
      "c81f7ffec79508bf78d0f2c67d8abe96",
      "14b8c62304f65a06653b9b35dfe12d97",
      "e0893310042511275ae04e5186ee5326",
      "b4f05903a6191093be719794417ac6fd",
      "2d7f75dcd73b9528c8396279ff09ff3a",
      "5a63cd1841e4ed470e4ca5ef845f2281",
      "610d899ca945fbead33287d4335a8b32",
      "6bafaad81fce37be46730187e78d8b11",
  },
  {
      // 32X16
      "290b23c9f5a1de7905bfa71a942da29b",
      "701e7b82593c66da5052fc4b6afd79ce",
      "4da828c5455cd246735a663fbb204989",
      "e3fbeaf234efece8dbd752b77226200c",
      "4d1d8c969f05155a7e7e84cf7aad021b",
      "c22e4877c2c946d5bdc0d542e29e70cf",
      "ffd86b234d65c2e1386a5b5b5c188a69",
      "50aaaa7d90e300b635ab18cdd73e189b",
      "a945dc7429df168e2169d81b58a15859",
      "66725070d7fad02dee78730ba0843e19",
      "33d873cb05d45df2af4ff59033833db7",
      "0dd783695b69271f65d56f5516fa6dc0",
      "8ac1ce815e7780500f842b0beb0bb980",
      "9fee2e2502b507f25bfad30a55b0b610",
      "4ced9c212ec6f9956e27f68a91b59fef",
      "4a7a0b93f138bb0863e4e465b01ec0b1",
  },
  {
      // 32X64
      "ad9cfc395a5c5644a21d958c7274ac14",
      "f29d6d03c143ddf96fef04c19f2c8333",
      "a8bdc852ef704dd4975c61893e8fbc3f",
      "7d0bd7dea26226741dbca9a97f27fa74",
      "45c27c5cca9a91b6ae8379feb0881c9f",
      "8a0b78df1e001b85c874d686eac4aa1b",
      "29d0a3c23a5d0d6a81e664b305f0d356",
      "840a14b5f883a786bd6fb71561e6dee5",
      "1167478a1908d0ce304bbd085b0c9bce",
      "397a4ab65842a71023b073ac95b945c1",
      "93160cafd18ee45cf9c4589ab2fced70",
      "a5ee4a18c10428cb1062533ed9c2fceb",
      "ce9fa75fac54a3f6c0cc3f2083b938f1",
      "c0dca10d88762c954af18dc9e3791a39",
      "61df229eddfccab913b8fda4bb02f9ac",
      "4f4df6bc8d50a5600b573f0e44d70e66",
  },
  {
      // 64X32
      "db9d82921fd88b24fdff6f849f2f9c87",
      "5ecc7fdc52d2f575ad4f2d0e9e6b1e11",
      "b4581311a0a73d95dfac7f8f44591032",
      "68bd283cfd1a125f6b2ee47cee874d36",
      "804179f05c032908a5e36077bb87c994",
      "fc5fd041a8ee779015394d0c066ee43c",
      "36f36b429d76f0f7ed041738b82bb124",
      "cdcfb9da92c9d8eda66c43dfebd3a0a2",
      "a7eae27932fc7419b0d9f1590ff97bec",
      "62ec8395cd043410971a213846a6f19c",
      "b701d37ce40499d3579d6790928fbe16",
      "3ed8b08ef3ff433295f190c28d4c2d37",
      "68f5579ccadfe9a1baafb158334a3db2",
      "fe237e45e215ab06d79046da9ad71e84",
      "9a8a938a6824551bf7d21b8fd1d70ea1",
      "eb7332f2017cd96882c76e7136aeaf53",
  },
  {
      // 4X16
      "7bafa307d507747b8132e7735b7f1c73",
      "e58bc2d8213a97d1fea9cfb73d7a9633",
      "435f8a8e8bbf14dbf2fe16b2be9e97aa",
      "1d0e767b68d84acbfb50b7a04e633836",
      "5f713bd7b324fe73bb7063e35ee14e5e",
      "0dac4e1fa3d59814202715468c01ed56",
      "6e70ee2fdea08901ce7dcd2d433dca52",
      "4fd15bcde944e8eb49efc2e6bd33d13a",
      "5c907ad981a3197cc8c70ebd4f81f858",
      "d4e9d92b0e639270ba13e47fb80d8e3e",
      "331091f344b211819feb8d59ba5cbfda",
      "79add77e52fdd022c5050668f1ad907f",
      "47709d1db4a330c7a8900f450e6fddd1",
      "258e0b930bb27db28f05da9cf7d1ee7c",
      "36cf030fbae767912593efea045bfff5",
      "248d7aceabb7499febae663fae41a920",
  },
  {
      // 16X4
      "04dde98e632670e393704742c89f9067",
      "8c72543f1664651ae1fa08e2ac0adb9b",
      "2354a2cdc2773aa2df8ab4010db1be39",
      "6300ad3221c26da39b10e0e6d87ee3be",
      "8ea30b661c6ba60b28d3167f19e449b8",
      "fb6c1e4ff101a371cede63c2955cdb7e",
      "fa4e37d472479e1f754a4a75eeaf1d27",
      "0c90ad3c84ff5301bbb167d7e1b75f47",
      "03831946677cb5c1f245b6b3885bc696",
      "522967fd65c8d0067acdd7214e77d9ef",
      "a4190b92b5ccaddb3d6c0186c8c31a3c",
      "159711a69ed5dc9acc7e953a8a199a52",
      "a517c06433d6d7927b16a72184a23e92",
      "393828be5d62ab6c48668bea5e2f801a",
      "b1e510c542013eb9d6fb188dea2ce90a",
      "569a8f2fe01679ca216535ecbcdccb62",
  },
  {
      // 8X32
      "9d541865c185ca7607852852613ac1fc",
      "b96be67f08c6b5fa5ebd3411299c2f7c",
      "75a2dcf50004b9d188849b048239767e",
      "429492ff415c9fd9b050d73b2ad500f8",
      "64b3606c1ccd036bd766bd5711392cf4",
      "cb59844a0f01660ac955bae3511f1100",
      "23940565bcaf8d5fde956fa9d0a4c56f",
      "e535e8f6fd3ad53043ab4d6b5cb58a2b",
      "29a49ef89eea9ff3ac69d0d4d1559633",
      "273c16dac32eb70633c94bd5bc127459",
      "043477dce84f7dd7b42b0a603a3d722c",
      "3ce5be0f55cecf2eb7fe9b301f56b866",
      "3e076155b7a70e8828618e3f33b51e3d",
      "ed2d1f597ab7c50beff690f737cf9726",
      "7909c6a26aaf20c59d996d3e5b5f9c29",
      "965798807240c98c6f7cc9b457ed0773",
  },
  {
      // 32X8
      "36f391aa31619eec1f4d9ee95ea454cc",
      "b82648f14eeba2527357cb50bc3223cb",
      "7a7b2adf429125e8bee9d1d00a66e13f",
      "4198e4d6ba503b7cc2d7e96bb845f661",
      "96c160d2ec1be9fe0cdea9682f14d257",
      "19a450bcebaa75afb4fc6bd1fd6434af",
      "26f3d938e6420f60d5aa15f286dcbff9",
      "64158828e2c80464906e8bd6adcb91ca",
      "5b699c879712a68a6c51875ad11cb603",
      "61e16aaae96229f4a56d85abfc977db6",
      "ee9684a3dfa64b13d67065dbe4a9fb60",
      "af7748d603604be574396edbb598ec36",
      "2bd2e35967d43d0ec1c6587a36f204d5",
      "49799a99aa4ccfbd989bee92a99422f1",
      "955530e99813812a74659edeac3f5475",
      "f0316b84e378a19cd11b19a6e40b2914",
  },
  {
      // 16X64
      "8cba1b70a0bde29e8ef235cedc5faa7d",
      "96d00ddc7537bf7f196006591b733b4e",
      "cbf69d5d157c9f3355a4757b1d6e3414",
      "3ac1f642019493dec1b737d7a3a1b4e5",
      "35f9ee300d7fa3c97338e81a6f21dcd4",
      "aae335442e77c8ebc280f16ea50ba9c7",
      "2c66c8dfea52e876db12dcfa3e20c70d",
      "703c6a8ec5db2278c11142530d7949bd",
      "5dc885a48b5d8e342db1a13b1578265e",
      "6d1deb79294dcc02a822af338dfa7cb7",
      "fbc9c1e5c210ef7cbdfc90f4621a7716",
      "e40ad09551d8791742b2270425d410c0",
      "a6140fdac2278644328be094d88731db",
      "2df93621b6ff100f7008432d509f4161",
      "c77bf5aee39e7ed4a3dd715f816f452a",
      "02109bd63557d90225c32a8f1338258e",
  },
  {
      // 64X16
      "a5e2f9fb685d5f4a048e9a96affd25a4",
      "1348f249690d9eefe09d9ad7ead2c801",
      "525da4b187acd81b1ff1116b60461141",
      "e99d072de858094c98b01bd4a6772634",
      "873bfa9dc24693f19721f7c8d527f7d3",
      "0acfc6507bd3468e9679efc127d6e4b9",
      "59afa6cbcd0735c674e65f8d40fc1ea9",
      "dd6bb2d2caa307e3a697ac4b7355b860",
      "b05b5723b1f174b14cf11766182fac08",
      "738fb2ce29a8c11a28bf423349366373",
      "580764bea9806b8bbc7590283cd2a6e7",
      "1fc5810b5f37e80bb231a725a04b4162",
      "57d03f8d079c7264854e22ac1157cfae",
      "6c2c4036f70c7d957a9399b5436c0774",
      "42b8e4a97b7f8416c72a5148c031c0b1",
      "a38a2c5f79993dfae8530e9e25800893",
  },
};

}  // namespace

#define HIGHBD_INTRA_PRED_TEST(arch, tx_size, dc, dc_left, dc_top, dc_128, v, \
                               h, d45e, d135, d117, d153, d207e, d63e, paeth, \
                               smooth, smooth_v, smooth_h)                    \
  TEST(arch, DISABLED_##TestHighbdIntraPred_##tx_size) {                      \
    static const AvxHighbdPredFunc aom_intra_pred[] = {                       \
      dc,   dc_left, dc_top, dc_128, v,     h,      d45e,     d135,           \
      d117, d153,    d207e,  d63e,   paeth, smooth, smooth_v, smooth_h        \
    };                                                                        \
    TestHighbdIntraPred(tx_size, aom_intra_pred, kHighbdSignatures[tx_size]); \
  }

// -----------------------------------------------------------------------------
// 4x4, 4x8, 4x16

HIGHBD_INTRA_PRED_TEST(
    C_1, TX_4X4, aom_highbd_dc_predictor_4x4_c,
    aom_highbd_dc_left_predictor_4x4_c, aom_highbd_dc_top_predictor_4x4_c,
    aom_highbd_dc_128_predictor_4x4_c, aom_highbd_v_predictor_4x4_c,
    aom_highbd_h_predictor_4x4_c, aom_highbd_d45e_predictor_4x4_c,
    aom_highbd_d135_predictor_4x4_c, aom_highbd_d117_predictor_4x4_c,
    aom_highbd_d153_predictor_4x4_c, aom_highbd_d207e_predictor_4x4_c,
    aom_highbd_d63e_predictor_4x4_c, aom_highbd_paeth_predictor_4x4_c,
    aom_highbd_smooth_predictor_4x4_c, aom_highbd_smooth_v_predictor_4x4_c,
    aom_highbd_smooth_h_predictor_4x4_c)

HIGHBD_INTRA_PRED_TEST(
    C_2, TX_4X8, aom_highbd_dc_predictor_4x8_c,
    aom_highbd_dc_left_predictor_4x8_c, aom_highbd_dc_top_predictor_4x8_c,
    aom_highbd_dc_128_predictor_4x8_c, aom_highbd_v_predictor_4x8_c,
    aom_highbd_h_predictor_4x8_c, aom_highbd_d45e_predictor_4x8_c,
    aom_highbd_d135_predictor_4x8_c, aom_highbd_d117_predictor_4x8_c,
    aom_highbd_d153_predictor_4x8_c, aom_highbd_d207e_predictor_4x8_c,
    aom_highbd_d63e_predictor_4x8_c, aom_highbd_paeth_predictor_4x8_c,
    aom_highbd_smooth_predictor_4x8_c, aom_highbd_smooth_v_predictor_4x8_c,
    aom_highbd_smooth_h_predictor_4x8_c)

HIGHBD_INTRA_PRED_TEST(
    C_3, TX_4X16, aom_highbd_dc_predictor_4x16_c,
    aom_highbd_dc_left_predictor_4x16_c, aom_highbd_dc_top_predictor_4x16_c,
    aom_highbd_dc_128_predictor_4x16_c, aom_highbd_v_predictor_4x16_c,
    aom_highbd_h_predictor_4x16_c, aom_highbd_d45e_predictor_4x16_c,
    aom_highbd_d135_predictor_4x16_c, aom_highbd_d117_predictor_4x16_c,
    aom_highbd_d153_predictor_4x16_c, aom_highbd_d207e_predictor_4x16_c,
    aom_highbd_d63e_predictor_4x16_c, aom_highbd_paeth_predictor_4x16_c,
    aom_highbd_smooth_predictor_4x16_c, aom_highbd_smooth_v_predictor_4x16_c,
    aom_highbd_smooth_h_predictor_4x16_c)

#if HAVE_SSE2
HIGHBD_INTRA_PRED_TEST(
    SSE2_1, TX_4X4, aom_highbd_dc_predictor_4x4_sse2,
    aom_highbd_dc_left_predictor_4x4_sse2, aom_highbd_dc_top_predictor_4x4_sse2,
    aom_highbd_dc_128_predictor_4x4_sse2, aom_highbd_v_predictor_4x4_sse2,
    aom_highbd_h_predictor_4x4_sse2, aom_highbd_d45e_predictor_4x4_sse2,
    aom_highbd_d135_predictor_4x4_sse2, aom_highbd_d117_predictor_4x4_sse2,
    aom_highbd_d153_predictor_4x4_sse2, aom_highbd_d207e_predictor_4x4_sse2,
    aom_highbd_d63e_predictor_4x4_sse2, NULL, NULL, NULL, NULL)

HIGHBD_INTRA_PRED_TEST(
    SSE2_2, TX_4X8, aom_highbd_dc_predictor_4x8_sse2,
    aom_highbd_dc_left_predictor_4x8_sse2, aom_highbd_dc_top_predictor_4x8_sse2,
    aom_highbd_dc_128_predictor_4x8_sse2, aom_highbd_v_predictor_4x8_sse2,
    aom_highbd_h_predictor_4x8_sse2, aom_highbd_d45e_predictor_4x8_sse2, NULL,
    NULL, NULL, aom_highbd_d207e_predictor_4x8_sse2,
    aom_highbd_d63e_predictor_4x8_sse2, NULL, NULL, NULL, NULL)
#endif

// -----------------------------------------------------------------------------
// 8x8, 8x4, 8x16, 8x32

HIGHBD_INTRA_PRED_TEST(
    C_1, TX_8X8, aom_highbd_dc_predictor_8x8_c,
    aom_highbd_dc_left_predictor_8x8_c, aom_highbd_dc_top_predictor_8x8_c,
    aom_highbd_dc_128_predictor_8x8_c, aom_highbd_v_predictor_8x8_c,
    aom_highbd_h_predictor_8x8_c, aom_highbd_d45e_predictor_8x8_c,
    aom_highbd_d135_predictor_8x8_c, aom_highbd_d117_predictor_8x8_c,
    aom_highbd_d153_predictor_8x8_c, aom_highbd_d207e_predictor_8x8_c,
    aom_highbd_d63e_predictor_8x8_c, aom_highbd_paeth_predictor_8x8_c,
    aom_highbd_smooth_predictor_8x8_c, aom_highbd_smooth_v_predictor_8x8_c,
    aom_highbd_smooth_h_predictor_8x8_c)

HIGHBD_INTRA_PRED_TEST(
    C_2, TX_8X4, aom_highbd_dc_predictor_8x4_c,
    aom_highbd_dc_left_predictor_8x4_c, aom_highbd_dc_top_predictor_8x4_c,
    aom_highbd_dc_128_predictor_8x4_c, aom_highbd_v_predictor_8x4_c,
    aom_highbd_h_predictor_8x4_c, aom_highbd_d45e_predictor_8x4_c,
    aom_highbd_d135_predictor_8x4_c, aom_highbd_d117_predictor_8x4_c,
    aom_highbd_d153_predictor_8x4_c, aom_highbd_d207e_predictor_8x4_c,
    aom_highbd_d63e_predictor_8x4_c, aom_highbd_paeth_predictor_8x4_c,
    aom_highbd_smooth_predictor_8x4_c, aom_highbd_smooth_v_predictor_8x4_c,
    aom_highbd_smooth_h_predictor_8x4_c)

HIGHBD_INTRA_PRED_TEST(
    C_3, TX_8X16, aom_highbd_dc_predictor_8x16_c,
    aom_highbd_dc_left_predictor_8x16_c, aom_highbd_dc_top_predictor_8x16_c,
    aom_highbd_dc_128_predictor_8x16_c, aom_highbd_v_predictor_8x16_c,
    aom_highbd_h_predictor_8x16_c, aom_highbd_d45e_predictor_8x16_c,
    aom_highbd_d135_predictor_8x16_c, aom_highbd_d117_predictor_8x16_c,
    aom_highbd_d153_predictor_8x16_c, aom_highbd_d207e_predictor_8x16_c,
    aom_highbd_d63e_predictor_8x16_c, aom_highbd_paeth_predictor_8x16_c,
    aom_highbd_smooth_predictor_8x16_c, aom_highbd_smooth_v_predictor_8x16_c,
    aom_highbd_smooth_h_predictor_8x16_c)

HIGHBD_INTRA_PRED_TEST(
    C_4, TX_8X32, aom_highbd_dc_predictor_8x32_c,
    aom_highbd_dc_left_predictor_8x32_c, aom_highbd_dc_top_predictor_8x32_c,
    aom_highbd_dc_128_predictor_8x32_c, aom_highbd_v_predictor_8x32_c,
    aom_highbd_h_predictor_8x32_c, aom_highbd_d45e_predictor_8x32_c,
    aom_highbd_d135_predictor_8x32_c, aom_highbd_d117_predictor_8x32_c,
    aom_highbd_d153_predictor_8x32_c, aom_highbd_d207e_predictor_8x32_c,
    aom_highbd_d63e_predictor_8x32_c, aom_highbd_paeth_predictor_8x32_c,
    aom_highbd_smooth_predictor_8x32_c, aom_highbd_smooth_v_predictor_8x32_c,
    aom_highbd_smooth_h_predictor_8x32_c)

#if HAVE_SSE2
HIGHBD_INTRA_PRED_TEST(
    SSE2_1, TX_8X8, aom_highbd_dc_predictor_8x8_sse2,
    aom_highbd_dc_left_predictor_8x8_sse2, aom_highbd_dc_top_predictor_8x8_sse2,
    aom_highbd_dc_128_predictor_8x8_sse2, aom_highbd_v_predictor_8x8_sse2,
    aom_highbd_h_predictor_8x8_sse2, aom_highbd_d45e_predictor_8x8_sse2, NULL,
    NULL, NULL, aom_highbd_d207e_predictor_8x8_sse2,
    aom_highbd_d63e_predictor_8x8_sse2, NULL, NULL, NULL, NULL)
HIGHBD_INTRA_PRED_TEST(
    SSE2_2, TX_8X4, aom_highbd_dc_predictor_8x4_sse2,
    aom_highbd_dc_left_predictor_8x4_sse2, aom_highbd_dc_top_predictor_8x4_sse2,
    aom_highbd_dc_128_predictor_8x4_sse2, aom_highbd_v_predictor_8x4_sse2,
    aom_highbd_h_predictor_8x4_sse2, aom_highbd_d45e_predictor_8x4_sse2, NULL,
    NULL, NULL, aom_highbd_d207e_predictor_8x4_sse2,
    aom_highbd_d63e_predictor_8x4_sse2, NULL, NULL, NULL, NULL)
HIGHBD_INTRA_PRED_TEST(SSE2_3, TX_8X16, aom_highbd_dc_predictor_8x16_sse2,
                       aom_highbd_dc_left_predictor_8x16_sse2,
                       aom_highbd_dc_top_predictor_8x16_sse2,
                       aom_highbd_dc_128_predictor_8x16_sse2,
                       aom_highbd_v_predictor_8x16_sse2,
                       aom_highbd_h_predictor_8x16_sse2,
                       aom_highbd_d45e_predictor_8x16_sse2, NULL, NULL, NULL,
                       aom_highbd_d207e_predictor_8x16_sse2,
                       aom_highbd_d63e_predictor_8x16_sse2, NULL, NULL, NULL,
                       NULL)
#endif

#if HAVE_SSSE3
HIGHBD_INTRA_PRED_TEST(SSSE3, TX_8X8, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       aom_highbd_d135_predictor_8x8_ssse3,
                       aom_highbd_d117_predictor_8x8_ssse3,
                       aom_highbd_d153_predictor_8x8_ssse3, NULL, NULL, NULL,
                       NULL, NULL, NULL)
#endif

// -----------------------------------------------------------------------------
// 16x16, 16x8, 16x32, 16x4, 16x64

HIGHBD_INTRA_PRED_TEST(
    C_1, TX_16X16, aom_highbd_dc_predictor_16x16_c,
    aom_highbd_dc_left_predictor_16x16_c, aom_highbd_dc_top_predictor_16x16_c,
    aom_highbd_dc_128_predictor_16x16_c, aom_highbd_v_predictor_16x16_c,
    aom_highbd_h_predictor_16x16_c, aom_highbd_d45e_predictor_16x16_c,
    aom_highbd_d135_predictor_16x16_c, aom_highbd_d117_predictor_16x16_c,
    aom_highbd_d153_predictor_16x16_c, aom_highbd_d207e_predictor_16x16_c,
    aom_highbd_d63e_predictor_16x16_c, aom_highbd_paeth_predictor_16x16_c,
    aom_highbd_smooth_predictor_16x16_c, aom_highbd_smooth_v_predictor_16x16_c,
    aom_highbd_smooth_h_predictor_16x16_c)

HIGHBD_INTRA_PRED_TEST(
    C_2, TX_16X8, aom_highbd_dc_predictor_16x8_c,
    aom_highbd_dc_left_predictor_16x8_c, aom_highbd_dc_top_predictor_16x8_c,
    aom_highbd_dc_128_predictor_16x8_c, aom_highbd_v_predictor_16x8_c,
    aom_highbd_h_predictor_16x8_c, aom_highbd_d45e_predictor_16x8_c,
    aom_highbd_d135_predictor_16x8_c, aom_highbd_d117_predictor_16x8_c,
    aom_highbd_d153_predictor_16x8_c, aom_highbd_d207e_predictor_16x8_c,
    aom_highbd_d63e_predictor_16x8_c, aom_highbd_paeth_predictor_16x8_c,
    aom_highbd_smooth_predictor_16x8_c, aom_highbd_smooth_v_predictor_16x8_c,
    aom_highbd_smooth_h_predictor_16x8_c)

HIGHBD_INTRA_PRED_TEST(
    C_3, TX_16X32, aom_highbd_dc_predictor_16x32_c,
    aom_highbd_dc_left_predictor_16x32_c, aom_highbd_dc_top_predictor_16x32_c,
    aom_highbd_dc_128_predictor_16x32_c, aom_highbd_v_predictor_16x32_c,
    aom_highbd_h_predictor_16x32_c, aom_highbd_d45e_predictor_16x32_c,
    aom_highbd_d135_predictor_16x32_c, aom_highbd_d117_predictor_16x32_c,
    aom_highbd_d153_predictor_16x32_c, aom_highbd_d207e_predictor_16x32_c,
    aom_highbd_d63e_predictor_16x32_c, aom_highbd_paeth_predictor_16x32_c,
    aom_highbd_smooth_predictor_16x32_c, aom_highbd_smooth_v_predictor_16x32_c,
    aom_highbd_smooth_h_predictor_16x32_c)

HIGHBD_INTRA_PRED_TEST(
    C_4, TX_16X4, aom_highbd_dc_predictor_16x4_c,
    aom_highbd_dc_left_predictor_16x4_c, aom_highbd_dc_top_predictor_16x4_c,
    aom_highbd_dc_128_predictor_16x4_c, aom_highbd_v_predictor_16x4_c,
    aom_highbd_h_predictor_16x4_c, aom_highbd_d45e_predictor_16x4_c,
    aom_highbd_d135_predictor_16x4_c, aom_highbd_d117_predictor_16x4_c,
    aom_highbd_d153_predictor_16x4_c, aom_highbd_d207e_predictor_16x4_c,
    aom_highbd_d63e_predictor_16x4_c, aom_highbd_paeth_predictor_16x4_c,
    aom_highbd_smooth_predictor_16x4_c, aom_highbd_smooth_v_predictor_16x4_c,
    aom_highbd_smooth_h_predictor_16x4_c)

HIGHBD_INTRA_PRED_TEST(
    C_5, TX_16X64, aom_highbd_dc_predictor_16x64_c,
    aom_highbd_dc_left_predictor_16x64_c, aom_highbd_dc_top_predictor_16x64_c,
    aom_highbd_dc_128_predictor_16x64_c, aom_highbd_v_predictor_16x64_c,
    aom_highbd_h_predictor_16x64_c, aom_highbd_d45e_predictor_16x64_c,
    aom_highbd_d135_predictor_16x64_c, aom_highbd_d117_predictor_16x64_c,
    aom_highbd_d153_predictor_16x64_c, aom_highbd_d207e_predictor_16x64_c,
    aom_highbd_d63e_predictor_16x64_c, aom_highbd_paeth_predictor_16x64_c,
    aom_highbd_smooth_predictor_16x64_c, aom_highbd_smooth_v_predictor_16x64_c,
    aom_highbd_smooth_h_predictor_16x64_c)

#if HAVE_SSE2
HIGHBD_INTRA_PRED_TEST(SSE2_1, TX_16X16, aom_highbd_dc_predictor_16x16_sse2,
                       aom_highbd_dc_left_predictor_16x16_sse2,
                       aom_highbd_dc_top_predictor_16x16_sse2,
                       aom_highbd_dc_128_predictor_16x16_sse2,
                       aom_highbd_v_predictor_16x16_sse2,
                       aom_highbd_h_predictor_16x16_sse2, NULL, NULL, NULL,
                       NULL, aom_highbd_d207e_predictor_16x16_sse2, NULL, NULL,
                       NULL, NULL, NULL)
HIGHBD_INTRA_PRED_TEST(SSE2_2, TX_16X8, aom_highbd_dc_predictor_16x8_sse2,
                       aom_highbd_dc_left_predictor_16x8_sse2,
                       aom_highbd_dc_top_predictor_16x8_sse2,
                       aom_highbd_dc_128_predictor_16x8_sse2,
                       aom_highbd_v_predictor_16x8_sse2,
                       aom_highbd_h_predictor_16x8_sse2, NULL, NULL, NULL, NULL,
                       aom_highbd_d207e_predictor_16x8_sse2, NULL, NULL, NULL,
                       NULL, NULL)
HIGHBD_INTRA_PRED_TEST(SSE2_3, TX_16X32, aom_highbd_dc_predictor_16x32_sse2,
                       aom_highbd_dc_left_predictor_16x32_sse2,
                       aom_highbd_dc_top_predictor_16x32_sse2,
                       aom_highbd_dc_128_predictor_16x32_sse2,
                       aom_highbd_v_predictor_16x32_sse2,
                       aom_highbd_h_predictor_16x32_sse2, NULL, NULL, NULL,
                       NULL, aom_highbd_d207e_predictor_16x32_sse2, NULL, NULL,
                       NULL, NULL, NULL)
#endif

#if HAVE_SSSE3
HIGHBD_INTRA_PRED_TEST(SSSE3_1, TX_16X16, NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, aom_highbd_d135_predictor_16x16_ssse3,
                       aom_highbd_d117_predictor_16x16_ssse3,
                       aom_highbd_d153_predictor_16x16_ssse3, NULL, NULL, NULL,
                       NULL, NULL, NULL)
#endif

#if HAVE_AVX2
HIGHBD_INTRA_PRED_TEST(AVX2_1, TX_16X16, NULL, NULL, NULL, NULL, NULL, NULL,
                       aom_highbd_d45e_predictor_16x16_avx2, NULL, NULL, NULL,
                       NULL, aom_highbd_d63e_predictor_16x16_avx2, NULL, NULL,
                       NULL, NULL)

HIGHBD_INTRA_PRED_TEST(AVX2_2, TX_16X8, NULL, NULL, NULL, NULL, NULL, NULL,
                       aom_highbd_d45e_predictor_16x8_avx2, NULL, NULL, NULL,
                       NULL, aom_highbd_d63e_predictor_16x8_avx2, NULL, NULL,
                       NULL, NULL)

HIGHBD_INTRA_PRED_TEST(AVX2_3, TX_16X32, NULL, NULL, NULL, NULL, NULL, NULL,
                       aom_highbd_d45e_predictor_16x32_avx2, NULL, NULL, NULL,
                       NULL, aom_highbd_d63e_predictor_16x32_avx2, NULL, NULL,
                       NULL, NULL)
#endif

// -----------------------------------------------------------------------------
// 32x32, 32x16, 32x64, 32x8

HIGHBD_INTRA_PRED_TEST(
    C_1, TX_32X32, aom_highbd_dc_predictor_32x32_c,
    aom_highbd_dc_left_predictor_32x32_c, aom_highbd_dc_top_predictor_32x32_c,
    aom_highbd_dc_128_predictor_32x32_c, aom_highbd_v_predictor_32x32_c,
    aom_highbd_h_predictor_32x32_c, aom_highbd_d45e_predictor_32x32_c,
    aom_highbd_d135_predictor_32x32_c, aom_highbd_d117_predictor_32x32_c,
    aom_highbd_d153_predictor_32x32_c, aom_highbd_d207e_predictor_32x32_c,
    aom_highbd_d63e_predictor_32x32_c, aom_highbd_paeth_predictor_32x32_c,
    aom_highbd_smooth_predictor_32x32_c, aom_highbd_smooth_v_predictor_32x32_c,
    aom_highbd_smooth_h_predictor_32x32_c)

HIGHBD_INTRA_PRED_TEST(
    C_2, TX_32X16, aom_highbd_dc_predictor_32x16_c,
    aom_highbd_dc_left_predictor_32x16_c, aom_highbd_dc_top_predictor_32x16_c,
    aom_highbd_dc_128_predictor_32x16_c, aom_highbd_v_predictor_32x16_c,
    aom_highbd_h_predictor_32x16_c, aom_highbd_d45e_predictor_32x16_c,
    aom_highbd_d135_predictor_32x16_c, aom_highbd_d117_predictor_32x16_c,
    aom_highbd_d153_predictor_32x16_c, aom_highbd_d207e_predictor_32x16_c,
    aom_highbd_d63e_predictor_32x16_c, aom_highbd_paeth_predictor_32x16_c,
    aom_highbd_smooth_predictor_32x16_c, aom_highbd_smooth_v_predictor_32x16_c,
    aom_highbd_smooth_h_predictor_32x16_c)

HIGHBD_INTRA_PRED_TEST(
    C_3, TX_32X64, aom_highbd_dc_predictor_32x64_c,
    aom_highbd_dc_left_predictor_32x64_c, aom_highbd_dc_top_predictor_32x64_c,
    aom_highbd_dc_128_predictor_32x64_c, aom_highbd_v_predictor_32x64_c,
    aom_highbd_h_predictor_32x64_c, aom_highbd_d45e_predictor_32x64_c,
    aom_highbd_d135_predictor_32x64_c, aom_highbd_d117_predictor_32x64_c,
    aom_highbd_d153_predictor_32x64_c, aom_highbd_d207e_predictor_32x64_c,
    aom_highbd_d63e_predictor_32x64_c, aom_highbd_paeth_predictor_32x64_c,
    aom_highbd_smooth_predictor_32x64_c, aom_highbd_smooth_v_predictor_32x64_c,
    aom_highbd_smooth_h_predictor_32x64_c)

HIGHBD_INTRA_PRED_TEST(
    C_4, TX_32X8, aom_highbd_dc_predictor_32x8_c,
    aom_highbd_dc_left_predictor_32x8_c, aom_highbd_dc_top_predictor_32x8_c,
    aom_highbd_dc_128_predictor_32x8_c, aom_highbd_v_predictor_32x8_c,
    aom_highbd_h_predictor_32x8_c, aom_highbd_d45e_predictor_32x8_c,
    aom_highbd_d135_predictor_32x8_c, aom_highbd_d117_predictor_32x8_c,
    aom_highbd_d153_predictor_32x8_c, aom_highbd_d207e_predictor_32x8_c,
    aom_highbd_d63e_predictor_32x8_c, aom_highbd_paeth_predictor_32x8_c,
    aom_highbd_smooth_predictor_32x8_c, aom_highbd_smooth_v_predictor_32x8_c,
    aom_highbd_smooth_h_predictor_32x8_c)

#if HAVE_SSE2
HIGHBD_INTRA_PRED_TEST(SSE2_1, TX_32X32, aom_highbd_dc_predictor_32x32_sse2,
                       aom_highbd_dc_left_predictor_32x32_sse2,
                       aom_highbd_dc_top_predictor_32x32_sse2,
                       aom_highbd_dc_128_predictor_32x32_sse2,
                       aom_highbd_v_predictor_32x32_sse2,
                       aom_highbd_h_predictor_32x32_sse2, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL)
HIGHBD_INTRA_PRED_TEST(SSE2_2, TX_32X16, aom_highbd_dc_predictor_32x16_sse2,
                       aom_highbd_dc_left_predictor_32x16_sse2,
                       aom_highbd_dc_top_predictor_32x16_sse2,
                       aom_highbd_dc_128_predictor_32x16_sse2,
                       aom_highbd_v_predictor_32x16_sse2,
                       aom_highbd_h_predictor_32x16_sse2, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL, NULL)
#endif

#if HAVE_SSSE3
HIGHBD_INTRA_PRED_TEST(SSSE3_1, TX_32X32, NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL, aom_highbd_d135_predictor_32x32_ssse3,
                       aom_highbd_d117_predictor_32x32_ssse3,
                       aom_highbd_d153_predictor_32x32_ssse3, NULL, NULL, NULL,
                       NULL, NULL, NULL)
#endif

#if HAVE_AVX2
HIGHBD_INTRA_PRED_TEST(AVX2_1, TX_32X32, NULL, NULL, NULL, NULL, NULL, NULL,
                       aom_highbd_d45e_predictor_32x32_avx2, NULL, NULL, NULL,
                       aom_highbd_d207e_predictor_32x32_avx2,
                       aom_highbd_d63e_predictor_32x32_avx2, NULL, NULL, NULL,
                       NULL)

HIGHBD_INTRA_PRED_TEST(AVX2_2, TX_32X16, NULL, NULL, NULL, NULL, NULL, NULL,
                       aom_highbd_d45e_predictor_32x16_avx2, NULL, NULL, NULL,
                       aom_highbd_d207e_predictor_32x16_avx2,
                       aom_highbd_d63e_predictor_32x16_avx2, NULL, NULL, NULL,
                       NULL)
#endif

// -----------------------------------------------------------------------------
// 64x64, 64x32, 64x16

HIGHBD_INTRA_PRED_TEST(
    C_1, TX_64X64, aom_highbd_dc_predictor_64x64_c,
    aom_highbd_dc_left_predictor_64x64_c, aom_highbd_dc_top_predictor_64x64_c,
    aom_highbd_dc_128_predictor_64x64_c, aom_highbd_v_predictor_64x64_c,
    aom_highbd_h_predictor_64x64_c, aom_highbd_d45e_predictor_64x64_c,
    aom_highbd_d135_predictor_64x64_c, aom_highbd_d117_predictor_64x64_c,
    aom_highbd_d153_predictor_64x64_c, aom_highbd_d207e_predictor_64x64_c,
    aom_highbd_d63e_predictor_64x64_c, aom_highbd_paeth_predictor_64x64_c,
    aom_highbd_smooth_predictor_64x64_c, aom_highbd_smooth_v_predictor_64x64_c,
    aom_highbd_smooth_h_predictor_64x64_c)

HIGHBD_INTRA_PRED_TEST(
    C_2, TX_64X32, aom_highbd_dc_predictor_64x32_c,
    aom_highbd_dc_left_predictor_64x32_c, aom_highbd_dc_top_predictor_64x32_c,
    aom_highbd_dc_128_predictor_64x32_c, aom_highbd_v_predictor_64x32_c,
    aom_highbd_h_predictor_64x32_c, aom_highbd_d45e_predictor_64x32_c,
    aom_highbd_d135_predictor_64x32_c, aom_highbd_d117_predictor_64x32_c,
    aom_highbd_d153_predictor_64x32_c, aom_highbd_d207e_predictor_64x32_c,
    aom_highbd_d63e_predictor_64x32_c, aom_highbd_paeth_predictor_64x32_c,
    aom_highbd_smooth_predictor_64x32_c, aom_highbd_smooth_v_predictor_64x32_c,
    aom_highbd_smooth_h_predictor_64x32_c)

HIGHBD_INTRA_PRED_TEST(
    C_3, TX_64X16, aom_highbd_dc_predictor_64x16_c,
    aom_highbd_dc_left_predictor_64x16_c, aom_highbd_dc_top_predictor_64x16_c,
    aom_highbd_dc_128_predictor_64x16_c, aom_highbd_v_predictor_64x16_c,
    aom_highbd_h_predictor_64x16_c, aom_highbd_d45e_predictor_64x16_c,
    aom_highbd_d135_predictor_64x16_c, aom_highbd_d117_predictor_64x16_c,
    aom_highbd_d153_predictor_64x16_c, aom_highbd_d207e_predictor_64x16_c,
    aom_highbd_d63e_predictor_64x16_c, aom_highbd_paeth_predictor_64x16_c,
    aom_highbd_smooth_predictor_64x16_c, aom_highbd_smooth_v_predictor_64x16_c,
    aom_highbd_smooth_h_predictor_64x16_c)

// -----------------------------------------------------------------------------

#include "test/test_libaom.cc"
