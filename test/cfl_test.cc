
#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "aom_ports/aom_timer.h"
#include "./av1_rtcd.h"
#include "test/util.h"
#include "test/acm_random.h"

using std::tr1::make_tuple;

using libaom_test::ACMRandom;

#define AVERAGE_ITERATIONS (100)

namespace {
typedef void (*subtract_fn)(int16_t *pred_buf_q3, int width, int height,
                            int16_t avg_q3);

typedef std::tr1::tuple<subtract_fn> CFL_Param;

class CFLAverageTest : public ::testing::TestWithParam<CFL_Param> {
 public:
  virtual ~CFLAverageTest() {}
  virtual void SetUp() { subtract = GET_PARAM(0); }

 protected:
  subtract_fn subtract;
};

class CFLAverageSpeedTest : public ::testing::TestWithParam<CFL_Param> {
 public:
  virtual ~CFLAverageSpeedTest() {}
  virtual void SetUp() {
    subtract = GET_PARAM(0);
    int16_t *pred_buf_q3 = pred_buf_q3_data;
    int k = 0;
    for (int j = 0; j < CFL_BUF_LINE; j++) {
      for (int i = 0; i < CFL_BUF_LINE; i++) {
        pred_buf_q3[i] = k++;
      }
      pred_buf_q3 += CFL_BUF_LINE;
    }
  }

 protected:
  subtract_fn subtract;
  int16_t pred_buf_q3_data[CFL_BUF_SQUARE];
};

static INLINE void setup(int16_t *pred_buf_q3) {
  int k = 0;
  for (int j = 0; j < CFL_BUF_LINE; j++) {
    for (int i = 0; i < CFL_BUF_LINE; i++) {
      pred_buf_q3[i] = k++;
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
}

static INLINE void test_average(subtract_fn subtract, int width, int height,
                                int16_t avg_q3) {
  int16_t pred_buf_q3_data[CFL_BUF_SQUARE];
  int16_t *pred_buf_q3 = pred_buf_q3_data;
  int k = 0;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      pred_buf_q3[i] = k++;
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
  k = -avg_q3;
  subtract(pred_buf_q3_data, width, height, avg_q3);
  pred_buf_q3 = pred_buf_q3_data;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      ASSERT_EQ(pred_buf_q3[i], k++);
    }
    pred_buf_q3 += CFL_BUF_LINE;
  }
}

void doTest(subtract_fn subtract) {
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  for (int k = 0; k < AVERAGE_ITERATIONS; k++) {
    test_average(subtract, 4, 4, rnd.Rand15Signed());
    test_average(subtract, 4, 8, rnd.Rand15Signed());
    test_average(subtract, 8, 4, rnd.Rand15Signed());
    test_average(subtract, 8, 8, rnd.Rand15Signed());
    test_average(subtract, 8, 16, rnd.Rand15Signed());
    test_average(subtract, 16, 8, rnd.Rand15Signed());
    test_average(subtract, 16, 16, rnd.Rand15Signed());
    test_average(subtract, 16, 32, rnd.Rand15Signed());
    test_average(subtract, 32, 16, rnd.Rand15Signed());
    test_average(subtract, 32, 32, rnd.Rand15Signed());
  }
}

void test_cfl_speed(subtract_fn subtract, subtract_fn ref_subtract,
                    int16_t *pred_buf_q3, int bsize, int iterations) {
  aom_usec_timer ref_timer;
  aom_usec_timer timer;

  setup(pred_buf_q3);
  aom_usec_timer_start(&ref_timer);
  for (int k = 0; k < iterations; k++) {
    ref_subtract(pred_buf_q3, bsize, bsize, k);
  }
  aom_usec_timer_mark(&ref_timer);
  int ref_elapsed_time = (int)aom_usec_timer_elapsed(&ref_timer);

  setup(pred_buf_q3);
  aom_usec_timer_start(&timer);
  for (int k = 0; k < iterations; k++) {
    subtract(pred_buf_q3, bsize, bsize, k);
  }
  aom_usec_timer_mark(&timer);
  int elapsed_time = (int)aom_usec_timer_elapsed(&timer);

#if 1
  std::cout.precision(2);
  std::cout << "[          ] " << bsize << "x" << bsize
            << ": C time = " << ref_elapsed_time
            << " us, SIMD time = " << elapsed_time << " us"
            << " (~" << ref_elapsed_time / (double)elapsed_time << "x) "
            << std::endl;
#endif

  EXPECT_GT(ref_elapsed_time, elapsed_time)
      << "Error: CDEFSpeedTest, SIMD slower than C." << std::endl
      << "C time: " << ref_elapsed_time << " us" << std::endl
      << "SIMD time: " << elapsed_time << " us" << std::endl;
}

TEST_P(CFLAverageSpeedTest, AverageSpeedTest4) {
  test_cfl_speed(subtract, &subtract_average_c, pred_buf_q3_data, 4,
                 UINT16_MAX);
}

TEST_P(CFLAverageSpeedTest, AverageSpeedTest8) {
  test_cfl_speed(subtract, &subtract_average_c, pred_buf_q3_data, 8,
                 UINT16_MAX);
}

TEST_P(CFLAverageSpeedTest, AverageSpeedTest16) {
  test_cfl_speed(subtract, &subtract_average_c, pred_buf_q3_data, 16,
                 UINT16_MAX);
}

TEST_P(CFLAverageSpeedTest, AverageSpeedTest32) {
  test_cfl_speed(subtract, &subtract_average_c, pred_buf_q3_data, 32,
                 UINT16_MAX);
}

TEST_P(CFLAverageTest, AverageSpeedTest32) { doTest(subtract); }

INSTANTIATE_TEST_CASE_P(C, CFLAverageTest,
                        ::testing::Values(make_tuple(&subtract_average_c)));
#if HAVE_SSE2
INSTANTIATE_TEST_CASE_P(SSE2, CFLAverageSpeedTest,
                        ::testing::Values(make_tuple(&subtract_average_sse2)));
#endif

}  // namespace
