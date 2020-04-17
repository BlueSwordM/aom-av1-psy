/*
 *  Copyright (c) 2019, Alliance for Open Media. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#include <stdlib.h>
#include <tuple>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

#include "config/aom_dsp_rtcd.h"

#include "test/acm_random.h"
#include "test/clear_system_state.h"
#include "test/register_state_check.h"
#include "test/util.h"

namespace {

using libaom_test::ACMRandom;

template <typename Pixel>
class AverageTestBase : public ::testing::Test {
 public:
  AverageTestBase(int width, int height)
      : width_(width), height_(height), source_data_(NULL), source_stride_(0),
        bit_depth_(8) {}

  virtual void TearDown() {
    aom_free(source_data_);
    source_data_ = NULL;
    libaom_test::ClearSystemState();
  }

 protected:
  // Handle blocks up to 4 blocks 64x64 with stride up to 128
  static const int kDataAlignment = 16;
  static const int kDataBlockSize = 64 * 128;

  virtual void SetUp() {
    source_data_ = static_cast<Pixel *>(
        aom_memalign(kDataAlignment, kDataBlockSize * sizeof(source_data_[0])));
    ASSERT_TRUE(source_data_ != NULL);
    source_stride_ = (width_ + 31) & ~31;
    bit_depth_ = 8;
    rnd_.Reset(ACMRandom::DeterministicSeed());
  }

  // Sum Pixels
  static unsigned int ReferenceAverage8x8(const Pixel *source, int pitch) {
    unsigned int average = 0;
    for (int h = 0; h < 8; ++h) {
      for (int w = 0; w < 8; ++w) average += source[h * pitch + w];
    }
    return (average + 32) >> 6;
  }

  static unsigned int ReferenceAverage4x4(const Pixel *source, int pitch) {
    unsigned int average = 0;
    for (int h = 0; h < 4; ++h) {
      for (int w = 0; w < 4; ++w) average += source[h * pitch + w];
    }
    return (average + 8) >> 4;
  }

  void FillConstant(Pixel fill_constant) {
    for (int i = 0; i < width_ * height_; ++i) {
      source_data_[i] = fill_constant;
    }
  }

  void FillRandom() {
    for (int i = 0; i < width_ * height_; ++i) {
      source_data_[i] = rnd_.Rand16() & ((1 << bit_depth_) - 1);
    }
  }

  int width_, height_;
  Pixel *source_data_;
  int source_stride_;
  int bit_depth_;

  ACMRandom rnd_;
};
typedef unsigned int (*AverageFunction)(const uint8_t *s, int pitch);

// Arguments: width, height, pitch, block size, avg function.
typedef std::tuple<int, int, int, int, AverageFunction> AvgFunc;

class AverageTest : public AverageTestBase<uint8_t>,
                    public ::testing::WithParamInterface<AvgFunc> {
 public:
  AverageTest() : AverageTestBase(GET_PARAM(0), GET_PARAM(1)) {}

 protected:
  void CheckAverages() {
    const int block_size = GET_PARAM(3);
    unsigned int expected = 0;
    if (block_size == 8) {
      expected =
          ReferenceAverage8x8(source_data_ + GET_PARAM(2), source_stride_);
    } else if (block_size == 4) {
      expected =
          ReferenceAverage4x4(source_data_ + GET_PARAM(2), source_stride_);
    }

    unsigned int actual;
    ASM_REGISTER_STATE_CHECK(
        actual = GET_PARAM(4)(source_data_ + GET_PARAM(2), source_stride_));

    EXPECT_EQ(expected, actual);
  }
};

TEST_P(AverageTest, MinValue) {
  FillConstant(0);
  CheckAverages();
}

TEST_P(AverageTest, MaxValue) {
  FillConstant(255);
  CheckAverages();
}

TEST_P(AverageTest, Random) {
  // The reference frame, but not the source frame, may be unaligned for
  // certain types of searches.
  for (int i = 0; i < 1000; i++) {
    FillRandom();
    CheckAverages();
  }
}

typedef void (*IntProRowFunc)(int16_t hbuf[16], uint8_t const *ref,
                              const int ref_stride, const int height);

// Params: height, asm function, c function.
typedef std::tuple<int, IntProRowFunc, IntProRowFunc> IntProRowParam;

class IntProRowTest : public AverageTestBase<uint8_t>,
                      public ::testing::WithParamInterface<IntProRowParam> {
 public:
  IntProRowTest()
      : AverageTestBase(16, GET_PARAM(0)), hbuf_asm_(NULL), hbuf_c_(NULL) {
    asm_func_ = GET_PARAM(1);
    c_func_ = GET_PARAM(2);
  }

 protected:
  virtual void SetUp() {
    source_data_ = static_cast<uint8_t *>(
        aom_memalign(kDataAlignment, kDataBlockSize * sizeof(source_data_[0])));
    ASSERT_TRUE(source_data_ != NULL);

    hbuf_asm_ = static_cast<int16_t *>(
        aom_memalign(kDataAlignment, sizeof(*hbuf_asm_) * 16));
    hbuf_c_ = static_cast<int16_t *>(
        aom_memalign(kDataAlignment, sizeof(*hbuf_c_) * 16));
  }

  virtual void TearDown() {
    aom_free(source_data_);
    source_data_ = NULL;
    aom_free(hbuf_c_);
    hbuf_c_ = NULL;
    aom_free(hbuf_asm_);
    hbuf_asm_ = NULL;
  }

  void RunComparison() {
    ASM_REGISTER_STATE_CHECK(c_func_(hbuf_c_, source_data_, 0, height_));
    ASM_REGISTER_STATE_CHECK(asm_func_(hbuf_asm_, source_data_, 0, height_));
    EXPECT_EQ(0, memcmp(hbuf_c_, hbuf_asm_, sizeof(*hbuf_c_) * 16))
        << "Output mismatch\n";
  }

  void RunSpeedTest() {
    const int numIter = 5000000;
    printf("Height = %d number of iteration is %d \n", height_, numIter);
    aom_usec_timer c_timer_;
    aom_usec_timer_start(&c_timer_);
    for (int i = 0; i < numIter; i++) {
      c_func_(hbuf_c_, source_data_, 0, height_);
    }
    aom_usec_timer_mark(&c_timer_);

    aom_usec_timer asm_timer_;
    aom_usec_timer_start(&asm_timer_);

    for (int i = 0; i < numIter; i++) {
      asm_func_(hbuf_asm_, source_data_, 0, height_);
    }
    aom_usec_timer_mark(&asm_timer_);

    const int c_sum_time = static_cast<int>(aom_usec_timer_elapsed(&c_timer_));
    const int asm_sum_time =
        static_cast<int>(aom_usec_timer_elapsed(&asm_timer_));

    printf("c_time = %d \t simd_time = %d \t Gain = %4.2f \n", c_sum_time,
           asm_sum_time,
           (static_cast<float>(c_sum_time) / static_cast<float>(asm_sum_time)));

    EXPECT_EQ(0, memcmp(hbuf_c_, hbuf_asm_, sizeof(*hbuf_c_) * 16))
        << "Output mismatch\n";
  }

 private:
  IntProRowFunc asm_func_;
  IntProRowFunc c_func_;
  int16_t *hbuf_asm_;
  int16_t *hbuf_c_;
};

typedef int16_t (*IntProColFunc)(uint8_t const *ref, const int width);

// Params: width, asm function, c function.
typedef std::tuple<int, IntProColFunc, IntProColFunc> IntProColParam;

class IntProColTest : public AverageTestBase<uint8_t>,
                      public ::testing::WithParamInterface<IntProColParam> {
 public:
  IntProColTest() : AverageTestBase(GET_PARAM(0), 1), sum_asm_(0), sum_c_(0) {
    asm_func_ = GET_PARAM(1);
    c_func_ = GET_PARAM(2);
  }

 protected:
  void RunComparison() {
    ASM_REGISTER_STATE_CHECK(sum_c_ = c_func_(source_data_, width_));
    ASM_REGISTER_STATE_CHECK(sum_asm_ = asm_func_(source_data_, width_));
    EXPECT_EQ(sum_c_, sum_asm_) << "Output mismatch";
  }
  void RunSpeedTest() {
    const int numIter = 5000000;
    printf("Width = %d number of iteration is %d \n", width_, numIter);
    aom_usec_timer c_timer_;
    aom_usec_timer_start(&c_timer_);
    for (int i = 0; i < numIter; i++) {
      sum_c_ = c_func_(source_data_, width_);
    }
    aom_usec_timer_mark(&c_timer_);

    aom_usec_timer asm_timer_;
    aom_usec_timer_start(&asm_timer_);

    for (int i = 0; i < numIter; i++) {
      sum_asm_ = asm_func_(source_data_, width_);
    }
    aom_usec_timer_mark(&asm_timer_);

    const int c_sum_time = static_cast<int>(aom_usec_timer_elapsed(&c_timer_));
    const int asm_sum_time =
        static_cast<int>(aom_usec_timer_elapsed(&asm_timer_));

    printf("c_time = %d \t simd_time = %d \t Gain = %4.2f \n", c_sum_time,
           asm_sum_time,
           (static_cast<float>(c_sum_time) / static_cast<float>(asm_sum_time)));

    EXPECT_EQ(sum_c_, sum_asm_) << "Output mismatch \n";
  }

 private:
  IntProColFunc asm_func_;
  IntProColFunc c_func_;
  int16_t sum_asm_;
  int16_t sum_c_;
};

TEST_P(IntProRowTest, MinValue) {
  FillConstant(0);
  RunComparison();
}

TEST_P(IntProRowTest, MaxValue) {
  FillConstant(255);
  RunComparison();
}

TEST_P(IntProRowTest, Random) {
  FillRandom();
  RunComparison();
}

TEST_P(IntProRowTest, DISABLED_Speed) {
  FillRandom();
  RunSpeedTest();
}

TEST_P(IntProColTest, MinValue) {
  FillConstant(0);
  RunComparison();
}

TEST_P(IntProColTest, MaxValue) {
  FillConstant(255);
  RunComparison();
}

TEST_P(IntProColTest, Random) {
  FillRandom();
  RunComparison();
}

TEST_P(IntProColTest, DISABLED_Speed) {
  FillRandom();
  RunSpeedTest();
}

using std::make_tuple;

INSTANTIATE_TEST_SUITE_P(
    C, AverageTest,
    ::testing::Values(make_tuple(16, 16, 1, 8, &aom_avg_8x8_c),
                      make_tuple(16, 16, 1, 4, &aom_avg_4x4_c)));

#if HAVE_SSE2
INSTANTIATE_TEST_SUITE_P(
    SSE2, AverageTest,
    ::testing::Values(make_tuple(16, 16, 0, 8, &aom_avg_8x8_sse2),
                      make_tuple(16, 16, 5, 8, &aom_avg_8x8_sse2),
                      make_tuple(32, 32, 15, 8, &aom_avg_8x8_sse2),
                      make_tuple(16, 16, 0, 4, &aom_avg_4x4_sse2),
                      make_tuple(16, 16, 5, 4, &aom_avg_4x4_sse2),
                      make_tuple(32, 32, 15, 4, &aom_avg_4x4_sse2)));

INSTANTIATE_TEST_SUITE_P(
    SSE2, IntProRowTest,
    ::testing::Values(make_tuple(16, &aom_int_pro_row_sse2, &aom_int_pro_row_c),
                      make_tuple(32, &aom_int_pro_row_sse2, &aom_int_pro_row_c),
                      make_tuple(64, &aom_int_pro_row_sse2, &aom_int_pro_row_c),
                      make_tuple(128, &aom_int_pro_row_sse2,
                                 &aom_int_pro_row_c)));

INSTANTIATE_TEST_SUITE_P(
    SSE2, IntProColTest,
    ::testing::Values(make_tuple(16, &aom_int_pro_col_sse2, &aom_int_pro_col_c),
                      make_tuple(32, &aom_int_pro_col_sse2, &aom_int_pro_col_c),
                      make_tuple(64, &aom_int_pro_col_sse2, &aom_int_pro_col_c),
                      make_tuple(128, &aom_int_pro_col_sse2,
                                 &aom_int_pro_col_c)));
#endif

#if HAVE_NEON
INSTANTIATE_TEST_SUITE_P(
    NEON, AverageTest,
    ::testing::Values(make_tuple(16, 16, 0, 8, &aom_avg_8x8_neon),
                      make_tuple(16, 16, 5, 8, &aom_avg_8x8_neon),
                      make_tuple(32, 32, 15, 8, &aom_avg_8x8_neon),
                      make_tuple(16, 16, 0, 4, &aom_avg_4x4_neon),
                      make_tuple(16, 16, 5, 4, &aom_avg_4x4_neon),
                      make_tuple(32, 32, 15, 4, &aom_avg_4x4_neon)));
INSTANTIATE_TEST_SUITE_P(
    NEON, IntProRowTest,
    ::testing::Values(make_tuple(16, &aom_int_pro_row_neon, &aom_int_pro_row_c),
                      make_tuple(32, &aom_int_pro_row_neon, &aom_int_pro_row_c),
                      make_tuple(64, &aom_int_pro_row_neon, &aom_int_pro_row_c),
                      make_tuple(128, &aom_int_pro_row_neon,
                                 &aom_int_pro_row_c)));

INSTANTIATE_TEST_SUITE_P(
    NEON, IntProColTest,
    ::testing::Values(make_tuple(16, &aom_int_pro_col_neon, &aom_int_pro_col_c),
                      make_tuple(32, &aom_int_pro_col_neon, &aom_int_pro_col_c),
                      make_tuple(64, &aom_int_pro_col_neon, &aom_int_pro_col_c),
                      make_tuple(128, &aom_int_pro_col_neon,
                                 &aom_int_pro_col_c)));
#endif

typedef int (*SatdFunc)(const tran_low_t *coeffs, int length);
typedef ::testing::tuple<int, SatdFunc, SatdFunc> SatdTestParam;
class SatdTest : public ::testing::Test,
                 public ::testing::WithParamInterface<SatdTestParam> {
 protected:
  virtual void SetUp() {
    satd_size_ = GET_PARAM(0);
    satd_func_ref_ = GET_PARAM(1);
    satd_func_simd_ = GET_PARAM(2);

    rnd_.Reset(ACMRandom::DeterministicSeed());
    src_ = reinterpret_cast<tran_low_t *>(
        aom_memalign(32, sizeof(*src_) * satd_size_));
    ASSERT_TRUE(src_ != NULL);
  }
  virtual void TearDown() {
    libaom_test::ClearSystemState();
    aom_free(src_);
  }
  void FillConstant(const tran_low_t val) {
    for (int i = 0; i < satd_size_; ++i) src_[i] = val;
  }
  void FillRandom() {
    for (int i = 0; i < satd_size_; ++i) {
      src_[i] = static_cast<int16_t>(rnd_.Rand16());
    }
  }
  void Check(int expected) {
    int total_ref;
    ASM_REGISTER_STATE_CHECK(total_ref = satd_func_ref_(src_, satd_size_));
    EXPECT_EQ(expected, total_ref);

    int total_simd;
    ASM_REGISTER_STATE_CHECK(total_simd = satd_func_simd_(src_, satd_size_));
    EXPECT_EQ(expected, total_simd);
  }
  void RunComparison() {
    int total_ref;
    ASM_REGISTER_STATE_CHECK(total_ref = satd_func_ref_(src_, satd_size_));

    int total_simd;
    ASM_REGISTER_STATE_CHECK(total_simd = satd_func_simd_(src_, satd_size_));

    EXPECT_EQ(total_ref, total_simd);
  }
  void RunSpeedTest() {
    const int numIter = 500000;
    printf("size = %d number of iteration is %d \n", satd_size_, numIter);

    int total_ref;
    aom_usec_timer c_timer_;
    aom_usec_timer_start(&c_timer_);
    for (int i = 0; i < numIter; i++) {
      total_ref = satd_func_ref_(src_, satd_size_);
    }
    aom_usec_timer_mark(&c_timer_);

    int total_simd;
    aom_usec_timer simd_timer_;
    aom_usec_timer_start(&simd_timer_);

    for (int i = 0; i < numIter; i++) {
      total_simd = satd_func_simd_(src_, satd_size_);
    }
    aom_usec_timer_mark(&simd_timer_);

    const int c_sum_time = static_cast<int>(aom_usec_timer_elapsed(&c_timer_));
    const int simd_sum_time =
        static_cast<int>(aom_usec_timer_elapsed(&simd_timer_));

    printf(
        "c_time = %d \t simd_time = %d \t Gain = %4.2f \n", c_sum_time,
        simd_sum_time,
        (static_cast<float>(c_sum_time) / static_cast<float>(simd_sum_time)));

    EXPECT_EQ(total_ref, total_simd) << "Output mismatch \n";
  }
  int satd_size_;

 private:
  tran_low_t *src_;
  SatdFunc satd_func_ref_;
  SatdFunc satd_func_simd_;
  ACMRandom rnd_;
};

TEST_P(SatdTest, MinValue) {
  const int kMin = -32640;
  const int expected = -kMin * satd_size_;
  FillConstant(kMin);
  Check(expected);
}
TEST_P(SatdTest, MaxValue) {
  const int kMax = 32640;
  const int expected = kMax * satd_size_;
  FillConstant(kMax);
  Check(expected);
}
TEST_P(SatdTest, Random) {
  int expected;
  switch (satd_size_) {
    case 16: expected = 205298; break;
    case 64: expected = 1113950; break;
    case 256: expected = 4268415; break;
    case 1024: expected = 16954082; break;
    default:
      FAIL() << "Invalid satd size (" << satd_size_
             << ") valid: 16/64/256/1024";
  }
  FillRandom();
  Check(expected);
}
TEST_P(SatdTest, Match) {
  FillRandom();
  RunComparison();
}
TEST_P(SatdTest, DISABLED_Speed) {
  FillRandom();
  RunSpeedTest();
}
#if HAVE_NEON
INSTANTIATE_TEST_SUITE_P(
    NEON, SatdTest,
    ::testing::Values(make_tuple(16, &aom_satd_c, &aom_satd_neon),
                      make_tuple(64, &aom_satd_c, &aom_satd_neon),
                      make_tuple(256, &aom_satd_c, &aom_satd_neon),
                      make_tuple(1024, &aom_satd_c, &aom_satd_neon)));
#endif

}  // namespace
