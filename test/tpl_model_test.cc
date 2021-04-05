/*
 * Copyright (c) 2021, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */

#include "av1/encoder/tpl_model.h"

#include <cstdlib>

#include "third_party/googletest/src/googletest/include/gtest/gtest.h"

namespace {

double laplace_prob(double q_step, double b, double zero_bin_ratio,
                    int qcoeff) {
  int abs_qcoeff = abs(qcoeff);
  double z0 = fmax(exp(-zero_bin_ratio / 2 * q_step / b), TPL_EPSILON);
  if (abs_qcoeff == 0) {
    double p0 = 1 - z0;
    return p0;
  } else {
    assert(abs_qcoeff > 0);
    double z = fmax(exp(-q_step / b), TPL_EPSILON);
    double p = z0 / 2 * (1 - z) * pow(z, abs_qcoeff - 1);
    return p;
  }
}

TEST(TplModelTest, TransformCoeffEntropyTest1) {
  // Check the consistency between av1_estimate_coeff_entropy() and
  // laplace_prob()
  double b = 1;
  double q_step = 1;
  double zero_bin_ratio = 2;
  for (int qcoeff = -256; qcoeff < 256; ++qcoeff) {
    double rate = av1_estimate_coeff_entropy(q_step, b, zero_bin_ratio, qcoeff);
    double prob = laplace_prob(q_step, b, zero_bin_ratio, qcoeff);
    double ref_rate = -log2(prob);
    EXPECT_DOUBLE_EQ(rate, ref_rate);
  }
}

TEST(TplModelTest, TransformCoeffEntropyTest2) {
  // Check the consistency between av1_estimate_coeff_entropy(), laplace_prob()
  // and av1_laplace_entropy()
  double b = 1;
  double q_step = 1;
  double zero_bin_ratio = 2;
  double est_expected_rate = 0;
  for (int qcoeff = -20; qcoeff < 20; ++qcoeff) {
    double rate = av1_estimate_coeff_entropy(q_step, b, zero_bin_ratio, qcoeff);
    double prob = laplace_prob(q_step, b, zero_bin_ratio, qcoeff);
    est_expected_rate += prob * rate;
  }
  double expected_rate = av1_laplace_entropy(q_step, b, zero_bin_ratio);
  EXPECT_NEAR(expected_rate, est_expected_rate, 0.001);
}

}  // namespace
