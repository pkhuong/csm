#include "csm.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <utility>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
std::pair<double, double> EstimateQuantile(double quantile, size_t min_count,
                                           double eps,
                                           std::vector<double> *observations) {
  const size_t lower_index =
      csm_quantile_index(observations->size(), quantile, -1, std::log(eps));
  const size_t upper_index =
      csm_quantile_index(observations->size(), quantile, 1, std::log(eps));

  double low_value = -std::numeric_limits<double>::max();
  double high_value = std::numeric_limits<double>::max();

  if (lower_index < observations->size()) {
    std::nth_element(observations->begin(), observations->begin() + lower_index,
                     observations->end());
    low_value = observations->at(lower_index);
  }

  if (upper_index < observations->size()) {
    std::nth_element(observations->begin(), observations->begin() + upper_index,
                     observations->end());
    high_value = observations->at(upper_index);
  }

  return std::make_pair(low_value, high_value);
}

// Runs the quantile estimation procedure for niter iterations at
// 1-eps confidence level.  Returns whether the actual quantile was in
// range at all iterations.
bool QuantileInRange(const double quantile, const double eps, size_t niter) {
  static const size_t kMinObservations = 32;
  const double expected = 100 * quantile;

  std::random_device dev;
  std::mt19937 rng(dev());

  std::uniform_real_distribution<double> dist(0, 100);
  std::vector<double> observations;

  observations.reserve(niter);
  for (size_t i = 0; i < niter; ++i) {
    observations.push_back(dist(rng));

    if (i >= kMinObservations) {
      double lo, hi;
      std::tie(lo, hi) =
          EstimateQuantile(quantile, kMinObservations, eps, &observations);

      if (expected < lo || expected > hi) {
        std::cout << "fail after " << i + 1 << ": " << lo << " " << hi << "\n";
        return false;
      }
    }
  }

  return true;
}

class QuantileTest : public testing::TestWithParam<double> {};

// Compute a 95%-level quantile for a distribution that happens to be
// uniform [0, 100).  Make sure the range does include the population
// quantile at least 95% of the time.
TEST_P(QuantileTest, IncludedInRange) {
  static const size_t kInnerIter = 20000;
  static const size_t kMaxIter = 10000;

  static const double eps = 2.5e-2;
  const double quantile = GetParam();

#ifndef NDEBUG
  {
    static bool once = true;
    if (once) {
      once = false;
      std::cout << "This test suite needs a few minutes in "
                   "optimized mode. "
                   "Debug mode may take for approximately ever."
                << std::endl;
    }
  }
#endif

  std::cout << "Quantile width (fraction of sample size) for " << quantile
            << " at eps < " << eps << " after " << kInnerIter << " iterations: "
            << (1.0 / kInnerIter) * csm_quantile_index(kInnerIter, quantile, -1,
                                                       std::log(eps)) -
                   quantile
            << ", "
            << (1.0 / kInnerIter) * csm_quantile_index(kInnerIter, quantile, 1,
                                                       std::log(eps)) -
                   quantile
            << std::endl;

  size_t successes = 0;
  for (size_t i = 0; i < kMaxIter; ++i) {
    if (QuantileInRange(quantile, eps, kInnerIter)) {
      ++successes;
    }

    const size_t total = i + 1;
    if ((total % 100) == 0) {
      std::cout << "Current: " << successes << " / " << total << " ("
                << 1.0 * successes / total << ")." << std::endl;
    }

    if (csm(total, 1.0 - eps, successes, std::log(1e-4), nullptr)) {
      std::cout << "Actual rate " << 1.0 * successes / total << "\n";
      EXPECT_LE(1.0 - 1.0 * successes / total, eps)
          << successes << " / " << total;
      return;
    }
  }

  EXPECT_TRUE(false) << "Only " << successes << " successes out of " << kMaxIter
                     << ".";
}

INSTANTIATE_TEST_SUITE_P(QuantileTest, QuantileTest,
                         testing::Values(0.05, 0.5, 0.75, 0.9, 0.99));
}  // namespace
