#include "csm.h"

#include <cmath>
#include <random>

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace {
// Generate and test against a Bernoulli with `p`, for `repeat`
// iterations.  Returns true iff CSM (mistakenly) rejects the null.
bool TestNegative(size_t repeat, double p, double eps) {
  std::random_device dev;
  std::mt19937 rng(dev());

  std::bernoulli_distribution var_success(p);

  size_t total = 0;
  size_t success = 0;
  for (size_t i = 0; i < repeat; ++i) {
    ++total;
    if (var_success(rng)) {
      ++success;
    }

    if (csm(total, p, success, std::log(eps), nullptr) != 0) {
      return true;
    }
  }

  return false;
}

// Generate and test for a Bernoulli with p = 2/3.
// Make sure our false positive rate is low enough.
TEST(CsmStat, Negative) {
  const double false_positive_rate = 1e-2;
  const size_t n = 10000;

  size_t failures = 0;
  for (size_t i = 0; i < n; ++i) {
    if (TestNegative(10000, 2.0 / 3, false_positive_rate)) {
      ++failures;
    }
  }

  // We truncate the search, so the atual false positive rate should
  // be much less than our parameters.
  EXPECT_LE(failures, false_positive_rate * n);
}

// Same test as above, but let CSM make the call.
TEST(CsmStat, NegativeRecursive) {
  const double false_positive_rate = 5e-2;
  const size_t n = 10000;

  size_t failures = 0;
  for (size_t i = 0; i < n; ++i) {
    if (TestNegative(10000, 2.0 / 3, false_positive_rate)) {
      ++failures;
    }

    if (csm(i + 1, false_positive_rate, failures, std::log(1e-3), nullptr) !=
        0) {
      EXPECT_LT(1.0 * failures / (i + 1), false_positive_rate)
          << "n=" << i + 1 << " fail=" << failures;
      return;
    }
  }

  EXPECT_TRUE(false) << "Too many iterations. " << failures
                     << " failures after " << n << " iterations.\n";
}

bool TestPositive(size_t repeat, double p, double eps) {
  std::random_device dev;
  std::mt19937 rng(dev());

  std::bernoulli_distribution var_success(p);

  const double test_p = (p < 0.5) ? p + 1e-1 : p - 1e-1;

  size_t total = 0;
  size_t success = 0;
  for (size_t i = 0; i < repeat; ++i) {
    ++total;
    if (var_success(rng)) {
      ++success;
    }

    if (csm(total, test_p, success, std::log(eps), nullptr) != 0) {
      return true;
    }
  }

  return false;
}

// Generate and test for different Bernoullis.  We expect a fair
// amount of true positives.
TEST(CsmStat, Positive) {
  const double false_positive_rate = 1e-1;
  const size_t n = 10000;

  size_t successes = 0;
  for (size_t i = 0; i < n; ++i) {
    if (TestPositive(10000, 2.0 / 3, false_positive_rate)) {
      ++successes;
    }
  }

  // The false positive rate is lax, so we should really hit this.
  EXPECT_GE(successes, 0.99 * n);
}
}  // namespace
