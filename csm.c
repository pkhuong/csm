#include "csm.h"

#ifdef TEST_CSM
/* The internal tests don't work without assert. */
# undef NDEBUG
#endif

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/**
 * Inline unit test hack.
 */
#ifdef TEST_CSM
struct test_case {
	const char *name;
	void (*fn)(void);
	struct test_case *next;
};

static struct test_case *csm_test_cases;

#define TEST(NAME)							\
	__attribute__((__unused__)) static void TEST_##NAME(void);	\
									\
	__attribute__((__constructor__)) static void PUSH_##NAME(void)	\
	{								\
		static struct test_case test_case = {			\
			.name = #NAME "...",				\
			.fn = TEST_##NAME				\
		};							\
									\
		test_case.next = csm_test_cases;			\
		csm_test_cases = &test_case;				\
	}								\
									\
	static void TEST_##NAME(void)
#else /* defined(TEST_CSM) */
#define TEST(NAME) __attribute__((unused)) static void TEST_##NAME(void)
#endif

#define NO_TEST(NAME) __attribute__((unused)) static void NO_TEST_##NAME(void)

#ifdef NDEBUG
# undef assert
# define assert(X) do {} while (0 && (X))
#endif

/**
 * Float frobbing utilities, mostly to inc/dec doubles by ULPs.
 *
 * next/prev do something a bit strange around zeros: the double
 * immediately preceding +0 is -0, and the double after -0 is +0.
 *
 * log only returns 0.0 when it's called on exactly 1, so this is
 * safe.  Same for log1p.  Neither ever returns a negative zero.
 *
 * For others, we have to make sure the negative zeros are never
 * generated, or only if exact.  That's mostly an issue for next,
 * which only advances by 1 ULP (prev is not currently used).
 *
 * We can check the code under default rounding rules (round to
 * nearest, break to even).
 *
 * For next(a + b), a + b = -0.0 iff a = b = -0.0 (unless rounding to
 * -infty).  That addition loses no bit, so `next` isn't even
 * necessary to get an upper bound.  0.0 is definitely a valid upper
 * bound for -0.0 or 0.0.
 *
 * We have some instance of next(a * b), but |ab| >> 0, or we can only
 * get a * b = 0 if it's an exact zero (e.g., a = 0).  incbeta is more
 * complicated, but the products that might cancel to zero are always
 * positive, so they won't generate negative zeros.
 *
 * Similarly, for next(+/-1 / x), |x| is finite, so |1 / x| >> 0.
 */

/* We don't need fmin's NaN handling. */
static inline double
min(double x, double y)
{
	return (x < y) ? x : y;
}

static inline double
max(double x, double y)
{
	return (x > y) ? x : y;
}

static inline uint64_t
float_bits(double x)
{
	uint64_t bits;
	uint64_t mask;

	memcpy(&bits, &x, sizeof(bits));
	/* extract the sign bit. */
	mask = (int64_t)bits >> 63;
	/*
	 * If negative, flip the significand bits to convert from
	 * sign-magnitude to 2's complement.
	 */
	return bits ^ (mask >> 1);
}

static inline double
bits_float(uint64_t bits)
{
	double ret;
	uint64_t mask;

	mask = (int64_t)bits >> 63;
	/* Undo the bit-flipping above. */
	bits ^= (mask >> 1);
	memcpy(&ret, &bits, sizeof(ret));
	return ret;
}

static inline double
next_k(double x, uint64_t delta)
{
	return bits_float(float_bits(x) + delta);
}

__attribute__((__unused__)) static inline double
next(double x)
{
	return next_k(x, 1);
}

static inline double
prev_k(double x, uint64_t delta)
{
	return bits_float(float_bits(x) - delta);
}

__attribute__((__unused__)) static inline double
prev(double x)
{
	return prev_k(x, 1);
}

TEST(float_bits)
{
	assert((int64_t)float_bits(1.0) > 0);
	assert((int64_t)float_bits(-1.0) < 0);

	assert(float_bits(0.0) == 0);
	assert(float_bits(-0.0) == -1ULL);

	assert(float_bits(1.0) == 4607182418800017408ULL);
	assert(float_bits(-2.5) == -4612811918334230529ULL);
	assert(-float_bits(M_PI) - 1 == float_bits(-M_PI));
}

TEST(bits_float)
{
	assert(bits_float(float_bits(1.0)) == 1.0);
	assert(bits_float(float_bits(-M_PI)) == -M_PI);

	assert(bits_float(0) == 0.0);
	assert(copysign(1.0, bits_float(0)) == 1.0);
	assert(bits_float(-1ULL) == -0.0);
	assert(copysign(1.0, bits_float(-1ULL)) == -1.0);

	assert(bits_float(4607182418800017408ULL) == 1.0);
	assert(bits_float(-4612811918334230529ULL) == -2.5);
	assert(bits_float(4607182418800017409ULL) == 1.0000000000000002);
}

TEST(next)
{
	assert(next(-0.0) == 0.0);
	assert(copysign(1.0, next(-0.0)) == 1.0);
	assert(next_k(-0.0, 2) > 0.0);

	assert(next(0.0) == 5e-324);
	assert(next_k(-0.0, 2) == 5e-324);

	assert(next(1.0) == 1.0000000000000002);
	assert(next(-1.0) == -0.9999999999999999);

	assert(-1.0 < next(-1.0));
	assert(next(-1.0) < next_k(-1.0, 2));
	assert(next_k(1.0, 2) < 1 + 1e-15);

	assert(0 < next(0.0));
	assert(next(0.0) < next_k(0.0, 2));
	assert(next_k(0.0, 2) < 1e-15);

	assert(1 < next(1.0));
	assert(next(1.0) < next_k(1.0, 2));
	assert(next_k(1.0, 2) < 1 + 1e-15);

	assert(M_PI == next_k(M_PI, 0));
	assert(next_k(-M_PI, 2) == next(next(-M_PI)));
}

TEST(prev)
{
	assert(prev(0.0) == -0.0);
	assert(copysign(1.0, prev(0.0)) == -1.0);

	assert(prev_k(0.0, 2) < -0.0);

	assert(prev(-0.0) == -5e-324);
	assert(prev(1.0) == 0.9999999999999999);
	assert(prev(-1.0) == -1.0000000000000002);

	assert(-1.0 > prev(-1.0));
	assert(prev(-1.0) > prev_k(-1.0, 2));
	assert(prev_k(-1.0, 2) > -1 - 1e-15);

	assert(-0.0 > prev(-0.0));
	assert(prev(-0.0) > prev_k(-0.0, 2));
	assert(prev_k(-0.0, 2) > -1e-15);

	assert(1.0 > prev(1.0));
	assert(prev(1.0) > prev_k(1.0, 2));
	assert(prev_k(1.0, 2) > 1 - 1e-15);

	assert(M_PI == prev_k(M_PI, 0));
	assert(prev_k(-M_PI, 2) == prev(prev(-M_PI)));
}

/**
 * Wrappers for log, log1p with one-sided errors.
 */

/* Assume libm is off by < 4 ULPs. */
static const uint64_t libm_error_limit = 4;

static inline double
log_up(double x)
{
	return next_k(log(x), libm_error_limit);
}

static inline double
log_down(double x)
{
	return prev_k(log(x), libm_error_limit);
}

static inline double
log1p_up(double x)
{
	return next_k(log1p(x), libm_error_limit);
}

static inline double
log1p_down(double x)
{
	return prev_k(log1p(x), libm_error_limit);
}

TEST(log_up)
{
	assert(log_up(1.0) > 0.0);
	assert(log_up(1.0) < 1e-10);

	assert(log_up(1e-10) > log(1e-10));
	assert(log_up(1e-10) - log(1e-10) < 5e-13);

	assert(log_up(4) > log(4));
	assert(log_up(4) - log(4) < 1e-10);

	assert(log_up(1.0) == 2e-323);
	assert(log_up(0.1) == -2.3025850929940437);
	assert(log_up(20.0) == 2.9957322735539926);
}

TEST(log_down)
{
	assert(log_down(1.0) < 0.0);
	assert(log_down(1.0) > -1e-10);

	assert(log_down(1e-10) < log(1e-10));
	assert(log_down(1e-10) - log(1e-10) > -5e-13);

	assert(log_down(1.0) == -1.5e-323);
	assert(log_down(0.1) == -2.3025850929940472);
	assert(log_down(20.0) == 2.995732273553989);
}

TEST(log1p_up)
{
	assert(log1p_up(0.0) > 0.0);
	assert(log1p_up(0.0) < 1e-16);

	assert(log1p_up(-1e-10) > log1p(-1e-10));
	assert(log1p_up(-1e-10) - log1p_up(-1e-10) < 1e-20);

	assert(log1p_up(0.5) > log1p(0.5));
	assert(log1p_up(0.5) - log1p(0.5) < 1e-10);

	assert(log1p_up(0.0) == 2e-323);
	assert(log1p_up(-0.1) == -0.10536051565782625);
	assert(log1p_up(1e-4) == 9.99950003333084e-05);
}

TEST(log1p_down)
{
	assert(log1p_down(0.0) < 0.0);
	assert(log1p_down(0.0) > -1e-16);

	assert(log1p_down(-1e-10) < log1p(-1e-10));
	assert(log1p_down(-1e-10) - log1p(-1e-10) > -1e-20);

	assert(log1p_down(0.5) < log1p(0.5));
	assert(log1p_down(0.5) - log1p(0.5) > -1e-10);

	assert(log1p_down(0.0) == -1.5e-323);
	assert(log1p_down(-0.1) == -0.10536051565782636);
	assert(log1p_down(1e-4) == 9.999500033330829e-05);
}

/**
 * Kahan-style summations, with one-sided error.
 *
 * Represent the accumulator as an evaluated sum of two doubles.  As
 * long as the compensation term is initially 0, the result is a safe
 * upper bound on the real value, and the two terms are
 * "non-overlapping."  For more details, see "Adaptive Precision
 * Floating-Point Arithmetic and Fast Robust Geometric Predicates",
 * Shewchuk, 1997; Technical report CMU-CS-96-140R / Discrete & Comp
 * Geom 18(3), October 1997.  Theorem 6 in particular.
 */

struct sum {
	double acc;
	double err;
};

static double
sum_update_up(struct sum *sum, double term, bool ordered)
{
	const double acc = sum->acc, err = sum->err;
	double a, b;

	if (ordered || fabs(acc) >= fabs(term)) {
		a = acc;
		b = term;
	} else {
		a = term;
		b = acc;
	}

	{
		double shifted = next(b + err);

		b = min(shifted, (err <= 0) ? b : shifted);
	}

	/* |a| >= |b|, or a == 0. */
	{
		const double x = a + b;
		const double b_virtual = x - a;
		const double y = b - b_virtual;

		sum->acc = x;
		sum->err = y;
		return x;
	}
}

static inline double
sum_update_finish(struct sum sum)
{
	double result = next(sum.acc + sum.err);

	return min(result, (sum.err > 0) ? result : sum.acc);
}

static double
sum_up(const double *values, size_t n)
{
	struct sum sum = { 0.0, 0.0 };

	for (size_t i = 0; i < n; i++) {
		sum_update_up(&sum, values[i], false);
	}

	return sum_update_finish(sum);
}

#define SUM(X) (sum_up(&((X)[0]), sizeof((X)) / sizeof(double)))

TEST(sum_update_up)
{
	{
		struct sum sum = { 1.0, 0.0 };

		sum_update_up(&sum, 0.5, true);
		assert(sum.acc == 1.5);
		assert(sum.err == 0.0);
	}

	{
		struct sum sum = { 3.0, 1e-16 };

		sum_update_up(&sum, M_PI, false);
		assert(sum.acc == 6.141592653589793);
		assert(sum.err == 4.440892098500626e-16);
	}

	{
		struct sum sum = { 10e10, 1e-8 };

		sum_update_up(&sum, M_PI * 1e-10, true);
		assert(sum.acc == 100000000000.0);
		assert(sum.err == 1.0314159265358981e-08);
	}

	{
		struct sum sum = { M_PI * 1e-10, 1e-8 };

		sum_update_up(&sum, 10e10, false);
		assert(sum.acc == 10e10);
		assert(sum.err == 1.0314159265358981e-08);
	}
}

TEST(sum_update_finish)
{
	assert(6.141592653589794 == sum_update_finish(
	    (struct sum) { 6.141592653589793, 4.440892098500626e-16 }));

	assert(10e10 == sum_update_finish(
	    (struct sum) { 10e10, -1.0314159265358981e-08 }));

	assert(1.0000000000000002e16 == sum_update_finish(
	    (struct sum) { 1e16, 1 }));
}

TEST(SUM)
{
	{
		double args[] = { 1.0, -2.0, 3.0, -4.0 };

		assert(SUM(args) == -2.0);
	}

	{
		double args[] = { 2.5, 1e16, -1e16 };

		assert(SUM(args) == 4.0);
	}

	{
		double args[] = { 1.0, 0.5, 0.25, 1.0 / 8 };

		assert(SUM(args) + 1.0 / 8 == 2.0);
	}
}

/**
 * Upper bound for log c(n, s).
 *
 * Use Robbins's "A Remark on Stirling's Formula," The American
 * Mathematical Monthly, Vol 62, No 1 (Jan 1955), pp 26-29.
 * http://www.jstor.org/stable/2308012.
 *
 *
 * \sqrt{2\pi} n^{n + 1/2} exp[-n + 1/(12n + 1)]
 * < n! <
 * \sqrt{2\pi} n^{n + 1/2} exp[-n + 1/(12n)]
 *
 * to upper bound log c(n, s) = log(n!) - log(s!) - log((n - s)!).
 */

/* Smallest double value > -log sqrt(2 pi).
 *
 * Should be bitwise identical to scalbn(-8277062471433908.0, -53).
 */
static const double minus_log_sqrt_2pi = -0.9189385332046727;

/**
 * Compute a conservative upper bound for log c(n, s).
 */
static double
robbins_log_choose(uint64_t n, uint64_t s)
{
	assert(n < 1UL << 49);
	assert(s < 1UL << 49);
	assert(s <= n);

	/* Check for easy cases. */
	if (s == 0 || s == n) {
		return 0.0;
	}

	if (s == 1 || s == n - 1) {
		return log_up(n);
	}

	/* Do the real work. */
	{
		uint64_t n_s = n - s;
		const double values[] = {
			minus_log_sqrt_2pi,
			next((n + 0.5) * log_up(n)),
			next(-(s + 0.5) * log_down(s)),
			next(-(n_s + 0.5) * log_down(n_s)),
			next(1.0 / (12 * n)),
			next(-1.0 / (12 * s + 1)),
			next(-1.0 / (12 * n_s + 1))
		};

		return SUM(values);
	}
}

TEST(minus_log_sqrt_2pi)
{
	assert(minus_log_sqrt_2pi == scalbn(-8277062471433908.0, -53));
}

TEST(robbins_log_choose)
{
	assert(robbins_log_choose(5, 5) == 0);
	assert(robbins_log_choose(1, 0) == 0);

	assert(exp(robbins_log_choose(10, 9)) > 10);
	assert(exp(robbins_log_choose(10, 9)) < 10 + 1e-10);

	{
		double expected = log(10.0 * 9 * 8 * 7 * 6 / (5 * 4 * 3 * 2));

		assert(robbins_log_choose(10, 5) > expected);
		assert(robbins_log_choose(10, 5) - expected < 1e-2);
	}

	assert(robbins_log_choose(4, 2) == 1.7944835223684492);
	assert(robbins_log_choose(10000, 100) == 556.7980123668373);
	assert(robbins_log_choose(10000, 8000) == 4999.416373646588);
}

/**
 * Confidence Sequence Method.
 *
 * See "A simple method for implementing Monte Carlo tests,"
 * Ding, Gandy, and Hahn, 2017 (https://arxiv.org/abs/1611.01675).
 *
 * Correctness is a corollary of Robbins's "Statistical Methods
 * Related to the Law of the Iterated Logarithm" (Robbins,
 * Ann. Math. Statist. Vol 41, No 5 (1970), pp 1397-1409.
 * https://projecteuclid.org/euclid.aoms/1177696786.
 *
 * Let { x_i : i \in |N } be a sequence of i.i.d. Bernoulli random
 * variables with success probability 0 < P(x_i = 1) = p < 1, for
 * all i.
 *
 * Further let S_n = \sum_{i=1}^n x_i, i.e., the number of successes
 * in the first n terms, and b(n, p, s) = c(n, s) p^s (1 - p)^{n - s}.
 *
 * The probability of any n > 0 satisfying
 *    b(n, p, S_n) < eps / (n + 1)),
 * for 0 < eps < 1, is less than eps.
 *
 * We can thus check whether the inequality above is ever satisfied,
 * and, when it is, decide that the stream of Bernoullis observed has
 * P(x_i = 1) != p.
 *
 * Ding, Gandy, and Hahn show that we can also expect the empirical
 * success rate S_n/n to be on the same side of the threshold p (alpha
 * in this implementation) as the real but unknown success rate
 * of the i.i.d. Bernoullis.
 */

/**
 * Given n trials and s successes, can we conclude that the success rate
 * differs from alpha with exp(log_eps) false positive rate?
 *
 * Output the current log confidence level in OUT_log_level if non-NULL.
 */
int
csm(uint64_t n, double alpha, uint64_t s, double log_eps,
    double *OUT_log_level)
{
	assert(n < 1UL << 40);
	assert(alpha > 0);
	assert(alpha < 1);
	assert(s <= n);
	assert(log_eps <= 0);

	{
		const double values[] = {
			log_up(n + 1),
			robbins_log_choose(n, s),
			next(s * log_up(alpha)),
			next((n - s) * log1p_up(-alpha))
		};
		double log_level = SUM(values);

		if (OUT_log_level != NULL) {
			*OUT_log_level = log_level;
		}

		return log_level < log_eps;
	}
}

TEST(csm)
{
	double level = 0.0;

	assert(csm(1, 0.5, 1, -1.0, NULL) == 0);

	assert(csm(1, 0.5, 0, -1.0, &level) == 0);
	assert(fabs(1.0 - level / 9.992007221626409e-16) < 1e-10);

	assert(csm(100, 0.9, 98, log(1e-5), NULL) == 0);
	assert(csm(100, 0.9, 98, log(0.2), NULL) != 0);

	assert(csm((size_t)1e9, 0.99, (size_t)(0.99 * 1e9), -1e-2, NULL) == 0);

	assert(csm(10000, 0.01, 50, log(1e-9), NULL) == 0);
	assert(csm(10000, 0.01, 50, log(1e-3), NULL) != 0);

	assert(csm(10, 0.05, 1, -10.0, &level) == 0);
	assert(fabs(1.0 - level / 1.243108442750477) < 1e-10);

	assert(csm(100000, 0.05, 100, -10.0, &level) != 0);
	assert(fabs(1.0 - level / -4624.756745998) < 1e-10);

	assert(csm(1000000, 0.99, 990596, -10.0, &level) == 0);
	assert(fabs(1.0 - level / -9.977077184266818) < 1e-10);

	assert(csm(1000000, 0.99, 990597, -10.0, &level) != 0);
	assert(fabs(1.0 - level / -10.039129993485403) < 1e-10);
}

/**
 * Quantile confidence interval.
 *
 * The q-th (e.g., 0.9) quantile is the value V such that q of the
 * observations are below (or equal to) `V`.  In practice, any sample
 * of `n` observations can't expect to find the qth quantile at
 * exactly `qn` (if only because that value will not be integral).
 * However, we can use `csm` to bound the likeliness of the qth
 * quantile being very far from `qn`.  For example, if the quantile
 * were actually at `h > qn`, that would mean that `h` observations
 * were below (or equal to) the actual quantile `V`, where the actual
 * occurrence rate should be `q`.
 *
 * If `V` had a measure of 0, we could use binary search to find the
 * most extreme indices that pass the CSM test.	 However, for discrete
 * distributions, we must go out one more, and find the least extreme
 * indices that *do not pass* CSM: `V` might be exactly equal to that
 * out-by-one-more index.  However, `V` cannot be even more extreme:
 * in that case, we'd have enough values strictly less than or
 * strictly greater than `V` to reject the null hypothesis (i.e., the
 * probability of that happening is low enough to satisfy `log_eps`).
 */

/**
 * Binary searches for the least index for the quantile that violates
 * the CSM hypothesis, i.e, the maximum index that is so low that the
 * CSM tells us we'll never (modulo `log_eps`) observe that.
 *
 * The value at that index might not actually violate the CSM CI (if
 * it's exactly equal to the quantile value), but any quantile value
 * strictly less than the one at that index definitely would.
 */
static uint64_t
quantile_index_lo(uint64_t n, double quantile, double log_eps)
{
	uint64_t lo = 0;
	uint64_t hi = floor(n * quantile);

	/* If 0 isn't extreme enough, we have too few observations. */
	if (lo >= hi ||
            csm(n, quantile, lo, log_eps, NULL) == 0) {
		return UINT64_MAX;
	}

	/*
	 * Invariant: csm rejects lo and accepts hi.
	 *
	 * If the quantile were lower than the value at lo, we'd only
	 * have observed "lo" values less than that of the quantile,
	 * and that's highly unlikely to happen.
	 */
	while (lo + 1 < hi) {
		const uint64_t mid = lo + (hi - lo) / 2;
		if (csm(n, quantile, mid, log_eps, NULL) == 0) {
			hi = mid;
		} else {
			lo = mid;
		}
	}

	return lo;
}

static uint64_t
quantile_index_hi(uint64_t n, double quantile, double log_eps)
{
	uint64_t lo = ceil(n * quantile);
	uint64_t hi = n;

	/*
	 * If n isn't extreme enough, we have too few observations.
	 */
	if (lo >= hi ||
            csm(n, quantile, hi, log_eps, NULL) == 0) {
		return UINT64_MAX;
	}

	/*
	 * Invariant: csm accepts lo and rejects hi.
	 *
	 * If the quantile were higher than the value at `hi - 1`,
	 * we'd have observed "hi" values less than that of the
	 * quantile, and that's highly unlikely to happen.
	 */
	while (lo + 1 < hi) {
		const uint64_t mid = lo + (hi - lo) / 2;
		if (csm(n, quantile, mid, log_eps, NULL) == 0) {
			lo = mid;
		} else {
			hi = mid;
		}
	}

	return hi - 1;
}

static const double log2_up = -0.6931471805599454;

uint64_t
csm_quantile_index(uint64_t n, double quantile, int direction, double log_eps)
{
	assert(quantile >= 0 && quantile <= 1 &&
	       "The quantile value must be a fraction in [0, 1].");
	if (quantile < 0) {
		quantile = 0.0;
	}
	if (quantile >= 1.0) {
		quantile = 1.0;
	}
        
        if (n <= 0) {
                return UINT64_MAX;
        }

	if (direction == 0) {
		if (quantile == 0.0) {
			return 0;
		}

		if (quantile == 1.0) {
			return n - 1;
		}

		const double estimate = quantile * n;
		const uint64_t floored = (uint64_t)estimate;
		/*
		 * If we're, e.g., taking the 10 percentile of 100,
		 * our estimate is exactly the 10th value, with index
		 * 10 - 1 = 9.
		 */
		if (floored == estimate) {
                        assert(floored > 0);
			return floored - 1;
		}

		return floored;
	}

	/* Find the lower end of the confidence interval. */
	if (direction < 0) {
		if (quantile == 0.0) {
			return UINT64_MAX;
		}

		if (quantile == 1.0) {
			return n - 1;
		}

		return quantile_index_lo(n, quantile, log_eps + log2_up);
	}


	/* Find the upper end of the confidence interval. */
	if (quantile == 0.0) {
		return 0;
	}

	if (quantile == 1.0) {
		return UINT64_MAX;
	}

	return quantile_index_hi(n, quantile, log_eps + log2_up);
}

/**
 * Beta confidence intervals.
 *
 * Approximate the CDF of the Beta(a, b), the regularised incomplete
 * Beta function I_x(a, b) with an upper bound based on the
 * hypergeometric representation
 *
 *   I_x(a, b) = [Gamma(a + b)/(Gamma(a) Gamma(b)) x^a (1 - x)^b/a] * \
 *		 sum_s=0^\infty [(a + b)_s / (a + 1)_s] x^s,
 * where
 *   (a + b)_0 = 1,
 *   (a + b)_1 = 1 * (a + b) = a + b
 *   (a + b)_s = 1 * (a + b) * (a + b + 1) * ... * (a + b + s - 1)
 * and, similarly for (a + 1)_s,
 *   (a + 1)_0 = 1,
 *   (a + 1)_1 = 1 * (a + 1) = a + 1,
 *   (a + 1)_s = 1 * (a + 1) * (a + 2) * ... * (a + s).
 *
 * The summands [(a + b)_s / (a + 1)_s] x^s can thus be reformulated
 * as
 *  \pi_s(a, b, x) := [(a + b)_s / (a + 1)_s] x^s
 *		   = \prod_i=1^s [(a + b - 1 + i) / (a + i)]x
 *		   = \prod_i=1^s [1 + (b - 1) / (a + i)]x.
 *
 * The parameters a and b are positive integers, so we can also
 * compute
 *   Gamma(a + b)/(Gamma(a) Gamma(b)) x^a (1 - x)^b/a
 * as
 *   c(a + b - 1, a) x^a (1 - x)^b.
 *
 * This is a product of very small and very large terms, so we'll
 * work on log space for that initial value.  Once it's computed, the
 * summands monotonically approach 0 from above, so we can use normal
 * arithmetic.	We can also easily overapproximate every intermediate
 * value, starting with Robbins's approximation for
 * log(c(n, s)) = log(c(a + b - 1, a)).
 *
 * This series of products representation lets us compute upper and
 * lower bounds for the tail of a partial sum, by over- and under-
 * approximating the tail with geometric series
 *   \pi_s(a, b, x) \sum_j=1^\infty x^j
 *  < \sum_j=1^\infty \pi_{s + j}(a, b, c) <
 *   \pi_s(a, b, x) \sum_j=1^\infty \pi_s(a, b, x)^j
 *
 * and thus
 *
 *   \pi_s(a, b, x) [1 / (1 - x) - 1]
 *  < \sum_j=1^\infty \pi_{s + j}(a, b, c) <
 *   \pi_s(a, b, x) [1 / (1 - \pi_s(a, b, x)) - 1].
 *
 * Given conservative comparisons between threshold and the limits for
 * our one-sided over-approximation of I_x(a, b), we can execute a
 * bisection search and invert the over-approximation. The result is a
 * conservative lower bound for the confidence interval on the actual
 * Beta CDF I_x(a, b).
 */

static double
incbeta(uint64_t a, uint64_t b, double x, double threshold, uint64_t limit)
{
	assert(0 < a);
	assert(a < 1UL << 44);
	assert(0 < b);
	assert(b < 1UL << 44);
	assert(x < (1.0 * a) / (a + b));

	if (limit == 0) {
		limit = 10 * (a + b + 1000);
	}

	{
		double b_1 = b - 1.0;
		double a_i;
		double product;
		struct sum sum;

		{
			double log_values[] = {
				robbins_log_choose(a + b - 1, a),
				next(a * log_up(x)),
				next(b * log1p_up(-x))
			};
			double log_initial = SUM(log_values);

			product = next_k(exp(log_initial), libm_error_limit);
			sum = (struct sum) { product, 0.0 };
		}

		a_i = a + 1.0;
		for (uint64_t i = 1; i <= limit; i++, a_i++) {
			double ratio = next(b_1 / a_i);
			double mult = min(next(x * next(1 + ratio)), 1.0);
			double old_acc = sum.acc;

			product = next(product * mult);
			/* Values are monotonically decreasing. */
			sum_update_up(&sum, product, true);
			/* Easy termination check. */
			if (sum.acc > threshold) {
				return sum.acc;
			}

			/* Look harder, but not too often. */
			if (old_acc == sum.acc || (i % 128) == 0) {
				double lhi = log_up(mult) - log1p_down(-mult);
				double tail_hi = product * exp(lhi);
				double llo = log_down(x) - log1p_up(-x);
				double tail_lo = product * exp(llo);
				double delta = (threshold - sum.acc) - sum.err;

				if (tail_lo > 2 * delta) {
					return max(sum.acc + tail_lo, threshold);
				}

				if (tail_hi < 0.5 * delta) {
					return sum.acc;
				}
			}
		}
	}

	return nan("");
}

/* Never look for more precision than eps. */
static const double beta_icdf_eps = 1e-10;

/* Consider giving up when we're this precise. */
static const double beta_icdf_goal = 1e-3;

static double
beta_icdf_lo(uint64_t a, uint64_t b, double alpha)
{
	assert(0 < a);
	assert(a < 1UL << 44);
	assert(0 < b);
	assert(b < 1UL << 44);
	assert(alpha < 0.5);

	if (alpha <= 0) {
		return 0;
	}

	/* Bisection search. */
	{
		double lo = 0;
		double hi = (1.0 * a) / (a + b);

		while (hi > max(beta_icdf_eps, lo + beta_icdf_eps * lo)) {
			double x = .5 * (hi + lo);
			bool close_enough = hi < (lo + lo * beta_icdf_goal);
			size_t iter_count = close_enough ? 1024 : 0;
			double px = incbeta(a, b, x, alpha, iter_count);

			if (isnan(px)) {
				if (close_enough) {
					break;
				}

				/* When in doubt, assume the worst. */
				hi = x;
			} else if (px < alpha) {
				lo = x;
			} else {
				hi = x;
			}
		}

		return lo;
	}
}

/**
 * Compute a conservative lower bound for an alpha-level confidence
 * interval for a Beta(a, b) distribution.  If upper > 0, compute a
 * conservative upper bound for 1 - alpha.
 *
 * If alpha > 0.05, compute a 0.05-level bound.
 */
double
beta_icdf(uint64_t a, uint64_t b, double alpha, int direction)
{
	assert(0 < a);
	assert(a < 1UL << 44);
	assert(0 < b);
	assert(b < 1UL << 44);
	assert(alpha >= 0);

	if (alpha <= 0) {
		return (direction > 0) ? 1.0 : 0.0;
	}

	alpha = min(alpha, 0.05);
	if (direction > 0) {
		double lo = beta_icdf_lo(b, a, alpha);

		return min(next(1 - lo), 1.0);
	}

	return beta_icdf_lo(a, b, alpha);
}

TEST(incbeta)
{
	double x;

	x = incbeta(5, 5, 0.001, 1.2558053968507e-13, 0);
	assert(x > 1.2558053968507e-13);
	assert(x - 1.2558053968507e-13 < 2e-15);

	x = incbeta(100, 1000000, 1e-5, 5.425166381479153e-63, 0);
	assert(x > 5.425166381479153e-63);
	assert(x - 5.425166381479153e-63 < 4e-69);

	x = incbeta(10000, 1, 0.999, 4.517334597704867e-5, 0);
	assert(x > 4.517334597704867e-5);
	assert(x - 4.517334597704867e-5 < 5e-17);

	assert(isnan(incbeta(10000, 1, 0.999, 4.517334597704867e-05, 10)));

	assert(incbeta(5, 5, 0.001, 1e-13, 0) > 1e-13);
	assert(incbeta(100, 1000000, 1e-5, 5.5e-63, 0) < 5.5e-63);

	assert(incbeta(5, 5, 0.1, 0.1, 0) == 0.0008914881911461997);
	assert(incbeta(5, 5, 0.001, 1.2558053968507e-13, 0) ==
	    1.2566030059287187e-13);
	assert(incbeta(100, 1000000, 1e-5, 5.425166381479153e-63, 0) ==
	    5.42516983189825e-63);
	assert(incbeta(10000, 1, 0.999, 4.517334597704867e-05, 0) ==
	    4.5173345977071525e-05);
	assert(incbeta(5, 5, 0.001, 1e-13, 0) == 1.2566030059287187e-13);
	assert(incbeta(100, 1000000, 1e-5, 5.5e-63, 0) ==
	    5.425170242086257e-63);
}

TEST(beta_icdf_lo)
{
	double err;

	assert(beta_icdf_lo(4, 4, 0.0) == 0.0);

	err = beta_icdf_lo(5, 5, 0.01) / 0.1709651054824590931911 - 1;
	assert(err > -1e-2);
	assert(err <= 0);

	err = beta_icdf_lo(10000, 10, 0.0001) / 0.9973853034151490539281 - 1;
	assert(err > -1e-4);
	assert(err <= 0);

	err = beta_icdf_lo(100, 1000000, 1e-8) / 5.36173569850883957903e-5 - 1;
	assert(err > -1e-6);
	assert(err <= 0);

	assert(beta_icdf_lo(5, 5, 0.01) == 0.1709400317777181);
	assert(beta_icdf_lo(10000, 10, 0.0001) == 0.9973546960851649);
	assert(beta_icdf_lo(100, 1000000, 1e-8) == 5.361735618111455e-05);
}

TEST(beta_icdf)
{
	double err;

	assert(beta_icdf(4, 4, 0.0, -1) == 0.0);
	assert(beta_icdf(4, 4, 0.0, 1) == 1.0);

	err = beta_icdf(5, 5, 0.01, -1) / 0.1709651054824590931911 - 1.0;
	assert(err > -1e-2);
	assert(err <= 0);

	err = beta_icdf(5, 5, 0.01, 1) / 0.8290348945175409068089 - 1.0;
	assert(err >= 0);
	assert(err < 1e-2);

	err = beta_icdf(10000, 10, 0.0001, -1) / 0.9973853034151490539281 - 1.0;
	assert(err > -1e-4);
	assert(err <= 1.0);

	err = beta_icdf(10000, 10, 0.0001, 1) / 0.9997803648233339942553 - 1.0;
	assert(err >= 0);
	assert(err < 1e-4);

	err = beta_icdf(100, 1000000, 1e-8, -1) / 5.36173569850883957903e-5 - 1;
	assert(err > -1e-6);
	assert(err <= 0);

	err = beta_icdf(100, 1000000, 1e-8, 1) / 1.666077240191560021563e-4 - 1;
	assert(err >= 0);
	assert(err < 2e-2);

	assert(beta_icdf(5, 5, 0.01, 0) == 0.1709400317777181);
	assert(beta_icdf(5, 5, 0.01, 1) == 0.829059968222282);

	assert(beta_icdf(10000, 10, 0.0001, 0) == 0.9973546960851649);
	assert(beta_icdf(10000, 10, 0.0001, 1) == 0.999780366628727);

	assert(beta_icdf(100, 1000000, 1e-8, 0) == 5.361735618111455e-05);
	assert(beta_icdf(100, 1000000, 1e-8, 1) == 0.0001686476860127684);
}

#ifdef TEST_CSM
#include <stdio.h>

int
main()
{
	struct test_case *stack = NULL;
	size_t count = 0;

	/* Reverse the test case stack (sometimes useful). */
	while (csm_test_cases != NULL) {
		struct test_case *current = csm_test_cases;
		struct test_case *next = current->next;

		/* pop current from cases. */
		csm_test_cases = next;

		/* push current to stack. */
		current->next = stack;
		stack = current;
		count++;
	}

	printf("csm.c: %zu test cases.\n", count);
	for (size_t i = 1; stack != NULL; i++, stack = stack->next) {
		printf("\tRunning %3zu / %3zu: %-20s ",
		    i, count, stack->name);
		stack->fn();
		printf("\tsuccess!\n");
	}

	printf("csm.c: successfully executed %zu test cases.\n", count);
	return 0;
}
#endif
