#ifndef CONFIDENCE_SEQUENCE_METHOD_H
#define CONFIDENCE_SEQUENCE_METHOD_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
 * Given n trials and s successes, can we conclude that the success
 * rate differs from `alpha` with a total false positive rate of at
 * most `exp(log_eps)`ñ?  The aggressiveness of each call is adjusted
 * such that the total false positive rate for an unbounded stream of
 * tests is at most `exp(log_eps)`.
 *
 * Output the current log confidence level in OUT_log_level if non-NULL.
 */
int
csm(uint64_t n, double alpha, uint64_t s, double log_eps, double *OUT_log_level);

/**
 * Returns a confidence half-interval on the rank of the `quantile` in
 * `n` sorted i.i.d. observations from a given distribution.  The
 * intervals form a confidence sequence: the aggressivenes of each
 * call is adjusted such that the total flase positive rate for an
 * unbounded stream of tests is at most `exp(log_eps)`.  This
 * half-interval already accounts for multiple hypothesis testing when
 * computing both ends of the confidence interval (i.e., the pair of
 * half-intervals include the real quantile with probability
 * `exp(log_eps)`), and is correctly extended for discrete
 * distributions.
 *
 * The `direction` argument determines which end of the confidence
 * interval is computed.  `direction < 0` returns the lower end of the
 * inclusive interval, `direction > 0` the upper end, and `direction
 * = 0` the empirical estimate.
 *
 * Returns UINT64_MAX if the quantile has likely not been observed
 * yet.
 */
uint64_t
csm_quantile_index(uint64_t n, double quantile, int direction, double log_eps);

/**
 * Compute a conservative lower bound for an alpha-level confidence
 * interval for a Beta(a, b) distribution.  If upper > 0, compute a
 * conservative upper bound for 1 - alpha.
 *
 * If alpha > 0.05, compute a 0.05-level bound.
 */
double
beta_icdf(uint64_t a, uint64_t b, double alpha, int direction);

#ifdef __cplusplus
}  /* extern "C" */
#endif
#endif /* !CONFIDENCE_SEQUENCE_METHOD_H */
