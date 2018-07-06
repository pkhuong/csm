(defpackage #:csm
  (:use #:cl)
  (:export #:csm #:beta-icdf
           #:csm-driver #:csm-power))

(in-package #:csm)

;;; Float frobbing utilities.  Increment/decrement double floats by a
;;; few ULPs to direct rounding either up or down.

#+sbcl
(progn
  (declaim (inline %float-bits %bits-float next prev))
  (defun %float-bits (x)
    "Convert a double float x to sign-extended sign/magnitude, and
     then to 2's complement."
    (declare (type double-float x))
    (let* ((hi (sb-kernel:double-float-high-bits x))
           (lo (sb-kernel:double-float-low-bits x))
           (word (+ (ash (ldb (byte 31 0) hi) 32) lo)))
      ;; hi is the high half of the 64 bit sign-magnitude
      ;; representation… in two's complement.  Extract the significand,
      ;; and then apply the sign bit.  We want to preserve signed zeros,
      ;; so return -1 - word instead of -word.
      ;;
      ;; (- -1 word) = (lognot word) = (logxor word -1).
      (logxor word (ash hi -32))))

  (defun %bits-float (bits)
    "Convert 2's complement to sign-extended sign/magnitude, then
     double float."
    (declare (type (signed-byte 64) bits))
    ;; convert back to sign-magnitude: if bits is negative, all but the
    ;; sign bit must be flipped again.
    (let ((bits (logxor bits
                        (ldb (byte 63 0) (ash bits -64)))))
      (sb-kernel:make-double-float (ash bits -32)
                                   (ldb (byte 32 0) bits))))

  (defun next (x &optional (delta 1))
    "Increment x by delta ULPs."
    (declare (type double-float x)
             (type unsigned-byte delta))
    (%bits-float (+ (%float-bits x) delta)))

  (defun prev (x &optional (delta 1))
    "Decrement x by delta ULPs."
    (declare (type double-float x)
             (type unsigned-byte delta))
    (%bits-float (- (%float-bits x) delta))))

#-sbcl
(progn
  (declaim (inline next prev))
  (defun next (x &optional (delta 1))
    "Increment x by delta ULPs. Very conservative for
     small (0/denormalised) values."
    (declare (type double-float x)
             (type unsigned-byte delta))
    (let* ((exponent (nth-value 1 (decode-float x)))
           (ulp (max (scale-float double-float-epsilon exponent)
                     least-positive-normalized-double-float)))
      (+ x (* delta ulp))))

  (defun prev (x &optional (delta 1))
    "Decrement x by delta ULPs. Very conservative for small values."
    (declare (type double-float x)
             (type unsigned-byte delta))
    (let* ((exponent (nth-value 1 (decode-float x)))
           (ulp (max (scale-float double-float-epsilon exponent)
                     least-positive-normalized-double-float)))
      (- x (* delta ulp)))))

;;; Wrappers for log/log1p with one-sided rounding errors.

(declaim (type (unsigned-byte 31) *libm-error-limit*))
(defvar *libm-error-limit* 4
  "Assume libm is off by less than 4 ULPs.")

(declaim (inline log-up log-down))
(defun log-up (x)
  "Conservative upper bound on log(x)."
  (declare (type double-float x))
  (next (log x) *libm-error-limit*))

(defun log-down (x)
  "Conservative lower bound on log(x)."
  (declare (type double-float x))
  (prev (log x) *libm-error-limit*))

#+sbcl
(progn
  (declaim (inline log1p-up log1p-down))
  (defun log1p-up (x)
    "Convervative upper bound on log(1 + x)."
    (declare (type double-float x))
    (next (sb-kernel:%log1p x) *libm-error-limit*))

  (defun log1p-down (x)
    "Conservative lower bound on log(1 + x)"
    (declare (type double-float x))
    (prev (sb-kernel:%log1p x) *libm-error-limit*)))

#-sbcl
(progn
  (declaim (ftype (function (double-float) (values double-float &optional))
                  log1p-up log1p-down))
  (defun log1p-up (x)
    "Convervative upper bound on log(1 + x)."
    (declare (type double-float x))
    (if (> (abs x) 1d-8)
        ;; Easy case: x is large.
        (next (log (+ 1d0 x)) *libm-error-limit*)
        ;; Harder case: -1d-8 <= x <= 1d-8
        ;; log(1 + x) = x - x^2 / 2 + x^3 / 3 - ....
        ;; x is a valid upper bound… and |x|^3 / 3 << .5d-16 |x|.
        ;;
        ;; For the positive case, assuming we're working with
        ;; IEEE doubles, that's less than the machine epsilon
        ;; for (* .5d0 x x).
        ;; For the negative case, the tail is all negative,
        ;; so any truncated series is an upper bound.
        (min x (next (- x (* .5d0 (prev (* x x))))))))

  (defun log1p-down (x)
    "Conservative lower bound on log(1 + x)."
    (declare (type double-float x))
    (if (> (abs x) 1d-8)
        (prev (log (+ 1d0 x)) *libm-error-limit*)
        ;; log(1 + x) = x - x^2 / 2 + x^3 / 3 ...
        ;; x - x^2 / 2 is a valid lower bound if x >= 0.
        ;;
        ;; For the negative case, we perform the same bound
        ;; analysis as for log1p-up, and see that we can
        ;; bound the tail by the epsilon for (* .5d0 x x).
        (prev (- x (* .5d0 (next (* x x))))))))

;;; Kahan-style summation.
;;;
;;; Represent the accumulator as an evaluated sum of two doubles.  As
;;; long as the compensation term is initially 0, the result is a safe
;;; upper bound on the real value, and the two terms are
;;; "non-overlapping."  For more details, see "Adaptive Precision
;;; Floating-Point Arithmetic and Fast Robust Geometric Predicates",
;;; Shewchuk, 1997; Technical report CMU-CS-96-140R / Discrete & Comp
;;; Geom 18(3), October 1997.  Theorem 6 in particular.

(declaim (inline sum-update-up sum-update-finish))
(defun sum-update-up (accumulator compensation term &optional ordered)
  "Given an evaluated sum
     (accumulator + compensation),
   return a new unevaluated sum for an upper bound on
     (accumulator + compensation + term).

   If ordered, assume
     term < accumulator,
   or
     accumulator = compensation = 0."
  (declare (type double-float accumulator compensation
                 term))
  (when (and (not ordered)
             (< (abs accumulator) (abs term)))
    (rotatef accumulator term))
  (let* ((rest-1 (next (+ compensation term))) ; safe upper bound on c + t
         (rest (if (<= compensation 0d0)       ; tighter, still safe.
                   (min term rest-1)
                   rest-1))
         ;; Perform a Dekker sum of accumulator + rest.  The result is
         ;; exact, so no need for next/prev here.
         ;;
         ;; Precondition: |accumulator| >= |rest| (or accumulator = 0).
         (a accumulator)
         (b rest)
         (x (+ a b))
         (b-virtual (- x a))     ; b-virtual = value really added to a
         (y (- b b-virtual)))
    (values x y)))

(defun sum-update-finish (accumulator compensation)
  "Return a conservative upper bound for accumulator + compensation.

   In theory, (+ accumulator compensation) is equal to accumulator.
   In practice, it doesn't hurt to do this right.  The second return
   value is the new compensation term (should never be positive)."
  (declare (type double-float accumulator compensation))
  (let* ((raw-sum (next (+ accumulator compensation)))
         (sum (if (> compensation 0d0)
                  raw-sum
                  ;; if compensation <= 0, acc is already an upper
                  ;; bound.
                  (min accumulator raw-sum)))
         (delta (- sum accumulator)))
    (assert (>= delta compensation))
    (values sum (- compensation delta))))

(declaim (ftype (function (&rest double-float)
                          (values double-float double-float &optional))
                sum-up))
(defun sum-up (&rest values)
  "Conservative upper bound for the sum of values, with a Kahan
   summation loop."
  (let ((acc 0d0)
        (err 0d0))
    (dolist (value values (sum-update-finish acc err))
      (setf (values acc err)
            (sum-update-up acc err value)))))

;;; Upper bound for log c(n, s).
;;;
;;; Use Robbins's "A Remark on Stirling's Formula," The American
;;; Mathematical Monthly, Vol 62, No 1 (Jan 1955), pp 26-29.
;;; http://www.jstor.org/stable/2308012.
;;;
;;;
;;; \sqrt{2\pi} n^{n + 1/2} exp[-n + 1/(12n + 1)]
;;; < n! <
;;; \sqrt{2\pi} n^{n + 1/2} exp[-n + 1/(12n)]
;;;
;;; to upper bound log c(n, s) = log(n!) - log(s!) - log((n - s)!).

(declaim (type double-float *minus-log-sqrt-2pi*))
(defvar *minus-log-sqrt-2pi* -0.9189385332046727d0
  "Smallest double precision value > -log sqrt(2pi).")

(declaim (ftype (function ((unsigned-byte 49) (unsigned-byte 49))
                          (values double-float double-float &optional))
                robbins-log-choose))
(defun robbins-log-choose (n s)
  "Compute a conservative upper bound on log c(n, s) based on
   Robbins's bounds for k!."
  (check-type n (unsigned-byte 49)) ;; ensure 53 bit arith is exact.
  (check-type s (unsigned-byte 49))
  (assert (<= 0 s n))
  ;; Handle easy cases, where c(n, s) is 1 or n.
  (when (or (= n s)
            (zerop s))
    (return-from robbins-log-choose (values 0d0 0d0)))
  (when (or (= s 1)
            (= s (1- n)))
    (return-from robbins-log-choose (values (log-up (float n 1d0))
                                            0d0)))
  (let* ((n (float n 1d0))
         (s (float s 1d0))
         (n-s (float (- n s) 1d0))
         (l1 (next (* (+ n .5d0) (log-up n)))) ; (+ n .5d0) is exact.
         (l2 (next (- (* (+ s .5d0) (log-down s)))))
         (l3 (next (- (* (+ n-s .5d0) (log-down n-s)))))
         (r1 (next (/ (* 12d0 n))))          ; (* 12d0 n) is exact.
         (r2 (next (- (/ (1+ (* 12d0 s)))))) ; also exact.
         (r3 (next (- (/ (1+ (* 12d0 n-s)))))))
    (sum-up *minus-log-sqrt-2pi*
            l1 l2 l3
            r1 r2 r3)))

;;; Confidence Sequence Method.
;;;
;;; See "A simple method for implementing Monte Carlo tests,"
;;; Ding, Gandy, and Hahn, 2017 (https://arxiv.org/abs/1611.01675).
;;;
;;; Correctness is a direct corollary of Robbins's "Statistical
;;; Methods Related to the Law of the Iterated Logarithm" (Robbins,
;;; Ann. Math. Statist. Vol 41, No 5 (1970), pp 1397-1409.
;;; https://projecteuclid.org/euclid.aoms/1177696786.
;;;
;;; Let { x_i : i \in |N } be a sequence of i.i.d. Bernoulli random
;;; variables with success probability 0 < P(x_i = 1) = p < 1, for
;;; all i.
;;;
;;; Further let S_n = \sum_{i=1}^n x_i, i.e., the number of successes
;;; in the first n terms, and b(n, p, s) = c(n, s) p^s (1 - p)^{n - s}.
;;;
;;; The probability of any n > 0 satisfying
;;;    b(n, p, S_n) < eps / (n + 1)),
;;; for 0 < eps < 1, is less than eps.
;;;
;;; We can thus check whether the inequality above is ever satisfied,
;;; and, when it is, decide that the stream of Bernoullis observed has
;;; P(x_i = 1) != p.
;;;
;;; Ding, Gandy, and Hahn show that we can also expect the empirical
;;; success rate S_n/n to be on the same side of the threshold p (alpha
;;; in this implementation) as the real but unknown success rate
;;; of the i.i.d. Bernoullis.

(declaim (ftype (function ((unsigned-byte 49)
                           (real (0) (1))
                           (unsigned-byte 49)
                           real)
                          (values boolean double-float &optional))
                csm))
(defun csm (n alpha s log-eps)
  "Given n trials and s sucesses, are we reasonably sure that the
  success rate is *not* alpha (with a false positive rate < exp(log-eps))?

  Answer that question with Ding, Gandy, and Hahn's confidence
  sequence method (CSM). The second return value is an estimate of the
  false positive target rate we would need to stop here.  This value
  should only be used for reporting; the target rate eps should always
  be fixed before starting the experiment."
  (check-type n (unsigned-byte 49))
  (check-type alpha (real (0) (1)))
  (check-type s (unsigned-byte 49))
  (check-type log-eps real)
  (assert (<= 0 s n))
  (let* ((log-choose (robbins-log-choose n s))
         (n (float n 1d0))
         (alpha (float alpha 1d0))
         (s (float s 1d0))
         (log-eps (float log-eps 1d0))
         (log-level (sum-up (log-up (1+ n))
                            log-choose
                            (next (* s (log-up alpha)))
                            (next (* (- n s) (log1p-up (- alpha)))))))
    (values (< log-level log-eps) log-level)))

;;; Beta confidence intervals.
;;;
;;; Approximate the CDF of the Beta(a, b), the regularised incomplete
;;; Beta function I_x(a, b) with an upper bound based on the
;;; hypergeometric representation
;;;
;;;   I_x(a, b) = [Gamma(a + b)/(Gamma(a) Gamma(b)) x^a (1 - x)^b/a] * \
;;;               sum_s=0^\infty [(a + b)_s / (a + 1)_s] x^s,
;;; where
;;;   (a + b)_0 = 1,
;;;   (a + b)_1 = 1 * (a + b) = a + b
;;;   (a + b)_s = 1 * (a + b) * (a + b + 1) * ... * (a + b + s - 1)
;;; and, similarly for (a + 1)_s,
;;;   (a + 1)_0 = 1,
;;;   (a + 1)_1 = 1 * (a + 1) = a + 1,
;;;   (a + 1)_s = 1 * (a + 1) * (a + 2) * ... * (a + s).
;;;
;;; The summands [(a + b)_s / (a + 1)_s] x^s can thus be reformulated
;;; as
;;;  \pi_s(a, b, x) := [(a + b)_s / (a + 1)_s] x^s
;;;                 = \prod_i=1^s [(a + b - 1 + i) / (a + i)]x
;;;                 = \prod_i=1^s [1 + (b - 1) / (a + i)]x.
;;;
;;; The parameters a and b are positive integers, so we can also
;;; compute
;;;   Gamma(a + b)/(Gamma(a) Gamma(b)) x^a (1 - x)^b/a
;;; as
;;;   c(a + b - 1, a) x^a (1 - x)^b.
;;;
;;; This is a product of very small and very large terms, so we'll
;;; work on log space for that initial value.  Once it's computed, the
;;; summands monotonically approach 0 from above, so we can use normal
;;; arithmetic.  We can also easily overapproximate every intermediate
;;; value, starting with Robbins's approximation for
;;; log(c(n, s)) = log(c(a + b - 1, a)).
;;;
;;; This series of products representation lets us compute upper and
;;; lower bounds for the tail of a partial sum, by over- and under-
;;; approximating the tail with geometric series
;;;   \pi_s(a, b, x) \sum_j=1^\infty x^j
;;;  < \sum_j=1^\infty \pi_{s + j}(a, b, c) <
;;;   \pi_s(a, b, x) \sum_j=1^\infty \pi_s(a, b, x)^j
;;;
;;; and thus
;;;
;;;   \pi_s(a, b, x) [1 / (1 - x) - 1]
;;;  < \sum_j=1^\infty \pi_{s + j}(a, b, c) <
;;;   \pi_s(a, b, x) [1 / (1 - \pi_s(a, b, x)) - 1].
;;;
;;; Given conservative comparisons between threshold and the limits
;;; for our one-sided over-approximation of I_x(a, b), we can execute
;;; a bisection search and invert the over-approximation.  The result
;;; is a conservative lower bound for the confidence interval on the
;;; actual Beta CDF I_x(a, b).

(deftype %count ()
  "Assume we only have at most 2^44 events (successes or failures)."
  `(integer 1 ,(ash 1 44)))

(declaim (ftype (function (%count
                           %count
                           (double-float (0d0) (1d0))
                           double-float
                           (or null (and unsigned-byte fixnum)))
                          (values (or null double-float) &optional))))
(defun %incbeta (a b x threshold limit)
  "Iteratively evaluate I_x(a, b) with a hypergeometric representation.

   Stop with an approximate value as soon as we know if I_x(a, b) is
   less than or greater than threshold; the approximate valus is on the
   same side of threshold.

   If the iteration limit is reached, return NIL."
  (declare (type %count a b)
           (type (double-float (0d0) (1d0)) x)
           (type double-float threshold)
           (type (or null (and unsigned-byte fixnum)) limit))
  ;; Correctness criterion for (quick) convergence.
  (assert (< x (/ (float a 1d0) (+ a b))))
  (let* ((limit (or limit (min most-positive-fixnum
                               (* 10 (+ a b 1000)))))
         (log-initial (sum-up (robbins-log-choose (+ a b -1) a)
                              (next (* a (log-up x)))
                              (next (* b (log1p-up (- x))))))
         (b-1 (float (1- b) 1d0))
         ;; Keep a running product for the summands.
         (product (next (exp log-initial) *libm-error-limit*))
         ;; Use round-up Kahan summation for the series.
         (acc product)
         (err 0d0))
    (declare (type (and unsigned-byte fixnum) limit)
             (type double-float product acc err))
    (labels ((maybe-terminate-slow (multiplicand)
               "Check for tighter termination conditions, by bounding
                the tail of the series."
               (declare (type double-float multiplicand))
               ;; Compute upper and lower bounds for the tail of the
               ;; series with the geometric series limit.  There is
               ;; some error here (if only because the loop is
               ;; conservative), so only return when we're strongly
               ;; bounded away from threshold.
               ;;
               ;; \sum_i=1^\infty r^i = [\sum_i=0^\infty r^i] - 1
               ;;                     = [1 / (1 - r)] - 1
               ;;                     = r / (1 - r)
               (let ((tail-hi (* product
                                 (exp (- (log-up multiplicand)
                                         (log1p-down (- multiplicand))))))
                     (tail-lo (* product
                                 (exp (- (log-down x)
                                         (log1p-up (- x))))))
                     (delta (- (- threshold acc) err)))
                 (when (> tail-lo (* 2d0 delta))
                   (return-from %incbeta (max (+ acc tail-lo)
                                              threshold)))
                 (when (< tail-hi (* .5d0 delta))
                   (return-from %incbeta acc))))
             (maybe-terminate (multiplicand slow)
               "Check for termination; if slow, perform tighter but
                slower checks as well."
               (declare (type double-float multiplicand))
               (when (> acc threshold)  ; |err| is < 1 ulp of acc, so
                 (return-from %incbeta acc)) ; acc + err > threshold.
               ;; Next tests are tighter, but slower.  Only do them from
               ;; time to time.
               (when slow
                 (maybe-terminate-slow multiplicand))))
      (loop for i from 1 upto (* 128 (ceiling limit 128)) do
           (let* ((ratio (next (/ b-1 (+ a i))))
                  (multiplicand (min (next (* x (next (1+ ratio))))
                                     1d0))
                  (old-acc acc))
             (setf product (next (* product multiplicand))
                   (values acc err) (sum-update-up acc err product))
             ;; Iteration count is rounded up to 128, so this will
             ;; always trigger on the last iteration.
             (maybe-terminate multiplicand (or (= old-acc acc)
                                               (zerop (mod i 128)))))))))

(declaim (type (double-float (0d0)) *beta-icdf-eps* *beta-icdf-goal*))
(defvar *beta-icdf-eps* 1d-10
  "Always stop when we have the inverse up to this precision.")

(defvar *beta-icdf-goal* 1d-3
  "Let the search stop as soon as the inverse is up to that precision.")

(declaim (ftype (function (%count %count (double-float * (0.5d0)))
                          (values double-float double-float &optional))
                %beta-icdf-lo))
(defun %beta-icdf-lo (a b alpha)
  "Maximise x s.t I_x(a, b) < alpha < 0.5. Assume x \in (0, a / (a + b)).

   We compute an interval that contains a conservative lower bound on x."
  (check-type a %count)
  (check-type b %count)
  (check-type alpha (double-float * (0.5d0)))
  (let* ((alpha (float alpha 1d0))
         (p (/ (float a 1d0) (+ a b)))
         (lo 0d0)
         (hi p))
    (when (<= alpha 0d0)
      (return-from %beta-icdf-lo (values 0d0 0d0)))
    (flet ((update (x)
             (assert (< 0 x p))
             (let* ((goodish (< hi (+ lo (* lo *beta-icdf-goal*))))
                    (limit (and goodish 1000))
                    (px (%incbeta a b x alpha limit)))
               (cond ((and (not px) goodish)
                      (return-from %beta-icdf-lo (values lo hi)))
                     ((and px (< px alpha))
                      (setf lo x))
                     (t
                      ;; Always safe to assume I_x(a, b) >= alpha.
                      (setf hi x))))))
      (loop while (and (> hi (+ lo (* lo *beta-icdf-eps*)))
                       (> hi *beta-icdf-eps*))
         do (update (* .5 (+ hi lo)))
         finally (return (values lo hi))))))

(declaim (ftype (function (%count %count (real 0)
                                  &optional t)
                          (values (double-float 0d0 1d0)
                                  (double-float 0d0 1d0)
                                  &optional))
                beta-icdf))
(defun beta-icdf (a b alpha &optional upper)
  "Compute a conservative lower bound for alpha-level confidence
   interval for a Beta(a, b) distribution.  If upper, compute
   a conservative upper bound for 1 - alpha.

   Assumes alpha is relatively small (at least 0.05); if that's not
   the case, force alpha to 0.05.

   Return the bound value, and an estimate of the (conservative)
   error on that bound."
  (check-type a (integer 1 #. (ash 1 44)))
  (check-type b (integer 1 #. (ash 1 44)))
  (check-type alpha (real 0))
  (when (<= alpha 0)
    (return-from beta-icdf (if upper
                               (values 1d0 0d0)
                               (values 0d0 0d0))))
  (let ((alpha (min (float alpha 1d0) 0.05d0)))
    (if upper
        (multiple-value-bind (lo hi)
            (%beta-icdf-lo b a alpha)
          (values (min (next (- 1 lo)) 1d0)
                  (min (next (- hi lo)) 1d0)))
        (multiple-value-bind (lo hi)
            (%beta-icdf-lo a b alpha)
          (values (max 0d0 lo)
                  (min (next (- hi lo)) 1d0))))))

;;; Basic utilities on top of csm and beta-icdf.

(declaim (ftype (function ((function (unsigned-byte) t)
                           (real (0) (1)) (real (0) (1))
                           &key (:alpha-hi (real (0) (1)))
                                (:min-count (or null real))
                                (:max-count (or null real))
                                (:stream t))
                          (values boolean boolean
                                  double-float unsigned-byte unsigned-byte
                                  double-float double-float))
                csm-driver report))

(defun csm-driver (generator alpha eps
                   &key (alpha-hi alpha) min-count max-count stream)
  "Run a CSM test, given generator(i) -> boolean.  Determines whether
   the success rate for generator is statistically different from
   alpha, with false positive rate < eps.

   Perform at least min-count iterations and at most max-count, if
   provided; log to stream if non-NIL.

   On termination (successful or not), the return values are:

   1. Is the estimated success rate probably different from alpha? (T or NIL)
   2. Is the estimated success rate probably different from alpha-hi? (T or NIL)
   3. The estimated success rate.
   4. The number of observations.
   5. The number of successes
   6. The lower end of 1 - bound-eps CI on the success rate
   7. The upper end of 1 - bound-eps CI on the success rate

   If the value in 1 or 2 is T, the odds that the estimated success
   rate and the actual success rate are on different sides of alpha or
   alpha-hi are less than eps."
  (check-type alpha (real (0) (1)))
  (check-type alpha-hi (real (0) (1)))
  (check-type eps (real (0) (1)))
  (check-type min-count (or null real))
  (check-type max-count (or null real))
  ;; Apply a Bonferroni correction.
  (setf eps (max 0d0 (prev (/ (float eps 1d0)
                              (if (= alpha alpha-hi) 2d0 3d0)))))
  (when (> alpha alpha-hi)
    (rotatef alpha alpha-hi))
  (labels ((bound (s n upper)
             (beta-icdf (1+ s) (1+ (- n s))
                        (* 0.5d0 eps)
                        upper))
           (log-out (s n next-print print-increment log-level
                       &aux (*read-default-float-format* 'double-float))
             (when (< n next-print)
               (return-from log-out (values next-print print-increment)))
             (format stream "~10,D ~,3E ~,3E ~,3E ~,3E~%"
                     n
                     (float (/ s n) 1d0)
                     (bound s n nil) (bound s n t)
                     (/ log-level (- (log 10d0))))
             (force-output stream)
             (let ((next-print (+ next-print print-increment))
                   (next-increment (* 10 print-increment)))
               (values next-print
                       (if (= next-print next-increment)
                           next-increment
                           print-increment))))
           (maybe-terminate (s n stop-alpha stop-hi)
             (when (and (or (not min-count)
                            (>= n min-count))
                        (or stop-alpha
                            stop-hi
                            (and max-count
                                 (>= n max-count))))
               (return-from csm-driver
                 (values stop-alpha stop-hi
                         (float (/ s n) 1d0) n s
                         (bound s n nil) (bound s n t))))))
    (loop
       with log-eps = (log-down eps)
       with next-print = 0
       with print-increment = 10
       with s = 0
       for n upfrom 1 do
         (incf s (if (funcall generator (1- n)) 1 0))
         (multiple-value-bind (stop-alpha log-level-alpha)
             (csm n alpha s log-eps)
           (multiple-value-bind (stop-hi log-level-hi)
               (csm n alpha-hi s log-eps)
             (when stream
               (setf (values next-print print-increment)
                     (log-out s n next-print print-increment
                              (min log-level-alpha log-level-hi))))
             (maybe-terminate s n stop-alpha stop-hi))))))

(defun report (generator alpha eps
               &key (alpha-hi alpha) min-count max-count (stream t))
  (setf alpha (float alpha 0d0)
        eps (float eps 0d0)
        alpha-hi (float alpha-hi 0d0))
  (when (> alpha alpha-hi)
    (rotatef alpha alpha-hi))
  (multiple-value-bind (stop stop-hi p n s ci-lo ci-hi)
      (csm-driver generator alpha eps
                  :alpha-hi alpha-hi
                  :min-count min-count
                  :max-count max-count
                  :stream stream)
    (multiple-value-prog1 (values stop stop-hi p n s ci-lo ci-hi)
      (let ((lower-bound nil)
            (upper-bound nil))
        (flet ((update-lb (x)
                 (when (or (null lower-bound)
                           (> x lower-bound))
                   (setf lower-bound x)))
               (update-ub (x)
                 (when (or (null upper-bound)
                           (< x upper-bound))
                   (setf upper-bound x))))
          (update-lb ci-lo)
          (update-ub ci-hi)
          (when stop
            (if (> p alpha) ;; assume p > alpha; alpha is a new lb.
                (update-lb alpha)
                (update-ub alpha)))
          (when stop-hi
            (if (> p alpha-hi)
                (update-lb alpha-hi)
                (update-ub alpha-hi)))
          (format stream "~
~D iterations, ~D successes (false positive rate < ~,6,,,,,'eG)~%~
success rate p ~~ ~,6,,,,,'eG~%~
confidence interval [~,6,,,,,'eG, ~,6,,,,,'eG]~%"
                  n s eps
                  p
                  ci-lo ci-hi)
          (cond ((and lower-bound upper-bound
                      (<= alpha lower-bound upper-bound alpha-hi))
                 (format stream "~,6,,,,,'eG < p < ~,6,,,,,'eG~%"
                         alpha alpha-hi))
                ((and lower-bound
                      (or (> lower-bound alpha)
                          (> lower-bound alpha-hi)))
                 (format stream "p > ~,6,,,,0,'eG~%"
                         (if (> lower-bound alpha-hi)
                             alpha-hi
                             alpha)))
                ((and upper-bound
                      (or (< upper-bound alpha)
                          (< upper-bound alpha-hi)))
                 (format stream "p < ~,6,,,,,'eG~%"
                         (if (< upper-bound alpha)
                             alpha
                             alpha-hi)))
                ((/= alpha alpha-hi)
                 (format stream "~
No statistically significant difference from ~,6,,,,,'eG or ~,6,,,,,'eG~%"
                         alpha alpha-hi))
                (t
                 (format stream "~
No statistically significant difference from ~,6,,,,,'eG~%"
                         alpha))))))))

(defun csm-power (p alpha max-count
                  &key (alpha-hi alpha) (eps 1d-5) (success-rate 0.99) stream)
  "Estimate the probability of successfully determining that p and
   alpha differ, given max-count iterations and a target false
   positive rate of eps.

   Attempts to determine if the probability is less than or greater than
   the success rate 0.99 by default (with a false positive rate for
   that outer approximation of 1d-9)."
  (let ((max-iter 0))
    (labels ((bernoulli (i)
               (declare (ignore i))
               (< (random 1d0) p))
             (generator (i)
               (declare (ignore i))
               (multiple-value-bind (success success-hi estimate n)
                   (csm-driver #'bernoulli alpha eps
                               :alpha-hi alpha-hi
                               :max-count max-count)
                 (when (or success success-hi)
                   (setf max-iter (max max-iter n)))
                 (let ((correct-alpha (if (< p alpha)
                                          (< estimate alpha)
                                          (> estimate alpha)))
                       (correct-hi (if (< p alpha-hi)
                                       (< estimate alpha-hi)
                                       (> estimate alpha-hi))))
                   (cond ((and success success-hi)
                          (and correct-alpha correct-hi))
                         (success
                          correct-alpha)
                         (success-hi
                          correct-hi)
                         (t
                          nil))))))
      (multiple-value-prog1
          (report #'generator success-rate 1d-9
                  :min-count 100
                  :stream stream)
        (format stream "max inner iteration count: ~D~%"
                (if (zerop max-iter)
                    max-count
                    max-iter))))))

(defun report-power (alpha alpha-hi max-count
              &key (eps 1d-5) (success-rate 0.99) stream)
  (csm-power (/ (+ alpha alpha-hi) 2)  alpha max-count
             :alpha-hi alpha-hi
             :eps eps
             :success-rate success-rate
             :stream stream))
