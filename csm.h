#ifndef CONFIDENCE_SEQUENCE_METHOD_H
#define CONFIDENCE_SEQUENCE_METHOD_H
#include <stdint.h>

int
csm(double *OUT_log_level, uint64_t n, double alpha, uint64_t s, double log_eps);

double
beta_icdf(uint64_t a, uint64_t b, double alpha, int direction);
#endif /* !CONFIDENCE_SEQUENCE_METHOD_H */
