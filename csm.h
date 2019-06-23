#ifndef CONFIDENCE_SEQUENCE_METHOD_H
#define CONFIDENCE_SEQUENCE_METHOD_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int
csm(double *OUT_log_level, uint64_t n, double alpha, uint64_t s, double log_eps);

double
beta_icdf(uint64_t a, uint64_t b, double alpha, int direction);
  
#ifdef __cplusplus
}  /* extern "C" */
#endif
#endif /* !CONFIDENCE_SEQUENCE_METHOD_H */
