#include <stdlib.h>

static int compare_fcn(const void * pa, const void * pb) {
  double *a = (double*) pa;
  double *b = (double*) pb;
  return *a - *b > 0;
}

void cqr_sort_array_(int * n, double * roots) {
  qsort(roots, *n, sizeof(double), compare_fcn);
}
