/*
 * This file requires libmpfr-dev, and libmpc-dev to be installed 
 * on the system. Compile from MATLAB with:
 *
 *  mex backward_error_pol_mp.c -lgmp -lmpfr -lmpc 
 */
 
#include <mex.h>
#include <mpfr.h>
#include <mpc.h>

void mexFunction(int nlhs, mxArray* plhs[],
		 int nrhs, const mxArray* prhs[])
{
  int n, d, i, j;
  double *mp = NULL, *mxr = NULL, *mxi = NULL, *mw = NULL;
  double *mwr = NULL, *mwi = NULL;

  mpfr_t *y = NULL;
  mpfr_t *p = NULL;
  mpc_t *x = NULL;
  mpc_t *w = NULL;

  mpfr_t tmp1, mpi, tmp2;
  mpc_t  ctmp1, ctmp2;

  if (nrhs != 3) {
    mexErrMsgIdAndTxt("backward_error_pol_mp:number_of_params",
		      "Wrong number of parameters");
		      
  }

  /* 
   * We expect that the parameters are: 
   *  - p: the cofficients of the polynomial, length n + 1
   *  - x: the roots computed by the method, length n
   *  - d: the number of digits to use
   */
  n = mxGetNumberOfElements(prhs[1]);
  d = (int) mxGetScalar(prhs[2]);

  /* We allocate data and copy it from MATLAB into the
   * multiprecision types. */
  y = malloc(sizeof(mpfr_t) * (n + 1));
  p = malloc(sizeof(mpfr_t) * (n + 1));
  x = malloc(sizeof(mpc_t) * n);
  w = malloc(sizeof(mpc_t) * (n + 1));

  mpfr_init2(tmp1, d);
  mpfr_init2(tmp2, d);
  mpfr_init2(mpi,  d); mpfr_const_pi(mpi, MPFR_RNDN);
  mpc_init2(ctmp1, d);
  mpc_init2(ctmp2, d);

  for (i = 0; i < n+1; i++) {
    mpfr_init2(y[i], d);
    mpfr_init2(p[i], d);
    mpc_init2(w[i], d);
    mpc_set_si(w[i], 1, MPFR_RNDN);
  }

  for (i = 0; i < n; i++) {
    mpc_init2(x[i], d);
  }

  /* Copy the data */
  mp = mxGetPr(prhs[0]);
  mxr = mxGetPr(prhs[1]);
  mxi = mxGetPi(prhs[1]);
  for (i = 0; i < n + 1; i++) {
    mpfr_set_d(p[i], mp[i], MPFR_RNDN);
  }
  for (i = 0; i < n; i++) {
    mpfr_set_d(mpc_realref(x[i]), mxr[i], MPFR_RNDN);
    mpfr_set_d(mpc_imagref(x[i]), mxi[i], MPFR_RNDN);    
  }  

  /* Create Chebyshev points */ 
  for (i = 0; i < n + 1; i++) {
    mpfr_set_si (tmp1, n - i, MPFR_RNDN);
    mpfr_div_si (tmp1, tmp1, n, MPFR_RNDN);
    mpfr_mul (tmp1, tmp1, mpi, MPFR_RNDN);
    mpfr_cos (y[i], tmp1, MPFR_RNDN);
  }

  /* Evaluate the polynomial */
  for (i = 0; i < n+1; i++) {
    for (j = 0; j < n; j++) {
      mpc_set_fr(ctmp1, y[i], MPFR_RNDN);
      mpc_sub(ctmp1, ctmp1, x[j], MPFR_RNDN);
      mpc_mul (w[i], w[i], ctmp1, MPFR_RNDN);
    }
  }

  /* Copy the vector in the output */
  plhs[0] = mxCreateDoubleMatrix(n+1, 1, mxCOMPLEX);
  mwr = mxGetPr(plhs[0]);
  mwi = mxGetPi(plhs[0]);

  /* Compute the norm of the vector */
  mpfr_set_si (tmp1, 0, MPFR_RNDN);
  for (i = 0; i < n + 1; i++) {
    mpc_abs(tmp2, w[i], MPFR_RNDN);
    mpfr_mul (tmp2, tmp2, tmp2, MPFR_RNDN);
    mpfr_add (tmp1, tmp1, tmp2, MPFR_RNDN);
  }
  mpfr_sqrt(tmp1, tmp1, MPFR_RNDN);

  for (i = 0; i < n+1; i++) {
    mpc_div_fr (w[i], w[i], tmp1, MPFR_RNDN);
    mwr[i] = mpfr_get_d(mpc_realref(w[i]), MPFR_RNDN);
    mwi[i] = mpfr_get_d(mpc_imagref(w[i]), MPFR_RNDN);
  }  

  /* Release all resources */
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  mpfr_clear(mpi);
  mpc_clear(ctmp1);
  mpc_clear(ctmp2);
  
  for (i = 0; i < n+1; i++) {
    mpfr_clear(y[i]);
    mpfr_clear(p[i]);
    mpc_clear(w[i]);
  }

  for (i = 0; i < n; i++) {
    mpc_clear(x[i]);
  }  

  free(y);
  free(p);
  free(x);
  free(w);

}
