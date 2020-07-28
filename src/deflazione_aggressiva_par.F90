! mex FFLAGS="$FFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" deflazione_aggressiva_par.F90 fastqr_par.f90 -llapack -lblas
#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  implicit none
  
  mwPointer plhs(*), prhs(*)
  mwSize mxGetM
  mwPointer mxGetPr, mxGetPi, mxCreateDoubleMatrix
  integer nrhs, nlhs, its
  character(256) :: buffer
  integer n
  complex*16, allocatable :: d(:), beta(:), u(:), v(:)
  double precision :: tmp

  if (nrhs .ne. 4) then
     call mexErrMsgTxt("The function fastqr_f requires 4 arguments")
  end if
  
  n = mxGetM(prhs(1))
  
  allocate(d(n), beta(n-1), u(n), v(n))
  
  call mxCopyPtrToComplex16(mxGetPr(prhs(1)), mxGetPi(prhs(1)), d, n)
  call mxCopyPtrToComplex16(mxGetPr(prhs(2)), mxGetPi(prhs(2)), beta, n-1)	
  call mxCopyPtrToComplex16(mxGetPr(prhs(3)), mxGetPi(prhs(3)), u, n)	
  call mxCopyPtrToComplex16(mxGetPr(prhs(4)), mxGetPi(prhs(4)), v, n)
  
  call aggressive_deflation(n, d, beta, u, v, 20,4)
  
  plhs(1) = mxCreateDoubleMatrix(n, 1, 1)

  
  call mxCopyComplex16ToPtr(d, mxGetPr(plhs(1)), mxGetPi(plhs(1)), n)
  
  deallocate(d, beta, u, v)
end subroutine mexFunction