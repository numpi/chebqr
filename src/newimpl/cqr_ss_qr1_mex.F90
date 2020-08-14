! mex cqr_ss_qr1_mex.F90 cqr_ss.f90 -llapack -lblas
#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  implicit none
  
  mwPointer plhs(*), prhs(*)
  mwSize mxGetM, mxGetN
  mwPointer :: mxGetPr, mxGetPi, mxCreateDoubleMatrix, mxCreateDoubleScalar, mxGetData
  integer :: nrhs, nlhs, its, info
  mwSize :: n
  complex*16, allocatable :: d(:), beta(:), u(:), v(:), Q(:,:)
  double precision :: tmp, mxGetScalar
  logical :: schurv = .false., aed = .true.
  integer :: nQ, ldQ
  integer*4 :: data

  if (nrhs .lt. 4 .or. nrhs .gt. 7) then
     call mexErrMsgTxt("The function cqr_ss_qr1_mex requires between 4 and 7 arguments")
  end if
  
  n = mxGetM(prhs(1))

  ! Option parsing, mainly for Schur vectors
  if (nrhs .ge. 5) then
     data = mxGetData(prhs(5))
     aed = data /= 0
  end if

  if (nrhs .ge. 6) then
     data = mxGetData(prhs(6))     
     schurv = data /= 0
  end if

  if (schurv) then
     nQ = mxGetN(prhs(7))
     ldQ = mxGetM(prhs(7))
     allocate(Q(ldQ,nQ))
     call mxCopyPtrToComplex16(mxGetPr(prhs(7)), mxGetPi(prhs(7)), Q, ldQ * nQ)
  end if
  
  allocate(d(n), beta(n-1), u(n), v(n))
  
  call mxCopyPtrToComplex16(mxGetPr(prhs(1)), mxGetPi(prhs(1)), d, n)
  call mxCopyPtrToComplex16(mxGetPr(prhs(2)), mxGetPi(prhs(2)), beta, n-1)	
  call mxCopyPtrToComplex16(mxGetPr(prhs(3)), mxGetPi(prhs(3)), u, n)	
  call mxCopyPtrToComplex16(mxGetPr(prhs(4)), mxGetPi(prhs(4)), v, n)
  
  call cqr_ss_qr1(n, d, beta, u, v, aed, schurv, Q, nQ, ldQ, its, info)
  
  plhs(1) = mxCreateDoubleMatrix(n, 1, 1)

  if (nlhs .ge. 2) then
    tmp = its
    plhs(2) = mxCreateDoubleScalar(tmp)
  end if
  
  call mxCopyComplex16ToPtr(d, mxGetPr(plhs(1)), mxGetPi(plhs(1)), n)

  if (schurv) then
     plhs(3) = mxCreateDoubleMatrix(ldQ, nQ, 1)
     call mxCopyComplex16ToPtr(Q, mxGetPr(plhs(3)), mxGetPi(plhs(3)), ldQ * nQ)
     deallocate(Q)
  end if
  
  deallocate(d, beta, u, v)
end subroutine mexFunction
