!
! This file is part of cqr_eig, a package for computing eigenvalues
! of hermitian-plus-rank-one matrices, with the special focuse of
! Chebyshev rootfinding.
!
! To compile it: 
!   mex cqr_eig.F90 single_shift.f90 -llapack -lblas

#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  implicit none
  
  mwPointer plhs(*), prhs(*)
  mwSize mxGetM, mxGetN, mxGetNumberOfElements
  logical :: mxIsDouble
  mwPointer mxGetPr, mxGetPi, mxCreateDoubleMatrix
  integer nrhs, nlhs, its
  character :: jobbw
  character(256) :: buffer
  integer n, i, k, np, omp_get_max_threads
  complex*16, allocatable :: d(:), beta(:), u(:), v(:)
  double precision :: tmp, mxGetScalar, bw

  ! Check arguments
  if (nrhs .lt. 4 .or. nrhs .gt. 5) then
     call mexErrMsgTxt("The function cqr_eig requires either 4 or 5 arguments")
  end if

  ! Argument 1, 2, 3, and 4 needs to be double matrices,
  ! of sizes n, n-1, n, and n, respectively.
  if (.not. (mxIsDouble(prhs(1)) .and. mxIsDouble(prhs(2)) .and. &
       mxIsDouble(prhs(3)) .and. mxIsDouble(prhs(4)))) then
     call mexErrMsgTxt("Arguments 1, ..., 4 are not double precision.")
  end if

  ! Check to be handling vectors, and not matrices
  do i = 1, 4
     if (mxGetM(prhs(i)) .ne. 1 .and. mxGetN(prhs(i)) .ne. 1) then
        call mexErrMsgTxt("One of the input arguments is not a vector")
     end if
  end do
  
  n = mxGetNumberOfElements(prhs(1))

  ! Check the right dimensions of the input
  if (mxGetNumberOfElements(prhs(2)) .ne. n - 1) then
     call mexErrMsgTxt("The vector beta has the wrong number of elements")
  end if

  if (mxGetNumberOfElements(prhs(3)) .ne. n) then
     call mexErrMsgTxt("The vector u has the wrong number of elements")
  end if

  if (mxGetNumberOfElements(prhs(4)) .ne. n) then
     call mexErrMsgTxt("The vector v has the wrong number of elements")
  end if

  ! Decide the number of processors to use. The user may have specified
  ! a choice, or otherwise we use the value provided by OpenMP.
  if (nrhs .eq. 5) then
     np = INT(mxGetScalar(prhs(5)))

     if (np .lt. 1) then
        write(buffer, *) 'Invalid number of processors specified, resetting to 1'
        call mexPrintf(buffer)
        np = 1
     end if
  else
     np = omp_get_max_threads()
  end if

  ! Compute the optimal size of the aggressive early deflation; this is related to
  ! the number of available processing threads. Note that we may set this too a
  ! large value for the problem at hand, but in that case cqr_eig() will adjust
  ! it automatically
  k = MAX(20, 6 * np)
  
  ! Prepare the allocated memory, copy the data
  allocate(d(n), beta(n-1), u(n), v(n))
  
  call mxCopyPtrToComplex16(mxGetPr(prhs(1)), mxGetPi(prhs(1)), d, n)
  call mxCopyPtrToComplex16(mxGetPr(prhs(2)), mxGetPi(prhs(2)), beta, n-1)	
  call mxCopyPtrToComplex16(mxGetPr(prhs(3)), mxGetPi(prhs(3)), u, n)	
  call mxCopyPtrToComplex16(mxGetPr(prhs(4)), mxGetPi(prhs(4)), v, n)

  if (nlhs .gt. 1) then
    jobbw = 'y'
  end if
  
  call cqr_single_eig(n, d, beta, u, v, k, np, bw, jobbw)

  ! Save the computed eigenvalues  
  plhs(1) = mxCreateDoubleMatrix(n, 1, 1)
  
  call mxCopyComplex16ToPtr(d, mxGetPr(plhs(1)), mxGetPi(plhs(1)), n)

  ! Save the upper bound for the backward error amplification
  if (nlhs .ge. 2) then
    plhs(2) = mxCreateDoubleMatrix(1, 1, 0)
    call mxCopyReal8ToPtr(bw, mxGetPr(plhs(2)), 1)
  endif
  
  deallocate(d, beta, u, v)
end subroutine mexFunction
