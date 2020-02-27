! Package ChebQR -- computing roots of smooth functions over [-1, 1]
!
! Authors:
!  Angelo Casulli <angelo.casulli@sns.it>
!  Leonardo Robol <leonardo.robol@unipi.it>
!

!
! SUBROUTINE CQR_ZEROS
!
! This subroutine takes as input a function handle and compute its
! roots over [-1, 1], approximating it with a polynomial at accuracy
! at least EPS.
!
! INPUT PARAMETERS:
!
!  F    Function that takes a DOUBLE PRECISION and returns the same type.
!       The function needs to implement the following interface, which need
!       to be included in the calling program to correctly pass the
!       function handle:
!
!       interface
!           function f(x)
!               double precision :: f, x
!           end function f
!        end interface
!
!  EPS   DOUBLE PRECISION, must be positive real. It gives the accuracy
!        at which F(X) should be approximated.
!
!  ZEROS DOUBLE_PRECISION, DIMENSION(:), ALLOCATABLE.
!        This variable will be allocated as a vector of length N, and will
!        contain the zeros of the given function in [-1, 1].
!
!  N     INTEGER
!        This variable is set as output, and contains the size of the vector
!        ZEROS, which is allocated by this function.
!
! Since this function allocates a vector, it needs to be called with an
! interface, which for conveniente is provided in the chebqr.h file. 
!
! Hence, one may call the function by first including the file in the
! subroutine of interest, using INCLUDE '/path/to/chebqr.h'.
!
subroutine cqr_zeros(f, eps, zeros, n)
  implicit none
  
  integer :: n, k = 25, j, nroots, i
  double precision :: eps, dr, di
  double precision, allocatable :: zeros(:), coeffs(:)

  ! This is set to true if double shift is enabled. As of now,
  ! root extraction is not yet implemented correctly for double
  ! shift, so this needs to be false. 
  logical :: ds = .true.

  ! Working variables for this routine
  double precision, allocatable :: rd(:), ru(:), rv(:), rbeta(:)
  complex(8), allocatable :: d(:), u(:), v(:), beta(:)
  double precision :: tstart, tend

  interface
     function f(x)
       double precision :: f, x
     end function f
  end interface

  interface
     subroutine cqr_interp(f, eps, coeffs, n)
       integer :: n
       double precision :: eps
       double precision, allocatable :: coeffs(:)
       interface
          function f(x)
            double precision :: f, x
          end function f
       end interface
     end subroutine cqr_interp
  end interface

  ! Construct the Chebyshev interpolant on [-1, 1] for
  ! the given function. This call allocates coeffs
  !call cpu_time(tstart)
  call cqr_interp(f, eps, coeffs, n)
  !call cpu_time(tend)
  ! print *, 'Interpolation time = ', tend - tstart

  ! Build the companion linearization in the Chebyshev basis
  ! and compute the eigenvalues.
  if (ds) then
     allocate(rd(n-1), rbeta(n-2), ru(n-1), rv(n-1))
     call cqr_build_colleague_ds(n, coeffs, rd, rbeta, ru, rv)
     call cqr_fastfastqr_ds(n-1, rd, rbeta, ru, rv, k)
  else
     allocate(u(1:n-1), v(1:n-1), d(1:n-1), beta(1:n-2))
     call cqr_build_colleague(n, coeffs, d, beta, u, v)
     call fastfastqr(n-1, d, beta, u, v, k)
  end if

  ! Count the number of real roots, and then actually extract them
  nroots = 0

  if (ds) then
     call cqr_extract_roots_ds(n-1, rd, rbeta, ru, rv, zeros, nroots, eps)
  else
     call cqr_extract_roots(n-1, d, beta, u, v, zeros, nroots, eps)
  end if
  
  if (nroots .gt. 0) then
     allocate(zeros(nroots))

     if (ds) then
        call cqr_extract_roots_ds(n-1, rd, rbeta, ru, rv, zeros, nroots, eps)
     else
        call cqr_extract_roots(n-1, d, beta, u, v, zeros, nroots, eps)
     end if
  end if

  n = nroots

  ! sort the roots
  call cqr_sort_array(n, zeros)

  deallocate(coeffs)
  
  if (ds) then     
     deallocate(rd,ru,rv,rbeta)
  else
     deallocate(d,u,v,beta)
  end if
     
end subroutine cqr_zeros

! Extract roots which are contained in [-1, 1], and real according
! to the given precision.
!
! INPUT PARAMETERS:
!
!  N    INTEGER, number of computed eigenvalues.
!
!  D    COMPLEX(*), DIMENSION(N), diagonal entries of the symmetric
!       part of the Schur form of the colleague matrix.
!
!  BETA COMPLEX(8), DIMENSION(N-1), subdiagonal entries of the symmetric
!       part of the Schur form of the colleague matrix.
!
!  U    COMPLEX(8), DIMENSION(N), left low-rank factor of the rank-1 correction.
!
!  V    COMPLEX(8), DIMENSION(N), right low-rank factor of the rank-1 correction.
!
!  ZEROS DOUBLE PRECISION, DIMENSION(K), on output, this vector contains the real roots
!       in the interval [-1, 1], according to the given precision tol. If the roots are
!       more than K, then only the first K are stored.
!
!  K    INTEGER, on input the size of the vector ZEROS. On output, the total number of
!       real roots inside [-1, 1], according to the given tolerance.
!
!  TOL  DOUBLE PRECISION, the tolerance used to determine if the complex roots are
!       in fact perturbed real roots. 
!
subroutine cqr_extract_roots(n, d, beta, u, v, zeros, k, tol)
  implicit none
  
  integer :: n, k, nroots, j
  complex(8) :: d(n), beta(n-1), u(n), v(n)
  double precision :: zeros(k), dr, di, tol

  nroots = 0
  do j = 1, n
     dr = realpart(d(j))
     di = imagpart(d(j))
     
     if ( (dr .ge. -1.d0) .and. &
          (dr .le.  1.d0) .and. &
          abs(di) .le. 1.0d-12) then
        nroots = nroots + 1

        if (nroots .le. k) then
           zeros(nroots) = dr
        end if
     end if
  end do

  k = nroots;
  
end subroutine cqr_extract_roots

! Extract roots which are contained in [-1, 1], and real according
! to the given precision.
!
! INPUT PARAMETERS:
!
!  N    INTEGER, number of computed eigenvalues.
!
!  D    DOUBLE PRECISION, DIMENSION(N), diagonal entries of the symmetric
!       part of the Schur form of the colleague matrix.
!
!  BETA DOUBLE PRECISION, DIMENSION(N-1), subdiagonal entries of the symmetric
!       part of the Schur form of the colleague matrix.
!
!  U    DOUBLE PRECISION, DIMENSION(N), left low-rank factor of the rank-1 correction.
!
!  V    DOUBLE PRECISION, DIMENSION(N), right low-rank factor of the rank-1 correction.
!
!  ZEROS DOUBLE PRECISION, DIMENSION(K), on output, this vector contains the real roots
!       in the interval [-1, 1], according to the given precision tol. If the roots are
!       more than K, then only the first K are stored.
!
!  K    INTEGER, on input the size of the vector ZEROS. On output, the total number of
!       real roots inside [-1, 1], according to the given tolerance.
!
!  TOL  DOUBLE PRECISION, the tolerance used to determine if the complex roots are
!       in fact perturbed real roots. 
!
subroutine cqr_extract_roots_ds(n, d, beta, u, v, zeros, k, tol)
  implicit none
  
  integer :: n, k, nroots, j, blk_sz
  double precision :: d(n), beta(n-1), u(n), v(n)
  double precision :: zeros(k), dr, di, tol

  nroots = 0
  j = 1

  do while (j .le. n)
     ! First, determine the size of the block: is it 1x1 or 2x2?
     blk_sz = 1
     if (j .lt. n) then
        if (beta(j) .ne. 0.d0) then
           blk_sz = 2
        end if
     end if

     ! Read real and imaginary part of the eigenvalues
     dr = d(j)

     if (blk_sz .eq. 2) then
        di = sqrt(-beta(j)*(beta(j)-u(j+1)*v(j)+v(j+1)*u(j)))
     else
        di = 0.d0
     end if

     if ( (dr .ge. -1.d0) .and. &
          (dr .le.  1.d0) .and. &
          abs(di) .le. 1.0d-12) then
        nroots = nroots + blk_sz

        if (nroots + blk_sz - 1 .le. k) then
           zeros(nroots : nroots + blk_sz - 1) = dr
        end if        
     end if

     j = j + blk_sz
  end do

  k = nroots;
  
end subroutine cqr_extract_roots_ds


subroutine cqr_build_colleague_ds(n, coeffs, d, beta, u, v)
  implicit none

  integer :: n, j
  double precision :: coeffs(n-1)
  double precision :: d(n-1), u(n-1), v(n-1), beta(n-2)

  ! Build the rank 1 correction
  u = 0.d0
  u(1) = 1

  do j = 1, n - 2
     v(j) = -coeffs(n-j) / coeffs(n)
  end do
  v(n-1) = -coeffs(1) / coeffs(n) * sqrt(2.d0)
  v = v / 2

  ! ... and the symmetric part
  d = 0.d0
  beta = .5d0
  beta(n-2) = 1.d0 / sqrt(2.d0)
end subroutine cqr_build_colleague_ds

subroutine cqr_build_colleague(n, coeffs, d, beta, u, v)
  implicit none

  integer :: n, j
  double precision :: coeffs(n-1)
  complex(8) :: d(n-1), u(n-1), v(n-1), beta(n-2)

  ! Build the rank 1 correction
  u = 0.d0
  u(1) = 1

  do j = 1, n - 2
     v(j) = -coeffs(n-j) / coeffs(n)
  end do
  v(n-1) = -coeffs(1) / coeffs(n) * sqrt(2.d0)
  v = v / 2

  ! ... and the symmetric part
  d = 0.d0
  beta = .5d0
  beta(n-2) = 1.d0 / sqrt(2.d0)
  
end subroutine cqr_build_colleague

! Interpolate a given function f over [-1, 1] as a Chebyshev polynomial,
! with accuracy eps. 
subroutine cqr_interp(f, eps, coeffs, n)
  implicit none
  
  double precision :: eps, y, xx, err, yy, tol, nrm
  integer :: n, j, maxdeg = 2**15
  integer*8 :: plan
  complex(8), allocatable :: fvals(:), fvals_tmp(:)
  double precision, allocatable :: coeffs(:)
  ! double precision, pointer :: coeffs_ptr(:)
  double precision :: PI = 4.0 * datan(1.d0)
  
  interface
     function f(x)
       double precision :: f, x
     end function f
  end interface

  ! Perform a preliminary evaluation of the function at MINDEG points
  n = 5

  allocate(fvals(n), coeffs(n))

  do j = 1, n
     xx = dcos(PI * (n-j) / (n - 1))
     fvals(j) = f(xx)
     tol = max(tol, abs(fvals(j)) * eps)
  end do

  call cqr_interp2(n, fvals, coeffs)

  do while (n .le. maxdeg)
     ! This vector will be used to store a temporary copy of the
     ! old evaluations at the Chebyshev points. 
     allocate(fvals_tmp(n))
     fvals_tmp = fvals

     deallocate(fvals); allocate(fvals(2*n - 1))
     call cqr_upsample(n, fvals_tmp, fvals)

     deallocate(fvals_tmp)
     
     ! Evaluate the function at the new Chebyshev points, and
     ! at the same time check the error on these new points of
     ! the previous interpolant
     err = 0.d0
     do j = 2, 2*n-1, 2
        ! j-th Chebyshev point
        xx = dcos(PI * (2*n-1-j) / (2*n - 2))
        ! print *, 'xx = ', xx
        fvals(j) = f(xx)
        tol = max(tol, abs(fvals(j)) * eps)
        
        ! FIXME: Check the interpolation error at the new points
        !call cqr_clenshaw(n, coeffs, xx, yy)
        !err = max(err, abs(yy - fvals(j)))
     end do

     ! Interpolate a polynomial of higher degree
     n = 2 * n - 1     
     deallocate(coeffs); allocate(coeffs(n))

     call cqr_interp2(n, fvals, coeffs)     

     err = 0.d0
     nrm = 0.d0
     
     do j = (n-1)/2, n
        err = err + abs(coeffs(j))
     end do
     do j = 1, n
        nrm = nrm + abs(coeffs(j))
     end do
     err = err / nrm    

     ! print *, 'Interpolation of degree =', n, ' error =', err, ' tol =', tol

     if (err .le. tol) then
        exit
     end if
     
  end do

  ! Make sure that the leading coefficient is not zero
  nrm = 0.d0
  do j = 1, n -1
     nrm = nrm + abs(coeffs(j))
  end do

  do while (abs(coeffs(n)) .lt. nrm * tol)
     n = n - 1
  end do
    
end subroutine cqr_interp

!
! Given a vector x of length n, set the components of y, which
! is expected to have length 2*n - 1, equal to:
!
!  y(2*i-1) = x(i),   i = 1, ..., n
!
subroutine cqr_upsample(n, x, y)
  integer :: n, i
  complex(8) :: x(n), y(2*n-1)

  do i = 1, n
     y(2*i-1) = x(i)
  end do
end subroutine cqr_upsample

!
! Compute the Chebyshev coefficients given the evaluations at the Chebyshev
! points, and the number of points.
!
subroutine cqr_interp2(n, fvals, c)
  integer :: n, j
  integer*8 :: plan
  complex(8) :: fvals(*)
  double precision :: c(*)
  complex(8), allocatable :: fft_in(:), fft_out(:)

  allocate(fft_in(2*n-2), fft_out(2*n-2))

  ! -1 and 0 are for FFTW_FORWARD and FFTW_ESTIMATE, respectively
  call fftw_f77_create_plan(plan, 2*n-2, -1, 0)

  do j = 1, n        
     fft_in(n-j+1) = fvals(j)
     if (j .lt. n) then
        fft_in(n+j-1) = fft_in(n-j+1)
     end if
  end do
  
  call fftw_f77_one(plan, fft_in, fft_out)
  call fftw_f77_destroy_plan(plan)

  fft_out = fft_out / (n-1)
  fft_out(1) = fft_out(1) / 2.d0
  fft_out(2*n-2) = fft_out(2*n-2) / 2.d0

  c(1:n) = real(fft_out(1:n))
  c(n) = c(n) / 2

  deallocate(fft_in, fft_out)
  
end subroutine cqr_interp2

subroutine cqr_clenshaw(n, coeffs, x, y)
  integer :: n, j
  double precision :: coeffs(*), x, y, y1, y2

  if (n .eq. 0) then
     y = 0.d0
  elseif (n .eq. 1) then
     y = coeffs(1) * 1.d0
  elseif (n .eq. 2) then
     y = coeffs(1) + coeffs(2) * x
  else
     y1 = coeffs(n)
     y2 = y1 * x + coeffs(n-1)

     do j = n - 2, 1, -1
        y = 2.d0 * x * y2 - y1 + coeffs(j) - coeffs(j+1) * x
        y1 = y2
        y2 = y
     end do
  end if
  
end subroutine cqr_clenshaw
