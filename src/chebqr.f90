! Package ChebQR -- computing roots of smooth functions over [-1, 1]
!
! Authors:
!  Angelo Casulli <angelo.casulli@sns.it>
!  Leonardo Robol <leonardo.robol@unipi.it>
!

subroutine cqr_zeros(f, eps, zeros, n)
  implicit none
  
  integer :: n, k = 25, j, nroots, i
  double precision :: eps, dr, di
  double precision, allocatable :: zeros(:), coeffs(:)

  ! Working variables for this routine
  complex(8), allocatable :: d(:), u(:), v(:), beta(:)

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
  call cqr_interp(f, eps, coeffs, n)

  !print *, 'n =', n
  !do j = 1, n
  !   print *, 'coeffs of deg', j-1, '=', coeffs(j)
  !end do
  
  ! Build the companion linearization in the Chebyshev basis
  ! and compute the eigenvalues.
  allocate(u(1:n-1), v(1:n-1), d(1:n-1), beta(1:n-2))

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

  ! Compute the eigenvalues
  call fastfastqr(n-1, d, beta, u, v, k)

  ! Extract the eigenvalues which are in [-1, 1]. First we count
  ! them and then allocate the array to store them
  nroots = 0
  do j = 1, n-1
     dr = realpart(d(j))
     di = imagpart(d(j))
     if ( (dr .ge. -1.d0) .and. &
          (dr .le.  1.d0) .and. &
          abs(di) .le. 1.0d-12) then
        nroots = nroots + 1
     end if
  end do

  ! print *, 'nroots = ', nroots

  if (nroots .gt. 0) then
     allocate(zeros(1:nroots))
     nroots = 0
     do j = 1, n-1
        dr = realpart(d(j))
        di = imagpart(d(j))
        if ( (dr .ge. -1.d0) .and. &
             (dr .le.  1.d0) .and. &
             abs(di) .le. 1.0d-12) then        
           nroots = nroots + 1
           zeros(nroots) = real(d(j))
        end if
     end do     
  end if

  n = nroots

  ! sort the roots
  call cqr_sort_array(n, zeros)

  deallocate(d,u,v,beta,coeffs)  
end subroutine cqr_zeros

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
        call cqr_clenshaw(n, coeffs, xx, yy)
        err = max(err, abs(yy - fvals(j)))
     end do

     ! print *, 'Interpolation of degree =', n, ' error =', err, ' tol =', tol

     if (err .le. tol) then
        exit
     end if
     
     ! Interpolate a polynomial of higher degree
     n = 2 * n - 1     
     deallocate(coeffs); allocate(coeffs(n))

     call cqr_interp2(n, fvals, coeffs)     
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

  call dfftw_plan_dft_1d(plan, 2*n-2, fft_in, fft_out, &
     -1, & ! -1: FFTW_FORWARD, 1 : FFTW_BACKWARD
     0   & ! 0: FFTW_ESTIMATE, 1: FFTW_MEASURE
  )

  do j = 1, n        
     fft_in(n-j+1) = fvals(j)
     if (j .lt. n) then
        fft_in(n+j-1) = fft_in(n-j+1)
     end if
  end do
  
  call dfftw_execute_dft(plan, fft_in, fft_out)
  call dfftw_destroy_plan(plan)

  fft_out = fft_out / (n-1)
  fft_out(1) = fft_out(1) / 2.d0
  fft_out(2*n-2) = fft_out(2*n-2) / 2.d0

  c(1:n) = real(fft_out(1:n))

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
     y2 = coeffs(n) * x

     do j = n - 1, 2, -1
        y = 2.d0 * x * y2 - y1 + coeffs(j) * x
        y1 = y2 + coeffs(j)
        y2 = y
     end do
     
  end if
  
end subroutine cqr_clenshaw