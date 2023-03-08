! Create the structured representation of the colleague matrix
! for the polynomial expressed in the Chebyshev basis:
!
!  p(x) = p(1) T_n + p(2) T_{n-1} + ... + p(n+1) T_0.
!
! INPUT PARAMETERS
!
! N    INTEGER. Degree of the polynomial.
!
! P    COMPLEX(8), DIMENSION(N+1). Coefficients of the polynomial
!      in the Chebyhshev basis, ordered from the highest to the
!      smallest degree.
!
! D    COMPLEX(8), DIMENSION(N). Diagonal of the colleague matrix,
!      returned as an output.
!
! BETA COMPLEX(8), DIMENSION(N-1). Subdiagonal of the colleague matrix,
!      returned as an output.
!
! U    COMPLEX(8), DIMENSION(N). Left factor of the rank-1 representation
!      of the upper diagonal part.
!
! V    COMPLEX(8), DIMENSION(N). Right factor of the rank-1 representation
!      of the upper diagonal part.
!
subroutine cqr_colleague_real(n, p, d, beta, u, v)
  implicit none

  integer :: n
  real(8) :: d(1:n), beta(1:n-1), u(1:n), v(1:n), p(1:n+1)

  d = 0.d0
  beta = .5d0
  beta(n-1) = 1.d0 / dsqrt(2.d0)

  v = -p(2:n+1) / p(1);
  v(n) = v(n) * dsqrt(2.d0)
  v = v / 2.d0

  u = 0.d0
  u(1) = 1

  d(1) = v(1)


end subroutine

! Create the structured representation of the colleague matrix
! for the polynomial expressed in the Chebyshev basis:
!
!  p(x) = p(1) T_n + p(2) T_{n-1} + ... + p(n+1) T_0.
!
! INPUT PARAMETERS
!
! N    INTEGER. Degree of the polynomial.
!
! P    COMPLEX(8), DIMENSION(N+1). Coefficients of the polynomial
!      in the Chebyhshev basis, ordered from the highest to the
!      smallest degree.
!
! D    COMPLEX(8), DIMENSION(N). Diagonal of the colleague matrix,
!      returned as an output.
!
! BETA COMPLEX(8), DIMENSION(N-1). Subdiagonal of the colleague matrix,
!      returned as an output.
!
! U    COMPLEX(8), DIMENSION(N). Left factor of the rank-1 representation
!      of the upper diagonal part.
!
! V    COMPLEX(8), DIMENSION(N). Right factor of the rank-1 representation
!      of the upper diagonal part.
!
subroutine cqr_colleague(n, p, d, beta, u, v)
  implicit none

  integer :: n
  complex(8) :: d(1:n), beta(1:n-1), u(1:n), v(1:n), p(1:n+1)

  d = dcmplx(0.d0, 0.d0)
  beta = dcmplx(.5d0, 0.d0)
  beta(n-1) = 1.d0 / dsqrt(2.d0)

  v = -p(2:n+1) / p(1);
  v(n) = v(n) * dsqrt(2.d0)
  v = v / 2.d0

  u = dcmplx(0.d0, 0.d0)
  u(1) = 1

  d(1) = v(1)


end subroutine

! Compute the roots of a polynomial p(x) of the form
!
!  p(x) = p(1) T_n + p(2) T_{n-1} + ... + p(n+1) T_0,
!
! where T_j are the Chebyshev polynomials of the first kind.
!
! INPUT PARAMETERS
!
! N   INTEGER. Degree of the polynomial.
!
! P   REAL(8), DIMENSION(N+1). Coefficients of the polynomial,
!     ordered from the highest to the lowest degree.
!
! X   COMPLEX(8), DIMENSION(N). Roots of the polynomial, returned
!     in output.
subroutine cqr_roots_real(n, p, x)
  implicit none

  integer :: n, i
  real(8) :: p(1:n+1), M(2,2), c, s, rt1r, rt1i, rt2r, rt2i
  complex(8) :: x(1:n)
  real(8), allocatable, dimension(:) :: d, beta, u, v
  integer, parameter :: k = 15

  allocate(d(n), beta(n-1), u(n), v(n))

  call cqr_colleague_real(n, p, d, beta, u, v)
  call cqr_double_eig(n, d, beta, u, v, 0, 'n')

  x = 0.d0
  i = 1

  do while (i .le. n)
    if (i .lt. n) then
      if (beta(i) .ne. 0) then
        M(1,1) = d(i)
        M(1,2) = beta(i) - u(i+1) * v(i) + u(i) * v(i+1)
        M(2,1) = beta(i)
        M(2,2) = d(i+1)
        call dlanv2(M(1,1), M(1,2), M(2,1), M(2,2), rt1r, rt1i, rt2r, rt2i, c, s)
        x(i) = dcmplx(rt1r, rt1i)
        x(i+1) = dcmplx(rt2r, rt2i)
        i = i + 2
      else
        x(i) = dcmplx(d(i), 0.d0)
        i = i + 1
      end if
    else
      x(i) = dcmplx(d(i), 0.d0)
      i = i + 1
    end if
  end do

  deallocate(d, beta, u, v)
end

! Compute the roots of a polynomial p(x) of the form
!
!  p(x) = p(1) T_n + p(2) T_{n-1} + ... + p(n+1) T_0,
!
! where T_j are the Chebyshev polynomials of the first kind.
!
! INPUT PARAMETERS
!
! N   INTEGER. Degree of the polynomial.
!
! P   COMPLEX(8), DIMENSION(N+1). Coefficients of the polynomial,
!     ordered from the highest to the lowest degree.
!
! X   COMPLEX(8), DIMENSION(N). Roots of the polynomial, returned
!     in output.
subroutine cqr_roots(n, p, x)
  implicit none

  integer :: n
  complex(8) :: p(1:n+1), x(1:n)
  complex(8), allocatable :: beta(:), u(:), v(:)
  integer, parameter :: k = 15

  allocate(beta(n-1), u(n), v(n))

  call cqr_colleague(n, p, x, beta, u, v)
  call cqr_single_eig(n, x, beta, u, v, k, 1, 0.d0, 'n', 0, 'n')

  deallocate(beta, u, v)

end


