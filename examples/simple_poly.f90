!
! Compile this file as follows:
!
!  gfortran -o simple_poly simple_poly.f90 ../src/single_eig.f90 ../src/utils.f90 -lblas -llapack
!
! or just run 'make'
program simple_poly
  implicit none

  integer, parameter :: n = 8, k = 15
  integer :: j
  complex(8) :: p(n+1), x(n)
  
  p = cmplx(0.d0, 0.d0)
  p(1) = cmplx(1.d0, 0.d0)
  
  call cqr_roots(n, p, x)
  
  print *, '> Roots of the monic Chebyshev polynomial of degree:', n
  print *, ''
  
  do j = 1, n
    print *, 'x[', j, ']', x(j)
  end do
  
  ! Do another test, with the polynomial with coefficients 
  ! of degree j equal to (j+1)
  do j = 1, n + 1
    p(j) = n-j+2 ! coeff of degree n-j+1
  end do
  
  call cqr_roots(n, p, x)
   
  print *, ''
  print *, '> Roots of the Chebyshev polynomial of degree', n
  print *, '> with coefficients of degree j equal to (j+1):'
  print *, ''
  
  do j = 1, n
    print *, 'x[', j, ']', x(j)
  end do
  
  ! Do another test, with the polynomial with coefficients 
  ! of degree j equal to (j+1) * (-1)^j
  p = (/ 1.d0, 0.d0, 8.d0, 0.d0, 28.d0, 0.d0, 56.0d0, 0.d0, 163.0d0 /)
  call cqr_roots(n, p, x)
   
  print *, ''
  print *, '> Roots of x^8 + 1'
  print *, ''
  
  do j = 1, n
    print *, 'x[', j, ']', x(j)
  end do
  

end program
