! Example 1: Computing the roots of a smooth function in [-1, 1],
! for which exact roots are known. 

function f(x)
  implicit none
  double precision :: f, x
  f = dsin(400 * x);
end function f

program example1
  implicit none

  integer j, n
  double precision, allocatable :: zeros(:)
  double precision :: y, PI = 4.d0 * atan(1.d0), tol = 1.d-12

  interface
     function f(x)
       double precision :: f, x
     end function f
  end interface

  ! This is needed to have the interface for calling cqr_zeros,
  ! which allocates the vector zeros(:). 
  INCLUDE '../chebqr.h'
  
  call cqr_zeros(f, tol, zeros, n)

  print '(A)', 'Function: sin(400*x) over [-1, 1]'
  print '(A, E8.3)', 'Approximation tolerance: ', tol
  print '(A, I4)', 'Number of expected roots: ', floor(2 * 400 / PI + 1)
  print '(A, I4)', 'Number of identified roots in [-1, 1]:', n
  print *, ''

  do j = 1, n
     print '(A, I3, A, F8.5, A, E10.3, A, F8.5, A, E10.3)', &
          'root #', j, ' = ', zeros(j), &
          ' | f(root) = ', f(zeros(j)), &
          ' | exact root = ', PI * (j - n/2 - 1) / 400, &
          ' | error = ', abs(zeros(j) - PI * (j - n/2 - 1) / 400)
  end do

  deallocate(zeros)

end program example1

