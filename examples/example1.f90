! Example 1: Computing the roots of a smooth function in [-1, 1]

function f(x)
  implicit none
  double precision :: f, x
  f = dsin(100 * x);
end function f

program example1
  implicit none

  integer j, n
  double precision, allocatable :: zeros(:)
  double precision :: y

  interface
     function f(x)
       double precision :: f, x
     end function f
  end interface

  INCLUDE '../chebqr.h'
  
  call cqr_zeros(f, 1.d-12, zeros, n)

  print *, 'Number of roots in [-1, 1]:', n

  do j = 1, n
     print *, 'root #', j, ' = ', zeros(j), 'f(root) =', f(zeros(j))
  end do

  deallocate(zeros)

end program example1

