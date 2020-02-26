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

    
interface
  subroutine cqr_zeros(f, eps, zeros, n)
      integer :: n
      double precision :: eps
      double precision, allocatable :: zeros(:)
      interface
         function f(x)
         double precision :: f, x
         end function f
      end interface
  end subroutine cqr_zeros
end interface
