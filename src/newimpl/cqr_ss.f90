recursive subroutine cqr_ss_qr1(n, d, b, u, v, &
     aed, schurv, Q, nQ, ldQ, its, info)
  implicit none

  integer :: n, its, imin, imax, last_deflation, deflated, info, nQ, ldQ
  complex(8) :: d(n), b(n-1), u(n), v(n), s, Q(*)
  logical :: schurv, aed
  character(len=256) :: buffer

  its = 0
  info = 0

  imin = 1
  imax = n

  last_deflation = 0

  do while (imax .gt. imin)
     ! Compute a Wilkinson shift for use in the next iteration,
     ! unless we have skipped a few iterations without deflations,
     ! then use a random shift.
     if (its .ge. last_deflation + 20) then
        call cqr_ss_random_shift(s)
        last_deflation = its
     else
        call cqr_ss_wilk_shift(d, b, u, v, imax, s)
     end if

     ! Introduce the bulge and chase it to the bottom
     call cqr_ss_chase(d, b, u, v, imin, imax, s, schurv, Q, nQ, ldQ)

     ! Check deflations at top and bottom of the matrix
     call cqr_ss_check_deflations(d, b, imin, imax, deflated)

     ! Check for middle deflations, and run a QR iteration on the subproblem
     if (mod(its, 5) .eq. 0) then
        call cqr_ss_check_middle_deflations(d, b, u, v, imin, imax, aed, &
             schurv, Q, nQ, ldQ, its, info)
     end if

     if (deflated .gt. 0) then
        last_deflation = its
     end if

     ! If we are doing too many iterations, return and set failure    
     if (its .ge. 10 * n + 100) then
        info = 1
     end if

     ! If an error code has been set somewhere, return
     if (info .ne. 0) then
        return
     end if

     its = its + 1
  end do  
  
end subroutine cqr_ss_qr1

subroutine cqr_ss_random_shift(s)
  implicit none

  complex(8) :: s
  double precision :: x, y

  call random_number(x)
  call random_number(y)

  s = cmplx(x, y)
  
end subroutine cqr_ss_random_shift

subroutine cqr_ss_wilk_shift(d, b, u, v, imax, s)
  implicit none

  integer :: imax
  complex(8) :: d(*), b(*), u(*), v(*), s, dA, trA, DD, l1, l2

  ! Compute determinant and trace of the trailing 2x2 matrix
  dA = d(imax-1) * d(imax) - &
       b(imax-1) * (conjg(b(imax-1))-conjg(u(imax))*v(imax-1) &
                    + u(imax-1)*conjg(v(imax)));
  trA = d(imax-1) + d(imax);
  DD = sqrt(trA**2 - 4*dA);

  ! FIXME: In principle, we may compute these in a slightly more
  ! stable manner -- not that is particularly relevant for the
  ! computation of the shift
  l1 = (trA + DD)/2
  l2 = (trA - DD)/2

  if ( abs(l1 - d(imax)) .lt. abs(l2 - d(imax)) ) then
     s = l1
  else
     s = l2
  end if

end subroutine cqr_ss_wilk_shift

subroutine cqr_ss_chase(d, b, u, v, imin, imax, s, schurv, Q, nQ, ldQ)
  implicit none

  integer :: imin, imax, i, nQ, ldQ
  complex(8) :: d(*), b(*), u(*), v(*), Q(*)
  complex(8) :: w1, w2, M(3,2), s
  double precision :: c
  logical :: schurv

  ! Handle trivial cases first
  if (imin .eq. imax) then
     return
  end if

  ! Determine the first rotation, that will introduce the bulge
  w1 = d(imin) - s
  w2 = b(imin)
  call zrotg(w1, w2, c, s)

  ! Construct the window that we use to work on the active part
  ! of the matrix.
  M(1,1) = d(imin)
  M(2,1) = b(imin)
  M(3,1) = 0
  M(1,2) = conjg(b(imin)) - conjg(u(imin+1)) * v(imin) + u(imin) * conjg(v(imin+1))
  M(2,2) = d(imin+1)

  ! We only need this row in case the matrix is at least 3x3
  if (imax .gt. imin + 1) then
     M(3,2) = b(imin+1)
  end if

  ! Update M with G from left and right
  call zrot(2, M(1,1), 3, M(2,1), 3, c, s)
  call zrot(3, M(1,1), 1, M(1,2), 1, c, conjg(s))

  ! Update the vectors u, v
  call zrot(1, u(imin), 1, u(imin+1), 1, c, s)
  call zrot(1, v(imin), 1, v(imin+1), 1, c, s)

  ! Update the Schur vectors, if requested
  call cqr_ss_update_schur_vectors(c, s, imin, schurv, Q, nQ, ldQ)

  d(imin) = M(1,1)

  do i = imin, imax - 2
     ! Annihilate the bulge, which is now stored in M
     w1 = M(2,1)
     w2 = M(3,1)
     call zrotg(w1, w2, c, s)

     ! Save the new subdiagonal entry
     b(i) = w1

     ! Update M, and then apply the rotations
     M(1,1) = M(2,2)
     M(2,1) = M(3,2)
     M(1,2) = conjg(M(2,1)) - conjg(u(i+2)) * v(i+1) + u(i+1) * conjg(v(i+2))
     M(2,2) = d(i+2)

     if (i .lt. imax - 2) then
        M(3,2) = b(i+2)
        M(3,1) = 0
     end if

     call zrot(2, M(1,1), 3, M(2,1), 3, c, s)
     call zrot(3, M(1,1), 1, M(1,2), 1, c, conjg(s))

     d(i+1) = M(1,1)

     call zrot(1, u(i+1), 1, u(i+2), 1, c, s)
     call zrot(1, v(i+1), 1, v(i+2), 1, c, s)

     ! Update the Schur vectors, if requested
     call cqr_ss_update_schur_vectors(c, s, i+1, schurv, Q, nQ, ldQ)
  end do

  ! Store the updated trailing entries before exiting
  d(imax) = M(2,2)
  b(imax-1) = M(2,1)
  
end subroutine cqr_ss_chase

subroutine cqr_ss_check_deflations(d, b, imin, imax, deflated)
  implicit none

  integer :: imin, imax, deflated
  complex(8) :: d(*), b(*)
  double precision :: eps, dlamch

  deflated = 0
  eps = dlamch('E')

  do while (&
       imin .lt. imax .and. &
       (abs(d(imin)) + abs(d(imin+1))) * eps .gt. abs(b(imin)))
     imin = imin + 1;
     deflated = deflated + 1
  end do

  do while (&
       imax .gt. imin .and. &
       (abs(d(imax-1)) + abs(d(imax))) * eps > abs(b(imax-1)))
     imax = imax - 1;
     deflated = deflated + 1
  end do
  
end subroutine cqr_ss_check_deflations

subroutine cqr_ss_check_middle_deflations(d, b, u, v, imin, imax, aed, &
     schurv, Q, nQ, ldQ, its, info)
  implicit none

  logical :: schurv, aed
  integer :: its, imin, imax, lits, i, info, nQ, ldQ
  complex(8) :: d(*), b(*), u(*), v(*), Q(*)
  double precision :: eps, dlamch, di, dim1

  di = abs(d(imin))
  eps = dlamch('E')

  do i = imin + 1, imax - 2
     ! This avoids computing absolute values twice
     dim1 = di
     di = abs(d(i))
     
     if (abs(b(i)) .lt. (di + dim1) * eps) then
        lits = 0
        call cqr_ss_qr1(i - imin + 1, d(imin), b(imin), u(imin), v(imin), &
             aed, schurv, Q(imin), nQ, ldQ, lits, info)

        if (info .ne. 0) then
           return
        end if
        
        its = its + lits
        imin = i + 1
     end if
  end do
  
end subroutine cqr_ss_check_middle_deflations

subroutine cqr_ss_update_schur_vectors(c, s, i, schurv, Q, nQ, ldQ)
  implicit none

  double precision :: c
  complex(8) :: s, Q(*)
  integer :: i, ldQ, nQ
  logical :: schurv

  if (schurv) then
     call zrot(nQ, Q(i), ldQ, Q(i+1), ldQ, c, s)
  end if
  
end subroutine cqr_ss_update_schur_vectors
