!--------------------------------------------------
!SUBROUTINE cqr_double_eig
!
! This subroutine computes the eigenvalues of a matrix which is the
! sum  of a hermitian and a rank one matrices, using a structured double shift
! QR algorithm with a structured aggressive early deflation.
! If JOBITER=y the number of iterations is stored.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix
!
! D    REAL(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA REAL(8), DIMENSION(N-1). Vector that contains the subdiagonal
!      entries of the matrix.
!
! U,V  REAL(8), DIMENSION(N). Vectors such that the rank one part of
!      the matrix is UV*.
!
! K    INTEGER. Number of QR steps performed before the aggressive
!      early deflation is applied.
!
! ITER  INTEGER. Number of iterations required for convergence.
!
! JOBITER CHARACTER.  JOBITER=y if the number of iterations is requested,
!		      JOBITER=n otherwise.
! If JOBITER=y the number of iterations is stored.
subroutine cqr_double_eig(n, d, beta, u, v, iter, jobiter)
  implicit none

  integer, intent(in) :: n
  character, intent(in) :: jobiter
  integer, intent(out) :: iter
  real(8), dimension(n), intent(inout) :: d, u, v
  real(8), dimension(n-1), intent(inout) :: beta

  if (jobiter .eq. 'y') then
    iter = 0
  end if

  call cqr_fastfastqr_ds(n, d, beta, u, v, iter, jobiter)

end subroutine

!SUBROUTINE cqr_check_deflation
!
! Check for top and bottom deflation in the matrix, and adjust the imin and 
! imax integers accordingly.
!
!  INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of
!      the matrix is UV*.
subroutine cqr_check_deflation(n, d, beta, u, v, imin, imax, count)
  implicit none

  integer, intent(in) :: n
  integer, intent(out) :: count
  integer, intent(inout) :: imin, imax
  real(8), intent(inout) :: d(n), beta(n-1), u(n), v(n)
  real(8) ::  eps, dlamch
  logical :: deflation_performed

  count = 0
  eps = dlamch('e')
  deflation_performed = .true.

  ! This cycle keeps going as long as at least one deflation is performed. 
  do while (deflation_performed)
    deflation_performed = .false.

    ! 1x1 top deflation
    if (imin .lt. imax) then
      if (abs(beta(imin)) .le. eps * (abs(d(imin)) + abs(d(imin+1)))) then        
        beta(imin) = 0.d0
        deflation_performed = .true.
        count = count + 1
        imin = imin + 1
      end if
    end if

    ! 2x2 top deflation
    if (imin + 1 .lt. imax) then
      if (abs(beta(imin+1)) .le. eps * (abs(d(imin+1)) + abs(d(imin+2)))) then        
        beta(imin+1) = 0.d0
        call cqr_fastqr6_ds_in(2, d(imin:imin+1), beta(imin:imin), u(imin:imin+1), v(imin:imin+1))
        imin = imin + 2
        deflation_performed = .true.
        count = count + 2
      end if
    end if

    ! 1x1 bottom deflation
    if ((imin .gt. 1) .and. (imin .lt. imax)) then
      if (abs(beta(imax-1)) .le. eps * (abs(d(imax-1)) + abs(d(imax)))) then
        beta(imax-1) = 0.d0
        deflation_performed = .true.
        count = count + 1
        imax = imax - 1        
      end if
    end if

    ! ! 2x2 bottom deflation
    if ((imin .gt. 2) .and. (imin + 1 .lt. imax)) then
      if (abs(beta(imax-2)) .le. eps * (abs(d(imax-2)) + abs(d(imax-1)))) then
        beta(imax-2) = 0.d0
        call cqr_fastqr6_ds_in(2, d(imax-1:imax), beta(imax-1:imax-1), u(imax-1:imax), v(imax-1:imax))
        deflation_performed = .true.
        count = count + 2        
        imax = imax - 2
      end if
    end if
  end do
  
end subroutine

!SUBROUTINE cqr_fastfastqr_ds
!
! This subroutine performs a double shift structured QR algorithm, without
! aggressive early deflation, to compute the eigenvalue of a matrix which is
! the sum of a hermitian and a rank one matrix.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of
!      the matrix is UV*.

recursive subroutine cqr_fastfastqr_ds(n,d,beta,u,v,iter,jobiter)
  implicit none
  integer, intent(in)  :: n
  integer :: imin, imax , cont, i, deflated_count
  integer, intent(inout) :: iter
  character, intent(in) :: jobiter
  real(8), dimension(n), intent(inout) :: d, u, v
  real(8), dimension(n-1), intent(inout) :: beta
  real(8) :: eps = 2.22e-16

  imax=n
  imin=1
  cont=0

  ! Try to deflate some eigenvalue
  call cqr_check_deflation(n, d, beta, u, v, imin, imax, deflated_count)

  do while (imax-imin .gt. 1)

    ! Compute a step of the double shift QR algorithm.
    call cqr_fastqr6_ds_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
    if (jobiter .eq. 'y') then
      iter = iter + 1
    end if

!Try to deflate some eigenvalue
    call cqr_check_deflation(n, d, beta, u, v, imin, imax, deflated_count)
    if (deflated_count .gt. 0) then
      cont = 0
    end if

    do i=imin+1,imax-1
! If a deflation occurs in the middle of the matrix,
! compute the eigenvalues of the smallest diagonal block,
! using a recursive structured QR algorithm.
      if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
        beta(i)=0
        if ( (i.le. (imax-imin)/2)) then
          call cqr_fastfastqr_ds(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i), iter, jobiter)
          do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
            beta(imin)=0
            imin = imin + 1
            cont=0
          end do
        else
          call cqr_fastfastqr_ds(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax), iter, jobiter)
          do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
            beta(imax-1)=0
            imax = imax - 1
            cont=0
          end do
        end if
      end if
    end do

    cont=cont+1
! If after some QR iteration there is not deflation, perform a structured
! QR step using a random shift vector.
    if (cont==10) then
      call fastqr7_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
      cont=0
    end if

  end do

  if (imax .gt. imin) then
    call cqr_fastqr6_ds_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
    if (jobiter .eq. 'y') then
      iter = iter + 1
    end if
  end if

end subroutine cqr_fastfastqr_ds
!------------------------------------
!SUBROUTINE cqr_fastqr6_ds_in
!
! This subroutine performs a step of the double shift structured QR algorithm,
! to compute the eigenvalue of a matrix which is the sum of a hermitian and a
! rank one matrix.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of
!      the matrix is UV*.

subroutine cqr_fastqr6_ds_in(n,d,beta,u,v)
  implicit none
  integer, intent(in)  :: n
  real(8), dimension(n), intent(inout) :: d, u,v
  real(8), dimension(n-1), intent(inout) :: beta
  real(8), dimension(n-1) :: gamm
  real(8), dimension(n-2) :: alph
  real(8), dimension(2) :: l
  real(8) :: rho
  real(8), dimension(4,3) :: R
  integer :: i
  real(8) :: S, C, ss, cc, r1,i1,r2,i2
  real(8):: z,zz,z2,zz2
  double precision :: dlamch
  real(8):: eps
  
  eps =dlamch('e')

  if (n>3) then
    gamm(1)=(beta(1)-v(1)*u(2))+u(1)*v(2)
    z=d(n-1)*d(n)-(beta(n-1)*((beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))
    rho=(d(n-1)+d(n))**2-4*z
    if (rho.ge.0) then
      rho=sqrt(rho)
      l(1)=(d(n-1)+d(n)+rho)/2
      l(2)=(d(n-1)+d(n)-rho)/2
      if (abs(l(1)-d(n))<abs(l(2)-d(n))) then
        rho=l(1);
      else
        rho=l(2);
      endif

      z=d(1)-rho
      z2=beta(1)
      call drotg(z,z2,C,S)


      R(1,1)=d(1)
      R(2,1)=beta(1)
      R(1,2)=gamm(1)
      R(2,2)=d(2)
      R(3,1)=0
      R(3,2)=beta(2)
      R(4,:)=0
      R(:,3)=0

      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(3, R(1,1), 1, R(1,2), 1, C, (S))

      d(1)=R(1,1)
      beta(1)=R(2,1)
      d(2)=R(2,2)
      beta(2)=R(3,2)

      call drot(1, u(1), 1, u(2), 1, C, S)
      call drot(1, v(1),1,v(2), 1, C, (S))


      do i=1,n-3


        gamm(i+1)=(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
        z=beta(i)
        z2=R(3,1)

        call drotg(z,z2,C,S)
        call drot(1, beta(i), 1, R(3,1), 1, C, S)


        R(1,1)=d(i+1)
        R(2,1)=beta(i+1)
        R(1,2)=gamm(i+1)
        R(2,2)=d(i+2)
        R(3,1)=0
        R(3,2)=beta(i+2)


        call drot(2, R(1,1), 4, R(2,1), 4, C, S)
        call drot(3, R(1,1), 1, R(1,2), 1, C, (S))


        d(i+1)=R(1,1)
        beta(i+1)=R(2,1)
        d(i+2)=R(2,2)
        beta(i+2)=R(3,2)

        call drot(1, u(i+1), 1, u(i+2), 1, C, S)
        call drot(1, v(i+1),1,v(i+2), 1, C, (S))
      enddo
      gamm(n-1)=(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
      z=beta(n-2)
      z2=R(3,1)
      call drotg(z,z2,C,S)
      call drot(1, beta(n-2), 1, R(3,1), 1, C, S)

      beta(n-2)=z
      R(1,1)=d(n-1)
      R(2,1)=beta(n-1)
      R(1,2)=gamm(n-1)
      R(2,2)=d(n)

      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(2, R(1,1), 1, R(1,2), 1, C, (S))

      d(n-1)=R(1,1)
      beta(n-1)=R(2,1)
      d(n)=R(2,2)

      call drot(1, u(n-1), 1, u(n), 1, C, S)
      call drot(1, v(n-1),1,v(n), 1, C, (S))

    else
      z=d(1)**2+gamm(1)*beta(1)-(d(n-1)+d(n))*d(1)+z
      zz=beta(1)*(d(1)+d(2)-d(n-1)-d(n))
      R(1,1)=zz
      zz2=beta(1)*beta(2)
      R(2,1)=zz2
      call drotg(zz,zz2,CC,SS)
      call drot(1, R(1,1), 1, R(2,1), 1, CC, SS)
      z2=R(1,1)
      call drotg(z,z2,C,S)

      gamm(2)=(beta(2)-v(2)*u(3))+u(2)*v(3)
      alph(1)=(-v(1)*u(3))+u(1)*v(3)

      R(1,1)=d(1)
      R(2,1)=beta(1)
      R(3,1)=0
      R(4,1)=0
      R(1,2)=gamm(1)
      R(2,2)=d(2)
      R(3,2)=beta(2)
      R(4,2)=0
      R(1,3)=alph(1)
      R(2,3)=gamm(2)
      R(3,3)=d(3)
      R(4,3)=beta(3)


      call drot(3, R(2,1), 4, R(3,1), 4, CC, SS)
      call drot(4, R(1,2), 1, R(1,3), 1, CC, (SS))
      call drot(3, R(1,1), 4, R(2,1), 4, C, S)
      call drot(4, R(1,1), 1, R(1,2), 1, C, (S))

      d(1)=R(1,1)
      beta(1)=R(2,1)
      d(2)=R(2,2)
      beta(2)=R(3,2)
      d(3)=R(3,3)
      beta(3)=R(4,3)


      call drot(1, u(2), 1, u(3), 1, CC, SS)
      call drot(1, v(2), 1, v(3), 1, CC, (SS))
      call drot(1, u(1), 1, u(2), 1, C, S)
      call drot(1, v(1),1,v(2), 1, C, (S))


      do i=1,n-4

        gamm(i+1)=(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
        gamm(i+2)=(beta(i+2)-v(i+2)*u(i+3))+u(i+2)*v(i+3)
        alph(i+1)=(R(4,2)-v(i+1)*u(i+3))+u(i+1)*v(i+3)
        zz=R(3,1)
        z=beta(i)
        zz2=R(4,1)
        call drotg(zz,zz2,CC,SS)
        call drot(1, R(3,1), 1, R(4,1), 1, CC, SS)
        z2=R(3,1)
        call drotg(z,z2,C,S)
        call drot(1, beta(i), 1, R(3,1), 1, C, S)

        R(1,1)=d(i+1)
        R(2,1)=beta(i+1)
        R(3,1)=R(4,2)
        R(4,1)=0
        R(1,2)=gamm(i+1)
        R(2,2)=d(i+2)
        R(3,2)=beta(i+2)
        R(4,2)=0
        R(1,3)=alph(i+1)
        R(2,3)=gamm(i+2)
        R(3,3)=d(i+3)
        R(4,3)=beta(i+3)


        call drot(3, R(2,1), 4, R(3,1), 4, CC, SS)
        call drot(4, R(1,2), 1, R(1,3), 1, CC, (SS))
        call drot(3, R(1,1), 4, R(2,1), 4, C, S)
        call drot(4, R(1,1), 1, R(1,2), 1, C, (S))


        d(i+1)=R(1,1)
        beta(i+1)=R(2,1)
        d(i+2)=R(2,2)
        beta(i+2)=R(3,2)
        d(i+3)=R(3,3)
        beta(i+3)=R(4,3)

        call drot(1, u(i+2), 1, u(i+3), 1, CC, SS)
        call drot(1, v(i+2), 1, v(i+3), 1, CC, (SS))
        call drot(1, u(i+1), 1, u(i+2), 1, C, S)
        call drot(1, v(i+1),1,v(i+2), 1, C, (S))
      enddo


      gamm(n-2)=(beta(n-2)-v(n-2)*u(n-1))+u(n-2)*v(n-1)
      gamm(n-1)=(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n)
      alph(n-2)=(R(4,2)-v(n-2)*u(n))+u(n-2)*v(n)
      zz=R(3,1)
      z=beta(n-3)
      zz2=R(4,1)
      call drotg(zz,zz2,CC,SS)
      call drot(1, R(3,1), 1, R(4,1), 1, CC, SS)
      z2=R(3,1)
      call drotg(z,z2,C,S)
      call drot(1, beta(n-3), 1, R(3,1), 1, C, S)

      R(1,1)=d(n-2)
      R(2,1)=beta(n-2)
      R(3,1)=R(4,2)
      R(1,2)=gamm(n-2)
      R(2,2)=d(n-1)
      R(3,2)=beta(n-1)
      R(1,3)=alph(n-2)
      R(2,3)=gamm(n-1)
      R(3,3)=d(n)


      call drot(3, R(2,1), 4, R(3,1), 4, CC, SS)
      call drot(3, R(1,2), 1, R(1,3), 1, CC, (SS))
      call drot(3, R(1,1), 4, R(2,1), 4, C, S)
      call drot(3, R(1,1), 1, R(1,2), 1, C, (S))


      d(n-2)=R(1,1)
      beta(n-2)=R(2,1)
      d(n-1)=R(2,2)
      beta(n-1)=R(3,2)
      d(n)=R(3,3)

      call drot(1, u(n-1), 1, u(n), 1, CC, SS)
      call drot(1, v(n-1), 1, v(n), 1, CC, (SS))
      call drot(1, u(n-2), 1, u(n-1), 1, C, S)
      call drot(1, v(n-2),1,v(n-1), 1, C, (S))

      gamm(n-1)=(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
      z=beta(n-2)
      z2=R(3,1)
      call drotg(z,z2,C,S)
      call drot(1, beta(n-2), 1, R(3,1), 1, C, S)

      R(1,1)=d(n-1)
      R(2,1)=beta(n-1)
      R(1,2)=gamm(n-1)
      R(2,2)=d(n)


      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(2, R(1,1), 1, R(1,2), 1, C,(S))



      d(n-1)=R(1,1)
      beta(n-1)=R(2,1)
      d(n)=R(2,2)



      call drot(1, u(n-1), 1, u(n), 1, C, S)
      call drot(1, v(n-1),1,v(n), 1, C, (S))




    end if
  else if (n==3) then
    gamm(1)=(beta(1)-v(1)*u(2))+u(1)*v(2)
    z=d(n-1)*d(n)-(beta(n-1)*((beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))
    rho=(d(n-1)+d(n))**2-4*z
    if (rho.ge.0) then
      rho=sqrt(rho)
      l(1)=(d(n-1)+d(n)+rho)/2
      l(2)=(d(n-1)+d(n)-rho)/2
      if (abs(l(1)-d(n))<abs(l(2)-d(n))) then
        rho=l(1);
      else
        rho=l(2);
      endif

      z=d(1)-rho
      z2=beta(1)
      call drotg(z,z2,C,S)

      R(1,1)=d(1)
      R(2,1)=beta(1)
      R(1,2)=gamm(1)
      R(2,2)=d(2)
      R(3,1)=0
      R(3,2)=beta(2)
      R(4,:)=0
      R(:,3)=0


      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(3, R(1,1), 1, R(1,2), 1, C, (S))

      d(1)=R(1,1)
      beta(1)=R(2,1)
      d(2)=R(2,2)
      beta(2)=R(3,2)

      call drot(1, u(1), 1, u(2), 1, C, S)
      call drot(1, v(1),1,v(2), 1, C, (S))


      gamm(n-1)=(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
      z=beta(n-2)
      z2=R(3,1)
      call drotg(z,z2,C,S)
      call drot(1, beta(n-2), 1, R(3,1), 1, C, S)


      R(1,1)=d(n-1)
      R(2,1)=beta(n-1)
      R(1,2)=gamm(n-1)
      R(2,2)=d(n)


      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(2, R(1,1), 1, R(1,2), 1, C, (S))



      d(n-1)=R(1,1)
      beta(n-1)=R(2,1)
      d(n)=R(2,2)



      call drot(1, u(n-1), 1, u(n), 1, C, S)
      call drot(1, v(n-1),1,v(n), 1, C, (S))

    else

      z=d(1)**2+gamm(1)*beta(1)-(d(n-1)+d(n))*d(1)+z
      zz=beta(1)*(d(1)+d(2)-d(n-1)-d(n))
      R(1,1)=zz
      zz2=beta(1)*beta(2)
      call drotg(zz,zz2,CC,SS)
      R(2,1)=beta(1)*beta(2)
      call drot(1, R(1,1), 1, R(2,1), 1, CC, SS)
      z2=R(1,1)
      call drotg(z,z2,C,S)

      gamm(2)=(beta(2)-v(2)*u(3))+u(2)*v(3)
      alph(1)=(-v(1)*u(3))+u(1)*v(3)

      R(1,1)=d(1)
      R(2,1)=beta(1)
      R(3,1)=0
      R(4,1)=0
      R(1,2)=gamm(1)
      R(2,2)=d(2)
      R(3,2)=beta(2)
      R(4,2)=0
      R(1,3)=alph(1)
      R(2,3)=gamm(2)
      R(3,3)=d(3)
      R(4,3)=0


      call drot(3, R(2,1), 4, R(3,1), 4, CC, SS)
      call drot(3, R(1,2), 1, R(1,3), 1, CC, (SS))
      call drot(3, R(1,1), 4, R(2,1), 4, C, S)
      call drot(3, R(1,1), 1, R(1,2), 1, C, (S))

      d(1)=R(1,1)
      beta(1)=R(2,1)
      d(2)=R(2,2)
      beta(2)=R(3,2)
      d(3)=R(3,3)


      call drot(1, u(2), 1, u(3), 1, CC, SS)
      call drot(1, v(2), 1, v(3), 1, CC, (SS))
      call drot(1, u(1), 1, u(2), 1, C, S)
      call drot(1, v(1),1,v(2), 1, C, (S))

      gamm(n-1)=(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
      z=beta(n-2)
      z2=R(3,1)
      call drotg(z,z2,C,S)
      call drot(1, beta(n-2), 1, R(3,1), 1, C, S)

      R(1,1)=d(n-1)
      R(2,1)=beta(n-1)
      R(1,2)=gamm(n-1)
      R(2,2)=d(n)


      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(2, R(1,1), 1, R(1,2), 1, C,(S))



      d(n-1)=R(1,1)
      beta(n-1)=R(2,1)
      d(n)=R(2,2)

      call drot(1, u(n-1), 1, u(n), 1, C, S)
      call drot(1, v(n-1),1,v(n), 1, C, (S))
    end if

  else if (n==2) then
    gamm(1)=(beta(1)-u(2)*v(1))+u(1)*v(2)

    call DLANV2(d(1),gamm(1),beta(1),d(2), r1,i1,r2,i2,C,S)
    call drot(1, u(1), 1, u(2), 1, C, S) ! forse qui va -S
    call drot(1, v(1),1,v(2), 1, C, S)   ! e anche qui

  endif
end subroutine cqr_fastqr6_ds_in

!-----------------------------------------------------------

!SUBROUTINE fastqr7_in
!
! This subroutine performs a step of the single shift structured QR algorithm,
! to compute the eigenvalue of a matrix which is the sum of a hermitian and a
! rank one matrix. The peculiarity of this subroutine is to use a random number
! as shift.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of
!      the matrix is UV*.

subroutine fastqr7_in(n,d,beta,u,v)

  implicit none
  integer, intent(in)  :: n
  real(8), dimension(n), intent(inout) :: d, u,v
  real(8), dimension(n-1), intent(inout) :: beta
  real(8), dimension(n-1) :: gamm
  real(8) :: rho
  real(8), dimension(4,3) :: R
  integer :: i
  real(8) :: S, C, r1,i1,r2,i2
  real(8):: z,z2
  real(8):: eps, dlamch

  eps = dlamch('e')

  if (n>2) then
    gamm(1)=(beta(1)-v(1)*u(2))+u(1)*v(2)
    ! Compute a random shift.
    call random_number(rho)

    ! Perform the structured QR step.
    z=d(1)-rho
    z2=beta(1)
    call drotg(z,z2,C,S)


    R(1,1)=d(1)
    R(2,1)=beta(1)
    R(1,2)=gamm(1)
    R(2,2)=d(2)
    R(3,1)=0
    R(3,2)=beta(2)
    R(4,:)=0
    R(:,3)=0

    call drot(2, R(1,1), 4, R(2,1), 4, C, S)
    call drot(3, R(1,1), 1, R(1,2), 1, C, (S))

    d(1)=R(1,1)
    beta(1)=R(2,1)
    d(2)=R(2,2)
    beta(2)=R(3,2)

    call drot(1, u(1), 1, u(2), 1, C, S)
    call drot(1, v(1),1,v(2), 1, C, (S))


    do i=1,n-3


      gamm(i+1)=(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
      z=beta(i)
      z2=R(3,1)

      call drotg(z,z2,C,S)
      call drot(1, beta(i), 1, R(3,1), 1, C, S)


      R(1,1)=d(i+1)
      R(2,1)=beta(i+1)
      R(1,2)=gamm(i+1)
      R(2,2)=d(i+2)
      R(3,1)=0
      R(3,2)=beta(i+2)


      call drot(2, R(1,1), 4, R(2,1), 4, C, S)
      call drot(3, R(1,1), 1, R(1,2), 1, C, (S))


      d(i+1)=R(1,1)
      beta(i+1)=R(2,1)
      d(i+2)=R(2,2)
      beta(i+2)=R(3,2)

      call drot(1, u(i+1), 1, u(i+2), 1, C, S)
      call drot(1, v(i+1),1,v(i+2), 1, C, (S))
    enddo
    gamm(n-1)=(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
    z=beta(n-2)
    z2=R(3,1)
    call drotg(z,z2,C,S)
    call drot(1, beta(n-2), 1, R(3,1), 1, C, S)

    beta(n-2)=z
    R(1,1)=d(n-1)
    R(2,1)=beta(n-1)
    R(1,2)=gamm(n-1)
    R(2,2)=d(n)


    call drot(2, R(1,1), 4, R(2,1), 4, C, S)
    call drot(2, R(1,1), 1, R(1,2), 1, C, (S))



    d(n-1)=R(1,1)
    beta(n-1)=R(2,1)
    d(n)=R(2,2)



    call drot(1, u(n-1), 1, u(n), 1, C, S)
    call drot(1, v(n-1),1,v(n), 1, C, (S))



  else if (n==2) then
    gamm(1)=(beta(1)-u(2)*v(1))+u(1)*v(2)

    call DLANV2(d(1),gamm(1),beta(1),d(2), r1,i1,r2,i2,C,S)
    call drot(1, u(1), 1, u(2), 1, C, S)  ! forse qui va -S
    call drot(1, v(1),1,v(2), 1, C, S)    ! e anche qui

  endif
end subroutine fastqr7_in
