!SUBROUTINE cqr_chasing
!
! This subroutine performs a structured bulge chasing. Moreover, if JOBH=y
! this subroutine updates the arbitrary vector h, using the Givens rotations
!created during the process.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! BULGE   COMPLEX(8). Value of the bluge.
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! JOBH CHARACTER.  JOBH=y if the arbirtary vector H has to be updated,
!		   JOB=N otherwise.
subroutine cqr_chasing(n,d,beta,u,v,bulge,h,jobh)
  implicit none
  
  integer, intent(in)  :: n
  complex(8), dimension(n), intent(inout) :: d, u,v,h
  character, intent(in):: jobh
  complex(8), dimension(n+1), intent(inout) :: beta
  complex(8) :: gamm
  complex(8), intent (inout) :: bulge
  complex(8), dimension(3,2) :: R
  complex(8) :: S, C, z
  integer :: i  


  do i=1,n-1
     gamm=conjg(beta(i+1)-v(i)*u(i+1))+u(i)*v(i+1)
     z=beta(i)

     call zrotg(z,bulge,C,S)

     beta(i)=z
     R(1,1)=d(i)
     R(2,1)=beta(i+1)
     R(1,2)=gamm
     R(2,2)=d(i+1)
     R(3,1)=0
     R(3,2)=beta(i+2)

     call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
     call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))

     d(i)=R(1,1)
     beta(i+1)=R(2,1)
     d(i+1)=R(2,2)
     beta(i+2)=R(3,2)

     call zrot(1, u(i), 1, u(i+1), 1, C, S)
     call zrot(1, v(i),1,v(i+1), 1, C, conjg(S))

     if (jobh.eq.'y') then	
	call zrot(1, h(i), 1, h(i+1), 1, C, S)
     end if

     d(i)=real(d(i)-u(i)*v(i))+(u(i)*v(i))
     d(i+1)=real(d(i+1)-u(i+1)*v(i+1))+(u(i+1)*v(i+1))

     bulge=R(3,1)
  end do
end subroutine cqr_chasing

!--------------------------------------------------
!SUBROUTINE cqr_single_eig 
!
! This subroutine computes the eigenvalues of a matrix which is the 
! sum  of a hermitian and a rank one matrix, using a structured single shift
! QR algorithm.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! K    INTEGER. Number of QR steps performed before the aggressive
!      early deflation is applied.
!
! NP    INTEGER. Number of cores available for parallelization.           
subroutine cqr_single_eig(n,d,beta,u,v,k,np)
  implicit none
  
  integer, intent(in)  :: n,np,k
  integer:: npnp,kk
  complex(8), dimension(n), intent(inout) :: d, u, v
  complex(8), dimension(n-1), intent(inout) :: beta

  npnp=np
  kk=k 
  if(n.lt.350)then
     ! Perform the structured QR algorithm without aggressive early
     ! deflation; note that we are passing a reference to u as the
     ! parameter h; since the last parameter is 'n', this won't be
     ! referenced. 
     call cqr_single_eig_small(n, d, beta, u, v, u, 'n')
  else
     ! Perform the structured QR algorithm with aggressive early
     ! deflation
     if (np.eq.1) then

	call cqr_single_eig_aed(n,d,beta,u,v,k) 
     else
	! Perform the structured QR algorithm with aggressive early
	! deflation in a parallel way.
	call cqr_single_eig_aed_par(n,d,beta,u,v,kk,npnp) 
     end if
  end if
end subroutine cqr_single_eig

!--------------------------------------------------
!SUBROUTINE cqr_single_eig_aed
!
! This subroutine computes the eigenvalues of a matrix which is the 
! sum  of a hermitian and a rank one matrices, using a structured single shift
! QR algorithm with a structured aggressive early deflation.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! K    INTEGER. Number of QR steps performed before the aggressive
!      early deflation is applied.

subroutine cqr_single_eig_aed(n,d,beta,u,v,k)
  implicit none
  
  integer, intent(in)  :: n,k
  integer :: imin, imax ,its, cont,i
  complex(8), dimension(n), intent(inout) :: d, u, v
  complex(8), dimension(n-1), intent(inout) :: beta
  complex(8), dimension(k) :: rho
  double precision :: eps, dlamch 
  real(8):: z
  real:: start, finish

  eps=dlamch('e')

  imax=n
  imin=1 
  cont=0
  rho=0

  ! Compute the first shift vector, of size k, using the Aggressive early deflation.
  call cqr_aggressive_deflation(n,d,beta,u,v,k,rho)  

  ! Try to do some deflation.
  do while ( imin.lt.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
     beta(imin)=0
     imin = imin + 1
     cont=0
  end do
  do while (imin .lt. imax .and. (beta(imax-1)==0))
     imax = imax - 1
     cont=0
  end do

  its=1
  do while (imax-imin .ge. 350)
     its=its+1
     ! Try to do some deflation.
     i=imin+1
     do while (i.ge.imin+1 .and. i.le.imax-1)
        if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
           beta(i)=0
           ! If a deflation occurs in the middle of the matrix, 
           ! compute the eigenvalues of the smallest diagonal block, 
           ! using a recursive structured QR algorithm. 
           if (i.le. (imax-imin)/2) then
              call cqr_single_eig_small(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i), u(imin:i), 'n')
              imin=i+1
           else
               call cqr_single_eig_small(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax), u(imin:i),'n')
              imax=i
           end if
        end if
        i=i+1
     end do

    
     ! Perform k steps of the structured QR algorithm using
     ! k shifts that are given as input, and return a shift vector of size k.
     call  cqr_multishift_sweep(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho)

     ! Try to do some deflation.
     do while ( imin.lt.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
        beta(imin)=0
        imin = imin + 1
        cont=0
     end do
     do while (imin.lt.imax .and. beta(imax-1)==0)
        imax = imax - 1
        cont=0
     end do

     cont=cont+1

     ! If after some QR iteration there is not delation, compute a random shift vector.
     if (cont==10) then
        do i=1,k
           call random_number(z)
           rho(i)=z
        end do
        ! Perform k seps of the structured QR algorithm using
        ! k shifts that are given as input, and return a shift vector of size k.
        call cqr_multishift_sweep(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho) 
        cont=0
     end if
  end do
  
  ! When the size of the matrix becames small, perform a structured QR algorithm
  ! without aggressive early deflation.
  call cqr_single_eig_small(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax), u(imin:imax), 'n')
  
end subroutine cqr_single_eig_aed

!----------------------------------------------------
!SUBROUTINE cqr_single_sweep 
!
! This subroutine performs a step of the single shift structured QR algorithm,
! where the shift is given as an input,  to compute the eigenvalue of a matrix 
!which is the sum of a hermitian and a  rank one matrix.  Moreover, if JOBH=y
! this subroutine updates the arbitrary vector h, using the Givens rotations
!created during the process.
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
!
! RHO  COMPLEX(8). Value of the shift.
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! JOBH CHARACTER.  JOBH=y if the arbirtary vector H has to be updated,
!		   JOB=N otherwise.
subroutine cqr_single_sweep(n,d,beta,u,v,rho,h,jobh)
  implicit none
  
  integer, intent(in)  :: n
  character, intent(in) ::jobh
  complex(8), dimension(n), intent(inout) :: d, u,v,h
  complex(8), dimension(n-1), intent(inout) :: beta
  complex(8) :: gamm, bulge
  complex(8), dimension(2) :: l
  complex(8), intent(in) :: rho
  complex(8), dimension(2,2) :: R
  complex(8) :: S, C
  complex(8):: z
  double precision :: eps, dlamch


  eps=dlamch('e')

  if (n>2) then
     ! Perform the structured QR step.
     call cqr_create_bulge_single(d(1:2),beta(1:2),u(1:2),v(1:2),rho,bulge,h(1:2),jobh)

     call cqr_chasing(n-2,d(2:n-1),beta,u(2:n-1),v(2:n-1),bulge,h(2:n-1),jobh) 

     call cqr_delete_bulge_single(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),bulge,h(n-1:n),jobh)
  
  
  else if (n==2) then

  	gamm=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
  	z=d(1)-rho
  	
  	R(1,1)=d(1)
  	R(2,1)=beta(1)
  	R(1,2)=gamm
  	R(2,2)=d(2)
  	
  	call zrotg(z,beta(1),C,S)

  	call zrot(2, R(1,1), 2, R(2,1), 2, C, S)
 	call zrot(2, R(1,1), 1, R(1,2), 1, C, conjg(S))

  	d(1)=R(1,1)
  	beta(1)=R(2,1)
  	d(2)=R(2,2)
  
  	call zrot(1, u(1), 1, u(2), 1, C, S)
  	call zrot(1, v(1),1,v(2), 1, C, conjg(S))

  	if (jobh.eq.'y') then	
     		call zrot(1, h(1), 1, h(2), 1, C, S)
  	end if

  	d(1)=real(d(1)-u(1)*v(1))+(u(1)*v(1))
  	d(2)=real(d(2)-u(2)*v(2))+(u(2)*v(2))
 end if
end subroutine cqr_single_sweep


!--------------------------------------------
!SUBROUTINE cqr_single_eig_small
!
! This subroutine performs a single shift structured QR algorithm, without 
! aggressive early deflation, to compute the eigenvalue of a matrix which is
! the sum of a hermitian and a rank one matrix. Moreover, if JOBH=y
! this subroutine computes the vector Q*H, where Q is the unitary matrix 
! obtaied from the QR decomposition of the input matrix, and H is a vector 
! given in input.
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
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! JOBH CHARACTER.  JOBH=y if the arbirtary vector H has to be updated,
!		   JOB=N otherwise.
recursive subroutine cqr_single_eig_small(n,d,beta,u,v,h,jobh)
  implicit none
  integer, intent(in)  :: n
  integer :: imin, imax ,cont,i
  character, intent(in):: jobh
  complex(8), dimension(n), intent(inout) :: d, u, v, h
  complex(8), dimension(n-1), intent(inout) :: beta
  complex(8), dimension(2) :: l
  complex(8) :: rho
  real(8) :: z
  double precision,  dimension(n)::w,ww
  double precision :: eps, dlamch

  eps=dlamch('e')
  
  !open(11, file = 'w.dat')
  !read(11, fmt = *) w
  !close(11)

  imax=n
  imin=1 
  cont=0
  !Try to deflate some eigenvalue
  do while ( imin.lt.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
     beta(imin)=0
     imin = imin + 1
  end do
  do while (imin.lt.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
     beta(imax-1)=0
     imax = imax - 1
  end do
  do while (imax-imin .gt. 0)
  

  ! do i=1,n
  !ww(i)=abs(v(i))/w(i)
  !end do
  
  !print*, norm2(ww)/n
  
     call cqr_wilkinson_shift(d(imax-1:imax),beta(imax-1),u(imax-1:imax),v(imax-1:imax),rho)

     call cqr_single_sweep(imax-imin+1,d(imin:imax),beta(imin:imax-1),u(imin:imax),v(imin:imax),rho,h(imin:imax),jobh)
    

     !Try to deflate some eigenvalue
     do while (imin.lt.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
        beta(imin)=0
        imin = imin + 1
        cont=0
     end do
     do while (imin.lt.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
        beta(imax-1)=0
        imax = imax - 1
        cont=0
     end do

     i=imin+1  
     do while (i.ge.imin+1 .and. i.le.imax-1)
        if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
           beta(i)=0
           ! If a deflation occurs in the middle of the matrix, 
           ! compute the eigenvalues of the smallest diagonal block, 
           ! using a recursive structured QR algorithm. 
           if (i.le. (imax-imin)/2) then
              call cqr_single_eig_small(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i),h(imin:i),jobh)
              imin=i+1
           else
              call cqr_single_eig_small(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax),h(i+1:imax),jobh)
              imax=i
           end if
        end if
        i=i+1
     end do

     cont=cont+1

     ! If after some QR iteration there is not delation, perform a structured
     ! QR step using a random shift vector.
     if (cont==10) then
	call random_number(z)
	rho=z
	call cqr_single_sweep(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,h(imin:imax),jobh)
        cont=0
     end if

  end do
end subroutine cqr_single_eig_small

!-----------------------------------------------

!SUBROUTINE cqr_wilkinson shift
!
! This subroutine return an eigenvalue of a (structured) 2x2 matrix for the 
! purpose of using it as a shift.
!
! INPUT PARAMETERS
!
! D    COMPLEX(8), DIMENSION(2). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(1). Vector that contains the subdiagonal 
!      entry of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(2). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! RHO  COMPLEX(8). Value of the output eigenvalue.


 subroutine cqr_wilkinson_shift(d,beta,u,v,rho)
 complex(8), dimension(2), intent(in):: d, u,v
 complex(8), dimension(1), intent(in):: beta
 complex(8), intent(out):: rho
 complex(8), dimension(2):: dd, uu, vv,l
 complex(8),dimension(1)::b
 double precision :: eps, dlamch

  eps=dlamch('e')
  

 rho=sqrt((d(1)+d(2))**2-4*(d(1)*d(2)-(beta(1)*(conjg(beta(1)-u(2)*v(1))+u(1)*v(2)))))
 l(1)=(d(1)+d(2)+rho)/2
 l(2)=(d(1)+d(2)-rho)/2
 dd=d
 uu=u
 vv=v
 b=beta
 do while (abs(b(1)).ne.0)
 	if (abs(l(1)-dd(2))<abs(l(2)-dd(2))) then
       		rho=l(1);
	else
       		rho=l(2);
	endif
	call cqr_single_sweep(2,dd,b,uu,vv,rho,uu,'n')
	if(abs(b(1)).lt.eps*(abs(dd(1))+abs(dd(2)))) then
		b(1)=0
	end if
 end do
 if (abs(dd(1)-d(2))<abs(dd(2)-d(2))) then
	rho=dd(1);
 else
       	rho=dd(2);
 endif
 end subroutine



!------------------------------------------------------
!SUBROUTINE cqr_multishift_sweep
!
! This subroutine performs h steps of the structured QR algorithm, using
! h shifts that are given as input, and returns a shift vector of size k.
!
! INPUT PARAMETERS
!
! NN    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! H   INTEGER. Number of input shifts.
!
! K   INTEGER. Number of output shifts, K>=H.
!
! RHRH   COMPLEX(8), DIMENSION(k). Vector that contains the shifts.
subroutine cqr_multishift_sweep(nn,d,beta,u,v,h,k, RHRH)
  implicit none
  
  integer, intent(in)  :: nn,k,h
  integer :: n
  complex(8), dimension(k), intent(inout) :: RHRH
  complex(8), dimension(nn), intent(inout) :: d, u,v
  complex(8), dimension(nn-1), intent(inout) :: beta
  complex(8), dimension(nn-1) :: gamm
  complex(8), dimension(2) :: l
  complex(8) :: rho
  complex(8), dimension(3,2) :: R
  integer :: i,p
  complex(8) :: S, C
  complex(8) :: z
  double precision :: eps, dlamch

  eps=dlamch('e')

  n=nn

  if (n>3/2*k) then
     ! Perform h steps of the structured QR algorithm.
     do p=1,h
        rho=RHRH(p)
        call cqr_single_sweep(n,d,beta,u,v,rho,u,'n')
     end do
     ! Try to do some deflation in the final part of the matrix.
     do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
        beta(n-1)=0
        n=n-1
     end do
     ! Perform a step of the aggressive early deflation and compute the new 
     ! shift vector.
     call cqr_aggressive_deflation(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RHRH)
  else
     ! If the size of the matrix is small, compute the egenvalues performing
     ! a structured QR algorithm without aggressive early deflation.
     call cqr_single_eig_small(n,d,beta,u,v,u,'n')
     RHRH=0
  end if
end subroutine cqr_multishift_sweep


!-------------------------------------------------------------
!SUBROUTINE cqr_aggressive_deflation
!
! This subroutine performs the structured aggressive early deflation.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! W   INTEGER. Number of output shifts.
!
! RHO   COMPLEX(8), DIMENSION(k). Vector that contains the output shifts.
recursive subroutine cqr_aggressive_deflation (n,d,beta,u,v,w,RHO)
  implicit none
  
  integer, intent(in)  :: n,w
  complex(8), dimension(n), intent(inout) :: d, u,v
  complex(8), dimension(n-1), intent(inout) :: beta
  complex(8), dimension(w*3/2) :: h,hatd
  complex(8), dimension (w), intent(out) :: RHO
  complex(8) :: S, C 
  integer :: f,i,j,K, p
  complex(8) :: l
  complex(8) :: z
  complex, dimension(2,2) :: G
  double precision :: eps, dlamch



  eps=dlamch('e')

  K=w*3/2
  if (n.gt.K) then
     h=0
     h(1)=beta(n-K)
     ! Compute the structured Schur form of the kxk trailing principal submatrix, updating the 
     ! vector h.

     call cqr_single_eig_small(K, d(n-K+1:n),beta(n-K+1:n-1),u(n-K+1:n),v(n-K+1:n),h(1:K),'y')

     i=K
     j=0
     ! Try to deflate some eigenvalue from the k+1xk+1 trailing principal submatrix
     do while ((i .gt. 0) .AND. (j .lt. i))

	if (abs(beta(n-K)).lt. abs(d(i+n-K))) then
           l=beta(n-K)
	else
           l=d(i+n-K)
	end if
        if (abs(h(i)) .lt. eps*abs(l)) then
           h(i)=0
           i=i-1
        else
           j=j+1
           do f=i-1,1,-1
              z=d(f+1+n-K)-d(f+n-K)
              call zrotg(z,u(f+n-K)*v(f+1+n-K)-conjg(u(f+1+n-K)*v(f+n-K)),C,S)
              l=u(f+n-K)
              u(f+n-K)=u(f+n-K+1)
              u(f+n-K+1)=l
              call zrot(1, u(f+n-K), 1, u(f+1+n-K), 1, C,S)
              l=h(f)
              h(f)=h(f+1)
              h(f+1)=l
              call zrot(1, h(f), 1, h(f+1), 1, C, S)
              l=v(f+n-K)
              v(f+n-K)=v(f+n-K+1)
              v(f+n-K+1)=l
              call zrot(1, v(f+n-K),1,v(f+1+n-K), 1, C, conjg(S))
              l=d(f+n-K)
              d(f+n-K)=d(f+n-K+1)
              d(f+n-K+1)=l
           end do
        end if
     end do
     ! Store the non deflated eigenvalues of the kxk trailing principal submatrix.
     hatd(1:i)=d(1+n-K:i+n-K)
     ! Bring back the matrix in Hessenberg form.
     call cqr_hessenberg_reduction(i,h(1:i),d(1+n-K:i+n-K),beta(n-K+1:n-K+i-1),u(1+n-K:i+n-K),v(1+n-K:i+n-K)) 
     beta(n-K)=h(1)
     if (i< w) then
        ! The stored eigenvalues are not enough to produce w shifts, hence perform
        ! a new aggressive deflation step. 
        call cqr_aggressive_deflation(n-K+i,d(1:n-K+i),beta(1:n-K+i-1),u(1:n-K+i),v(1:n-K+i),w,RHO)
     else
  	! Store the smallest (in magnitude) w elements of hatd as shifts.
	call cqr_sort(i,hatd(1:i))
	RHO=hatd(1:w)
     end if
  else
     ! If the size of the matrix is small, compute the egenvalues performing
     ! a structured QR algorithm without aggressive early deflation.
     call cqr_single_eig_small(n,d,beta,u,v,u,'n')
     RHO=0
  end if

end subroutine cqr_aggressive_deflation


!--------------------------------------------------------
!SUBROUTINE cqr_hessenberg_reduction  
!
! This subroutine takes in input an arbitrary vector and the structured
! representation of a triangular matrix that is the sum of a Hermitian
! and a rank one matrix. The subroutine annihilates all the components 
! of the arbitrary vector, except the first, by means of Givens rotations.
! Moreover the subroutine recuces in Hessenberg form the triangular matrix 
! transformed by similarity with the Givens rotations used before.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix. 
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary vector.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
subroutine cqr_hessenberg_reduction(n,h,d,beta,u,v) 
  implicit none

  integer, intent(in)  :: n
  complex(8), dimension(n), intent(inout) :: d, u,v
  complex(8), dimension(n-1), intent(out) :: beta
  complex(8), dimension(n), intent(inout) :: h
  complex(8), dimension(3,2) :: R
  integer :: i,j,p
  complex(8) :: C,S
  complex(8), dimension(2,2) :: G
  complex(8) :: z,gamm
  if (n>1) then
     ! Annihilate the last component of the arbitrary vector.
     z=h(n-1)
     call zrotg(z,h(n),C,S)
     h(n-1)=z
     h(n)=0
     ! Transform the triangular matrix by means the previous 
     ! Givens rotation.
     gamm=conjg(-v(n-1)*u(n))+u(n-1)*v(n)
     R(1,1)=d(n-1)
     R(2,2)=d(n)
     R(1,2)=gamm
     R(2,1)=0
     R(3,1)=0
     R(3,2)=0
     call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
     call zrot(2, R(1,1), 1, R(1,2), 1, C, conjg(S))
     call zrot(1, u(n-1), 1, u(n), 1, C, S)
     call zrot(1, v(n-1), 1, v(n), 1, C, conjg(S))
     d(n-1)=R(1,1)
     d(n)=R(2,2)
     beta(n-1)=R(2,1)
     do i=n-2,1,-1
        ! Annihilate the last nonzero entry of the arbitrary vector.
	z=h(i)
	call zrotg(z,h(i+1),C,S)
        h(i)=z
        h(i+1)=0
        ! Transform the triangular matrix by means the previous 
        ! Givens rotation.
	gamm=conjg(-v(i)*u(i+1))+u(i)*v(i+1)
	R(1,1)=d(i)
	R(2,2)=d(i+1)
	R(1,2)=gamm
	R(2,1)=0
	R(3,1)=0
	R(3,2)=beta(i+1)
	call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
	call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))
	call zrot(1, u(i), 1, u(i+1), 1, C, S)
	call zrot(1, v(i), 1, v(i+1), 1, C, conjg(S))
	d(i)=R(1,1)
	d(i+1)=R(2,2)
	beta(i)=R(2,1)
	beta(i+1)=R(3,2)

        ! Chase the bulge to bring back the matrix in Hessenberg form.
        
        call cqr_chasing(n-i-1,d(i+1:n-1),beta(i:n-1),u(i+1:n-1),v(i+1:n-1),R(3,1),u(i+1:n-1),'n')
        
	call cqr_delete_bulge_single(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),R(3,1),u(n-1:n),'n')

     end do
  end if
end subroutine cqr_hessenberg_reduction


!------------------------------------------------------------------------------
!SUBROUTINE cqr_sort
! This subroutine sorts an input vector in decreasing order with respect
! to the magnitute of its entries.  
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input vector. 
!
! A    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
subroutine cqr_sort(n,a) 
  implicit none 
  integer, intent(in) :: n
  complex(8), dimension(n), intent(inout) :: a
  complex(8) :: b
  integer :: i
  logical :: keepgoing


  keepgoing=.true.
  do while (keepgoing)
     keepgoing=.false.
     do i=1,n-1
	if (abs(a(i))>abs(a(i+1))) then 
           b=a(i+1)
           a(i+1)=a(i)
           a(i)=b
           keepgoing=.true.
	end if
     end do
  end do
end subroutine cqr_sort

!--------------------------------------------------------------------
!SUBROUTINE cqr_create_bulge_single 
! This subroutine creates the first bulge.  Moreover, if JOBH=y
! this subroutine updates the arbitrary vector h, using the Givens rotations
!created during the process.
!
!
! INPUT PARAMETERS
!
! D      COMPLEX(8), DIMENSION(2). Vector that contains the first 2 
!        diagonal entries of the matrix. 
!
! BETA   COMPLEX(8), DIMENSION(2). Vector that contains the first 2
!        subdiagonal entries of the matrix.
!
! U,V    COMPLEX(8), DIMENSION(N). Vectors that contain the first 2
!        entries of the rank-one vectors.
!
! RHO    COMPLEX(8). Value of the shift.
!
! BULGE  COMPLEX(8). Value of the bulge.
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! JOBH CHARACTER.  JOBH=y if the arbirtary vector H has to be updated,
!		   JOB=N otherwise.
subroutine cqr_create_bulge_single(d,beta,u,v,rho,bulge,h,jobh)
  implicit none
  complex(8), dimension(2), intent(inout) :: d,beta,u,v,h
  character, intent(in) :: jobh
  complex(8), intent(in):: rho
  complex(8), intent(out):: bulge
  complex(8):: z, gamm, S, C
  complex(8), dimension(3,2) :: R

  gamm=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
  z=d(1)-rho
  
  R(1,1)=d(1)
  R(2,1)=beta(1)
  R(1,2)=gamm
  R(2,2)=d(2)
  R(3,1)=0
  R(3,2)=beta(2)
  
  call zrotg(z,beta(1),C,S)


  call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
  call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))

  d(1)=R(1,1)
  beta(1)=R(2,1)
  d(2)=R(2,2)
  beta(2)=R(3,2)

  call zrot(1, u(1), 1, u(2), 1, C, S)
  call zrot(1, v(1),1,v(2), 1, C, conjg(S))

  if (jobh.eq.'y') then	
     call zrot(1, h(1), 1, h(2), 1, C, S)
  end if

  d(1)=real(d(1)-u(1)*v(1))+(u(1)*v(1))
  d(2)=real(d(2)-u(2)*v(2))+(u(2)*v(2))

  bulge=R(3,1)

end subroutine cqr_create_bulge_single

!-------------------------------------------------------------------
!SUBROUTINE cqr_delete_bulge_single 
! This subroutine deletes a bulge that is in the last row of the matrix. 
! Moreover, if JOBH=y this subroutine updates the arbitrary vector h, 
! using the Givens rotations created during the process.
!
!
! INPUT PARAMETERS
!
! D      COMPLEX(8), DIMENSION(2). Vector that contains the last 2 
!        diagonal entries of the matrix. 
!
! BETA   COMPLEX(8), DIMENSION(2). Vector that contains the last 2
!        subdiagonal entries of the matrix.
!
! U,V    COMPLEX(8), DIMENSION(N). Vectors that contain the last 2
!        entries of the rank-one vectors.
!
! BULGE  COMPLEX(8). Value of the bulge.
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! JOBH CHARACTER.  JOBH=y if the arbirtary vector H has to be updated,
!		   JOB=N otherwise.
subroutine cqr_delete_bulge_single (d,beta,u,v,bulge,h,jobh)
  implicit none
  complex(8), dimension(2), intent(inout) :: d,beta,u,v,h
  character, intent(in)::jobh
  complex(8), intent(in):: bulge
  complex(8):: z, gamm, C, S
  complex(8), dimension(2,2) :: R

  gamm=conjg(beta(2)-v(1)*u(2))+u(1)*v(2);
  z=beta(1)
  call zrotg(z,bulge,C,S)

  beta(1)=z
  R(1,1)=d(1)
  R(2,1)=beta(2)
  R(1,2)=gamm
  R(2,2)=d(2)

  call zrot(2, R(1,1), 2, R(2,1), 2, C, S)
  call zrot(2, R(1,1), 1, R(1,2), 1, C, conjg(S))

  d(1)=R(1,1)
  beta(2)=R(2,1)
  d(2)=R(2,2)

  call zrot(1, u(1), 1, u(2), 1, C, S)
  call zrot(1, v(1),1,v(2), 1, C, conjg(S))

  if (jobh .eq. 'y') then	
     call zrot(1, h(1), 1, h(2), 1, C, S)
  end if

  d(1)=real(d(1)-u(1)*v(1))+(u(1)*v(1))
  d(2)=real(d(2)-u(2)*v(2))+(u(2)*v(2))

end subroutine cqr_delete_bulge_single

!-----------------------------------------------------------------------------------------------------------------------------------------
!PARALLELIZATION
!-----------------------------------------------------------------------------------------------------------------------------------------

!SUBROUTINE cqr_single_sweep_par
!
! This subroutine performs in parallel k steps of the single shift structured 
! QR algorithm, using k shifts that are given as input.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! RH   COMPLEX(8), DIMENSION(K). Vector that contains the shifts.
!
! K    INTEGER. Number of input shifts.
!
! NP   INTEGER. Number of cores available for parallelization. 
subroutine cqr_single_sweep_par(n,d,beta,u,v,rh,k,np)
  implicit none
  
  integer, intent(in)  :: n,k,np
  integer :: siz
  complex(8), dimension(k), intent(inout) :: RH
  complex(8), dimension(n), intent(inout) :: d, u,v
  complex(8), dimension(n-1), intent(inout) :: beta
  complex(8), dimension(n-1) :: gamm
  complex(8), dimension(2) :: l
  complex(8) :: rho
  complex(8), dimension(3,2) :: R
  integer :: i,p,q,ii,iii, qq(np)
  complex(8) :: S, C
  complex(8) :: z
  double precision :: eps, dlamch
  complex(8), dimension(np):: bulge 
  real :: start, finish
  logical:: check=.false.

  eps=dlamch('e')

  siz=(n-3)/np

  if (n>k*3/2) then
     do p=1,np-1
	if (modulo(p,10).eq.0)then
           check=.true.
	else
           check=.false.
	end if
 
        rho=RH(p)

        call cqr_create_bulge_single(d(1:2),beta(1:2),u(1:2),v(1:2),rho,bulge(1),u(1:2),'n')

        !$OMP PARALLEL DO
        do q=1,p 
           qq(q)=(q-1)*(siz)+1
           call cqr_chasing(siz-1, d(qq(q)+1:qq(q)+siz-1), beta(qq(q):qq(q)+siz-1), u(qq(q)+1:qq(q)+siz-1), &
                v(qq(q)+1:qq(q)+siz-1),bulge(q),u(qq(q)+1:qq(q)+siz-1),'n')

           if (check) then
              do i=qq(q),qq(q)+siz-3
                 if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                    beta(i)=0
                 end if
              end do
           end if
        end do
        !$OMP END PARALLEL DO

        do q=p,1,-1
           call cqr_chasing(3, d(qq(q)+siz-1:qq(q)+siz+1), beta(qq(q)+siz-2:qq(q)+siz+1), u(qq(q)+siz-1:qq(q)+siz+1), &
                v(qq(q)+siz-1:qq(q)+siz+1),bulge(q),u(qq(q)+siz-1:qq(q)+siz+1),'n')
           bulge(q+1)=bulge(q)

           if (check) then
              do i=qq(q)+siz-2,qq(q)+siz-1
                 if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                    beta(i)=0
                 end if
              end do
           end if

        end do

     end do


     do p=np,k

	if (modulo(p,10).eq.0)then
           check=.true.
	else
           check=.false.
	end if

        rho=RH(p)
        call cqr_create_bulge_single(d(1:2),beta(1:2),u(1:2),v(1:2),rho,bulge(1),u(1:2),'n')

        !$OMP PARALLEL DO		
        do q=1,np 		
           qq(q)=(q-1)*(siz)+1

           call cqr_chasing(siz-1, d(qq(q)+1:qq(q)+siz-1), beta(qq(q):qq(q)+siz-1), u(qq(q)+1:qq(q)+siz-1), &
		v(qq(q)+1:qq(q)+siz-1),bulge(q), u(qq(q)+1:qq(q)+siz-1),'n')


           if (check) then
              do i=qq(q),qq(q)+siz-3
                 if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                    beta(i)=0
                 end if
              end do
           end if

        end do
        !$OMP END PARALLEL DO

        do q=np,1,-1
           call cqr_chasing(3, d(qq(q)+siz-1:qq(q)+siz+1), beta(qq(q)+siz-2:qq(q)+siz+1), u(qq(q)+siz-1:qq(q)+siz+1), &
                v(qq(q)+siz-1:qq(q)+siz+1),bulge(q),u(qq(q)+siz-1:qq(q)+siz+1),'n')

           if (check) then
              do i=qq(q)+siz-2,qq(q)+siz-1
                 if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                    beta(i)=0
                 end if
              end do
           end if
        end do

        call cqr_chasing(n-(qq(np)+siz+1),d(qq(np)+siz+1:n-1),beta(qq(np)+siz:n-1),u(qq(np)+siz+1:n-1),v(qq(np)+siz+1:n-1), &
             bulge(np),u(qq(np)+siz+1:n-1),'n')

        if (check) then
           do i=qq(np)+siz,n-1
              if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                 beta(i)=0
              end if
           end do
        end if

        call cqr_delete_bulge_single(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),bulge(np),u(n-1:n),'n')

        do q=np-1,1,-1
           bulge(q+1)=bulge(q)
        end do
     end do

     do p=2,np
	if (modulo(p,10).eq.0)then
           check=.true.
	else
           check=.false.
	end if
 
        !$OMP PARALLEL DO
        do q=p,np 
           qq(q)=(q-1)*(siz)+1
           call cqr_chasing(siz-1, d(qq(q)+1:qq(q)+siz-1), beta(qq(q):qq(q)+siz-1), u(qq(q)+1:qq(q)+siz-1), &
                v(qq(q)+1:qq(q)+siz-1),bulge(q),u(qq(q)+1:qq(q)+siz-1),'n')

           if (check) then
              do i=qq(q),qq(q)+siz-3
                 if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                    beta(i)=0
                 end if
              end do
           end if
        end do
        !$OMP END PARALLEL DO

        do q=p,np
           call cqr_chasing(3, d(qq(q)+siz-1:qq(q)+siz+1), beta(qq(q)+siz-2:qq(q)+siz+1), u(qq(q)+siz-1:qq(q)+siz+1), &
                v(qq(q)+siz-1:qq(q)+siz+1),bulge(q),u(qq(q)+siz-1:qq(q)+siz+1),'n')

           if (check) then
              do i=qq(q)+siz-2,qq(q)+siz-1
                 if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                    beta(i)=0
                 end if
              end do
           end if

        end do
        call cqr_chasing(n-(qq(np)+siz+1),d(qq(np)+siz+1:n-1),beta(qq(np)+siz:n-1),u(qq(np)+siz+1:n-1),&
             v(qq(np)+siz+1:n-1),bulge(np),u(qq(np)+siz+1:n-1),'n')

        if (check) then
           do i=qq(np)+siz,n-1
              if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
                 beta(i)=0
              end if
           end do
        end if

        call cqr_delete_bulge_single(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),bulge(np),u(n-1:n),'n')

        do q=np-1,p,-1
           bulge(q+1)=bulge(q)
        end do
     end do
  else
     call cqr_single_sweep(n,d,beta,u,v,rho,u,'n')
  end if
end subroutine cqr_single_sweep_par

!---------------------------------------------------------
!SUBROUTINE cqr_multishift_sweep_par
!
! This subroutine performs in a parallel way h steps of the structured 
! QR algorithm, using h shifts that are given as input, and returns a 
! shift vector of size k.
!
! INPUT PARAMETERS
!
! NN    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! H   INTEGER. Number of input shifts.
!
! K   INTEGER. Number of output shifts, K>=H.
!
! RHRH   COMPLEX(8), DIMENSION(k). Vector that contains the shifts.
!
! NP    INTEGER. Number of cores available for parallelization. 
subroutine cqr_multishift_sweep_par(nn,d,beta,u,v,h,k, RHRH,np)
  implicit none
  
  integer, intent(in)  :: nn,k,h,np
  integer :: n
  complex(8), dimension(k), intent(inout) :: RHRH
  complex(8), dimension(nn), intent(inout) :: d, u,v
  complex(8), dimension(nn-1), intent(inout) :: beta
  complex(8), dimension(nn-1) :: gamm
  complex(8), dimension(2) :: l
  complex(8) :: rho
  complex(8), dimension(3,2) :: R
  integer :: i,p
  complex(8) :: S, C
  complex(8) :: z
  double precision :: eps, dlamch

  eps=dlamch('e')

  n=nn

  if (n>3/2*k) then
     ! Perform h steps of the structured QR algorithm.
     call cqr_single_sweep_par(n,d,beta,u,v,RHRH,h,np)

     ! Try to do some deflation in the final part of the matrix.
     do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
        beta(n-1)=0
        n=n-1
     end do
     
     ! Perform a step of the aggressive early deflation and compute the new 
     ! shift vector.
     call cqr_aggressive_deflation(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RHRH)
  else
     ! If the size of the matrix is small, compute the egenvalues performing
     ! a structured QR algorithm without aggressive early deflation.
     call cqr_single_eig_small(n,d,beta,u,v,u,'n')
     
     RHRH=0
  end if
end subroutine cqr_multishift_sweep_par

!-------------------------------------------------------------------------------------------
!SUBROUTINE cqr_single_eig_aed_par
!
! This subroutine computes the eigenvalues of a matrix which is the 
! sum  of a hermitian and a rank one matrices, using a structured single shift
! QR algorithm with a structured aggressive early deflation run in a parallel 
! way.
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix 
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.
!
! K    INTEGER. Number of QR steps performed before the aggressive
!      early deflation is applied.
!
! NP    INTEGER. Number of cores available for parallelization. 
subroutine cqr_single_eig_aed_par(n,d,beta,u,v,k,np)
  implicit none
  integer, intent(in)  :: n
  integer, intent(inout) :: np,k
  integer :: imin, imax ,its, cont,i
  complex(8), dimension(n), intent(inout) :: d, u, v
  complex(8), dimension(n-1), intent(inout) :: beta
  complex(8), dimension(k) :: rho
  double precision :: eps, dlamch
  real(8):: z
  real:: finish, start

  imax=n
  imin=1 
  cont=0
  rho=0

  eps=dlamch('e')

  ! Compute the first shift vector, of size k, using the Aggressive early deflation.
  call cqr_aggressive_deflation(n,d,beta,u,v,k,rho)

  ! Try to do some deflation.
  do while ( imin.lt.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
     beta(imin)=0
     imin = imin + 1
     cont=0
  end do
  do while (imin .lt. imax .and. beta(imax-1)==0 )
     imax = imax - 1
     cont=0
  end do
  its=1
  do while (imax-imin .ge. 350)
     do while (imax-imin .ge. 64*np)

	its=its+1
 
        ! Try to do some middle deflation.
     i=imin+1
     do while (i.ge.imin+1 .and. i.le.imax-1)
        if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
           beta(i)=0
           ! If a deflation occurs in the middle of the matrix, 
           ! compute the eigenvalues of the smallest diagonal block, 
           ! using a recursive structured QR algorithm. 
           if (i.le. (imax-imin)/2) then
              call cqr_single_eig_small(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i),u(imin:i),'n')
              imin=i+1
           else
               call cqr_single_eig_small(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax),u(i+1:imax),'n')
              imax=i
           end if
        end if
        i=i+1
     end do
 
        ! Perform k steps of the structured QR algorithm using
        ! k shifts that are given as input, and return a shift vector of size k.
        call  cqr_multishift_sweep_par(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho,np) 

        ! Try to do some deflation.
	do while (imin.lt.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
           beta(imin)=0
           imin = imin + 1
           cont=0
	end do
	do while (imin.lt.imax .and. beta(imax-1)==0)
           imax = imax - 1
           cont=0
	end do


        ! If after some QR iteration there is not delation, compute a random shift vector.
	if (cont==10) then
           do i=1,k
              call random_number(z)
              rho(i)=z
           end do
           ! Perform k seps of the structured QR algorithm using
           ! k shifts that are given as input, and return a shift vector of size k.
           call cqr_multishift_sweep_par(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho,np) 
           cont=0
	end if
     end do
     np=np/2
     k=np*6
  end do
  
  ! When the size of the matrix becames small, perform a structured QR algorithm
  ! without aggressive early deflation.
  call cqr_single_eig_small(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),u(imin:imax),'n')
end subroutine cqr_single_eig_aed_par

