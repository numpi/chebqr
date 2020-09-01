!SUBROUTINE chasing
!
! This subroutine performs a structured bulge chasing.
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



subroutine chasing(n,d,beta,u,v,bulge)
implicit none
integer, intent(in)  :: n
complex(8), dimension(n), intent(inout) :: d, u,v
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
	
	d(i)=real(d(i)-u(i)*v(i))+(u(i)*v(i))
	d(i+1)=real(d(i+1)-u(i+1)*v(i+1))+(u(i+1)*v(i+1))

	bulge=R(3,1)
end do
end subroutine chasing

!--------------------------------------------------
!SUBROUTINE fastfastqr
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
! T    INTEGER. Per facilitare gli esperimenti ho inserito questo parametro:
!      se t=0 viene eseguito in sequenziale, altrimenti in parallelo.           magari poi modificare la cosa in modo piu formale


subroutine fastfastqr(n,d,beta,u,v,k,np)
implicit none
integer, intent(in)  :: n,k,np
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta

 
if(n.lt.350)then

	! Perform the structured QR algorithm without aggressive early
	! deflation
	call fastqr6(n,d,beta,u,v)
else
	! Perform the structured QR algorithm with aggressive early
	! deflation
	if (np.eq.1) then
	
	call aggressive_deflation(n,d,beta,u,v,k)
	else
	
	call aggressive_deflation_par(n,d,beta,u,v,k,np)
	end if
end if
end subroutine

!--------------------------------------------------

!SUBROUTINE aggressive_deflation
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

subroutine aggressive_deflation(n,d,beta,u,v,k)
implicit none
integer, intent(in)  :: n,k
integer :: imin, imax ,its, cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(k) :: rho
double precision :: eps = 2.22e-16
real(8):: z
real:: start, finish

imax=n
imin=1 
cont=0
rho=0
! Compute the first shift vector, of size 2, using the Wilkinson shift.
!!rho(2)=sqrt((d(n-1)+d(n))**2-4*(d(n-1)*d(n)-beta(n-1)*conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))
!!rho(1)=(d(n-1)+d(n)+rho(2))/2
!!rho(2)=(d(n-1)+d(n)-rho(2))/2
! Perform 2 seps of the structured QR algorithm using
! 2 shifts that are given as input, and return a shift vector of size k.
!!call  fastqr12_in(n,d,beta,u,v,2,k,rho)

! Compute the first shift vector, of size k, using the Aggressive early deflation.
call aggressive_deflation_in(n,d,beta,u,v,k,rho)

! Try to do some deflation.
do while ( imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
	beta(imin)=0
	imin = imin + 1
	cont=0
end do
do while (imin .le. imax .and. (beta(imax-1)==0))
	imax = imax - 1
	cont=0
end do


its=1
do while (imax-imin .ge. 350)
	its=its+1
	! Try to do some deflation.
	do i=imin+1,imax-1
		if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
			beta(i)=0
				! If a deflation occurs in the middle of the matrix, 
				! compute the eigenvalues of the smallest diagonal block, 
				! using a structured QR algorithm without aggressive
				! early deflation. 
				if (i.le. (imax-imin)/2) then
				call fastqr6(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
					do while (imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
					beta(imin)=0
					imin = imin + 1
					cont=0
				end do
			else
				call fastqr6(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
					do while ( imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
					beta(imax-1)=0
					imax = imax - 1
					cont=0
				end do
			end if
		end if
	end do
		! Perform k steps of the structured QR algorithm using
                ! k shifts that are given as input, and return a shift vector of size k.
              !  call cpu_time(start)
        call  fastqr12_in(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho) 
        !call cpu_time(finish)
        !print'("time3=",f6.3," seconds")', finish-start
		! Try to do some deflation.
	do while ( imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
		beta(imin)=0
		imin = imin + 1
		cont=0
	end do
	do while (imin.le.imax .and. beta(imax-1)==0)
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
                call fastqr12_in(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho) 
		cont=0
	end if
end do
! When the size of the matrix becames small, perform a structured QR algorithm
! without aggressive early deflation.
call fastqr6(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
end subroutine aggressive_deflation

!----------------------------------------------------

!SUBROUTINE fastqr_ss_in
!
! This subroutine performs a step of the single shift structured QR algorithm,
! where the shift is given as an input,  to compute the eigenvalue of a matrix 
!which is the sum of a hermitian and a  rank one matrix. 
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

subroutine fastqr_ss_in(n,d,beta,u,v,rho)
implicit none
integer, intent(in)  :: n
complex(8), dimension(n), intent(inout) :: d, u,v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8) :: gamm, bulge
complex(8), dimension(2) :: l
complex(8), intent(in) :: rho
complex(8), dimension(2,2) :: A
integer :: i
complex(8) :: S, C
complex(8):: z
double precision :: eps = 2.22e-16
real:: finish, start

if (n>2) then
	! Perform the structured QR step.
	
	call create_bulge(d(1:2),beta(1:2),u(1:2),v(1:2),rho,bulge)
	
	call chasing(n-2,d(2:n-1),beta,u(2:n-1),v(2:n-1),bulge)
	
	call delete_bulge(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),bulge)

else
    if (n==2) then 
	gamm=conjg(beta(1)-u(2)*v(1))+u(1)*v(2)
	
	A(1,1)=d(1)
	A(2,1)=beta(1)
	A(1,2)=gamm
	A(2,2)=d(2)
	if (A(2,1).NE.0) then	
	gamm=sqrt((A(1,1)+A(2,2))**2-4*(A(1,1)*A(2,2)-A(2,1)*A(1,2)))
	l(1)=(A(1,1)+A(2,2)+gamm)/2
	l(2)=(A(1,1)+A(2,2)-gamm)/2
	if (abs(l(1)-d(2))<abs(l(2)-d(2))) then
	   d(1)=l(1);
	else
    		d(1)=l(2);
	endif
	z=A(1,1)-d(1)
	call zrotg(z,beta(1),C,S)
	call zrot(2, A(1,1), 2, A(2,1), 2 , C, S)
	call zrot(2, A(1,1), 1, A(1,2), 1, C, conjg(S))
	call zrot(1, u(1), 1, u(2), 1, C, S)
	call zrot(1, v(1),1,v(2), 1, C, conjg(S))
	d(1)=A(1,1)
	d(2)=A(2,2)
        beta(1)=0
	end if
    endif
endif
end subroutine fastqr_ss_in
!--------------------------------------------
!SUBROUTINE fastqr6
!
! This subroutine performs a single shift structured QR algorithm, without 
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

recursive subroutine fastqr6(n,d,beta,u,v)
implicit none
integer, intent(in)  :: n
integer :: imin, imax ,cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(2) :: l
complex(8) :: rho
real(8) :: z
double precision :: eps = 2.22e-16

imax=n
imin=1 
cont=0

!Try to deflate some eigenvalue
do while ( imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
beta(imin)=0
imin = imin + 1
end do
do while (imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
beta(imax-1)=0
imax = imax - 1
end do
do while (imax-imin .gt. 0)

! Compute a step of the QR algorithm.
	rho=sqrt((d(imax-1)+d(imax))**2-4*(d(imax-1)*d(imax)-(beta(imax-1)*(conjg(beta(imax-1)-u(imax)*v(imax-1))+u(imax-1)*v(imax)))))
	l(1)=(d(imax-1)+d(imax)+rho)/2
	l(2)=(d(imax-1)+d(imax)-rho)/2
	if (abs(l(1)-d(imax))<abs(l(2)-d(imax))) then
	    rho=l(1);
	else
    		rho=l(2);
	endif
	call fastqr_ss_in(imax-imin+1,d(imin:imax),beta(imin:imax-1),u(imin:imax),v(imin:imax),rho)

!Try to deflate some eigenvalue
do while (imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
beta(imin)=0
imin = imin + 1
cont=0
end do
do while (imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
beta(imax-1)=0
imax = imax - 1
cont=0
end do

do i=imin+1,imax-1
if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
beta(i)=0
! If a deflation occurs in the middle of the matrix, 
! compute the eigenvalues of the smallest diagonal block, 
! using a recursive structured QR algorithm. 
if (i.le. (imax-imin)/2) then
call fastqr6(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
do while (imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
beta(imin)=0
imin = imin + 1
cont=0
end do
else
call fastqr6(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
do while (imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
beta(imax-1)=0
imax = imax - 1
cont=0
end do
end if
end if
end do




cont=cont+1

! If after some QR iteration there is not delation, perform a structured
! QR step using a random shift vector.
if (cont==10) then
	call random_number(z)
	rho=z
	call fastqr_ss_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho)
cont=0
end if

end do

end subroutine fastqr6
!------------------------------------------------------

!SUBROUTINE fastqr12_in
!
! This subroutine performs h steps of the structured QR algorithm, using
! 2 shifts that are given as input, and returns a shift vector of size k.
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
! RH   COMPLEX(8), DIMENSION(H). Vector that contains the input shifts.
!
! H   INTEGER. Number of input shifts.
!
! K   INTEGER. Number of output shifts, K>=H.
!
! RHRH   COMPLEX(8), DIMENSION(k). Vector that contains the output shifts.

subroutine fastqr12_in(nn,d,beta,u,v,h,k, RHRH)
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
double precision :: eps = 2.22e-16

n=nn

if (n>3/2*k) then
! Perform h steps of the structured QR algorithm.
	do p=1,h
		rho=RHRH(p)
		call fastqr_ss_in(n,d,beta,u,v,rho)
	end do
	! Try to do some deflation in the final part of the matrix.
	do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
                beta(n-1)=0
                 n=n-1
        end do
        ! Perform a step of the aggressive early deflation and compute the new 
	! shift vector.
	call aggressive_deflation_in(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RHRH)
else
	! If the size of the matrix is small, compute the egenvalues performing
	! a structured QR algorithm without aggressive early deflation.
	call fastqr6(n,d,beta,u,v)
	RHRH=0
end if
end subroutine fastqr12_in

!-------------------------------------------------------------

recursive subroutine aggressive_deflation_in (n,d,beta,u,v,w,RHO)
implicit none
integer, intent(in)  :: n,w
complex(8), dimension(n), intent(inout) :: d, u,v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(w*3/2) :: h
complex(8), dimension (w), intent(out) :: RHO
complex(8), dimension(:), allocatable :: hatd
complex(8) :: S, C 
integer :: f,i,j,K, p
complex(8) :: l
complex(8) :: z
complex, dimension(2,2) :: G
double precision :: eps = 2.22e-16

K=w*3/2
if (n>K) then
    h=0
    h(1)=beta(n-K)
    ! Compute the structured Schur form of the kxk trailing principal submatrix, updating the 
    ! vector h.
  
    call fastqr11(K, h(1:K), d(n-K+1:n),beta(n-K+1:n-1),u(n-K+1:n),v(n-K+1:n))
 
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
allocate(hatd(i))
! Store the non deflated eigenvalues of the kxk trailing principal submatrix.
hatd=d(1+n-K:i+n-K)
! Bring back the matrix in Hessenberg form.
   call Hessenberg_reduction(i,h(1:i),d(1+n-K:i+n-K),beta(n-K+1:n-K+i-1),u(1+n-K:i+n-K),v(1+n-K:i+n-K))
	beta(n-K)=h(1)
   if (i< w) then
    ! The stored eigenvalues are not enough to produce w shifts, hence perform
    ! a new aggressive deflation step. 
       call aggressive_deflation_in(n-K+i,d(1:n-K+i),beta(1:n-K+i-1),u(1:n-K+i),v(1:n-K+i),w,RHO)
  else
  	! Store the smallest (in magnitude) w elements of hatd as shifts.
	call sort_1(i,hatd)
	RHO=hatd(1:w)
	deallocate(hatd)
end if
else
! If the size of the matrix is small, compute the egenvalues performing
! a structured QR algorithm without aggressive early deflation.
call fastqr6(n,d,beta,u,v)
RHO=0
end if



end subroutine aggressive_deflation_in

!--------------------------------------------------------

!SUBROUTINE Hessenberg_reduction
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



subroutine Hessenberg_reduction(n,h,d,beta,u,v) 
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
	do j=i,n-3

    		gamm=conjg(beta(j+1)-v(j+1)*u(j+2))+u(j+1)*v(j+2)
		z=beta(j)

	    	call zrotg(z,R(3,1),C,S)
	
    		beta(j)=z
		R(1,1)=d(j+1)
		R(2,1)=beta(j+1)
		R(1,2)=gamm
		R(2,2)=d(j+2)
    		R(3,1)=0
		R(3,2)=beta(j+2)
	
	
		call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
		call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))
	
	
    		d(j+1)=R(1,1)
    		beta(j+1)=R(2,1)
    		d(j+2)=R(2,2)
   		beta(j+2)=R(3,2)
	
		call zrot(1, u(j+1), 1, u(j+2), 1, C, S)
		call zrot(1, v(j+1),1,v(j+2), 1, C, conjg(S))
	end do
	
    	gamm=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,R(3,1),C,S)
	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm
	R(2,2)=d(n)


	call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
	call zrot(2, R(1,1), 1, R(1,2), 1, C, conjg(S))
	


    	d(n-1)=R(1,1)
    	beta(n-1)=R(2,1)
    	d(n)=R(2,2)

  	

	call zrot(1, u(n-1), 1, u(n), 1, C, S)
	call zrot(1, v(n-1),1,v(n), 1, C, conjg(S))	

	d(n-1)=real(d(n-1)-u(n-1)*v(n-1))+(u(n-1)*v(n-1))
    	d(n)=real(d(n)-u(n)*v(n))+(u(n)*v(n))
    	
end do
end if
end subroutine Hessenberg_reduction
!---------------------------------------------------------
!SUBROUTINE fastqr11
!
! This subroutine performs a single shift structured QR algorithm, without 
! aggressive early deflation, to compute the eigenvalue of a matrix which is
! the sum of a hermitian and a rank one matrix. Moreover this subroutine
! computes the vector Q*H, where Q* is the unitary matrix of the Schur 
! decomposition, and H is a vector given in input. 
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix. 
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.

recursive subroutine fastqr11(n,h,d,beta,u,v) 
implicit none
integer, intent(in)  :: n
integer :: imin, imax, p,i
complex(8), dimension(n), intent(inout) :: d, u, v, h
complex(8), dimension(n-1), intent(inout) :: beta
double precision :: eps = 2.22e-16



imax=n
imin=1 
do while (imax-imin .gt. 0)
	! Compute a step of the QR algorithm updating h.
	
	call fastqr11_in(imax-imin+1, h(imin:imax),d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
	
	!Try to deflate some eigenvalue
	do while (imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
		beta(imin)=0
		imin = imin + 1
	end do
	do while (imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
		beta(imax-1)=0
		imax = imax - 1
	end do
	do i=imin+1,imax-1
if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
beta(i)=0
! If a deflation occurs in the middle of the matrix, 
! compute the eigenvalues of the smallest diagonal block, 
! using a recursive structured QR algorithm. 
if (i.le. (imax-imin)/2) then
call fastqr11(i-imin+1,h(imin:i), d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
do while (imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
beta(imin)=0
imin = imin + 1
end do
else
call fastqr11(imax-i,h(i+1:imax), d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
do while (imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
beta(imax-1)=0
imax = imax - 1
end do
end if
end if
end do

end do


end subroutine fastqr11


!--------------------------------------------------------------
!SUBROUTINE fastqr11_in
!
! This subroutine performs a step of the single shift structured QR algorithm, 
! to compute the eigenvalue of a matrix which is the sum of a hermitian and a 
! rank one matrix. Moreover this subroutine computes the vector Q*H, where Q 
! is the unitary matrix obtaied from the QR decomposition of the input matrix, 
! and H is a vector given in input. 
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input matrix. 
!
! H    COMPLEX(8), DIMENSION(N). Arbitrary input vector.
!
! D    COMPLEX(8), DIMENSION(N). Vector that contains the diagonal entries
!      of the matrix.
!
! BETA COMPLEX(8), DIMENSION(N-1). Vector that contains the subdiagonal 
!      entries of the matrix.
!
! U,V  COMPLEX(8), DIMENSION(N). Vectors such that the rank one part of 
!      the matrix is UV*.

subroutine fastqr11_in(n,h,d,beta,u,v)

implicit none
integer, intent(in)  :: n
complex(8), dimension(n), intent(inout) :: d, u,v,h
complex(8), dimension(n-1), intent(inout) :: beta
complex(8) :: gamm
complex(8), dimension(2) :: l
complex(8) :: rho
complex(8), dimension(3,2) :: R
complex(8), dimension(2,2) :: A
integer :: i ,p
complex(8) :: S, C
complex(8):: z
double precision :: eps = 2.22e-16



if (n>2) then
	! Compute the Wilkinson shift.
	gamm=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
	rho=sqrt((d(n-1)+d(n))**2-4*(d(n-1)*d(n)-(beta(n-1)*(conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))))
	l(1)=(d(n-1)+d(n)+rho)/2
	l(2)=(d(n-1)+d(n)-rho)/2
	if (abs(l(1)-d(n))<abs(l(2)-d(n))) then
	    rho=l(1);
	else
    		rho=l(2);
	endif
	! Perform the structured QR step.
	z=d(1)-rho
	call zrotg(z,beta(1),C,S)

	R(1,1)=d(1)
	R(2,1)=beta(1)
	R(1,2)=gamm
	R(2,2)=d(2)
	R(3,1)=0
	R(3,2)=beta(2)
	
	call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
	call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))

	d(1)=R(1,1)
	beta(1)=R(2,1)
	d(2)=R(2,2)
	beta(2)=R(3,2)

	call zrot(1, u(1), 1, u(2), 1, C, S)
	call zrot(1, h(1), 1, h(2), 1, C, S)
	call zrot(1, v(1),1,v(2), 1, C, conjg(S))	

	d(1)=real(d(1)-u(1)*v(1))+(u(1)*v(1))
	d(2)=real(d(2)-u(2)*v(2))+(u(2)*v(2))
	
	
	do i=1,n-3
    		gamm=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
	z=beta(i)
    	call zrotg(z,R(3,1),C,S)

    	beta(i)=z
	R(1,1)=d(i+1)
	R(2,1)=beta(i+1)
	R(1,2)=gamm
	R(2,2)=d(i+2)
    	R(3,1)=0
	R(3,2)=beta(i+2)

	call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
	call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))

    	d(i+1)=R(1,1)
    	beta(i+1)=R(2,1)
    	d(i+2)=R(2,2)
   	beta(i+2)=R(3,2)
	
	call zrot(1, u(i+1), 1, u(i+2), 1, C, S)
	call zrot(1, h(i+1), 1, h(i+2), 1, C, S)
	call zrot(1, v(i+1),1,v(i+2), 1, C, conjg(S))

	enddo
    	gamm=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,R(3,1),C,S)

	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm
	R(2,2)=d(n)
	call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
	call zrot(2, R(1,1), 1, R(1,2), 1, C, conjg(S))
	

    	d(n-1)=R(1,1)
    	beta(n-1)=R(2,1)
    	d(n)=R(2,2)
  	

	call zrot(1, u(n-1), 1, u(n), 1, C, S)
	call zrot(1, h(n-1), 1, h(n), 1, C, S)
	call zrot(1, v(n-1),1,v(n), 1, C, conjg(S))	

	d(n-1)=real(d(n-1)-u(n-1)*v(n-1))+(u(n-1)*v(n-1))
    	d(n)=real(d(n)-u(n)*v(n))+(u(n)*v(n))

else
    if (n==2) then 
	gamm=conjg(beta(1)-u(2)*v(1))+u(1)*v(2)
	
	A(1,1)=d(1)
	A(2,1)=beta(1)
	A(1,2)=gamm
	A(2,2)=d(2)
if (A(2,1).NE.0) then	
	rho=sqrt((A(1,1)+A(2,2))**2-4*(A(1,1)*A(2,2)-A(2,1)*A(1,2)))
	l(1)=(A(1,1)+A(2,2)+rho)/2
	l(2)=(A(1,1)+A(2,2)-rho)/2
	if (abs(l(1)-d(2))<abs(l(2)-d(2))) then
	   d(1)=l(1);
	else
    		d(1)=l(2);
	endif
	z=A(1,1)-d(1)
	call zrotg(z,beta(1),C,S)
	call zrot(2, A(1,1), 2, A(2,1), 2 , C, S)
	call zrot(2, A(1,1), 1, A(1,2), 1, C, conjg(S))
	call zrot(1, u(1), 1, u(2), 1, C, S)
	call zrot(1, h(1), 1, h(2), 1, C, S)
	call zrot(1, v(1),1,v(2), 1, C, conjg(S))
	d(1)=A(1,1)
	d(2)=A(2,2)
        beta(1)=0
	end if
    endif
endif
end subroutine fastqr11_in
!------------------------------------------------------------------------------
!SUBROUTINE sort_1
! This subroutine sorts an input vector in decreasing order with respect
! to the magnitute of its entries.  
!
! INPUT PARAMETERS
!
! N    INTEGER. Size of the input vector. 
!
! A    COMPLEX(8), DIMENSION(N). Arbitrary input vector.

subroutine sort_1(n,a) 
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
end subroutine sort_1

!--------------------------------------------------------------------

!SUBROUTINE create_bulge
! This subroutine creates the first bulge.
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

subroutine create_bulge(d,beta,u,v,rho,bulge)
implicit none
complex(8), dimension(2), intent(inout) :: d,beta,u,v
complex(8), intent(in):: rho
complex(8), intent(out):: bulge
complex(8):: z, gamm, S, C
complex(8), dimension(3,2) :: R

gamm=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
z=d(1)-rho
call zrotg(z,beta(1),C,S)

R(1,1)=d(1)
R(2,1)=beta(1)
R(1,2)=gamm
R(2,2)=d(2)
R(3,1)=0
R(3,2)=beta(2)
	
call zrot(2, R(1,1), 3, R(2,1), 3, C, S)
call zrot(3, R(1,1), 1, R(1,2), 1, C, conjg(S))

d(1)=R(1,1)
beta(1)=R(2,1)
d(2)=R(2,2)
beta(2)=R(3,2)

call zrot(1, u(1), 1, u(2), 1, C, S)
call zrot(1, v(1),1,v(2), 1, C, conjg(S))

d(1)=real(d(1)-u(1)*v(1))+(u(1)*v(1))
d(2)=real(d(2)-u(2)*v(2))+(u(2)*v(2))

bulge=R(3,1)

end subroutine create_bulge

!-------------------------------------------------------------------

!SUBROUTINE delete_bulge
! This subroutine deletes a bulge that is in the last row of the matrix. 
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

subroutine delete_bulge (d,beta,u,v,bulge)
implicit none
complex(8), dimension(2), intent(inout) :: d,beta,u,v
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

	d(1)=real(d(1)-u(1)*v(1))+(u(1)*v(1))
    	d(2)=real(d(2)-u(2)*v(2))+(u(2)*v(2))

end subroutine delete_bulge

!-----------------------------------------------------------------------------------------------------------------------------------------

!PARALLELIZATION

!-----------------------------------------------------------------------------------------------------------------------------------------

!SUBROUTINE fastqr_ss_in_par
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
! NP   INTEGER. Number of processors.

subroutine fastqr_ss_in_par(n,d,beta,u,v,rh,k,np)
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
double precision :: eps = 2.22e-16
complex(8), dimension(np):: bulge 
real :: start, finish

siz=(n-3)/np

if (n>k*3/2) then

	do p=1,np-1
		rho=RH(p)
		
		call create_bulge(d(1:2),beta(1:2),u(1:2),v(1:2),rho,bulge(1))
					
		!$OMP PARALLEL DO
		do q=1,p 
		
			qq(q)=(q-1)*(siz)+1
			call chasing(siz-1, d(qq(q)+1:qq(q)+siz-1), beta(qq(q):qq(q)+siz-1), u(qq(q)+1:qq(q)+siz-1), &
			v(qq(q)+1:qq(q)+siz-1),bulge(q))
		end do
               	!$OMP END PARALLEL DO
               	
		do q=p,1,-1
			call chasing(3, d(qq(q)+siz-1:qq(q)+siz+1), beta(qq(q)+siz-2:qq(q)+siz+1), u(qq(q)+siz-1:qq(q)+siz+1), &
			v(qq(q)+siz-1:qq(q)+siz+1),bulge(q))
			bulge(q+1)=bulge(q)
		end do  
	end do
	
	do p=np,k
		rho=RH(p)
		call create_bulge(d(1:2),beta(1:2),u(1:2),v(1:2),rho,bulge(1))
		
		!$OMP PARALLEL DO
		do q=1,np 
		
			qq(q)=(q-1)*(siz)+1
			call chasing(siz-1, d(qq(q)+1:qq(q)+siz-1), beta(qq(q):qq(q)+siz-1), u(qq(q)+1:qq(q)+siz-1), &
			v(qq(q)+1:qq(q)+siz-1),bulge(q))
		end do
               	!$OMP END PARALLEL DO
               	   	
               	do q=np,1,-1
			call chasing(3, d(qq(q)+siz-1:qq(q)+siz+1), beta(qq(q)+siz-2:qq(q)+siz+1), u(qq(q)+siz-1:qq(q)+siz+1), &
			v(qq(q)+siz-1:qq(q)+siz+1),bulge(q))
		end do 
	
		call chasing(n-(qq(np)+siz+1),d(qq(np)+siz+1:n-1),beta(qq(np)+siz:n-1),u(qq(np)+siz+1:n-1),v(qq(np)+siz+1:n-1),bulge(np))
	
		call delete_bulge(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),bulge(np))
	
		do q=np-1,1,-1
			bulge(q+1)=bulge(q)
		end do
	end do

	do p=2,np
		!$OMP PARALLEL DO
		do q=p,np 
		
			qq(q)=(q-1)*(siz)+1
			call chasing(siz-1, d(qq(q)+1:qq(q)+siz-1), beta(qq(q):qq(q)+siz-1), u(qq(q)+1:qq(q)+siz-1), &
			v(qq(q)+1:qq(q)+siz-1),bulge(q))
		end do
               	!$OMP END PARALLEL DO
               	   	
               	do q=np,p,-1
               	
			call chasing(3, d(qq(q)+siz-1:qq(q)+siz+1), beta(qq(q)+siz-2:qq(q)+siz+1), u(qq(q)+siz-1:qq(q)+siz+1), &
			v(qq(q)+siz-1:qq(q)+siz+1),bulge(q))
		end do 
	
		call chasing(n-(qq(np)+siz+1),d(qq(np)+siz+1:n-1),beta(qq(np)+siz:n-1),u(qq(np)+siz+1:n-1),&
		v(qq(np)+siz+1:n-1),bulge(np))
		
		call delete_bulge(d(n-1:n),beta(n-2:n-1),u(n-1:n),v(n-1:n),bulge(np))
	
		do q=np-1,p,-1
			bulge(q+1)=bulge(q)
		end do
	end do
else
	call fastqr_ss_in(n,d,beta,u,v,rho)
end if
end subroutine fastqr_ss_in_par
!---------------------------------------------------------

!SUBROUTINE fastqr12_in_par
!
! This subroutine performs h steps of the structured QR algorithm, using
! 2 shifts that are given as input, and returns a shift vector of size k.
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
! RH   COMPLEX(8), DIMENSION(H). Vector that contains the input shifts.
!
! H   INTEGER. Number of input shifts.
!
! K   INTEGER. Number of output shifts, K>=H.
!
! RHRH   COMPLEX(8), DIMENSION(k). Vector that contains the output shifts.

subroutine fastqr12_in_par(nn,d,beta,u,v,h,k, RHRH,np)
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
double precision :: eps = 2.22e-16

n=nn

if (n>3/2*k) then
! Perform h steps of the structured QR algorithm.
	call fastqr_ss_in_par(n,d,beta,u,v,RHRH,h,np)

	! Try to do some deflation in the final part of the matrix.
	do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
                beta(n-1)=0
                 n=n-1
        end do
        ! Perform a step of the aggressive early deflation and compute the new 
	! shift vector.
	
	call aggressive_deflation_in(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RHRH)
	
else
	! If the size of the matrix is small, compute the egenvalues performing
	! a structured QR algorithm without aggressive early deflation.
	call fastqr6(n,d,beta,u,v)
	RHRH=0
end if
end subroutine fastqr12_in_par

!-------------------------------------------------------------------------------------------

!SUBROUTINE aggressive_deflation_par
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
!
! NP   INTEGER. Number of cores used for parallelization.

subroutine aggressive_deflation_par(n,d,beta,u,v,k,np)
implicit none
integer, intent(in)  :: n,k,np
integer :: imin, imax ,its, cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(k) :: rho
double precision :: eps = 2.22e-16
real(8):: z
real:: finish, start

imax=n
imin=1 
cont=0
rho=0

! Compute the first shift vector, of size k, using the Aggressive early deflation.
call aggressive_deflation_in(n,d,beta,u,v,k,rho)

! Try to do some deflation.
do while ( imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
	beta(imin)=0
	imin = imin + 1
	cont=0
end do
do while (imin .le. imax .and. beta(imax-1)==0 )
	imax = imax - 1
	cont=0
end do

its=1
do while (imax-imin .ge. 350)
	its=its+1
	! Try to do some middle deflation.
	do i=imin+1,imax-1
		if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
			beta(i)=0
				! If a deflation occurs in the middle of the matrix, 
				! compute the eigenvalues of the smallest diagonal block, 
				! using a structured QR algorithm without aggressive
				! early deflation. 
				if (i.le. (imax-imin)/2) then
				call fastqr6(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
					do while ( imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
					beta(imin)=0
					imin = imin + 1
					cont=0
				end do
			else
				call fastqr6(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
					do while (imin.le.imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
					beta(imax-1)=0
					imax = imax - 1
					cont=0
				end do
			end if
		end if
	end do
		! Perform k steps of the structured QR algorithm using
                ! k shifts that are given as input, and return a shift vector of size k.
         !       call cpu_time(start)
        call  fastqr12_in_par(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho,np) 
        !call cpu_time(finish)
        !print'("time4=",f6.3," seconds")', finish-start
		! Try to do some deflation.
	do while (imin.le.imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
		beta(imin)=0
		imin = imin + 1
		cont=0
	end do
	do while (imin.le.imax .and. beta(imax-1)==0)
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
                call fastqr12_in_par(imax-imin+1,d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),k,k, rho,np) 
		cont=0
	end if
end do
! When the size of the matrix becames small, perform a structured QR algorithm
! without aggressive early deflation.
call fastqr6(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
end subroutine aggressive_deflation_par

