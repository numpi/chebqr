
!-------------------------------------
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
!NP    INTEGER Number of used processors.


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
	call aggressive_deflation(n,d,beta,u,v,k,np)
end if
end subroutine

!------------------------------------------------
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

subroutine aggressive_deflation(n,d,beta,u,v,k,np)
implicit none
integer, intent(in)  :: n,k,np
integer :: imin, imax ,its, cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(2) :: rhorho
complex(8), dimension(k) :: rho
double precision :: eps = 2.22e-16
real(8):: z
character(len=256) :: buffer

imax=n
imin=1 
cont=0

if (imin .ge. imax) then
   return
end if

! Compute the first shift vector, of size 2, using the Wilkinson shift.
rhorho(2)=sqrt((d(n-1)+d(n))**2-4*(d(n-1)*d(n)-beta(n-1)*conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))
rhorho(1)=(d(n-1)+d(n)+rhorho(2))/2
rhorho(2)=(d(n-1)+d(n)-rhorho(2))/2
! Perform 2 seps of the structured QR algorithm using
! 2 shifts that are given as input, and return a shift vector of size k.
call fastqr12_in(n,d,beta,u,v,rhorho,k,rho)
! Try to do some deflation.
do while (imin .lt. imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
	beta(imin)=0
	imin = imin + 1
	cont=0
end do
do while (imax .gt. 1 .and. beta(imax-1)==0)
	imax = imax - 1
	cont=0
end do


its=1
do while (imax-imin .ge. 35)
   its=its+1
   ! Try to do some deflation.
   do i=1,0 ! imin+1,imax-1
      if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
         beta(i)=0
         write(buffer, *) 'deflation at', i, '\n'
         call mexPrintf(buffer)
     
         ! If a deflation occurs in the middle of the matrix, 
         ! compute the eigenvalues of the smallest diagonal block, 
         ! using a structured QR algorithm without aggressive
         ! early deflation. 
         if (i.le. (imax-imin)/2) then
            call fastqr6(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
            do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
               beta(imin)=0
               imin = imin + 1
               cont=0
            end do
         else
            call fastqr6(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
            do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
               beta(imax-1)=0
               imax = imax - 1
               cont=0
            end do
         end if
      end if
   end do
   ! Perform k steps of the structured QR algorithm using
   ! k shifts that are given as input, and return a shift vector of size k.
   ! If the size of the martix is large enough perform it in parallel. 
   if (imax-imin .ge. 5) then
      call fastqr10_in_par(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,k,np)
   else	
      call fastqr10_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,k)
   end if
   ! Try to do some deflation.
   do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
      beta(imin)=0
      imin = imin + 1
      cont=0
   end do
   do while (beta(imax-1)==0.and.imin.le.imax)
      imax = imax - 1
      cont=0
   end do
   
   cont=cont+1
   
   ! If after some QR iteration there is not deflation, compute a random shift vector.
   if (cont==10) then
      do i=1,k
         call random_number(z)
         rho(i)=z
      end do
      ! Perform k seps of the structured QR algorithm using
      ! k shifts that are given as input, and return a shift vector of size k.
      ! If the size of the martix is large enough perform it in parallel. 
      if (imax-imin .ge. 5) then
         call fastqr10_in_par(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,k,np)
      else	
         call fastqr10_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,k)
      end if
      cont=0
   end if
end do
! When the size of the matrix becames small, perform a structured QR algorithm
! without aggressive early deflation.
if (imax .gt. imin) then
   call fastqr6(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
end if

end subroutine aggressive_deflation

!---------------------------------------------------


!SUBROUTINE fastqr10_in
!
! This subroutine performs k steps of the structured QR algorithm, using
! k shifts that are given as input, and returns a shift vector of size k.
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
! RH   COMPLEX(8), DIMENSION(K). Vector that contains the shifts.
!
! K    INTEGER. Number of output shifts.

subroutine fastqr10_in(nn,d,beta,u,v, RH,k)
implicit none
integer, intent(in)  :: nn,k
integer :: n
complex(8), dimension(k), intent(inout) :: RH
complex(8), dimension(nn), intent(inout) :: d, u,v
complex(8), dimension(nn-1), intent(inout) :: beta
complex(8), dimension(nn-1) :: gamm
complex(8), dimension(2) :: l
complex(8) :: rho
complex(8), dimension(3,2) :: R
integer :: i,p,q
complex(8) :: S, C
complex(8) :: z
double precision :: eps = 2.22e-16




n=nn
if (n>k*3/2) then
! Perform k steps of the structured QR algorithm.
	do p=1,k
		rho=RH(p)
		gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
		z=d(1)-rho
		call zrotg(z,beta(1),C,S)

		R(1,1)=d(1)
		R(2,1)=beta(1)
		R(1,2)=gamm(1)
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

		do i=1,n-3


		if(abs(R(3,1))<eps*(abs(beta(i)))) then
		R(3,1)=0
		end if
	
	    	gamm(i+1)=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
		z=beta(i)
	    	call zrotg(z,R(3,1),C,S)

	    	beta(i)=z
		R(1,1)=d(i+1)
		R(2,1)=beta(i+1)
		R(1,2)=gamm(i+1)
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
		call zrot(1, v(i+1),1,v(i+2), 1, C, conjg(S))

		enddo
    	gamm(n-1)=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,R(3,1),C,S)

	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm(n-1)
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
! Try to do some deflation in the final part of the matrix.
	do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
                beta(n-1)=0
                 n=n-1
         end do
! Perform a step of the aggressive early deflation and compute the new 
! shift vector.
  	call aggressive_deflation_in(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RH)
else
! If the size of the matrix is small, compute the egenvalues performing
! a structured QR algorithm without aggressive early deflation.
	call fastqr6(n,d,beta,u,v)
	RH=0
end if
end subroutine fastqr10_in
!------------------------------------------------------
!SUBROUTINE fastqr10_in_par
!
! This subroutine performs in parallel k steps of the structured QR algorithm, using
! k shifts that are given as input, and returns a shift vector of size k.
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
! RH   COMPLEX(8), DIMENSION(K). Vector that contains the shifts.
!
! K    INTEGER. Number of output shifts.
!
! NP   INTEGER. Number of processors.

subroutine fastqr10_in_par(nn,d,beta,u,v, RH,k,np)
implicit none
integer, intent(in)  :: nn,k,np
integer :: n, siz
complex(8), dimension(k), intent(inout) :: RH
complex(8), dimension(nn), intent(inout) :: d, u,v
complex(8), dimension(nn-1), intent(inout) :: beta
complex(8), dimension(nn-1) :: gamm
complex(8), dimension(2) :: l
complex(8) :: rho
complex(8), dimension(3,2) :: R
integer :: i,p,q,ii,iii, qq(np)
complex(8) :: S, C
complex(8) :: z
double precision :: eps = 2.22e-16
complex(8), dimension(np):: bulge 
real(8), dimension(np):: t1,t2



n=nn
siz=(n-4)/np

if (n>k*3/2) then
	do ii=1,k/np
		do p=1,np-1
		rho=RH(p+np*(ii-1))
		gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
		z=d(1)-rho
		call zrotg(z,beta(1),C,S)

		R(1,1)=d(1)
		R(2,1)=beta(1)
		R(1,2)=gamm(1)
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

		bulge(p)=R(3,1)
	
		!$OMP PARALLEL DO
		do q=1,p 
			qq(q)=(p-q)*(siz-1)
			call chasing(siz-1, d(qq(q)+2:qq(q)+siz), beta(qq(q)+1:qq(q)+siz), u(qq(q)+2:qq(q)+siz), v(qq(q)+2:qq(q)+siz),bulge(q))
		end do
        	!$OMP END PARALLEL DO
		do q=1,p
			qq(p)=(p-q)*(siz-1)
			call chasing(2,d(qq(q)+siz:qq(q)+siz+1),beta(qq(q)+siz-1:qq(q)+siz+1),u(qq(q)+siz:qq(q)+siz+1), &
        		v(qq(q)+siz:qq(q)+siz+1),bulge(q))
		end do
	end do
	p=np
	rho=RH(p+np*(ii-1))
	gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
	z=d(1)-rho
	call zrotg(z,beta(1),C,S)

	R(1,1)=d(1)
	R(2,1)=beta(1)
	R(1,2)=gamm(1)
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

	bulge(p)=R(3,1)

	do iii=1,1-np+(n-2)/(siz-1)-5

        !$OMP PARALLEL DO
	do p=1,np
	qq(p) = (np-p+iii-1)*(siz-1)
	call chasing(siz-1, d(qq(p)+2:qq(p)+siz), beta(qq(p)+1:qq(p)+siz), u(qq(p)+2:qq(p)+siz), v(qq(p)+2:qq(p)+siz),bulge(p))
	end do
        !$OMP END PARALLEL DO
	do p=1,np
	qq(p)=(np-p+iii-1)*(siz-1)
	call chasing(2,d(qq(p)+siz:qq(p)+siz+1),beta(qq(p)+siz-1:qq(p)+siz+1),u(qq(p)+siz:qq(p)+siz+1),v(qq(p)+siz:qq(p)+siz+1),bulge(p))
	end do
	end do
	iii=iii-1
	do q=1,np
	i=(np-1+iii)*(siz-1)+1
	call chasing(n-i-1,d(i+1:n-1),beta(i:n-1),u(i+1:n-1),v(i+1:n-1),bulge(1))

    	gamm(n-1)=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,bulge(1),C,S)

	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm(n-1)
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

	
	!$OMP PARALLEL DO
	do p=q+1,np 
	qq(p)=(np-p+q+iii-1)*(siz-1)
	call chasing(siz-1, d(qq(p)+2:qq(p)+siz), beta(qq(p)+1:qq(p)+siz), u(qq(p)+2:qq(p)+siz), v(qq(p)+2:qq(p)+siz),bulge(p-q+1))
	end do
        !$OMP END PARALLEL DO
	do p=q+1,np
	qq(p)=(np-p+q+iii-1)*(siz-1)
	call chasing(2,d(qq(p)+siz:qq(p)+siz+1),beta(qq(p)+siz-1:qq(p)+siz+1),u(qq(p)+siz:qq(p)+siz+1), &
                     v(qq(p)+siz:qq(p)+siz+1),bulge(p-q+1))
	end do
	do p=q+1,np
	bulge(p-q)=bulge(p-q+1)
	end do
	bulge(np-q+1)=0
	end do
end do
	do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
                beta(n-1)=0
                 n=n-1
         end do
  	call aggressive_deflation_in(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RH)
else
	call fastqr6(n,d,beta,u,v)
	RH=0
end if
end subroutine fastqr10_in_par
!------------------------------------------------------
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
complex(8), dimension(n) :: gamm
complex(8), intent (inout) :: bulge
complex(8), dimension(3,2) :: R
complex(8) :: S, C, z
integer :: i


do i=1,n-1
	gamm(i)=conjg(beta(i+1)-v(i)*u(i+1))+u(i)*v(i+1)
	z=beta(i)
    	call zrotg(z,bulge,C,S)

    	beta(i)=z
	R(1,1)=d(i)
	R(2,1)=beta(i+1)
	R(1,2)=gamm(i)
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

!------------------------------------------------------


!SUBROUTINE fastqr12_in
!
! This subroutine performs 2 steps of the structured QR algorithm, using
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
! RH   COMPLEX(8), DIMENSION(2). Vector that contains the input shifts.
!
! K    INTEGER. Number of output shifts.
!
! RHRH   COMPLEX(8), DIMENSION(k). Vector that contains the output shifts.

subroutine fastqr12_in(nn,d,beta,u,v, RH,k, RHRH)
implicit none
integer, intent(in)  :: nn,k
integer :: n
complex(8), dimension(2), intent(in) :: RH
complex(8), dimension(k), intent(out) :: RHRH
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
! Perform 2 steps of the structured QR algorithm.
	do p=1,2
		rho=RH(p)
		gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
		z=d(1)-rho
		call zrotg(z,beta(1),C,S)

		R(1,1)=d(1)
		R(2,1)=beta(1)
		R(1,2)=gamm(1)
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
	
		do i=1,n-3
    			gamm(i+1)=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
			z=beta(i)
 	   		call zrotg(z,R(3,1),C,S)
	
    			beta(i)=z
			R(1,1)=d(i+1)
			R(2,1)=beta(i+1)
			R(1,2)=gamm(i+1)
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
			call zrot(1, v(i+1),1,v(i+2), 1, C, conjg(S))
	
		enddo
    		gamm(n-1)=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
		z=beta(n-2)
		call zrotg(z,R(3,1),C,S)

		beta(n-2)=z
		R(1,1)=d(n-1)
		R(2,1)=beta(n-1)
		R(1,2)=gamm(n-1)
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


!-----------------------------------------
!SUBROUTINE aggressive_deflation_in
!
! This subroutine performs a single step of the structured aggressive early 
! deflation and produce a vector of shifts that can be used in subsequent
! steps of QR algorithm.
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
! W    INTEGER. Number of output shifts.
!
! RHO   COMPLEX(8), DIMENSION(w). Vector that contains the output shifts.

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
end if
else
! If the size of the matrix is small, compute the egenvalues performing
! a structured QR algorithm without aggressive early deflation.
call fastqr6(n,d,beta,u,v)
RHO=0
end if

end subroutine aggressive_deflation_in
!--------------------------------------------------
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

subroutine fastqr11(n,h,d,beta,u,v) 
implicit none
integer, intent(in)  :: n
integer :: imin, imax, p
complex(8), dimension(n), intent(inout) :: d, u, v, h
complex(8), dimension(n-1), intent(inout) :: beta
double precision :: eps = 2.22e-16



imax=n
imin=1 
do while (imax-imin .gt. 0)
	! Compute a step of the QR algorithm updating h.
	call fastqr11_in(imax-imin+1, h(imin:imax),d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
	!Try to deflate some eigenvalue
	do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
		beta(imin)=0
		imin = imin + 1
	end do
	do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
		beta(imax-1)=0
		imax = imax - 1
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
complex(8), dimension(n-1) :: gamm
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
	gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
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
	R(1,2)=gamm(1)
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
    		gamm(i+1)=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
	z=beta(i)
    	call zrotg(z,R(3,1),C,S)

    	beta(i)=z
	R(1,1)=d(i+1)
	R(2,1)=beta(i+1)
	R(1,2)=gamm(i+1)
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
    	gamm(n-1)=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,R(3,1),C,S)

	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm(n-1)
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
	gamm(1)=conjg(beta(1)-u(2)*v(1))+u(1)*v(2)
	
	A(1,1)=d(1)
	A(2,1)=beta(1)
	A(1,2)=gamm(1)
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

subroutine fastqr6(n,d,beta,u,v)
implicit none
integer, intent(in)  :: n
integer :: imin, imax ,cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
double precision :: eps = 2.22e-16

imax=n
imin=1 
cont=0

!Try to deflate some eigenvalue
do while (imin .lt. imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
   beta(imin)=0
   imin = imin + 1
end do

do while (imin .lt. imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
   beta(imax-1)=0
   imax = imax - 1
end do

do while (imax - imin .gt. 0)

! Compute a step of the QR algorithm.
call fastqr6_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))


!Try to deflate some eigenvalue
do while (imin .lt. imax .and. abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))))
   beta(imin)=0
   imin = imin + 1
   cont = 0
end do

do while (imin .lt. imax .and. abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))))
   beta(imax-1)=0
   imax = imax - 1
   cont = 0
end do

do i=imin+1,imax-1
if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
beta(i)=0
! If a deflation occurs in the middle of the matrix, 
! compute the eigenvalues of the smallest diagonal block, 
! using a recursive structured QR algorithm. 
if (i.le. (imax-imin)/2) then
call fastqr8(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
cont=0
end do
else
call fastqr8(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
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
 call fastqr7_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
cont=0
end if

end do

end subroutine fastqr6
!------------------------------------
!SUBROUTINE fastqr6_in
!
! This subroutine performs a step of the single shift structured QR algorithm, 
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

subroutine fastqr6_in(n,d,beta,u,v)
implicit none
integer, intent(in)  :: n
complex(8), dimension(n), intent(inout) :: d, u,v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(n-1) :: gamm
complex(8), dimension(2) :: l
complex(8) :: rho
complex(8), dimension(3,2) :: R
complex(8), dimension(2,2) :: A
integer :: i
complex(8) :: S, C
complex(8):: z
double precision :: eps = 2.22e-16

if (n>2) then
	! Compute the Wilkinson shift.
	gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
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
	R(1,2)=gamm(1)
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



	do i=1,n-3

	!if(abs(R(3,1))<eps*(abs(beta(i)))) then
	!R(3,1)=0
	!end if

    	gamm(i+1)=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
	z=beta(i)

    	call zrotg(z,R(3,1),C,S)

    	beta(i)=z
	R(1,1)=d(i+1)
	R(2,1)=beta(i+1)
	R(1,2)=gamm(i+1)
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
	call zrot(1, v(i+1),1,v(i+2), 1, C, conjg(S))
	enddo
    	gamm(n-1)=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,R(3,1),C,S)

	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm(n-1)
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

else
    if (n==2) then 
	gamm(1)=conjg(beta(1)-u(2)*v(1))+u(1)*v(2)
	
	A(1,1)=d(1)
	A(2,1)=beta(1)
	A(1,2)=gamm(1)
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
	call zrot(1, v(1),1,v(2), 1, C, conjg(S))
	d(1)=A(1,1)
	d(2)=A(2,2)
        beta(1)=0
	end if
    endif
endif
end subroutine fastqr6_in
!-----------------------------------
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
complex(8), dimension(n), intent(inout) :: d, u,v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(n-1) :: gamm
complex(8), dimension(2) :: l
real(8) :: rho
complex(8), dimension(3,2) :: R
complex(8), dimension(2,2) :: A
integer :: i
complex(8) :: S, C
complex(8):: z
double precision :: eps = 2.22e-16

if (n>2) then
	gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
	! Compute a random shift.
call random_number(rho)

	! Perform the structured QR step.
	z=d(1)-rho
	call zrotg(z,beta(1),C,S)

	R(1,1)=d(1)
	R(2,1)=beta(1)
	R(1,2)=gamm(1)
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



	do i=1,n-3


    	gamm(i+1)=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
	z=beta(i)
    	call zrotg(z,R(3,1),C,S)

    	beta(i)=z
	R(1,1)=d(i+1)
	R(2,1)=beta(i+1)
	R(1,2)=gamm(i+1)
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
	call zrot(1, v(i+1),1,v(i+2), 1, C, conjg(S))
	enddo
    	gamm(n-1)=conjg(beta(n-1)-v(n-1)*u(n))+u(n-1)*v(n);
	z=beta(n-2)
	call zrotg(z,R(3,1),C,S)

	beta(n-2)=z
	R(1,1)=d(n-1)
	R(2,1)=beta(n-1)
	R(1,2)=gamm(n-1)
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



else
    if (n==2) then 
	gamm(1)=conjg(beta(1)-u(2)*v(1))+u(1)*v(2)
	
	A(1,1)=d(1)
	A(2,1)=beta(1)
	A(1,2)=gamm(1)
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
	call zrot(1, v(1),1,v(2), 1, C, conjg(S))
	d(1)=A(1,1)
	d(2)=A(2,2)
        beta(1)=0
	end if
    endif
endif
end subroutine fastqr7_in
!-------------------------------------------

!SUBROUTINE fastqr8
!
! This subroutine performs  a recursive single shift structured QR algorithm, without 
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


recursive subroutine fastqr8(n,d,beta,u,v)
implicit none
integer, intent(in)  :: n
integer :: imin, imax ,cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
double precision :: eps = 2.22e-16

imax=n
imin=1 
cont=0
! Try to deflate some eigenvalue
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
end do
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
beta(imax-1)=0
imax = imax - 1
end do
do while (imax-imin .gt. 0)

! Perform a step of the structured QR algorithm.
call fastqr6_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))

! Try to deflate some eigenvalue.
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
cont=0
end do
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
beta(imax-1)=0
imax = imax - 1
cont=0
end do

do i=imin+1,imax-1
if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
beta(i)=0
! If a deflation occurs in the middle of the matrix, recursively
! compute the eigenvalues of the smallest diagonal block.
if (i.le. (imax-imin)/2) then
call fastqr8(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
cont=0
end do
else
call fastqr8(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
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
 call fastqr7_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
cont=0
end if

end do

end subroutine fastqr8

