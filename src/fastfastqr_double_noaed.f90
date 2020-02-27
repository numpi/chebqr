subroutine cqr_fastfastqr_ds(n,d,beta,u,v,k)
implicit none
integer, intent(in)  :: n,k
real(8), dimension(n), intent(inout) :: d, u, v
real(8), dimension(n-1), intent(inout) :: beta


!if(n.lt.350)then
call cqr_fastqr6_ds(n,d,beta,u,v)
!else
!call aggressive_deflation(n,d,beta,u,v,k)
!end if
end subroutine cqr_fastfastqr_ds

!----------------------------------------------------
recursive subroutine cqr_fastqr6_ds(n,d,beta,u,v)
implicit none
integer, intent(in)  :: n
integer :: imin, imax ,cont,i
real(8), dimension(n), intent(inout) :: d, u, v
real(8), dimension(n-1), intent(inout) :: beta
real(8) :: eps = 2.22e-16

imax=n
imin=1 
cont=0
do while ((abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax).or. &
abs(beta(imin+1))<eps*(abs(d(imin+1))+abs(d(imin+2))).and. imin+1.le.imax)

if (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1)))) then
beta(imin)=0
imin = imin + 1
else
beta(imin+1)=0
call cqr_fastqr6_ds_in(2, d(imin:imin+1), beta(imin:imin), u(imin:imin+1), v(imin:imin+1))
imin = imin + 2
end if
end do

do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax.or. &
abs(beta(imax-2))<eps*(abs(d(imax-2))+abs(d(imax-1))).and.imin.le.imax-1)


if (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax)))) then
beta(imax-1)=0
imax = imax - 1
else
beta(imax-2)=0
call cqr_fastqr6_ds_in(2, d(imax-1:imax), beta(imax-1:imax-1), u(imax-1:imax), v(imax-1:imax))
imax = imax - 2
end if
end do

do while (imax-imin .gt. 1)

call cqr_fastqr6_ds_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
print*, imax-imin+1

do while ((abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax).or. &
abs(beta(imin+1))<eps*(abs(d(imin+1))+abs(d(imin+2))).and. imin+1.le.imax)


if (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1)))) then
beta(imin)=0
imin = imin + 1
cont=0
else
beta(imin+1)=0
call cqr_fastqr6_ds_in(2, d(imin:imin+1), beta(imin:imin), u(imin:imin+1), v(imin:imin+1))
imin = imin + 2
cont=0
end if
end do
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax.or. &
abs(beta(imax-2))<eps*(abs(d(imax-2))+abs(d(imax-1))).and.imin.le.imax-1)

if (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax)))) then
beta(imax-1)=0
imax = imax - 1
cont=0
else
beta(imax-2)=0
call cqr_fastqr6_ds_in(2, d(imax-1:imax), beta(imax-1:imax-1), u(imax-1:imax), v(imax-1:imax))
imax = imax - 2
cont=0
end if
end do

do i=imin+1,imax-1
if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
beta(i)=0
if (i.le. (imax-imin)/2) then
call cqr_fastqr6_ds(i-imin+1, d(imin:i), beta(imin:i-1), u(imin:i), v(imin:i))
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
cont=0
end do
else
call cqr_fastqr6_ds(imax-i, d(i+1:imax), beta(i+1:imax-1), u(i+1:imax), v(i+1:imax))
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
beta(imax-1)=0
imax = imax - 1
cont=0
end do
end if
end if
end do




cont=cont+1
!if (cont==10) then
! call fastqr7_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
!cont=0
!print*, 'cont=',cont
!end if

end do
call cqr_fastqr6_ds_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))

end subroutine cqr_fastqr6_ds
!------------------------------------
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
real(8), dimension(2,2) :: A
integer :: i
real(8) :: S, C, ss, cc, r1,i1,r2,i2
real(8):: z,zz,z2,zz2
real(8):: eps = 2.22e-16


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
	call drot(1, u(1), 1, u(2), 1, C, S) !forse qui va -S
	call drot(1, v(1),1,v(2), 1, C, (S))	!e anche qui
	
endif
end subroutine cqr_fastqr6_ds_in
