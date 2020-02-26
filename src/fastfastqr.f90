


!-------------------------------------
subroutine fastfastqr(n,d,beta,u,v,k)
implicit none
integer, intent(in)  :: n,k
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta


if(n.lt.350)then
call fastqr6(n,d,beta,u,v)
else
call aggressive_deflation(n,d,beta,u,v,k)
end if
end subroutine

!--------------------------------------------------

subroutine aggressive_deflation(n,d,beta,u,v,k)
implicit none
integer, intent(in)  :: n,k
integer :: imin, imax ,its, cont,i
complex(8), dimension(n), intent(inout) :: d, u, v
complex(8), dimension(n-1), intent(inout) :: beta
complex(8), dimension(2) :: rhorho
complex(8), dimension(k) :: rho
double precision :: eps = 2.22e-16
real(8):: z

imax=n
imin=1 
cont=0
rhorho(2)=sqrt((d(n-1)+d(n))**2-4*(d(n-1)*d(n)-beta(n-1)*conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))
rhorho(1)=(d(n-1)+d(n)+rhorho(2))/2
rhorho(2)=(d(n-1)+d(n)-rhorho(2))/2
call fastqr12_in(n,d,beta,u,v,rhorho,k,rho)
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
cont=0
end do
do while (beta(imax-1)==0 .and. imin .le. imax)
imax = imax - 1
cont=0
end do


its=1
do while (imax-imin .ge. 350)
its=its+1
do i=imin+1,imax-1
if (abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
beta(i)=0
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
!print*, imax-imin
call fastqr10_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,k)
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
cont=0
end do
do while (beta(imax-1)==0.and.imin.le.imax)
imax = imax - 1
cont=0
!print*, imax-imin
end do

cont=cont+1
if (cont==10) then
do i=1,k
call random_number(z)
rho(i)=z
end do
call fastqr10_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax),rho,k)
cont=0
!print*, 'cont=',cont
!print*, 'rho=',rho
end if

end do
call fastqr6(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
!print*, its
end subroutine aggressive_deflation

!----------------------------------------------------
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
	do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
                beta(n-1)=0
                 n=n-1
         end do
  	call aggressive_deflation_in(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RH)
else
	call fastqr6(n,d,beta,u,v)
	RH=0
end if
end subroutine fastqr10_in
!------------------------------------------------------

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
	do while (n>1 .AND. abs(beta(n-1))<eps*(abs(d(n-1))+abs(d(n))))
                beta(n-1)=0
                 n=n-1
         end do
	call aggressive_deflation_in(n,d(1:n),beta(1:n-1),u(1:n),v(1:n),k,RHRH)
else
	call fastqr6(n,d,beta,u,v)
	RHRH=0
end if
end subroutine fastqr12_in

!-----------------------------------------
recursive subroutine aggressive_deflation_in (n,d,beta,u,v,w,RHO)
implicit none
integer, intent(in)  :: n,w
complex(8), dimension(n), intent(inout) :: d, u,v
!complex(8), dimension(1,n), intent(out) :: v
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
    call fastqr11(K, h(1:K), d(n-K+1:n),beta(n-K+1:n-1),u(n-K+1:n),v(n-K+1:n))
    i=K
    j=0
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
		!G(1,1)=C
		!G(2,1)=-conjg(S)
		!G(1,2)=S
		!G(2,2)=C
		l=u(f+n-K)
		u(f+n-K)=u(f+n-K+1)
		u(f+n-K+1)=l
		!u(f+n-K:f+1+n-K)=matmul(transpose(conjg(G)),u(f+n-K:f+1+n-K))
 		call zrot(1, u(f+n-K), 1, u(f+1+n-K), 1, C,S)
		l=h(f)
		h(f)=h(f+1)
		h(f+1)=l
		call zrot(1, h(f), 1, h(f+1), 1, C, S)
		!h(f:f+1)=matmul(transpose(conjg(G)),h(f:f+1))
		l=v(f+n-K)
		v(f+n-K)=v(f+n-K+1)
		v(f+n-K+1)=l
		!v(f+n-K:f+1+n-K)=matmul(v(f+n-K:f+1+n-K),G)
		call zrot(1, v(f+n-K),1,v(f+1+n-K), 1, C, conjg(S))
		!call zrot(1, d(f+n-K), 1, d(f+1+n-K), 1, 0, 1)
		l=d(f+n-K)
		d(f+n-K)=d(f+n-K+1)
		d(f+n-K+1)=l
            !u(f+n-K:f+1+n-K)=matmul(transpose(conjg(G)),u([f+1+n-K,f+n-K]))
            !h(f:f+1)=matmul(transpose(conjg(G)),h([f+1,f]))
           ! v(f+n-K:f+1+n-K)=matmul(v([f+1+n-K,f+n-K]),G)
           ! d([f+n-K,f+1+n-K])=d([f+1+n-K,f+n-K])
           end do
       end if
   end do
allocate(hatd(i))
hatd=d(1+n-K:i+n-K)
   call Hessenberg_piena(i,h(1:i),d(1+n-K:i+n-K),beta(n-K+1:n-K+i-1),u(1+n-K:i+n-K),v(1+n-K:i+n-K))
	beta(n-K)=h(1)
   if (i< w) then
       call aggressive_deflation_in(n-K+i,d(1:n-K+i),beta(1:n-K+i-1),u(1:n-K+i),v(1:n-K+i),w,RHO)
  else
	call sort_imbecille(i,hatd)
	RHO=hatd(1:w)
end if
else
call fastqr6(n,d,beta,u,v)
RHO=0
end if



end subroutine aggressive_deflation_in
!--------------------------------------------------
!function [d,beta,u,v]=Hessenberg_piena(h,d,u,v)
!n=size(d,1);
!if n>1
!T=diag(d);
!T=T+triu(u*v,1)-(tril(u*v,-1))';
!for i=n-1:-1:1
!    [G,y]=planerot(h(i:i+1));
!    h(i:i+1)=y;
!    T(i:i+1,:)=G*T(i:i+1,:);
!    T(:,i:i+1)=T(:,i:i+1)*G';
!    u(i:i+1)=G*u(i:i+1);
!    v(i:i+1)=v(i:i+1)*G';
!end
!for i=1:n-2
!    for j=n-1:-1:i+1
!        G=planerot(T(j:j+1,i));
!        T(j:j+1,:)=G*T(j:j+1,:);
!        T(:,j:j+1)=T(:,j:j+1)*G';
!        u(j:j+1)=G*u(j:j+1);
!        v(j:j+1)=v(j:j+1)*G';
!    end
!end
!d=diag(T);
!beta=diag(T,-1);
!else
!    beta=0;
!end
!-------------------------------------------
subroutine Hessenberg_piena(n,h,d,beta,u,v) 
implicit none
integer, intent(in)  :: n
complex(8), dimension(n), intent(inout) :: d, u,v
!complex(8), dimension(1,n), intent(out) :: v
complex(8), dimension(n-1), intent(out) :: beta
complex(8), dimension(n), intent(inout) :: h
integer :: i,j,p
complex(8), dimension(n,n) :: T
complex(8) :: C,S
complex(8), dimension(2,2) :: G
complex(8) :: z

if (n>1) then
	T=0
	do i=1,n
	   T(i,i)=d(i)	
	   do j=i+1,n
	   	T(i,j)=u(i)*v(j)-conjg(u(j)*v(i))
	   end do
	end do
do i=n-1,1,-1
	call zrotg(h(i),h(i+1),C,S)
	h(i+1)=0
	call zrot(n-i+1, T(i,i), n, T(i+1,i), n, C, S)
	call zrot(n, T(1,i), 1, T(1,i+1), 1, C, conjg(S))
	call zrot(1, u(i), 1, u(i+1), 1, C, S)
	call zrot(1, v(i), 1, v(i+1), 1, C, conjg(S))	
	!G=reshape([C , -conjg(S), S, C], shape(G))
	!h([i,i+1])=matmul(G,h([i,i+1]))
	!T(i:i+1,:)=matmul(G,T(i:i+1,:))
    	!T(:,i:i+1)=matmul(T(:,i:i+1),(transpose(conjg(G))))
   	!u(i:i+1)=matmul(G,u(i:i+1))
    	!v(i:i+1)=matmul(v(i:i+1),(transpose(conjg(G))))
end do
do i=1,n-2
	do j=n-1,i+1,-1
	call zrotg(T(j,i),T(j+1,i),C,S)
	T(j+1,i)=0
	call zrot(n-i, T(j,i+1), n, T(j+1,i+1), n, C, S)
	call zrot(n, T(1,j), 1, T(1,j+1), 1, C, conjg(S))
	call zrot(1, u(j), 1, u(j+1), 1, C, S)
	call zrot(1, v(j), 1, v(j+1), 1, C, conjg(S))
	!G=reshape([C , -conjg(S), S, C], shape(G))
	!T(j:j+1,:)=matmul(G,T(j:j+1,:))
        !T(:,j:j+1)=matmul(T(:,j:j+1),(transpose(conjg(G))))
        !u(j:j+1)=matmul(G,u(j:j+1))
        !v(j:j+1)=matmul(v(j:j+1),(transpose(conjg(G))))
    end do
end do
do i=1,n-1
d(i)=T(i,i)
beta(i)=T(i+1,i)
end do
d(n)=T(n,n)
end if
end subroutine Hessenberg_piena
!---------------------------------------------------------
!function [h,d, beta, u , v]=fastQR11(d,beta,u,v,h) %fastqr6 con vettore h che si aggiorna per aggressive deflation
!n=size(d,1);
!beta(n)=0;
!while norm(beta) ~= 0
!     c=0;
!    for l=1:n
!        if beta(l)==0
!            [d(c+1:l), beta(c+1:l-1), u(c+1:l) , v(c+1:l),h(c+1:l)]=fastQR11_in(d(c+1:l), beta(c+1:l-1), u(c+1:l) , v(c+1:l),h(c+1:l));
!            c=l;
!        end
!    end
!end

!-----------------------------------------------------------------------
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
call fastqr11_in(imax-imin+1, h(imin:imax),d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
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
	gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
	rho=sqrt((d(n-1)+d(n))**2-4*(d(n-1)*d(n)-(beta(n-1)*(conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))))
	l(1)=(d(n-1)+d(n)+rho)/2
	l(2)=(d(n-1)+d(n)-rho)/2
	if (abs(l(1)-d(n))<abs(l(2)-d(n))) then
	    rho=l(1);
	else
    		rho=l(2);
	endif
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
	!A(1,1)=d(1)
	!A(2,1)=beta(1)
	!A(1,2)=conjg(beta(1)-u(2)*v(1))+u(1)*v(2)
   	!A(2,2)=d(2)
	
	!lwork=16
	!call zgees ('V','N',.TRUE.,2, A ,2,0,d(1:2), sch , 2, work, lwork, rwork,.TRUE. , info)

	!beta(1)=0
	!l=u(1:2)
	!u(1)=conjg(sch(1,1))*l(1)+conjg(sch(2,1))*l(2)
	!u(2)=conjg(sch(1,2))*l(1)+conjg(sch(2,2))*l(2)
	!l=h(1:2)
	!h(1)=conjg(sch(1,1))*l(1)+conjg(sch(2,1))*l(2)
	!h(2)=conjg(sch(1,2))*l(1)+conjg(sch(2,2))*l(2)
	!l=v(1:2)
	!v(1)=sch(1,1)*l(1)+sch(2,1)*l(2)
	!v(2)=sch(1,2)*l(1)+sch(2,2)*l(2)
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

!--------------------------------------------------
subroutine sort_imbecille(n,a) 
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
end subroutine sort_imbecille
!--------------------------------------------
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
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
end do
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
beta(imax-1)=0
imax = imax - 1
end do
do while (imax-imin .gt. 0)



call fastqr6_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))

!print*, imin, 'imin'
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
if (cont==10) then
 call fastqr7_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
cont=0
!print*, 'cont=',cont
end if

end do

end subroutine fastqr6
!------------------------------------
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

	gamm(1)=conjg(beta(1)-v(1)*u(2))+u(1)*v(2)
	rho=sqrt((d(n-1)+d(n))**2-4*(d(n-1)*d(n)-(beta(n-1)*(conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n)))))
	l(1)=(d(n-1)+d(n)+rho)/2
	l(2)=(d(n-1)+d(n)-rho)/2
	if (abs(l(1)-d(n))<abs(l(2)-d(n))) then
	    rho=l(1);
	else
    		rho=l(2);
	endif

if (n==65)then
!print*, min(abs(beta(1)),abs(beta(2)),abs(beta(3)),abs(beta(4)),abs(beta(5)),abs(beta(6)),abs(beta(7)),abs(beta(8)),abs(beta(9)),abs(beta(10)),abs(beta(11)))
!print*,beta
end if

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

	!if(abs(beta(i))<eps*(abs(d(i))+abs(d(i+1)))) then
	!beta(i)=0
	!end if
	if(abs(R(3,1))<eps*(abs(beta(i)))) then
	R(3,1)=0
	end if

    	gamm(i+1)=conjg(beta(i+1)-v(i+1)*u(i+2))+u(i+1)*v(i+2)
	z=beta(i)
!print*,'R=', R(3,1)
!print*, 'z=',z
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
!if (isnan(abs(R(1,1)))) then
!print*,'1', beta(i), z
!print*,',cs', c, s
!print*, d(i+1),gamm(i+1)
!print*, beta(i+1), d(i+2)
!stop
!end if
!print*, R(1,1),R(1,2)
!print*, R(2,1),R(2,2)
!print*, R(3,1), R(3,2)


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
!if (abs(rho) .le. 100) then
!print*, d(n)
!print*, d(n-1), (conjg(beta(n-1)-u(n)*v(n-1))+u(n-1)*v(n))
!print*, beta(n-1), d(n)
!end if



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
call random_number(rho)

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
do while (abs(beta(imin))<eps*(abs(d(imin))+abs(d(imin+1))).and. imin.le.imax)
beta(imin)=0
imin = imin + 1
end do
do while (abs(beta(imax-1))<eps*(abs(d(imax-1))+abs(d(imax))).and.imin.le.imax)
beta(imax-1)=0
imax = imax - 1
end do
do while (imax-imin .gt. 0)

!print*, imax-imin, 
call fastqr6_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
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
if (cont==10) then
 call fastqr7_in(imax-imin+1, d(imin:imax), beta(imin:imax-1), u(imin:imax), v(imin:imax))
cont=0
!print*, 'cont=',cont
end if

end do

end subroutine fastqr8

