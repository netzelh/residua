MODULE pik
use tools
implicit none

CONTAINS

subroutine pik2(fs,n,f2,a2,sn2,r,czy_p2)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n !liczba wierszy
	real(kind=8) :: f2,a2,sn2
	real :: r ! rozdzielczosc
	real :: d=1.0027 ! doba
	integer :: i,j,k
	logical :: czy_p2
			

	call maximum(fs,n,f2,a2,dble(1.98),dble(2.02),k)
	call sn(fs,n,f2,0.2,sn2)
	if (sn2.gt.4.0) then
		czy_p2=.true.
		do j=1,9
			call aliasy_roczne(fs,n,f2+(j-1)*d,r,2.0)
		end do 
		call aliasy_dobowe (fs,n,f2,r)
	end if

end subroutine


subroutine szukaj_px(fs,n,fx,ax,snx,f,r,czy_mod,czy_al)
implicit none

	real (kind=8) :: fx,ax,snx,f
	real :: r
	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n,k
	logical :: czy_mod, czy_al
	
	real (kind=8) :: fmin,fmax
	

	fmin=f/0.64
	fmax=f/0.58

	call maximum(fs,n,fx,ax,fmin,fmax,k)
	call sn(fs,n,fx,0.2,snx)
	if (snx.gt.4.0) then
		czy_mod=.true.
		if (fs(3,k).eq.2.0) czy_al=.true.
	end if

end subroutine

subroutine pik_gl(fs,n,fp,ap,snp,f,r,czy_zo,czy_bl)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	real (kind=8) :: fp,ap,snp,f,fmin,fmax
	real :: r
	integer :: n,k
	logical :: czy_zo, czy_bl

	fmin=f-0.2
	fmax=f+0.2

	do
		call maximum(fs,n,fp,ap,fmin,fmax,k)
		call sn(fs,n,fp,0.2,snp)
		if (snp.gt.4.0) then
			if (abs(fp-f).lt.r) then
			!write (*,*) fp,snp,'Zmiana okresu',1/(f-fp)
				czy_zo=.true.
				call aliasy_dobowe	(fs,n,fp,r)
				call aliasy_roczne (fs,n,fp,r,3.0)			
			else
			!	write(*,*) fp,snp,'Blazko? Osobna analiza',1/(f-fp)
				!write (pom2,*) 1/(f-fp)
				!pom = trim(bl)
				!bl = trim(pom)//trim(pom2)
				if (fs(3,k).ne.3.0) then
					czy_bl=.true.
					!write (pom2,*) 1/(f-fp)
					!pom = trim(bl)
					!bl = trim(pom)//trim(pom2)
				end if
				call aliasy_roczne (fs,n,fp,r,3.0) 
				call aliasy_dobowe (fs,n,fp,r)
				fmax=max(f+2*abs(f-fp)+r,f+0.2)
				fmin=min(f-2*abs(f-fp)-r,f-0.2)
			end if
		else
			exit
		end if
		
	end do

end subroutine

subroutine trend(fs,n,ft,at,snt,r,czy_trend)
implicit none

	integer :: k
	real (kind=8) :: at,snt
	logical :: czy_trend
	
	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n !liczba wierszy
	real(kind=8) :: ft
	real :: r ! rozdzielczosc
	real :: d=1.0027 ! doba
	integer :: i,j
			

	call maximum(fs,n,ft,at,dble(0.0),dble(0.003),k)
	call sn(fs,n,ft,0.2,snt)
	if (snt.gt.4.0) then
		czy_trend=.true.
		do j=1,9
			call aliasy_roczne(fs,n,ft+(j-1)*d,r,2.0)
		end do 
		call aliasy_dobowe (fs,n,ft,r)
	end if

end subroutine


end module pik
