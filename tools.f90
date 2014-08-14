MODULE tools
implicit none

CONTAINS

subroutine aliasy_roczne (fs,n,freq,r,s)	!wszystkie czestosci w odleglosci [freq,freq+almax] sa zaznaczone jako niepewne (2.0) 
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n,i
	real(kind=8) :: freq
	real(kind=8),parameter :: almin=0.00274, almax=0.00278	
	real :: r !rozdzielczosc
	real :: s !liczba wpisywana zamiast 1.0
	logical :: czy_al
	
	
	do i=1,n
		if ((abs(fs(1,i)-freq).gt.0).and.(abs(fs(1,i)-freq).lt.(almax+r)).and.fs(3,i).eq.1.0) then
			fs(3,i)=s
		end if
	end do

end subroutine

subroutine aliasy_dobowe(fs,n,freq,r)

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n,i,j
	real(kind=8) :: freq
	real,parameter :: d=1.0027
	real :: r

	do j=1,9
		do i=1,n
			if((fs(1,i).ge.(freq+(j-1)*d-r)).and.(fs(1,i).le.(freq+(j-1)*d+r))) then
				fs(3,i)=0.0
			end if
		end do
	end do

end subroutine

subroutine maximum (fs,n,freqmax,ampmax,fmin,fmax,k)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	real(kind=8) :: freqmax,ampmax
	integer :: n,k,i
	real(kind=8)::fmin,fmax
	
	freqmax=-1.0
	ampmax=-1.0
	do i=1,n
		if ((fs(1,i).le.fmax).and.(fs(1,i).ge.fmin).and.(fs(3,i).ne.0.0)) then
			if (fs(2,i).gt.ampmax) then
				ampmax=fs(2,i)
				freqmax=fs(1,i)
				k=i
			end if
		end if
	end do
	
!	write (*,*) 'Znaleziono maksimum o amplitudzie: ',ampmax,' dla czestotliwosci: ',freqmax

end subroutine

subroutine sn(fs,n,freq,box,snr)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n
	real(kind=8) :: freq,amp,snr
	real :: box
	
	integer :: i,j
	real(kind=8) :: fmin_out,fmin_in,fmax_out,fmax_in
	real(kind=8) :: signal,noise
	

	
	fmin_out=max(freq-dble(box/2.0),dble(0.0))
	fmin_in=freq-dble(box/40.0)
	fmax_in=freq+dble(box/40.0)
	fmax_out=min(freq+dble(box/2.0),dble(10.0))
	
	j=0
	noise=dble(0)
	signal=dble(0)

	do i=1,n
		if(((fs(1,i).ge.fmin_out).and.(fs(1,i).lt.fmin_in)).or.((fs(1,i).gt.fmax_in).and.(fs(1,i).le.fmax_out))) then
			noise=noise+fs(2,i)
			j=j+1
		end if
		if((fs(1,i).ge.fmin_in).and.(fs(1,i).le.fmax_in)) then
			signal=max(signal,fs(2,i))
		end if
	end do
	noise=noise/dble(j)
	
	snr=signal/noise
	
	!write (*,*) freq,snr

end subroutine sn

end module tools
