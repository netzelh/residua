program res3bl !wywolanie programu: res <nazwa_pliku_z_residuami> <okres> [przedzial czasowy danych]
implicit none

	character (len=30) :: plik
	character (len=30) :: period_ch
	character (len=500) :: sys
	character (len=3) :: flaga=" " !T lub ZO lub TZO
	character (len=200) :: bl=" "
	
	integer :: fn=10 !file number
	integer :: n !liczba wierszy
	integer :: i,k
	integer :: win=20 ! numer pliku z oknem spektralnym
	integer :: trendwin=30 !numer pliku z oknem spektralnym trendu
	
	real(kind=8), dimension (:,:), allocatable :: fs

	real(kind=8) :: freqmax,ampmax,snr
	real(kind=8) :: fmin,fmax
	real(kind=8) :: period,f
	real(kind=8) :: ft,at,snt !trend
	real(kind=8) :: fp,ap,snp !w pobliżu piku głównego
	real(kind=8) :: f2,a2,sn2 !efekt instrumentalny w 2
!	real(kind=8) :: f1,a1,sn1 !efekt instrumentalny w 1
	real(kind=8) :: fp2,ap2,snp2 !w pobliżu 2f
	real(kind=8) :: fx,ax,snx !w szukanym zakresie, szukana czestosc X
	
	real :: r ! rozdzielczosc
	
	logical :: czy_zo=.false., czy_trend=.false., czy_bl=.false., czy_mod=.false., czy_al=.false., czy_p2=.false.



call get_arg (plik,period,f,r)				!wczytywanie argumentow wywolania
call get_file(fn,fs,n) 					!wczytywanie pliku res do tablicy

call trend	(fs,n,ft,r)		!max dla malych czestotliwosci f<0.003
call pik2 (fs,n,f2,r)		!pik w okolicach 2
call pik_gl (bl)			!max w poblizu piku glownego - sprawdzanie czy efekt Blazki (odrzucenie gwiazdy do innej analizy)
call szukaj_px 				!max w szukanym zakresie

!fragment do testowania
!open (50,file='ppp.txt')
!do i=1,n
!	write (50,*) fs(1,i),fs(2,i),fs(3,i)
!end do
!close (50)
!~~~~~~~~~~~~~~~~~~
	
deallocate(fs)
close(fn)

call flagi(bl)

CONTAINS

subroutine pik2(f2,n,ft,r)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n !liczba wierszy
	real(kind=8) :: f2
	real :: r ! rozdzielczosc
	real :: d=1.0027 ! doba
	integer :: i,j
			

	call maximum(fs,n,f2,a2,dble(1.98),dble(2.02),k)
	call sn(fs,n,f2,0.2,sn2)
	if (snt.gt.4.0) then
		czy_p2=.true.
		do j=1,9
			call aliasy_roczne(fs,n,f2+(j-1)*d,r,2.0)
		end do 
		call aliasy_dobowe (fs,n,f2,r)
	end if

end subroutine

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

subroutine szukaj_px

	fmin=f/0.64
	fmax=f/0.58

	call maximum(fs,n,fx,ax,fmin,fmax,k)
	call sn(fs,n,fx,0.2,snx)
	if (snx.gt.4.0) then
		czy_mod=.true.
		if (fs(3,k).eq.2.0) czy_al=.true.
	end if

end subroutine

subroutine pik_gl(bl)
implicit none

	character (len=*) :: bl
	character (len=200) :: pom,pom2

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
					write (pom2,*) 1/(f-fp)
					pom = trim(bl)
					bl = trim(pom)//trim(pom2)
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

subroutine flagi(bl)
implicit none

	character (len=500) :: line,pom
	character (len=20) :: sx_ch
	character (len=*) :: bl
	real :: sx
	
	sx=f/fx
	
	write (sx_ch,*) sx
	
	line=plik(0:len_trim(plik)-4)//" "//trim(period_ch)//" "//trim(sx_ch)
	
	if (czy_bl) then
		pom=line
		line=trim(pom)//" BL "
	end if
			
	if (czy_trend) then
		pom=line
		line=trim(pom)//" T "
	end if

	if (czy_zo) then
		pom=line
		line=trim(pom)//" ZO "
	end if
	
	if (czy_p2) then
		pom=line
		line=trim(pom)//" P2 "
	end if
	
	if (czy_mod) then
		pom=line
		line=trim(pom)//" MOD "
	end if
	
	if (czy_al) then
		pom=line
		line=trim(pom)//" AL "
	end if

	write (*,*) 

	pom=line
	line=trim(pom)//trim(bl)

	!write (*,*) trim(line)
	!write (*,*) trim(blazko)
	
	call system('echo " '//trim(line)//' " >> rrc2.txt')

end subroutine

subroutine trend(fs,n,ft,r)
implicit none

	integer :: k
	
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

subroutine get_arg(plik,period,f,r)
implicit none

	character (len=30) :: plik
	character (len=15) :: t_ch
	integer :: na
	integer :: u
	real :: t,r
	real, parameter :: df=5
	real(kind=8) :: period,f
	
	na=iargc()
	
	select case(na)
		case (2)
			r=0.001
		case (3)
			call getarg(3,t_ch)
			read (t_ch,*) t
			r=2.0/t
	end select
	
	call getarg(1,plik)
	call getarg(2,period_ch)
	read(period_ch,*) period	
	f=dble(1.0)/period

end subroutine

subroutine get_file(fn,fs,n)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: fn !file number
	integer :: n !liczba wierszy

	open(fn,file=plik)
	write (*,*) 'Plik: ',plik
	n=ile_wierszy(fn)
	allocate (fs(3,n)) !1-freq,2-amp,3-zero/jeden
	call getp(fs,fn,n)

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

subroutine getp(fs,fn,n)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: fn,n
	integer :: i
	
	do i=1,n
		read(fn,*) fs(1,i),fs(2,i)
		fs(3,i)=1.0
	end do
	
!	write (*,*) 'OK, wczytano dane do programu'
	
end subroutine

integer function ile_wierszy(fn)
implicit none

	integer :: fn
	real :: a,b !zmienne pomocnicze
	
	ile_wierszy=0
	
	do
		read (fn,*,end=100) a,b
		ile_wierszy=ile_wierszy+1
	end do

	100 rewind(fn)
	
	!write (*,*) 'Wierszy do wczytania: ',ile_wierszy
		
	return

end function

subroutine maximum (fs,n,freqmax,ampmax,fmin,fmax,k)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	real(kind=8) :: freqmax,ampmax
	integer :: n,k
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



end program res3bl
