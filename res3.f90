program res3bl !wywolanie programu: res <nazwa_pliku_z_residuami> <okres> [przedzial czasowy danych]
use io
use pik
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
call get_file(plik,fn,fs,n) 					!wczytywanie pliku res do tablicy

call trend	(fs,n,ft,at,snt,r,czy_trend)		!max dla malych czestotliwosci f<0.003
call pik2 (fs,n,f2,a2,sn2,r,czy_p2)		!pik w okolicach 2
call pik_gl (fs,n,fp,ap,snp,f,r,czy_zo,czy_bl)			!max w poblizu piku glownego - sprawdzanie czy efekt Blazki (odrzucenie gwiazdy do innej analizy)
call szukaj_px (fs,n,fx,ax,snx,f,r,czy_mod,czy_al) 				!max w szukanym zakresie

!fragment do testowania
!open (50,file='ppp.txt')
!do i=1,n
!	write (50,*) fs(1,i),fs(2,i),fs(3,i)
!end do
!close (50)
!~~~~~~~~~~~~~~~~~~
	
deallocate(fs)
close(fn)

call flagi(f,fx,plik,period_ch,czy_al,czy_bl,czy_mod,czy_p2,czy_trend,czy_zo)

end program res3bl
