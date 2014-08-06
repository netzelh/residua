program res3 !wywolanie programu: res <nazwa_pliku_z_residuami> <okres>
implicit none

	character (len=30) :: plik
	character (len=30) :: period_ch
	character (len=500) :: line,sys
	character (len=3) :: flaga=" " !T lub ZO lub TZO
	
	integer :: fn=10 !file name
	integer :: n !liczba wierszy
	integer :: i,k
	integer :: win=20 ! numer pliku z oknem spektralnym
	integer :: trendwin=30 !numer pliku z oknem spektralnym trendu
	integer :: u=50 !bedzie trzeba policzyc to u z rozdzielczosci transformaty
	
	real(kind=8), dimension (:,:), allocatable :: fs

	real(kind=8) :: freqmax,ampmax,snr
	real(kind=8) :: fmin,fmax
	real(kind=8) :: period,f
	real(kind=8) :: ft,at,snt !trend
	real(kind=8) :: fp,ap,snp !w pobliżu piku głównego
!	real(kind=8) :: f2,a2,sn2 !efekt instrumentalny w 2
!	real(kind=8) :: f1,a1,sn1 !efekt instrumentalny w 1
	real(kind=8) :: fp2,ap2,snp2 !w pobliżu 2f
	real(kind=8) :: fx,ax,snx !w szukanym zakresie, szukana czestosc X :)
	
	real :: e=0.001 ! odrzucenie (trzeba poprawic na rozdzielczosc transformaty, t biore z .dat)
	
	logical :: czy_zo=.false., czy_trend=.false., czy_bl=.false., czy_f2=.false.

!flagi beda ustawiane poprzez zmienne logiczne czy_*


! WCZYTYWANIE ARGUMENTOW
call getarg(1,plik)
call getarg(2,period_ch)
read(period_ch,*) period	
f=dble(1.0)/period

! WCZYTYWANIE PLIKU Z RESIDUAMI
open(fn,file=plik)
write (*,*) 'Plik: ',plik
n=ile_wierszy(fn)
allocate (fs(3,n)) !1-freq,2-amp,3-zero/jeden
call getp(fs,fn,n)

!max dla malych czestotliwosci f<0.003
call maximum(fs,n,ft,at,dble(0.0),dble(0.003),k)
call sn(fs,n,ft,at,0.2,snt)
!jesli trend istnieje to wywalanie z analizy pikow z okna spektralnego dla 0.0
!piki:
!0.0000227
!1.0027527
!2.0027346
!3.0080538
!4.0110110
!5.0110383
!6.0000049
!7.0191783
!7.9945178
if (snt.gt.4.0) then
	czy_trend=.true.
	write (*,*) 'jest trend'
	do i=1,n
		if ((fs(1,i).gt.0.0000227-e).and.(fs(1,i).lt.0.0000227+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.1.0027527-e).and.(fs(1,i).lt.1.0027527+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.2.0027346-e).and.(fs(1,i).lt.2.0027346+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.3.0080538-e).and.(fs(1,i).lt.3.0080538+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.4.0110110-e).and.(fs(1,i).lt.4.0110110+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.5.0110383-e).and.(fs(1,i).lt.5.0110383+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.6.0000049-e).and.(fs(1,i).lt.6.0000049+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.7.0191783-e).and.(fs(1,i).lt.7.0191783+e)) then
			fs(3,i)=0.0
		else if ((fs(1,i).gt.7.9945178-e).and.(fs(1,i).lt.7.9945178+e)) then
			fs(3,i)=0.0
		end if
	end do
end if

!max dla 1 (trzeba testowac czy to jest alias trendu czy jest wyzsza amp, a moze nie bo sie wyklucza nie wazne przez ktore z tych wywalimy czestoci?)

!max dla 2
!to powyzsze wyjdzie w praniu czy trzeba dodac - moze jakis sladowy trend wystarczy zeby to sie usunelo i jednak pojawia sie zawsze z malym trendem


!max w poblizu piku glownego - sprawdzanie czy efekt Blazki
!max w poblizu piku glownego
fmin=f-0.2
fmax=f+0.2

do
	call maximum(fs,n,fp,ap,fmin,fmax,k)
	call sn(fs,n,fp,ap,0.2,snp)

	if (snp.gt.4.0) then !sprawdzanie istotnosci maksimum
		if (abs(fp-f).lt.e) then !sprawdzanie czy zmiana okresu
			write (*,*) fp,snp,'Zmiana okresu',k,1/(f-fp)
			fs(3,k)=0
			do i=1,u
				if((k-i).ge.0) fs(3,k-i)=0
				if((k+i).le.n) fs(3,k+i)=0
			end do
			czy_zo=.true.
		else
			write(*,*) fp,snp,'Blazko?',1/(f-fp)
		!brakuje mi tu testu czy moge to przyjac za Blazko, ale po prostu wywalam do osobnej analizy tj musze przeskoczyc (goto?) z ustawieniem dobrej flagi
		czy_bl=.true.
		
			fs(3,k)=0
			do i=1,u
				if((k-i).ge.0) fs(3,k-i)=0
				if((k+i).le.n) fs(3,k+i)=0
			end do
			fmax=max(f+2*abs(f-fp)+e,f+0.2)
			fmin=min(f-2*abs(f-fp)-e,f-0.2)
		end if
	else
	!	if(czy_cos.eqv..false.) call system(' echo  " '//trim(plik)//' " >> test-nic.txt')
		write (*,*) 'Brak wyraznego sygnalu w poblizu glownej czestosci'
		!to znaczy ze mozna przeszukiwac obszar :)
		exit
	end if
	
end do

!i tu musi byc wyskoczenie do zapisywania gwiazdy jesli blazko! tzn do dalszej analizy jesli jest podejzana o Blazko bo moga wpadac aliasy modu 0.61

!odrzucenie aliasow jest jest zmiana okresu
if (czy_zo.eqv..true.) then
	
end if

! przeszukiwanie aliasow 2f ma sens wtedy, gdy jest zmiennosc modu podstawowego raczej - wyjdzie w praniu albo zapytac

!max w szukanym zakresie
fmin=f/0.64
fmax=f/0.58

call maximum(fs,n,fx,ax,fmin,fmax,k)
call sn(fs,n,fx,ax,0.2,snx)
write (*,*) ' fx: ',fx,'snx: ',snx, 'px/p=',f/fx
if (snx.gt.4.0) then
	write (*,*) 'JEST MOD!'
else
	write (*,*) 'nic'
end if
	
	
!write (line,'(a25,f11.7,6f20.16,a5)') plik, period,f1,sn1,f2,sn2,f3,sn3,flaga
!write (*,*) line
!call system('echo " '//line//' ">>residua.txt')

deallocate(fs)
close(fn)

CONTAINS

subroutine sn(fs,n,freq,amp,box,snr)
implicit none

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: n
	real(kind=8) :: freq,amp,snr
	real :: box
	
	integer :: i,j
	real(kind=8) :: fmin_out,fmin_in,fmax_out,fmax_in
	real(kind=8) :: signal,noise
	
	signal=amp
	
	fmin_out=max(freq-dble(box/2.0),dble(0.0))
	fmin_in=freq-dble(box/40.0)
	fmax_in=freq+dble(box/40.0)
	fmax_out=min(freq+dble(box/2.0),dble(10.0))
	
	j=0
	noise=0
	do i=1,n
		if(((fs(1,i).ge.fmin_out).and.(fs(1,i).le.fmin_in)).or.((fs(1,i).ge.fmax_in).and.(fs(1,i).le.fmax_out))) then
			noise=noise+fs(2,i)
			j=j+1
		end if
	end do
	noise=noise/dble(j)
	
	snr=signal/noise

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
	
	write (*,*) 'OK, wczytano dane do programu'
	
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
		if ((fs(1,i).le.fmax).and.(fs(1,i).ge.fmin).and.(fs(3,i).eq.1.0)) then
			if (fs(2,i).gt.ampmax) then
				ampmax=fs(2,i)
				freqmax=fs(1,i)
				k=i
			end if
		end if
	end do
	
!	write (*,*) 'Znaleziono maksimum o amplitudzie: ',ampmax,' dla czestotliwosci: ',freqmax

end subroutine



end program res3
