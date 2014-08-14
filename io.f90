MODULE io
implicit none

CONTAINS

subroutine flagi(f,fx,plik,period_ch,czy_al,czy_bl,czy_mod,czy_p2,czy_trend,czy_zo)
implicit none

	character (len=*) :: plik,period_ch
	real (kind=8) :: f,fx
	logical :: czy_al,czy_bl,czy_mod,czy_p2,czy_trend,czy_zo

	character (len=500) :: line,pom
	character (len=20) :: sx_ch
	!character (len=*) :: bl
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

	!write (*,*) 

	!pom=line
	!line=trim(pom)//trim(bl)

	!write (*,*) trim(line)
	!write (*,*) trim(blazko)
	
	call system('echo " '//trim(line)//' " >> rrc2.txt')

end subroutine

subroutine get_arg(plik,period,f,r)
implicit none

	character (len=30) :: plik
	character (len=30) :: period_ch
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

subroutine get_file(plik,fn,fs,n)
implicit none

	character (len=*) :: plik

	real(kind=8), dimension (:,:), allocatable :: fs
	integer :: fn !file number
	integer :: n !liczba wierszy

	open(fn,file=plik)
	write (*,*) 'Plik: ',plik
	n=ile_wierszy(fn)
	allocate (fs(3,n)) !1-freq,2-amp,3-zero/jeden
	call getp(fs,fn,n)

end subroutine

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

end module io
