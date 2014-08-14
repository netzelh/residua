res3: res3.o io.o tools.o pik.o
	gfortran res3.o io.o tools.o pik.o -o res3
	rm *.o *.mod
	
tools.o: tools.f90
	gfortran -c tools.f90
pik.o: pik.f90 tools.o
	gfortran -c pik.f90
io.o: io.f90
	gfortran -c io.f90
res3.o: res3.f90 io.o pik.o
	gfortran -c res3.f90
	
clean:
	rm *.o *.mod
