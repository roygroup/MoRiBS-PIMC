
#FC=ifort
FC=~pnroy/Dev/bin/gfortran

#LDFLAGS=-L/opt/intel/mkl/lib/intel64 -llapack

symrho.x:	symrho.o
		$(FC) -o $@ symrho.o 

.f.o:
	${FC} -O3 -c $<

clean:
	rm -f *.o symrho.x
