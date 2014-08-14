FC = gfortran
FCFLAGS = -c -O2# \
		  -ggdb \
		  -Wall -Wsurprising -Wextra -Wunderflow -pedantic \
		  -fcheck=all \
		  -fbacktrace \
		  -std=f2003 \
		  
FLFLAGS = #-ggdb -fbacktrace

main: main.o boxflow_solvers.o
	$(FC) $(FLFLAGS) -o main *.o

main.o: main.f03 boxflow_solvers.mod
	$(FC) $(FCFLAGS) main.f03

boxflow_solvers.o: boxflow_solvers.f03
	$(FC) $(FCFLAGS) boxflow_solvers.f03 

boxflow_solvers.mod: boxflow_solvers.f03
	$(FC) $(FCFLAGS) boxflow_solvers.f03

run: main
	./main

clean:
	rm *.o *.mod *.txt fort.6

cleanall:
	rm main *.o *.mod *.txt fort.6


