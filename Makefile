src = src/
bin = bin/

all:
	gfortran -c $(src)subroutines.F90 $(src)psi_energies.F90 $(src)Var_Dif.f90
	gfortran $(src)QMC.f90 -o $(bin)qmc.x $(src)subroutines.o $(src)psi_energies.o $(src)Var_Dif.o
tests:
	gfortran -DTEST	$(src)psi_energies.F90 -o $(bin)test_en.x
	gfortran -DTEST	$(src)subroutines.F90 -o $(bin)test_drift.x
