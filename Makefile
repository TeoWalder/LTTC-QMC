all:
	gfortran -c subroutines.F90 psi_energies.F90 Var_Dif.f90
	gfortran QMC.f90 -o qmc.x subroutines.o psi_energies.o Var_Dif.o
	rm *.o
tests:
	gfortran -DTEST	psi_energies.F90 -o test_en.x
	gfortran -DTEST	subroutines.F90 -o test_drift.x
