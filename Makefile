all:
	gfortran -c utils.f90 energies_psi.f90 Var_Diff.F90
	gfortran QMC.F90 -o qmc.x utils.o energies_psi.o Var_Diff.o
	rm *.o
