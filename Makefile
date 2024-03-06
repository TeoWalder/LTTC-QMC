all:
	gfortran -c utils.F90 psi_energies.F90 Var_Dif.F90
	gfortran QMC.F90 -o qmc.x utils.o psi_energies.o Var_Dif.o
	rm *.o
test:
	gfortran -DTEST	psi_energies.F90 -o test_en.x
