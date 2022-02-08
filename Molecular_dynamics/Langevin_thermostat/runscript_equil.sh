clear 
mpicc -g md_LJ_TS.c -lgsl -lgslcblas -lm -O3 -o md_ljts.exe
dmax=0.2


#### Long runs
mkdir -p data_equil

rm para_LJequil.par
echo "Natoms		256" >> para_LJequil.par
echo "ens         0" >> para_LJequil.par
echo "eps         1.0" >> para_LJequil.par
echo "sigma       1.0" >> para_LJequil.par
echo "rc          2.5" >> para_LJequil.par
echo "dt		0.001" >> para_LJequil.par
echo "EqStep		0" >> para_LJequil.par
echo "TotStep		40000" >> para_LJequil.par
echo "Dump	    1000" >> para_LJequil.par
echo "Thermo		10" >> para_LJequil.par
echo "Seed		34592" >> para_LJequil.par
echo "Init		0" >> para_LJequil.par
echo "Label       Conf_" >> para_LJequil.par
echo "startCry    1" >> para_LJequil.par
echo "gamma 0.1" >> para_LJequil.par


rm source_LJ_equil.dat

# echo "0.6     0.0    4.0" >> source_LJ_equil.dat
echo "0.8     0.0    2.0" >> source_LJ_equil.dat
# echo "0.6     0.0    1.1" >> source_LJ_equil.dat


#mpirun -np 1 valgrind --leak-check=yes --track-origins=yes ./md_ljts.exe para_LJequil.par source_LJequil.dat dataequil/
mpirun -np 1 ./md_ljts.exe para_LJequil.par source_LJ_equil.dat data_equil/



