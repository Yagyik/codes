rm *~
rm */*~
rm */*/*~
clear 
g++ AnaCode_wT.c++ -L/Data/.gsl-2.2.1/lib/ -lgsl -lgslcblas -lm -DFFTW3 -lfftw3 -o a_wT.out
## debug version
##g++ -g AnaCode.c++ -lgsl -lgslcblas -lm -DFFTW3 -lfftw3


## to run ./a.out anaFile.dat /path/to/input/file/folder/ output_folder/
## to run debug valgrind --leak-check=yes --track-origins=yes ./a.out <args>
