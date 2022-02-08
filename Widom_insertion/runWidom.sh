#!/bin/bash

gcc Widom.c -lm -O3 -o Widom.out
N=256
T=1.32
mkdir -p out_Widom_N${N}_T${T}
for rho in 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.9
do
rm paraWidom.par
echo "Step		200" >> paraWidom.par
echo "Ato			${N}" >> paraWidom.par
echo "FileName	Conf__N${N}_rho${rho}_T${T}-LAMMPSdump_" >> paraWidom.par
echo "StartAt		10000" >> paraWidom.par
echo "StopAt		40001" >> paraWidom.par
echo "Insert        1000" >> paraWidom.par
echo "eps         1.0" >> paraWidom.par
echo "sigma       1.0" >> paraWidom.par
echo "rc          2.5" >> paraWidom.par
echo "Temp        ${T}" >> paraWidom.par

./Widom.out paraWidom.par ../data_MSD_long/ out_Widom_N${N}_T${T}/Widom_thermo_rho${rho}.dat

done
