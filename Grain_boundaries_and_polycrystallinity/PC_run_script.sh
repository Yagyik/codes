mkdir -p anaFiles

for T in 1055 #1070
do
tt=`echo "scale=4; $T/25173" | bc`
mkdir -p out_T_${T}
for q6 in 0.05 0.04
do

mkdir -p out_T_${T}/q6_${q6}

for rhok in K10000-rho0.485K0 K0-rho0.485K0 #rho0.455K10000 rho0.465K10000 rho0.475K10000 rho0.485K10000 rho0.485K0 rho0.47K10000
#for rhok in rho0.47K10000 rho0.475K10000 rho0.485K10000 rho0.485K0
do

rm anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
# make the anafile

echo "Step		10000" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Ato			512" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "DeltaT		0.3830" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "FileName	Conf_q6${q6}-${rhok}_T0${tt}-bias50-LAMMPSdump_" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "StartAt		20000000" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "StopAt		50000001" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Gr			0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Gr_limit	10" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Msd			0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Msd_limit	1800" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Isf			0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Isf_limit	7" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "IsfK_value	3.50" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "OrPar		0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "OrPar_limit	10" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "OrPar_r;	1.4" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "OrPar_bin;	100" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Sq			0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Sq_limit	20" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Sq_max		40" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Sq_bin		400" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Coher		0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Coher_limit	101" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Coher_value	3.50" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "G5r			0	" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "G5r_limit	600" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "G5r_max		4.0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "ClSize		1" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Cl_limit	3000" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Cl_r		1.4" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Cl_solid	3" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "Cl_cut		1.8" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "CoordNO		0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "cn_limit 	6" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "CoordLim	1.43" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "interface	1" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "int_limit	3000" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "int_tol		0.5" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "opla		0" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "opla_lim	1" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "opla_bin	100" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "laqnbin		500" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "opla_r		1.4" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat
echo "opla_ra		1.8" >> anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat

# make the directory for the set

mkdir -p out_T_${T}/q6_${q6}/${rhok}
# copy the thermofile there for future reference
cp ../data_512_${T}_q6-${q6}/Conf_q6${q6}-${rhok}_T0${tt}-bias50-thermo.dat out_T_${T}/q6_${q6}/${rhok}/

#compile and run
sh compileMulti_nT.sh
#valgrind --leak-check=yes --track-origins=yes ./aMulti_nT.out anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat ../data_512_${T}_q6-${q6}/ out_T_${T}/q6_${q6}/${rhok}/
./aMulti_nT.out anaFiles/ana_q6pc-${q6}-${rhok}-T${T}.dat ../data_512_${T}_q6-${q6}/ out_T_${T}/q6_${q6}/${rhok}/
done
done

done
