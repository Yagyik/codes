# script to run the 2d free energy code that calculates the free energy two ways
# the first way is to stitch bDG(n) by minimising the error between adjacent windows (at the edges) and shifting the P(n,rho) along with P(n)
# here then bDG(n,rho)=-ln(P(n,rho)) - this is a special case of the wham code where we find the shift required to stitch the free energy 
# rather than find the self-consistent solution (FE) for each window that stitches to together.
# the second way is to use the 2dFE method bDG(n,rho) = bDG(n) - ln(P(n,rho)/P(n))

mkdir -p plotFiles
rm plotFiles/plot_rhoPE_histo.gnu
rm plotFiles/plot_2DFE_WHAM.gnu


echo "reset" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set terminal epslatex size 4cm,4cm" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set term eps" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set output output" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set pm3d at sb" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set samples 1000" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set isosamples 100" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set palette rgb 33,13,10 # colour map for the plot" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set tics scale 0.5" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set view map scale 0.8 #vertical view" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set contour base" >> plotFiles/plot_2DFE_WHAM.gnu
# echo "set cntrlabel format '%2.0g' font ',7'" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set cntrparam levels discrete 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 # contour min, spacing and max" >> plotFiles/plot_2DFE_WHAM.gnu
#echo "set cntrlabel onecolor" >> plotFiles/plot_2DFE_WHAM.gnu
echo "unset surface" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set label label at graph 0.75,0.12 center front tc default font \"Helvetica,20\" " >> plotFiles/plot_2DFE_WHAM.gnu
echo "show label" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set label \"C\" at screen 0.175,0.85 center front tc default font \"Helvetica,20\" " >> plotFiles/plot_2DFE_WHAM.gnu
echo "set xlabel \"{/Helvetica=20 n}\" offset 0,0 ; set xrange [0:14] ; set xtics 2 ; set mxtics 1" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set ylabel \"{/Helvetica=20 {/Symbol r} [g/cc]}\" offset -0.5 ; set yrange [2.25:2.55] ; set ytics 0.05 ; set mytics 3" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set cblabel \"{/Helvetica=20 {/Symbol bD}G(n,{/Symbol r}) \" offset 1.5; set cbrange [0:14] ; set cbtics 2.0 ; set mcbtics 1" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set size ratio 0.9" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set lmargin 0.01" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set rmargin 0.01" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set tmargin 0.01" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set bmargin 0.01" >> plotFiles/plot_2DFE_WHAM.gnu
echo "set datafile missing '9999.000000' # rows with this z value are masked/ignored" >> plotFiles/plot_2DFE_WHAM.gnu
#echo "set datafile missing 'inf' # rows with this z value are masked/ignored" >> plotFiles/plot_2DFE_WHAM.gnu
echo "splot datafile u 1:2:3 title \"\" w l lw 3"  >> plotFiles/plot_2DFE_WHAM.gnu


echo "reset" >> plotFiles/plot_rhoPE_histo.gnu
echo "set terminal epslatex size 4cm,4cm" >> plotFiles/plot_rhoPE_histo.gnu
echo "set term eps" >> plotFiles/plot_rhoPE_histo.gnu
echo "set output output" >> plotFiles/plot_rhoPE_histo.gnu
echo "set pm3d at sb" >> plotFiles/plot_rhoPE_histo.gnu
echo "set samples 1000" >> plotFiles/plot_rhoPE_histo.gnu
echo "set isosamples 100" >> plotFiles/plot_rhoPE_histo.gnu
echo "set palette rgb 33,13,10 # colour map for the plot" >> plotFiles/plot_rhoPE_histo.gnu
echo "set tics scale 0.5" >> plotFiles/plot_rhoPE_histo.gnu
echo "set view map scale 0.8 #vertical view" >> plotFiles/plot_rhoPE_histo.gnu
echo "set contour surface" >> plotFiles/plot_rhoPE_histo.gnu
echo "set cntrparam levels incremental -1.0,1.0,8.0 # contour min, spacing and max" >> plotFiles/plot_rhoPE_histo.gnu
echo "set cntrlabel onecolor" >> plotFiles/plot_rhoPE_histo.gnu
echo "unset surface" >> plotFiles/plot_rhoPE_histo.gnu
echo "set label label at graph 0.75,0.3 center front tc default font \"Helvetica,20\" " >> plotFiles/plot_rhoPE_histo.gnu
echo "show label" >> plotFiles/plot_rhoPE_histo.gnu
echo "set label label2 at screen 0.2,0.85 center front tc default font \"Helvetica,20\" " >> plotFiles/plot_rhoPE_histo.gnu
echo "set xlabel \"{/Helvetica=20 {/Symbol r} [g/cc]\" offset 0,0 ; set xrange [2.3:2.55] ; set xtics 0.05 ; set mxtics 3" >> plotFiles/plot_rhoPE_histo.gnu
echo "set ylabel \"{/Helvetica=20 PE [{/Symbol e}/N]}\" offset -0.5 ; set yrange [-1.9:-1.8] ; set ytics 0.02 ; set mytics 1" >> plotFiles/plot_rhoPE_histo.gnu
echo "set cblabel \"{/Helvetica=20 -ln [P({/Symbol r},PE)] \" offset 1.5; set cbrange [0:14] ; set cbtics 2.0 ; set mcbtics 1" >> plotFiles/plot_rhoPE_histo.gnu
echo "set size ratio 1.1" >> plotFiles/plot_rhoPE_histo.gnu
echo "set lmargin 0.01" >> plotFiles/plot_rhoPE_histo.gnu
echo "set rmargin 0.01" >> plotFiles/plot_rhoPE_histo.gnu
echo "set tmargin 0.01" >> plotFiles/plot_rhoPE_histo.gnu
echo "set bmargin 0.01" >> plotFiles/plot_rhoPE_histo.gnu
echo "set datafile missing '9999.000000' # rows with this z value are masked/ignored" >> plotFiles/plot_rhoPE_histo.gnu
#echo "set datafile missing 'inf' # rows with this z value are masked/ignored" >> plotFiles/plot_rhoPE_histo.gnu
echo "splot datafile u 1:2:3 title \"\" w l"  >> plotFiles/plot_rhoPE_histo.gnu


#!/bin/bash
conda activate base
rm make_unbiased_histo.exe
rm make_unbiased_histo_WHAM.exe
rm make_unbiased_histo_2DWHAM.exe
gcc -g make_unbiased_histo.c -lgsl -lgslcblas -lm -o make_unbiased_histo.exe
gcc -g make_unbiased_histo_WHAM.c -lgsl -lgslcblas -lm -o make_unbiased_histo_WHAM.exe
gcc -g make_unbiased_histo_2DWHAM.c -lgsl -lgslcblas -lm -o make_unbiased_histo_2DWHAM.exe


# temp=0.0395 #0.0387 0.0389 0.0391 0.0395 0.0399
# P=02
lim=3
n0=2

# temp=0.0387
P=00
##########################################
N=512
krho=4000 # N=512
bias=50 #N=512
skip=80000000 #N=512
last=100000000
path=data_TrhoPT_P${P}

# N=1000
# krho=6000 # N=1000
# bias=50 #N=1000
# skip=20000000 #N=1000
# last=27000000
# path=data_1k_TrhoPT_P${P}


# N=1500
# krho=8000 # N=1500
# bias=30 #N=1500
# skip=10000000 #N=1500
# last=15000000
# path=data_1k5_TrhoPT_P${P}
######################################################
# N=2000
# krho=10000 #N=2000
# bias=20 #N=2000
# skip=20000000 # N=2000
# last=30000000
# path=data_2k_TrhoPT_P${P}

##############################################
# P=04
# temp=0.0363 # 0.0355 0.0357 0.0359 0.0363
# 
# krho=2000 # N=512
# bias=50 #N=512
# skip=70000000 #N=512
# last=90000000
# path=data_512_rhoPT_P${P}

# lenAC=$(echo "scale=0; ${last}/500-${skip}/500" | bc)
# echo ${temp} ${lenAC}
# for temp in 0.0389 0.0391 0.0395
# do
# mkdir -p ${path}/unbiased_histo
# mkdir -p ${path}/AC_files/${temp}
# mkdir -p ${path}/PT_Track_plots/${temp}
# mkdir -p ${path}/PT_Track_plots/${temp}/nSp_max
# mkdir -p ${path}/PT_Track_plots/${temp}/rho
# mkdir -p ${path}/PT_Track_plots/${temp}/Tem
# mkdir -p ${path}/PT_stats/${temp}
# rm ${path}/AC_files/${temp}/IntegAC.dat
# 
# for rank in 2 5 8 11 14 17
# do
# for rho in 0.454 0.462 0.47 0.478 0.486 0.494
# do 
# python3 auto_corr_script.py ${path} ${rank} ${temp} ${rho} ${krho} ${bias} ${lenAC} ${skip}
# done # rho
# echo "AC"
# echo ${path} ${tempP} ${nTemp} ${nnTemp}
# done #rank
# 
# done #temp

numTem=4
tempArr="0.0387 0.0389 0.0391 0.0395"

rm ParaFiles/PT_para.par
echo "totSize 144" >> ParaFiles/PT_para.par
echo "num_nSp_max 6" >> ParaFiles/PT_para.par
echo "d_nSp_max 3" >> ParaFiles/PT_para.par
echo "nSp_max_0 2" >> ParaFiles/PT_para.par
echo "num_rho 6" >> ParaFiles/PT_para.par
echo "d_rho 0.008" >> ParaFiles/PT_para.par
echo "rho_0 0.454" >> ParaFiles/PT_para.par
echo "num_Tem ${numTem}" >> ParaFiles/PT_para.par
for i in ${tempArr}
do
echo "tem ${i}" >> ParaFiles/PT_para.par
done


python3 auto_corr_PT_track_Masterscript.py ${path} ParaFiles/PT_para.par ${temp} ${lenAC}



# for rho0 in 0.454 0.462 0.47 0.478 0.486 0.494
# do
# 
# mkdir -p ParaFiles
# rm ParaFiles/parafile.par
# 
# ## first append thermo files
# 
# 
# 
# ## later add a temp or n0 loop
# 
# 
# echo "rhomin    0.45" >> ParaFiles/parafile.par
# echo "rhomax    0.52" >> ParaFiles/parafile.par
# echo "nrhobin    60" >> ParaFiles/parafile.par
# echo "n0  ${n0}" >> ParaFiles/parafile.par
# echo "rho0  ${rho0}" >> ParaFiles/parafile.par
# echo "krho ${krho}" >> ParaFiles/parafile.par
# echo "temp ${temp}" >> ParaFiles/parafile.par
# echo "bias ${bias}" >> ParaFiles/parafile.par
# echo "skip ${skip}" >> ParaFiles/parafile.par
# echo "last ${last}" >> ParaFiles/parafile.par
# 
# 
# 
# 
# 
# for n0 in 2 5 8 11 14 17
# do
# rm ${path}/thermo_compile/Conf_nSp${n0}-KXX-rho${rho0}K${krho}_T${temp}-bias${bias}-thermo.dat
# paste <(awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${path}/thermo_compile/Conf_nSp${n0}-K0.01-rho${rho0}K${krho}_T${temp}-bias${bias}-thermo.dat) >> ${path}/thermo_compile/Conf_nSp${n0}-KXX-rho${rho0}K${krho}_T${temp}-bias${bias}-thermo.dat
# # paste <(awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' ${path}/thermo_compile/Conf_nSp${n0}-K0-rho${rho0}K${krho}_T${temp}-bias${bias}-thermo.dat) >> ${path}/thermo_compile/Conf_nSp${n0}-KXX-rho${rho0}K${krho}_T${temp}-bias${bias}-thermo.dat
# done

# valgrind --leak-check=yes --track-origins=yes ./make_unbiased_histo.exe ParaFiles/parafile.par ${path} ${lim}
./make_unbiased_histo.exe ParaFiles/parafile.par ${path} ${lim}


done

n0=2
rho0=0.454
drho0=0.008
drhow=0.010
drw=`printf "%0.3f\n" ${drhow}`
for temp in 0.0389 0.0391 0.0395
do

tempP=`printf "%0.4f\n" ${temp}`
nTemp=$(echo "scale=0; 1+25173*${tempP}" | bc)
nnTemp=`printf "%d\n" ${nTemp}`
rm ParaFiles/paraWHAM.par
echo "rhomin    0.45" >> ParaFiles/paraWHAM.par
echo "rhomax    0.52" >> ParaFiles/paraWHAM.par
echo "nrhobin    60" >> ParaFiles/paraWHAM.par
echo "n0  ${n0}" >> ParaFiles/paraWHAM.par
echo "rho0  ${rho0}" >> ParaFiles/paraWHAM.par
echo "drho0  ${drho0}" >> ParaFiles/paraWHAM.par
echo "drhow  ${drhow}" >> ParaFiles/paraWHAM.par
echo "krho ${krho}" >> ParaFiles/paraWHAM.par
echo "PEmin -1.9" >> ParaFiles/paraWHAM.par
echo "PEmax -1.8" >> ParaFiles/paraWHAM.par
echo "nPEbin 40" >> ParaFiles/paraWHAM.par
echo "nAto ${N}" >> ParaFiles/paraWHAM.par
echo "temp ${temp}" >> ParaFiles/paraWHAM.par
echo "bias ${bias}" >> ParaFiles/paraWHAM.par
echo "skip ${skip}" >> ParaFiles/paraWHAM.par
echo "last ${last}" >> ParaFiles/paraWHAM.par
echo "nWin 6" >> ParaFiles/paraWHAM.par
echo "tol 0.0001" >> ParaFiles/paraWHAM.par

if [ ${temp} = 0.0389 ]
then
label2=F
fi

if [ ${temp} = 0.0391 ]
then
label2=C
fi

if [ ${temp} = 0.0395 ]
then
label2=D
fi

# valgrind --leak-check=yes --track-origins=yes ./make_unbiased_histo_WHAM.exe ParaFiles/paraWHAM.par ${path} ${lim}
./make_unbiased_histo_WHAM.exe ParaFiles/paraWHAM.par ${path} ${lim}
gnuplot -e "datafile='${path}/whamPrhoPE-drw${drw}-${tempP}-lim${lim}.dat' ; output='${path}/WHAMrhoPE_histo-drw${drw}-${nnTemp}-lim${lim}.eps' ; label='T=${nnTemp}K'; label2='${label2}'" plotFiles/plot_rhoPE_histo.gnu


rm ParaFiles/paraWHAM2D.par
echo "rhomin    0.45" >> ParaFiles/paraWHAM2D.par
echo "rhomax    0.52" >> ParaFiles/paraWHAM2D.par
echo "nrhobin    35" >> ParaFiles/paraWHAM2D.par
echo "nrholines    100" >> ParaFiles/paraWHAM2D.par
echo "ncbin  20" >> ParaFiles/paraWHAM2D.par
echo "n0  ${n0}" >> ParaFiles/paraWHAM2D.par
echo "dn0  3" >> ParaFiles/paraWHAM2D.par
echo "dnw  2" >> ParaFiles/paraWHAM2D.par
echo "rho0  ${rho0}" >> ParaFiles/paraWHAM2D.par
echo "drho0  ${drho0}" >> ParaFiles/paraWHAM2D.par
echo "drhow  ${drhow}" >> ParaFiles/paraWHAM2D.par
echo "krho ${krho}" >> ParaFiles/paraWHAM2D.par
echo "temp ${temp}" >> ParaFiles/paraWHAM2D.par
echo "bias ${bias}" >> ParaFiles/paraWHAM2D.par
echo "nWin 36" >> ParaFiles/paraWHAM2D.par
echo "rnWin 6" >> ParaFiles/paraWHAM2D.par
echo "nnWin 6" >> ParaFiles/paraWHAM2D.par
echo "tol 0.00001" >> ParaFiles/paraWHAM2D.par


label2=A
# valgrind --leak-check=yes --track-origins=yes ./make_unbiased_histo_WHAM.exe ParaFiles/paraWHAM.par ${path} ${lim}
./make_unbiased_histo_2DWHAM.exe ParaFiles/paraWHAM2D.par ${path} ${lim}
gnuplot -e "datafile='${path}/whamPnrho-drw${drw}-${tempP}-lim${lim}.dat' ; output='${path}/WHAM2DFEnrho-drw${drw}-${nnTemp}-lim${lim}.eps' ; label='T=${nnTemp}K'; label2='${label2}'" plotFiles/plot_2DFE_WHAM.gnu




done


# for i in 0 # 02 #04
# do
# inp=data_2D_TPT_16_P${i}
# mkdir -p ${inp}/AC_files
# mkdir -p ${inp}/PT_Track_plots/
# # python3 auto_corr_PT_track_script.py ${inp}
# 
# if [ $i = 0 ]
# then
# mkdir -p ${inp}/AC_files
# 
# for temp in 0.0419 0.042 0.0423 0.0425 0.0429 #0.042
# do
# echo ${temp}
# tempP=`printf "%0.4f\n" ${temp}`
# nTemp=$(echo "scale=0; 1+25173*${tempP}" | bc)
# nnTemp=`printf "%d\n" ${nTemp}`
# echo ${inp} ${tempP} ${nTemp} ${nnTemp}
# mkdir -p ${inp}/AC_files/${temp}
# 
# # for rank in {1..16}
# # do
# # #python3 auto_corr_script.py ${inp} ${rank} ${temp} 
# # #echo "AC"
# # echo ${inp} ${tempP} ${nTemp} ${nnTemp}
# # done #rank
# ./2d_twoWays.out ParaFiles/parafile.par ${inp} ${temp} Pn Pnrho
# echo ${tempP}
# gnuplot -e "datafile='${inp}/Pnrho-${tempP}.dat' ;  output='${inp}/2DFE_Pnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_FE2D_nmaxrho.gnu
# gnuplot -e "datafile='${inp}/2df1dPnrho-${tempP}.dat' ;  output='${inp}/2DFE_2df1dPnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_2DFE_WHAM.gnu
# gnuplot -e "datafile='${inp}/raw2df1dPnrho-${tempP}.dat' ;  output='${inp}/raw2DFE_2df1dPnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_2DFE_WHAM.gnu
# 
# done # temp
# fi
# 
# 
# 
# 
# if [ $i = 02 ]
# then
# 
# for temp in 0.0383 0.0387 0.0391 0.0395 0.0399
# do
# tempP=`printf "%0.4f\n" ${temp}`
# nTemp=$(echo "scale=0; 1+25173*${tempP}" | bc)
# nnTemp=`printf "%d\n" ${nTemp}`
# echo ${inp} ${tempP} ${nTemp} ${nnTemp}
# mkdir -p ${inp}/AC_files/${tempP}
# # for rank in {1..16}
# # do
# # #python3 auto_corr_script.py ${inp} ${rank} ${temp} 
# # #echo "AC"
# # echo ${inp} ${tempP} ${nTemp} ${nnTemp}
# # done #rank
# ./2d_twoWays.out ParaFiles/parafile.par ${inp} ${temp} Pn Pnrho
# 
# gnuplot -e "datafile='${inp}/Pnrho-${temp}.dat' ;  output='${inp}/2DFE_Pnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_FE2D_nmaxrho.gnu
# gnuplot -e "datafile='${inp}/2df1dPnrho-${temp}.dat' ;  output='${inp}/2DFE_2df1dPnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_2DFE_WHAM.gnu
# gnuplot -e "datafile='${inp}/raw2df1dPnrho-${temp}.dat' ;  output='${inp}/raw2DFE_2df1dPnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_2DFE_WHAM.gnu
# done # temp
# fi
# 
# if [ $i = 04 ]
# then
# 
# for temp in 0.0351 0.0355 0.0359 0.0363 0.0367
# do
# tempP=`printf "%0.4f\n" ${temp}`
# nTemp=$(echo "scale=0; 1+25173*${tempP}" | bc)
# nnTemp=`printf "%d\n" ${nTemp}`
# echo ${inp} ${tempP} ${nTemp} ${nnTemp}
# mkdir -p ${inp}/AC_files/${tempP}
# # for rank in {1..16}
# # do
# # python3 auto_corr_script.py ${inp} ${rank} ${temp} 
# # # echo "AC"
# # # echo ${inp} ${tempP} ${nTemp} ${nnTemp}
# # done #rank
# ./2d_twoWays.out ParaFiles/parafile.par ${inp} ${temp} Pn Pnrho
# 
# gnuplot -e "datafile='${inp}/Pnrho-${temp}.dat' ;  output='${inp}/2DFE_Pnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_FE2D_nmaxrho.gnu
# gnuplot -e "datafile='${inp}/2df1dPnrho-${temp}.dat' ;  output='${inp}/2DFE_2df1dPnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_2DFE_WHAM.gnu
# gnuplot -e "datafile='${inp}/raw2df1dPnrho-${temp}.dat' ;  output='${inp}/raw2DFE_2df1dPnrho-${nnTemp}.eps'; label='T=${nnTemp}K'" plotFiles/plot_2DFE_WHAM.gnu
# done # temp
# fi
# 
# 
# done #i
# 
conda deactivate
