import numpy
import matplotlib.pyplot as plt
import math
import time
import sys





#### characteristic of the file(s)

length = int(sys.argv[7])
hLength = int(length/2)
step = 500
filerho=float(sys.argv[4])
nfiles=16
rank=int(sys.argv[2])
temp=float(sys.argv[3])
path=sys.argv[1]
########
#filename='data_test_1070_12_4/thermo_compile/Conf_nSp1-K0.01-rho0.456K0_T0.0425-bias40-thermo.dat')
#filename=sys.argv[1]+'/thermo_compile/Conf_nSp'+sys.argv[2]+'-K0.01-rho'+str(filerho)+'K0_T'+sys.argv[3]+'-bias40-thermo.dat'
filename=sys.argv[1]+'/thermo_compile/Conf_nSp'+sys.argv[2]+'-K0.01-rho'+str(filerho)+'K'+sys.argv[5]+'_T'+sys.argv[3]+'-bias'+sys.argv[6]+'-thermo.dat'
#filename=sys.argv[1]+'/thermo_compile/Conf_nSp'+sys.argv[2]+'-K0.01-rho'+str(filerho)+'K0_T'+sys.argv[3]+'-bias20-thermo.dat'

#filename='1k_h60_1055_4/thermo_compile/Conf_Hard60_nSp0-K0.01-rho'+str(filerho)+'K40000_T0.0419-bias50-thermo.dat'
print(filename)
infile=open(filename,'r')
lines=infile.readlines()
infile.close()
n_skip=int(sys.argv[8])
n_skip=int(n_skip/step)
print(length,n_skip)

l_skip=1
spam=0
n=0
val=numpy.zeros([3],numpy.float64)
en=numpy.zeros([length],numpy.float64)
rho=numpy.zeros([length],numpy.float64)
q6=numpy.zeros([length],numpy.float64)
for line in lines:
     if line[0] != '#' and line[0] != '@' and spam<length:
            if n >= n_skip and n%l_skip==0 :
                tokens = line.split()
                val[0] = float(tokens[1])
                val[1] = float(tokens[3])
                val[2] = float(tokens[4])
                en[spam]=val[0]
                rho[spam]=val[1]
                q6[spam]=val[2]
                spam+=1
            n+=1


print(q6[-1],rho[-1],en[-1])
print(q6[length-1],rho[length-1],en[length-1])
print(q6[0],rho[0],en[0])
#var_en=en.var()
x=numpy.arange(0,hLength*step,step)            
en=en-numpy.mean(en)
rho=rho-numpy.mean(rho)
q6=q6-numpy.mean(q6)
norm_en=numpy.sum(en**2)
norm_rho=numpy.sum(rho**2)
norm_q6=numpy.sum(q6**2)
result_en=numpy.correlate(en,en,mode='same')/norm_en
result_rho=numpy.correlate(rho,rho,mode='same')/norm_rho
result_q6=numpy.correlate(q6,q6,mode='same')/norm_q6
print(result_en.size,result_en.size/2)
res_final_en =result_en[int(result_en.size/2):]
res_final_rho =result_rho[int(result_rho.size/2):]
res_final_q6 =result_q6[int(result_q6.size/2):]


### log bin our data
#range is typically from 10^3 to 10^7 (3-7). make 200 bins
dlogt=(10.0-2.0)/300.0
mint=math.log(step,10)
nc=numpy.zeros(300,numpy.float64)
c_en=numpy.zeros(300,numpy.float64)
c_rho=numpy.zeros(300,numpy.float64)
c_q6=numpy.zeros(300,numpy.float64)
for i in range(hLength):
    if i >0:
        log10t=math.log(i*step,10)
        spam = int(log10t/dlogt)
    else:
        spam = 0
    if spam<300:
        nc[spam]+=1.0
        c_en[spam]+=res_final_en[i]
        c_rho[spam]+=res_final_rho[i]
        c_q6[spam]+=res_final_q6[i]

f=open(sys.argv[1]+"/AC_files/"+sys.argv[3]+"/AC_logbin_rho"+sys.argv[4]+"rank_"+sys.argv[2]+".dat","w+") # put proper output file name with temp and fblamb here
f.write("%8d %4.4f %4.4f %4.4f\n" % (0,1.0,1.0,1.0))
for i in range(300):
    if nc[i]>0:
        f.write("%8d %4.4f %4.4f %4.4f\n" % (math.pow(10.0,((i+0.5)*dlogt)),c_en[i]/nc[i],c_rho[i]/nc[i],c_q6[i]/nc[i]))
f.close()

f=open(sys.argv[1]+"/AC_files/"+sys.argv[3]+"/AC_rho"+sys.argv[4]+"rank_"+sys.argv[2]+".dat","w+") # put proper output file name with temp and fblamb here
for i in range(hLength):
     f.write("%8d %4.4f %4.4f %4.4f\n" % (i*step,res_final_en[i],res_final_rho[i],res_final_q6[i]))
f.close()
print("length is ",len(result_en[int(result_en.size/2):]))

plt.plot(x,res_final_en,lw=2,label='PE autocorr')
plt.plot(x,res_final_rho,lw=2,label='density autocorr')
plt.plot(x,res_final_q6,lw=2,label='Q6 autocorr')

plt.xscale('log')
plt.legend()
#plt.show(block=False)
#time.sleep(5)
#plt.pause(5)
plt.savefig(sys.argv[1]+'/AC_files/'+sys.argv[3]+'/AC_rho'+sys.argv[4]+'rank_'+sys.argv[2]+'.png')
plt.close()
    
#### integrated auto-correlation for density and q6
grho=0.0
gq6=0.0
for i in range(hLength):
    grho+=(1-i/hLength)*step*res_final_rho[i]
    gq6+=(1-i/hLength)*step*res_final_q6[i]

f=open(sys.argv[1]+"/AC_files/"+sys.argv[3]+"/IntegAC.dat","a")
#if(int(sys.argv[2])==2):
f.write("%s %s %4.4f %4.4f\n" % (sys.argv[4],sys.argv[2],grho,gq6))

f.close()
            
            
        
