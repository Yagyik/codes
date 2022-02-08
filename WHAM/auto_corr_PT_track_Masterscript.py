import numpy
import matplotlib.pyplot as plt
import time
import math
import sys




#### characteristic of the file(s)

length = int(sys.argv[4])
hLength = int(length/2)
step = 500
filerho=0.456
path=sys.argv[1]


# parameters
filename=sys.argv[2]
parafile=open(filename,'r')

line=parafile.readline()
tokens = line.split()
totSize=int(tokens[1])
line=parafile.readline()
tokens = line.split()
num_nSp_max=int(tokens[1])
line=parafile.readline()
tokens = line.split()
d_nSp_max=int(tokens[1])
line=parafile.readline()
tokens = line.split()
nSp_max_0=int(tokens[1])
line=parafile.readline()
tokens = line.split()
num_rho=int(tokens[1])
line=parafile.readline()
tokens = line.split()
d_rho=float(tokens[1])
line=parafile.readline()
tokens = line.split()
rho_0=float(tokens[1])
line=parafile.readline()
tokens = line.split()
num_Tem=int(tokens[1])

tems=numpy.zeros([num_Tem],numpy.float64)
for i in range(num_Tem):
    line=parafile.readline()
    tokens = line.split()
    tems[i]=float(tokens[1])
line=parafile.readline()
tokens = line.split()
plot_nSp_max_ac=int(tokens[1])
line=parafile.readline()
tokens = line.split()
plot_Tem_ac=int(tokens[1])
line=parafile.readline()
tokens = line.split()
plot_rho_ac=int(tokens[1])
parafile.close()

        

##totSize=int(sys.argv[2])
##num_nSp_max=int(sys.argv[3])
##d_nSp_max=float(sys.argv[4])
##nSp_max_0=float(sys.argv[5])
##num_rho=int(sys.argv[6])
##d_rho=float(sys.argv[7])
##rho_0=float(sys.argv[8])
##num_Tem=int(sys.argv[9])
##d_Tem=float(sys.argv[10])
##Tem_0=float(sys.argv[11])
skip=200
actual_length=0
outlay=numpy.zeros([totSize,3],numpy.float64)

nSp_max_home=numpy.zeros([totSize,2],numpy.float64)
rho_home=numpy.zeros([totSize,2],numpy.float64)
Tem_home=numpy.zeros([totSize,2],numpy.float64)

val=numpy.zeros([2],numpy.float64)
time=numpy.zeros([length],numpy.float64)
window=numpy.zeros([length],numpy.float64)
print("%d %d %d %d\n" % (num_Tem,num_rho,num_nSp_max,num_Tem*num_rho*num_nSp_max))
for i in range(num_Tem):
    for j in range(num_rho):
        for k in range(num_nSp_max):
            index=i*num_nSp_max*num_rho+j*num_nSp_max+k
            outlay[index,0]=nSp_max_0+k*d_nSp_max
            outlay[index,1]=rho_0+j*d_rho
            outlay[index,2]=tems[i]
            nSp_max_home[index,0]=index%num_nSp_max
            rho_home[index,0]=int(index/num_nSp_max)
            Tem_home[index,0]=int(index/(num_nSp_max*num_rho))
            nSp_max_home[index,1]=outlay[index,0]
            rho_home[index,1]=outlay[index,1]
            Tem_home[index,1]=outlay[index,2]
            print("%d ind %d %d %d n %4.4f %d %4.4f rho %4.4f tem\n" % (index,nSp_max_0,k,d_nSp_max,rho_0,j,d_rho,tems[i])) 

########


nv_nSp_max_home=numpy.zeros([totSize],numpy.float64)
nv_nSp_max_top=numpy.zeros([totSize],numpy.float64)
nv_nSp_max_bot=numpy.zeros([totSize],numpy.float64)
exc_len_nSp_max=numpy.zeros([totSize],numpy.float64)
num_exc_nSp_max=numpy.zeros([totSize],numpy.float64)
max_exc_nSp_max=numpy.zeros([totSize],numpy.float64)


if num_nSp_max>1 and plot_nSp_max_ac==1:
    for rank in range(totSize):
        filename=sys.argv[1]+'/PT_Track_compile/nSp_max/PT_Track_nSp_max_rank'+str(rank)+'.dat'
        infile=open(filename,'r')
        lines=infile.readlines()
        infile.close()
        if len(lines)==0:
            continue
        spam=0
        n=0
        val=numpy.zeros([2],numpy.float64)
        time=numpy.zeros([length],numpy.float64)
        window=numpy.zeros([length],numpy.float64)
        
        old_val=float(lines[0].split()[0]) ## start point
        start=float(lines[0].split()[0])-skip
        oldtime_nSp_max=start
        old_val=old_val-start
        print("start for rank %8d is %4.4f\n" % (rank,start))
        old_val1=nSp_max_home[rank]
        actual_length=0
        for line in lines:
            if actual_length>=300000:
                break
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                val[0] = float(tokens[0])
                val[1] = float(tokens[1])
                
                if val[1] == nSp_max_home[rank,1]: ## back home
                    nv_nSp_max_home[rank]+=1
                    exc_len_nSp_max[rank]+=val[0] - oldtime_nSp_max
                    num_exc_nSp_max[rank]+=1
                    
##                    print("%4.4f %4.4f %4.4f %4.4f\n" % (val[0],oldtime_nSp_max,val[0]-oldtime_nSp_max,max_exc_nSp_max[rank]))
                    if val[0] - oldtime_nSp_max > max_exc_nSp_max[rank]:
                        max_exc_nSp_max[rank]=val[0] - oldtime_nSp_max
                    oldtime_nSp_max=val[0]
                if val[1] == nSp_max_0:
                    nv_nSp_max_bot[rank]+=1
                if val[1] == nSp_max_0+(num_nSp_max-1)*d_nSp_max:
                   nv_nSp_max_top[rank]+=1

                
                spam=int((val[0] -start)/skip)   
                if spam<length:
                    if spam-old_val>1:                
                        for j in range(old_val,spam,1):
                            window[j]=old_val1
                    window[int(spam)]=val[1]
                    old_val=spam
                    old_val1=val[1]
                actual_length+=1
        #var_en=en.var()
        window_new=numpy.zeros([actual_length])
        window_new=window[:actual_length]-numpy.mean(window[:actual_length])
        norm_window_new=numpy.sum(window_new**2)

        result_window=numpy.correlate(window_new,window_new,mode='same')/norm_window_new
        print(result_window)
        print(window[:actual_length])
        print(result_window.size/2,actual_length)
        res_final_window =result_window[int(actual_length/2):]
        
        ### log bin our data
        #range is typically from 10^2 to 10^7 (3-7). make 200 bins
        dlogt=(8.0-2.0)/300.0
        mint=math.log(skip,10)
        nc=numpy.zeros(300,numpy.float64)
        c_pt=numpy.zeros(300,numpy.float64)
        
        for i in range(int(0.5*actual_length)):
            if i >0:
                log10t=math.log(i*skip,10)
                spam = int(log10t/dlogt)
            else:
                spam = 0
            if spam<300:
                nc[spam]+=1.0
                c_pt[spam]+=res_final_window[i]
                

        f=open(sys.argv[1]+"/PT_Track_plots/"+sys.argv[3]+"/nSp_max/track_ac_logbin_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+"_"+str(rank)+".dat","w+") # put proper output file name with temp and fblamb here
        f.write("%8d %4.4f\n" % (0,1.0))
        for i in range(300):
            if nc[i]>0:
                f.write("%8d %4.4f\n" % (math.pow(10.0,((i+0.5)*dlogt)),c_pt[i]/nc[i]))
        f.close()

        x=numpy.arange(0,0.5*actual_length*skip,skip)
        f=open(sys.argv[1]+"/PT_Track_plots/"+sys.argv[3]+"/nSp_max/track_ac_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+"_"+str(rank)+".dat","w+") # put proper output file name with temp and fblamb here
        for i in range(int(0.5*actual_length)):
             f.write("%8d %4.4f\n" % (i*skip,res_final_window[i]))
        f.close()
        print("length is ",len(result_window[int(result_window.size/2):]))

        plt.plot(x,res_final_window,lw=2,label='PT track autocorr')


        plt.xscale('log')
        plt.legend()
        #plt.show(block=False)
        #time.sleep(5)
        #plt.pause(5)
        plt.savefig(sys.argv[1]+'/PT_Track_plots/'+sys.argv[3]+'/nSp_max/track_ac_n'+str(outlay[rank][0])+'_r'+str(outlay[rank][1])+'_t'+str(outlay[rank][2])+'_'+str(rank)+'.png')
        plt.close()

        
            
    f=open(sys.argv[1]+"/PT_stats/"+sys.argv[3]+"/nSp_max_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+".dat","w+")
    for i in range(totSize):
        f.write("%8d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n" % (i,nSp_max_home[i,1],rho_home[i,1],\
                                                                                 Tem_home[i,1],nv_nSp_max_bot[i],\
                                                                                 nv_nSp_max_top[i],nv_nSp_max_home[i],\
                                                                     exc_len_nSp_max[i]/(num_exc_nSp_max[i]+1)\
                                                                        ,max_exc_nSp_max[i]));

    f.close
    print("printed PT stats nSp_max")

######################################
####### Do for temp###################


nv_Tem_home=numpy.zeros([totSize],numpy.float64)
nv_Tem_top=numpy.zeros([totSize],numpy.float64)
nv_Tem_bot=numpy.zeros([totSize],numpy.float64)
exc_len_Tem=numpy.zeros([totSize],numpy.float64)
num_exc_Tem=numpy.zeros([totSize],numpy.float64)
max_exc_Tem=numpy.zeros([totSize],numpy.float64)

if num_Tem>1 and plot_Tem_ac==1:
    for rank in range(totSize):
        filename=sys.argv[1]+'/PT_Track_compile/Tem/PT_Track_Tem_rank'+str(rank)+'.dat'
        infile=open(filename,'r')
        lines=infile.readlines()
        infile.close()
        if len(lines)==0:
            continue
        spam=0
        n=0
        val=numpy.zeros([2],numpy.float64)
        time=numpy.zeros([length],numpy.float64)
        window=numpy.zeros([length],numpy.float64)
        old_val=float(lines[0].split()[0]) ## start point
        start=float(lines[0].split()[0])-skip
        oldtime_Tem=start
        old_val=old_val-start
        print("start for rank %d is %4.4f\n" % (rank,start))
        old_val1=Tem_home[rank]
        actual_length=0
        for line in lines:
            if actual_length>=300000:
                break
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                val[0] = float(tokens[0])
                val[1] = float(tokens[1])
                
                if val[1] == Tem_home[rank,1]: ## back home
                    nv_Tem_home[rank]+=1
                    exc_len_Tem[rank]+=val[0] - oldtime_Tem
                    num_exc_Tem[rank]+=1
                    
                    if val[0] - oldtime_Tem > max_exc_Tem[rank]:
                        max_exc_Tem[rank]=val[0] - oldtime_Tem
                    oldtime_Tem=val[0]
                if val[1] == tems[0]:
                    nv_Tem_bot[rank]+=1
                if val[1] == tems[num_Tem-1]:
                   nv_Tem_top[rank]+=1

                
                spam=int((val[0] -start)/skip) 
                if spam<length:
                    if spam-old_val>1:                
                        for j in range(old_val,spam,1):
                            window[j]=old_val1
                    window[int(spam)]=val[1]
                    old_val=spam
                    old_val1=val[1]
                actual_length+=1
        #var_en=en.var()
        window_new=numpy.zeros([actual_length])
        window_new=window[:actual_length]-numpy.mean(window[:actual_length])
        norm_window_new=numpy.sum(window_new**2)

        result_window=numpy.correlate(window_new,window_new,mode='same')/norm_window_new
        print(result_window)
        print(window[:actual_length])
        print(result_window.size/2,actual_length)
        res_final_window =result_window[int(actual_length/2):]
        
        
        ### log bin our data
        #range is typically from 10^2 to 10^7 (3-7). make 200 bins
        dlogt=(8.0-2.0)/300.0
        mint=math.log(skip,10)
        nc=numpy.zeros(300,numpy.float64)
        c_pt=numpy.zeros(300,numpy.float64)
        
        for i in range(int(0.5*actual_length)):
            if i >0:
                log10t=math.log(i*skip,10)
                spam = int(log10t/dlogt)
            else:
                spam = 0
            print(spam,i,res_final_window[i])
            if spam<300:
                nc[spam]+=1.0
                c_pt[spam]+=res_final_window[i]
                

        f=open(sys.argv[1]+"/PT_Track_plots/"+sys.argv[3]+"/Tem/track_ac_logbin_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+"_"+str(rank)+".dat","w+") # put proper output file name with temp and fblamb here
        f.write("%8d %4.4f\n" % (0,1.0))
        for i in range(300):
            if nc[i]>0:
                f.write("%8d %4.4f\n" % (math.pow(10.0,((i+0.5)*dlogt)),c_pt[i]/nc[i]))
        f.close()



        x=numpy.arange(0,0.5*actual_length*skip,skip)
        f=open(sys.argv[1]+"/PT_Track_plots/"+sys.argv[3]+"/Tem/track_ac_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+"_"+str(rank)+".dat","w+") # put proper output file name with temp and fblamb here
        for i in range(int(0.5*actual_length)):
             f.write("%8d %4.4f\n" % (i*skip,res_final_window[i]))
        f.close()
        print("length is ",len(result_window[int(result_window.size/2):]))

        plt.plot(x,res_final_window,lw=2,label='PT track autocorr')


        plt.xscale('log')
        plt.legend()
        #plt.show(block=False)
        #time.sleep(5)
        #plt.pause(5)
        plt.savefig(sys.argv[1]+'/PT_Track_plots/'+sys.argv[3]+'/Tem/track_ac_n'+str(outlay[rank][0])+'_r'+str(outlay[rank][1])+'_t'+str(outlay[rank][2])+'_'+str(rank)+'.png')
        plt.close()

        
            
    f=open(sys.argv[1]+"/PT_stats/"+sys.argv[3]+"/Tem.dat","w+")
    for i in range(totSize):
        f.write("%8d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n" % (i,nSp_max_home[i,1],rho_home[i,1],Tem_home[i,1],\
                                                                     nv_Tem_bot[i],nv_Tem_top[i],nv_Tem_home[i],\
                                                                     exc_len_Tem[i]/(num_exc_Tem[i]+1),\
                                                                           max_exc_Tem[i]));

    f.close
    print("printed PT stats nSp_max")

    ######################################
####### Do for temp###################

nv_rho_home=numpy.zeros([totSize],numpy.float64)
nv_rho_top=numpy.zeros([totSize],numpy.float64)
nv_rho_bot=numpy.zeros([totSize],numpy.float64)
exc_len_rho=numpy.zeros([totSize],numpy.float64)
num_exc_rho=numpy.zeros([totSize],numpy.float64)
max_exc_rho=numpy.zeros([totSize],numpy.float64)

if num_rho>1 and plot_rho_ac==1:
    for rank in range(totSize):
        filename=sys.argv[1]+'/PT_Track_compile/rho/PT_Track_rho_rank'+str(rank)+'.dat'
        infile=open(filename,'r')
        lines=infile.readlines()
        infile.close()
        if lines[0] == '\n':
            continue
        spam=0
        n=0
        val=numpy.zeros([2],numpy.float64)
        time=numpy.zeros([length],numpy.float64)
        window=numpy.zeros([length],numpy.float64)
        old_val=float(lines[0].split()[0]) ## start point
        start=float(lines[0].split()[0])-skip
        oldtime_rho=start
        old_val=old_val-start
        print("start for rank %d is %4.4f\n" % (rank,start))
        old_val1=rho_home[rank]
        actual_length=0
        for line in lines:
            if actual_length>=300000:
                break
            if line[0] != '#' and line[0] != '@':
                tokens = line.split()
                val[0] = float(tokens[0])
                val[1] = float(tokens[1])
                
                if val[1] == rho_home[rank,1]: ## back home
                    nv_rho_home[rank]+=1
                    exc_len_rho[rank]+=val[0] - oldtime_rho
                    num_exc_rho[rank]+=1
                    
                    if val[0] - oldtime_rho > max_exc_rho[rank]:
                        max_exc_rho[rank]=val[0] - oldtime_rho
                    oldtime_rho=val[0]
                if val[1] == rho_0:
                    nv_rho_bot[rank]+=1
                if val[1] == rho_0+(num_rho-1)*d_rho:
                   nv_rho_top[rank]+=1

                
                spam=int((val[0] -start)/skip) 
                if spam<length:
                    if spam-old_val>1:                
                        for j in range(old_val,spam,1):
                            window[j]=old_val1
                    window[int(spam)]=val[1]
                    old_val=spam
                    old_val1=val[1]
                actual_length+=1
        #var_en=en.var()
        window_new=numpy.zeros([actual_length])
        window_new=window[:actual_length]-numpy.mean(window[:actual_length])
        norm_window_new=numpy.sum(window_new**2)

        result_window=numpy.correlate(window_new,window_new,mode='same')/norm_window_new
        print(result_window)
        print(window[:actual_length])
        print(result_window.size/2,actual_length)
        res_final_window =result_window[int(actual_length/2):]
        
        ### log bin our data
        #range is typically from 10^2 to 10^7 (3-7). make 200 bins
        dlogt=(8.0-2.0)/300.0
        mint=math.log(skip,10)
        nc=numpy.zeros(300,numpy.float64)
        c_pt=numpy.zeros(300,numpy.float64)
        
        for i in range(int(0.5*actual_length)):
            if i >0:
                log10t=math.log(i*skip,10)
                spam = int(log10t/dlogt)
            else:
                spam = 0
            if spam<300:
                nc[spam]+=1.0
                c_pt[spam]+=res_final_window[i]
                

        f=open(sys.argv[1]+"/PT_Track_plots/"+sys.argv[3]+"/rho/track_ac_logbin_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+"_"+str(rank)+".dat","w+") # put proper output file name with temp and fblamb here
        f.write("%8d %4.4f\n" % (0,1.0))
        for i in range(300):
            if nc[i]>0:
                f.write("%8d %4.4f\n" % (math.pow(10.0,((i+0.5)*dlogt)),c_pt[i]/nc[i]))
        f.close()


        x=numpy.arange(0,0.5*actual_length*skip,skip)
        f=open(sys.argv[1]+"/PT_Track_plots/"+sys.argv[3]+"/rho/track_ac_n"+str(outlay[rank][0])+"_r"+str(outlay[rank][1])+"_t"+str(outlay[rank][2])+"_"+str(rank)+".dat","w+") # put proper output file name with temp and fblamb here
        for i in range(int(0.5*actual_length)):
             f.write("%8d %4.4f\n" % (i*skip,res_final_window[i]))
        f.close()
        print("length is ",len(result_window[int(result_window.size/2):]))

        plt.plot(x,res_final_window,lw=2,label='PT track autocorr')


        plt.xscale('log')
        plt.legend()
        #plt.show(block=False)
        #time.sleep(5)
        #plt.pause(5)
        plt.savefig(sys.argv[1]+'/PT_Track_plots/'+sys.argv[3]+'/rho/track_ac_n'+str(outlay[rank][0])+'_r'+str(outlay[rank][1])+'_t'+str(outlay[rank][2])+'_'+str(rank)+'.png')
        plt.close()

        
            
    f=open(sys.argv[1]+"/PT_stats/"+sys.argv[3]+"/rho.dat","w+")
    for i in range(totSize):
        f.write("%8d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n" % (i,nSp_max_home[i,1],rho_home[i,1],Tem_home[i,1],\
                                                                      nv_rho_bot[i],nv_rho_top[i],nv_rho_home[i],\
                                                                      exc_len_rho[i]/(num_exc_rho[i]+1)\
                                                                           ,max_exc_rho[i]));

    f.close
    print("printed PT stats rho")

