#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gsl/gsl_math.h>
int main(int argc, char *argv[])
{
	// this program reads a histo2d file and rebins it.
    // first it reads the file - needs to know bounds and number of bins
    // then it groups it along the rho bins to rebin/ coarse-grain it
    

    // We calculate bDG(n) from the read P(n) and P(n,rho) using the shifted stitching employed for the 1D case.
    // The argument is as follows - from each sim, we get eq P(n) and P(n,rho) for n0L<=n<=n0U. We can change the normalisation constant
    // for each window independently, so long as it preserves the relative P within the window. 
    // we identify the shift so that we can compare results from adjacent windows iteratively to get the full free energy
    
	FILE *fp, *fpPara, *fp2, *fp3;
	char buffer[700]="";
	char buffer2[500], buffer3[500],buffer4[500],string[100];
	int i,j,k,l;
    double rhomin,rhomax,nmid;
    double pe,press,rho,q6,nSp,nSp_max,nl4,nl4_max,nl5,nl5_max,rho0;
    int mcstep;
    //double *rho_histo_unbias;
	
    int n0,bias,krho,skip,last;
    
    double accum_ebw,ebw;
    
    double temp;
    
    // WHAM vars
    
    double *totH,*rawH,*numPoints,*F,*oldF,*Pub,*grho,*gq6,*err;
    double **Hi,**rawHi,**embW;
    double *rPEtotH,*rPErawH,*rPEF,*rPEoldF,*rPEPub,*rPEerr;
    double **rPEHi,**rPErawHi,**rPEembW;
    int nWin,rnWin,nnWin,nAto;
    int nwid,rwid;
    double tol;
    double rho0L,rho0U,drho0,drhow;
    double n0L,n0U,dn0,dnw;
    int spamrho,nrhobin,nrholines,ncbin,spamn,spamindex;
    double spamP;

    int r=0;
    int rr=0;
    double spamgrho=0.0;
    double spamgq6=0.0;
    double crho=0.0;
    double cn=0.0;
	fpPara=fopen(argv[1],"r");
	
	if(fscanf(fpPara, "%s",string)>0); 				// highest x
	if(fscanf(fpPara, "%s",string)>0);
	rhomin=atof(string);
	if(fscanf(fpPara, "%s",string)>0); 				// lowest y
	if(fscanf(fpPara, "%s",string)>0);
	rhomax=atof(string);
	if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nrhobin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nrholines=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	ncbin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	n0=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	dn0=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	dnw=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	rho0=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	drho0=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	drhow=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	krho=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	temp=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	bias=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nWin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	rnWin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nnWin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	tol=atof(string);
    
    
    printf("temp %f\n",temp);
    double drho=1.0*(rhomax-rhomin)/nrhobin;
    
    fclose(fpPara);
    int lim=atoi(argv[3]);
    totH=calloc(nrhobin*ncbin,sizeof(double *));
    rawH=calloc(nrhobin*ncbin,sizeof(double *));
    err=calloc(nrhobin*ncbin,sizeof(double *));
    Pub=calloc(nrhobin*ncbin,sizeof(double *));
    numPoints=calloc(nWin,sizeof(double *));
    F=calloc(nWin,sizeof(double *));
    oldF=calloc(nWin,sizeof(double *));
    grho=calloc(nWin,sizeof(double *));
    gq6=calloc(nWin,sizeof(double *));
    
    
    Hi=calloc(nWin,sizeof(double *));
    rawHi=calloc(nWin,sizeof(double *));
    embW=calloc(nWin,sizeof(double *));

    for(i=0;i<nWin;i++)
    {
        Hi[i]=calloc(nrhobin*ncbin,sizeof(double *));
        rawHi[i]=calloc(nrhobin*ncbin,sizeof(double *));
        embW[i]=calloc(nrhobin*ncbin,sizeof(double *));
        
    }
   
    
    for(r=0;r<nWin;r++)
    {
        for(i=0;i<nrhobin*ncbin;i++)
            embW[r][i]=0.0;
    }
    
    sprintf(buffer,"%s",argv[2]);
    printf("overwritten buffer %s\n",buffer);
    if(((int)(temp*10000))%10 == 0)
        sprintf(buffer,"%s/AC_files/%0.3f/IntegAC.dat",argv[2],temp);
    else
        sprintf(buffer,"%s/AC_files/%0.4f/IntegAC.dat",argv[2],temp);
    
    fp=fopen(buffer,"r");
    if(fp==NULL)
    {
        printf("could not open IntegAC file %s\n",buffer);
    }
    for(rr=0;rr<nWin;rr++)
    {
        fscanf(fp,"%lf %lf %lf %lf\n",&rho,&nSp_max,&spamgrho,&spamgq6);
        grho[rr]=spamgrho;
        gq6[rr]=spamgq6;
    }
    fclose(fp);
    
    for(r=0;r<nWin;r++)
    {
        // first read and put in integrated autocorrelation time
        
        rwid=r%rnWin;
        nwid=(int)r/rnWin;
        cn=n0+nwid*dn0;
        n0L=cn-dnw;
        n0U=cn+dnw;
        
        crho=rho0+rwid*drho0;
        if(cn==2)
        {
            rho0L=crho-drhow;
            rho0U=crho+drhow;
        }
        else
        {
            rho0L=crho-drhow-0.002;
            rho0U=crho+drhow+0.002;
        }
        
        printf("%d wid %d r %d n %f %f %f rho %f %f %f\n",r,rwid,nwid,rho0L,crho,rho0U,n0L,cn,n0U);
        if(((int)(crho*1000))%10 == 0) // two decimal rho0
            if(((int)(temp*10000))%10 == 0)
                sprintf(buffer,"%s/histo2D2_compile/Conf_nSp%d-K0.01-rho%0.2fK%d_T%0.3f-bias%d-histo_2d2.dat",argv[2],(int)cn,crho,krho,temp,bias);
            else
                sprintf(buffer,"%s/histo2D2_compile/Conf_nSp%d-K0.01-rho%0.2fK%d_T%0.4f-bias%d-histo_2d2.dat",argv[2],(int)cn,crho,krho,temp,bias);
        else
            if(((int)(temp*10000))%10 == 0)
                sprintf(buffer,"%s/histo2D2_compile/Conf_nSp%d-K0.01-rho%0.3fK%d_T%0.3f-bias%d-histo_2d2.dat",argv[2],(int)cn,crho,krho,temp,bias);
            else
                sprintf(buffer,"%s/histo2D2_compile/Conf_nSp%d-K0.01-rho%0.3fK%d_T%0.4f-bias%d-histo_2d2.dat",argv[2],(int)cn,crho,krho,temp,bias);
            
        printf("%s\n",buffer);
        fp=fopen(buffer,"r");
        if(fp==NULL)
        {
            printf("could not open histo file %s\n",buffer);
        }
        
        accum_ebw=0.0;
        ebw=0.0;
        int old_mcstep=0;
        while(!(feof(fp)))
        {
            
            for(i=n0L;i<=n0U;i++)
            {
                for(j=0;j<nrholines;j++) // read all rho lines
                {
                    fscanf(fp,"%lf %lf %lf\n",&nSp_max,&rho,&spamP);
//                     spamn=(int)nSp_max;
                    if(rho>=rhomin && rho<rhomax)
                    {
                        spamrho=(int)((rho-rhomin)/drho);
//                         printf("%f %d rho spamrho %d %f\n",rho,spamrho,i,nSp_max);
                        spamindex=i*nrhobin+spamrho;
                    
                        if(rho0L<=rho && rho <=rho0U && spamrho<nrhobin && spamindex <nrhobin*ncbin)
                        {
                            if(spamP>=0.0)
                            {
                                Hi[r][spamindex]+=spamP;
                                totH[spamindex]+=spamP;
                                numPoints[r]+=spamP;
                                embW[r][spamindex]=1.0;
                            }
                            else
                            {
                                if(spamindex<nrhobin*ncbin)
                                    embW[r][spamindex]=0.0;
                            }
                        }
                    }
                }
            }
                    
        }
        if(fp!=NULL)
        {   
//             printf("tried to close %s\n",buffer);
            fclose(fp);
//             printf("managed to close %s\n",buffer);
        }
        else
        printf("%s is closed somehow\n",buffer);
        printf("done reading from  %s\n",buffer);
    }
    
    
    
    // initialise window FE/shift value and old window shift/ F value
    // wham iterate to convergence
    int wham=1;
    int iter=0;
    double conv=0.0;
    double oldconv=0.0;
    double denom=0.0;
    double tShift=0.0;
    double accum_Pub=0.0;
    for(i=0;i<nWin;i++)
    {
        F[i]=temp/(numPoints[r]);
//         for(r=0;r<nrhobin;r++)
//         {
//           embW[i][r]=1.0;   
//         }
    }
    while(wham==1)
    {
        // calc unbiased prob
        iter+=1;
//         printf("WHAM iter %d\n",iter);
        accum_Pub=0.0;
        for(r=0;r<nrhobin*ncbin;r++)
        {
            
        
//             printf("tot points %d %f\n",r,totH[r]);
            // calc denominator
            denom=0;
            err[r]=0.0;
            Pub[r]=0.0;
//             Pub[r] = totH[r]; // numerator
            for(i=0;i<nWin;i++)
            {
//                  denom+=numPoints[i]*embW[i][r]*exp(F[i]);
                if(grho[i]>0.0)
                {
                    Pub[r]+=Hi[i][r]/grho[i];
//                     printf("%d %d %f %f %f\n",i,r,rhomin+(r+0.5)*drho,Hi[i][r],grho[i]);
//                     printf("denom %d %f %f %f %f %f\n",i,numPoints[i],embW[i][r],temp,F[i],exp(F[i]/temp));
                    denom+=numPoints[i]*embW[i][r]*exp(F[i])/grho[i]; // no division by temp because we have been calculating beta F[i]
                    
                }
                else
                    printf("integ ac time for %d coming 0 or neg %f %f- something vei wong\n",i,grho[i],gq6[i]);
//                 err[r]+= rawHi[i][r];
                
                err[r]+=numPoints[i]*Hi[i][r]/grho[i];
                
            }
//             err[r]=totH[r];
//             printf("%d r %f %f\n",r,Pub[r],denom); 
            if(denom>0.0)
            Pub[r]/=denom;
            
            accum_Pub+=Pub[r];

            err[r]=sqrt(1.0/err[r]);
//             printf("unbiased prob %d %f %f %f\n",r,rhomin+drho*(r+0.5),Pub[r],err[r]);
        }
        
        //now find window FE/shift value
        tShift=0.0;
        for(i=0;i<nWin;i++)
        {
            F[i]=0.0;
            for(r=0;r<nrhobin*ncbin;r++)
            {
                if(isnan(fabs(Pub[r]))!=1 && isinf(fabs(embW[i][r]))!=1)
                {
                    F[i]+=Pub[r]*embW[i][r];
//                     printf("%d iter %d i %d n %f pubn %f embw %f F[i]\n",iter,i,r,Pub[r],embW[i][r],F[i]);
                }

            }
            
            if(F[i]>0.0)
                F[i]= -log(F[i]); // F[i] is actually beta F[i] 
            else
//                 printf("%f F[i] issue\n",F[i]);
//             printf("%d win %f and shift\n",i,F[i]);
            tShift+=F[i];
        }
        // check for convergence of F[i]
        
        conv=0.0;
        for(i=0;i<nWin;i++)
        {
//             printf("%d win %f %f\n",i,F[i],oldF[i]);
            conv+=fabs(F[i]-oldF[i]);
            
        }
        conv/=1.0*nWin;
//         printf("error is %d %f temp %f\n",iter,conv,temp);
//         if(conv<=tol || tShift==0.0 || fabs(conv-oldconv)<=tol)
        if(conv<=tol)
        {
            wham=0;
        }
        else
        {
            for(i=0;i<nWin;i++)
            {
                oldF[i]=F[i];
            }
            oldconv=conv;
        }
        if(iter==5000)
        wham=0;
    }
    printf("\nfinal shifts %f\n",temp);
    for(i=0;i<nWin;i++)
    {
        rwid=i%rnWin;
        nwid=(int)i/rnWin;
        crho=rho0+rwid*drho0;
        cn=n0+nwid*dn0;
        printf("%d %f %f %f\n",i,cn,crho,F[i]);
    }
    double logmin=1000.0;
    for(r=0;r<nrhobin*ncbin;r++)
    {
        if(-log(Pub[r]/accum_Pub) < logmin)
            logmin=-log(Pub[r]/accum_Pub);
    }
    
    sprintf(buffer3,"%s/whamPnrho-drw%0.3f-%0.4f-lim%d.dat",argv[2],drhow,temp,lim);
    fp2=fopen(buffer3,"w");
    double ranshift=0.0;
    if(temp==0.0395)
        ranshift=1.0;
    if(temp==0.0391)
        ranshift=2.0;
    if(temp==0.0389)
        ranshift=2.8;
    for(r=0;r<ncbin;r++)
    {
        for(rr=0;rr<nrhobin;rr++)
        {
            spamindex=r*nrhobin+rr;
//             printf("%d %d %f %f\n",rr,spamindex,rhomin+(rr+0.5)*drho,drho);
            if(r<5)
            {
            if(Pub[spamindex]>0.0)
            {
                if(isinf(err[spamindex])!=1)
                    fprintf(fp2,"%f %f %f %f %f\n",1.0*r,5.0571*(rhomin+(rr+0.5)*drho),-log(Pub[spamindex]/accum_Pub)-logmin,err[spamindex],Pub[spamindex]);
                else
                    fprintf(fp2,"%f %f %f %f %f\n",1.0*r,5.0571*(rhomin+(rr+0.5)*drho),-log(Pub[spamindex]/accum_Pub)-logmin,0.0,Pub[spamindex]);
            }
            else
            {
                fprintf(fp2,"%f %f %f %f %f\n",1.0*r,5.0571*(rhomin+(rr+0.5)*drho),9999.000000,0.0,0.0);
            }
            }
            if(r>=5)
            {
            if(Pub[spamindex]>0.0)
            {
                if(isinf(err[spamindex])!=1)
                    fprintf(fp2,"%f %f %f %f %f\n",1.0*r,5.0571*(rhomin+(rr+0.5)*drho),-log(Pub[spamindex]/accum_Pub)-logmin-ranshift-fabs((r-9)/ranshift),err[spamindex],Pub[spamindex]);
                else
                    fprintf(fp2,"%f %f %f %f %f\n",1.0*r,5.0571*(rhomin+(rr+0.5)*drho),-log(Pub[spamindex]/accum_Pub)-logmin-ranshift-fabs((r-9)/ranshift),0.0,Pub[spamindex]);
            }
            else
            {
                fprintf(fp2,"%f %f %f %f %f\n",1.0*r,5.0571*(rhomin+(rr+0.5)*drho),9999.000000,0.0,0.0);
            }    
            }
        }
        fprintf(fp2,"\n");
    }
	fclose(fp2);
     
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
