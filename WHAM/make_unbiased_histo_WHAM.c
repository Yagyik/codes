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
    int nWin,nAto;
    double tol;
    double rho0L,rho0U,drho0,drhow;
    int spamrho,nrhobin;
    double PEmin,PEmax;
    int spamPE,nPEbin,spamindex;
    
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
	n0=atoi(string);
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
	krho=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest x
	if(fscanf(fpPara, "%s",string)>0);
	PEmin=atof(string);
	if(fscanf(fpPara, "%s",string)>0); 				// lowest y
	if(fscanf(fpPara, "%s",string)>0);
	PEmax=atof(string);
	if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nPEbin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nAto=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	temp=atof(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	bias=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	skip=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	last=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	nWin=atoi(string);
    if(fscanf(fpPara, "%s",string)>0); 				// highest y
	if(fscanf(fpPara, "%s",string)>0);
	tol=atof(string);
    
    
    printf("temp %f\n",temp);
    double drho=1.0*(rhomax-rhomin)/nrhobin;
    double dPE=(PEmax-PEmin)/(1.0*nPEbin); //
    
    fclose(fpPara);
    int lim=atoi(argv[3]);
    totH=calloc(nrhobin,sizeof(double *));
    rawH=calloc(nrhobin,sizeof(double *));
    err=calloc(nrhobin,sizeof(double *));
    Pub=calloc(nrhobin,sizeof(double *));
    rPEtotH=calloc(nrhobin*nPEbin,sizeof(double *));
    rPErawH=calloc(nrhobin*nPEbin,sizeof(double *));
    rPEerr=calloc(nrhobin*nPEbin,sizeof(double *));
    rPEPub=calloc(nrhobin*nPEbin,sizeof(double *));
    numPoints=calloc(nWin,sizeof(double *));
    F=calloc(nWin,sizeof(double *));
    oldF=calloc(nWin,sizeof(double *));
    rPEF=calloc(nWin,sizeof(double *));
    rPEoldF=calloc(nWin,sizeof(double *));
    grho=calloc(nWin,sizeof(double *));
    gq6=calloc(nWin,sizeof(double *));
    
    
    Hi=calloc(nWin,sizeof(double *));
    rawHi=calloc(nWin,sizeof(double *));
    embW=calloc(nWin,sizeof(double *));
    rPEHi=calloc(nWin,sizeof(double *));
    rPErawHi=calloc(nWin,sizeof(double *));
    rPEembW=calloc(nWin,sizeof(double *));
    for(i=0;i<nWin;i++)
    {
        Hi[i]=calloc(nrhobin,sizeof(double *));
        rawHi[i]=calloc(nrhobin,sizeof(double *));
        embW[i]=calloc(nrhobin,sizeof(double *));
        rPEHi[i]=calloc(nrhobin*nPEbin,sizeof(double *));
        rPErawHi[i]=calloc(nrhobin*nPEbin,sizeof(double *));
        rPEembW[i]=calloc(nrhobin*nPEbin,sizeof(double *));
        
    }
   
    int r=0;
    int rr=0;
    double spamgrho=0.0;
    double spamgq6=0.0;
    double crho=0.0;
     for(r=0;r<nWin;r++)
    {
        for(i=0;i<nrhobin;i++)
            embW[r][i]=0.0;
    }
    for(r=0;r<nWin;r++)
    {
        // first read and put in integrated autocorrelation time
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
        crho=rho0+r*drho0;
        rho0L=crho-drhow;
        rho0U=crho+drhow;
        if(((int)(crho*1000))%10 == 0) // two decimal rho0
            if(((int)(temp*10000))%10 == 0)
                sprintf(buffer,"%s/thermo_compile/Conf_nSp%d-K0.01-rho%0.2fK%d_T%0.3f-bias%d-thermo.dat",argv[2],n0,crho,krho,temp,bias);
            else
                sprintf(buffer,"%s/thermo_compile/Conf_nSp%d-K0.01-rho%0.2fK%d_T%0.4f-bias%d-thermo.dat",argv[2],n0,crho,krho,temp,bias);
        else
            if(((int)(temp*10000))%10 == 0)
                sprintf(buffer,"%s/thermo_compile/Conf_nSp%d-K0.01-rho%0.3fK%d_T%0.3f-bias%d-thermo.dat",argv[2],n0,crho,krho,temp,bias);
            else
                sprintf(buffer,"%s/thermo_compile/Conf_nSp%d-K0.01-rho%0.3fK%d_T%0.4f-bias%d-thermo.dat",argv[2],n0,crho,krho,temp,bias);
            
        printf("%s\n",buffer);
        fp=fopen(buffer,"r");
        if(fp==NULL)
        {
            printf("could not open thermo file %s\n",buffer);
        }
        
        accum_ebw=0.0;
        ebw=0.0;
        int old_mcstep=0;
        while(!(feof(fp)))
        {
            
            fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&mcstep,&pe,&press,&rho,&q6,&nSp,&nSp_max,&nl4,&nl4_max,&nl5,&nl5_max,&ebw);
            if(mcstep==old_mcstep || mcstep>=last)
                break;
            spamrho=(int)((rho-rhomin)/drho);
            spamPE=(int)((pe/(1.0*nAto)-PEmin)/dPE);
            spamindex=spamrho*nPEbin+spamPE;
            if(mcstep>skip && nSp_max <=lim && mcstep <=last)
            {
//                 printf("%d %f %f %f %f\n",mcstep,rho,rho0,nSp_max,ebw);
//                 printf("%d %f %f %f ind and rho bounds\n",spamrho,rho,rho0L,rho0U);
                if(rho0L<=rho && rho <=rho0U && spamrho<nrhobin && spamPE >=0 && spamPE<nPEbin && spamindex <nrhobin*nPEbin)
                {
//                     printf("updating histo %d win rho %f L %f U %f spam %d ebw %f\n",r,rho,rho0L,rho0U,spamrho,ebw);
                    if(ebw>=0.0)
                    {
                        Hi[r][spamrho]+=ebw;
                        rawHi[r][spamrho]+=1.0;
                        totH[spamrho]+=ebw;
                        rawH[spamrho]+=1.0;
                        numPoints[r]+=ebw;
                        embW[r][spamrho]=1.0;
                        
                        rPEHi[r][spamindex]+=ebw;
                        rPErawHi[r][spamindex]+=1.0;
                        rPEtotH[spamindex]+=ebw;
                        rPErawH[spamindex]+=1.0;
                        rPEembW[r][spamindex]=1.0;
                    }
                }
                else
                {   
                    if(spamrho<nrhobin)
                    embW[r][spamrho]=0.0;
                    if(spamindex<nrhobin*nPEbin)
                    rPEembW[r][spamindex]=0.0;
                }
                accum_ebw+=ebw;
            }
            
            old_mcstep=mcstep;
        }
        if(fp!=NULL)
        {   
            printf("tried to close %s\n",buffer);
            fclose(fp);
            printf("managed to close %s\n",buffer);
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
    double accum_rPEPub=0.0;
    for(i=0;i<nWin;i++)
    {
        F[i]=temp/(1.0*i+1);
//         for(r=0;r<nrhobin;r++)
//         {
//           embW[i][r]=1.0;   
//         }
    }
    while(wham==1)
    {
        // calc unbiased prob
        iter+=1;
        printf("WHAM iter %d\n",iter);
        for(r=0;r<nrhobin;r++)
        {
            Pub[r] = totH[r]; // numerator
        
//             printf("tot points %d %f %f\n",r,rhomin+(r+0.5)*drho,totH[r]);
            // calc denominator
            denom=0;
            err[r]=0.0;
//             Pub[r]=0.0;
            for(i=0;i<nWin;i++)
            {
//                 if(grho[i]>0.0)
//                 {
//                     Pub[r]+=Hi[i][r]/grho[i];
// //                     printf("%d %d %f %f %f\n",i,r,rhomin+(r+0.5)*drho,Hi[i][r],grho[i]);
// //                     printf("denom %d %f %f %f %f %f\n",i,numPoints[i],embW[i][r],temp,F[i],exp(F[i]/temp));
//                     denom+=numPoints[i]*embW[i][r]*exp(F[i])/grho[i]; // no division by temp because we have been calculating beta F[i]
//                 }
//                 else
//                     printf("integ ac time for %d coming 0 or neg %f %f- something vei wong\n",i,grho[i],gq6[i]);
                denom+=numPoints[i]*embW[i][r]*exp(F[i]); // no division by temp because we have been calculating   
//                 err[r]+= rawHi[i][r];
                err[r]+=numPoints[i]*Hi[i][r]/grho[i];
                
            }
//             printf("%d %f r rho %f %f err %f\n",r,rhomin+(r+0.5)*drho,Pub[r],denom,err[r]); 
            if(denom>0.0)
            Pub[r]/=denom;

            err[r]=sqrt(1.0/err[r]);
//             printf("unbiased prob %d %f %f %f\n",r,rhomin+drho*(r+0.5),Pub[r],err[r]);
        }
        
        //now find window FE/shift value
        tShift=0.0;
        for(i=0;i<nWin;i++)
        {
            F[i]=0.0;
            for(r=0;r<nrhobin;r++)
            {
                if(isnan(fabs(Pub[r]))!=1 && isinf(fabs(embW[i][r]))!=1)
                {
                    F[i]+=Pub[r]*embW[i][r];
//                     printf("%d iter %d i %d n %f pubn %f embw %f F[i]\n",iter,i,r,Pub[r],embW[i][r],F[i]);
                }
//                 
//                 if(i==nWin-1)
//                     printf("%d r %f %f P %f emb %f F \n",r,rhomin+(r+0.5)*drho,Pub[r],embW[i][r],F[i]);
            }
            
            
//             n0L=i+1-1;
//             n0U=i+1+1;
//             for(n=n0L;n<=n0U;n++)
//             {
//                 if(Pub[n]>0.0)
//                 F[i]+=Pub[n]*embW[i][n];
//                 
//             }
            if(F[i]>0.0)
                F[i]= -log(F[i]); // F[i] is actually beta F[i] 
            else
                printf("%f F[i] issue\n",F[i]);
//             printf("%d win %f and shift\n",i,F[i]);
            tShift+=F[i];
        }
        // check for convergence of F[i]
        
        conv=0.0;
        for(i=0;i<nWin;i++)
        {
            printf("%d win %f %f\n",i,F[i],oldF[i]);
            conv+=fabs(F[i]-oldF[i]);
            
        }
        conv/=1.0*nWin;
        printf("error is %f temp %f\n",conv,temp);
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
        if(iter==500)
        wham=0;
    }
	
	
	// found the unbiased prob - print it
	sprintf(buffer3,"%s/whamPrho-drw%0.3f-%0.4f-lim%d.dat",argv[2],drhow,temp,lim);
    fp2=fopen(buffer3,"w");
    for(r=0;r<nrhobin;r++)
    {
        if(Pub[r]>0.0)
            if(isinf(err[r])!=1)
                fprintf(fp2,"%f %f %f %f\n",5.0571*(rhomin+(r+0.5)*drho),-log(Pub[r]),err[r],Pub[r]);
            else
                fprintf(fp2,"%f %f %f %f\n",5.0571*(rhomin+(r+0.5)*drho),-log(Pub[r]),0.0,Pub[r]);
    }
	fclose(fp2);
    
    
    wham=1;
    iter=0;
    conv=0.0;
    oldconv=0.0;
    denom=0.0;
    tShift=0.0;
    accum_rPEPub=0.0;
    for(i=0;i<nWin;i++)
    {
        rPEF[i]=temp/(1.0*i+1);
//         for(r=0;r<nrhobin;r++)
//         {
//           embW[i][r]=1.0;   
//         }
    }
    while(wham==1)
    {
        // calc unbiased prob
        iter+=1;
        printf("WHAM iter %d\n",iter);
        accum_rPEPub=0.0;
        for(r=0;r<nrhobin*nPEbin;r++)
        {
//             Pub[r] = totH[r]; // numerator
        
//             printf("tot points %d %f %f\n",r,rhomin+(r+0.5)*drho,totH[r]);
            // calc denominator
            denom=0;
            rPEerr[r]=0.0;
            rPEPub[r]=0.0;
            rPEPub[r]=rPEtotH[r];
            for(i=0;i<nWin;i++)
            {/*
                if(grho[i]>0.0)
                {
                    rPEPub[r]+=rPEHi[i][r]/grho[i];
//                     printf("%d %d %f %f %f\n",i,r,rhomin+(r+0.5)*drho,Hi[i][r],grho[i]);
//                     printf("denom %d %f %f %f %f %f\n",i,numPoints[i],embW[i][r],temp,F[i],exp(F[i]/temp));
                    denom+=numPoints[i]*rPEembW[i][r]*exp(rPEF[i])/grho[i]; // no division by temp because we have been calculating beta F[i]
                }
                else
                    printf("integ ac time for %d coming 0 or neg %f %f- something vei wong\n",i,grho[i],gq6[i]);*/
//                 err[r]+= rawHi[i][r];
                denom+=numPoints[i]*rPEembW[i][r]*exp(rPEF[i]); // no division by temp because we have been calculating beta F[i]
                rPEerr[r]+=numPoints[i]*rPEHi[i][r]/grho[i];
                
            }
//             printf("%d %f r rho %f %f err %f\n",r,rhomin+(r+0.5)*drho,Pub[r],denom,err[r]); 
            if(denom>0.0)
            rPEPub[r]/=denom;
            
            accum_rPEPub+=rPEPub[r];

            rPEerr[r]=sqrt(1.0/err[r]);
//             printf("unbiased prob %d %f %f %f\n",r,rhomin+drho*(r+0.5),Pub[r],err[r]);
        }
        
        //now find window FE/shift value
        tShift=0.0;
        for(i=0;i<nWin;i++)
        {
            rPEF[i]=0.0;
            for(r=0;r<nrhobin*nPEbin;r++)
            {
                if(isnan(fabs(rPEPub[r]))!=1 && isinf(fabs(rPEembW[i][r]))!=1)
                {
                    rPEF[i]+=rPEPub[r]*rPEembW[i][r];
//                     printf("%d iter %d i %d n %f pubn %f embw %f F[i]\n",iter,i,r,Pub[r],embW[i][r],F[i]);
                }

            }
            
            if(rPEF[i]>0.0)
                rPEF[i]= -log(rPEF[i]); // F[i] is actually beta F[i] 
            else
                printf("%f F[i] issue\n",rPEF[i]);
//             printf("%d win %f and shift\n",i,F[i]);
            tShift+=rPEF[i];
        }
        // check for convergence of F[i]
        
        conv=0.0;
        for(i=0;i<nWin;i++)
        {
            printf("%d win %f %f\n",i,rPEF[i],rPEoldF[i]);
            conv+=fabs(rPEF[i]-rPEoldF[i]);
            
        }
        conv/=1.0*nWin;
        printf("error is %f temp %f\n",conv,temp);
//         if(conv<=tol || tShift==0.0 || fabs(conv-oldconv)<=tol)
        if(conv<=tol)
        {
            wham=0;
        }
        else
        {
            for(i=0;i<nWin;i++)
            {
                rPEoldF[i]=rPEF[i];
            }
            oldconv=conv;
        }
        if(iter==500)
        wham=0;
    }
    double logmin=1000.0;
    for(r=0;r<nrhobin*nPEbin;r++)
    {
        if(-log(rPEPub[r]/accum_rPEPub) < logmin)
            logmin=-log(rPEPub[r]/accum_rPEPub);
    }
    
    sprintf(buffer3,"%s/whamPrhoPE-drw%0.3f-%0.4f-lim%d.dat",argv[2],drhow,temp,lim);
    fp2=fopen(buffer3,"w");
    for(r=0;r<nrhobin;r++)
    {
        for(rr=0;rr<nPEbin;rr++)
        {
            spamindex=r*nPEbin+rr;
            if(rPEPub[spamindex]>0.0)
            {
                if(isinf(err[spamindex])!=1)
                    fprintf(fp2,"%f %f %f %f %f\n",5.0571*(rhomin+(r+0.5)*drho),PEmin+(rr+0.5)*dPE,-log(rPEPub[spamindex]/accum_rPEPub)-logmin+0.6,rPEerr[spamindex],rPEPub[spamindex]);
                else
                    fprintf(fp2,"%f %f %f %f %f\n",5.0571*(rhomin+(r+0.5)*drho),PEmin+(rr+0.5)*dPE,-log(rPEPub[spamindex]/accum_rPEPub)-logmin+0.6,0.0,rPEPub[spamindex]);
            }
            else
            {
                fprintf(fp2,"%f %f %f %f %f\n",5.0571*(rhomin+(r+0.5)*drho),PEmin+(rr+0.5)*dPE,9999.000000,0.0,0.0);
            }
        }
        fprintf(fp2,"\n");
    }
	fclose(fp2);
     
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
