//#include"YgCoordNo.h"
struct CLSIZE
{
	int counter;          
	int fileCounter;
	int *oldclustSeb,**nucleo,*clustSeb,*csizeSeb;
	int **C,*neigh;
	double **q3mR,**q3mI,*r,*q3,**q6mR,**q6mI,*q6;
    double *qqdistro3,*qqdistro4,*qqdistro5,*qqdistro6, *q3NCDistro1, *q3NCDistro2, *qqNeig;
    double *qq6distro3,*qq6distro4,*qq6distro5,*qq6distro6,*q6NCDistro1, *q6NCDistro2;
	int **ClTime;	// Max cluster size e numero di solid particles,
	double **coordHisto, **CryNeigh, **bondHisto;	
	//	and histogram of coord nos. double because we are normalising histograms, also crystal specific g(r), histogram of bonds
	int **coordList;
	int maxCoord;
	int numfield; // complicated shit for normalising a la g(r)
	
} clsize;

void Init_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	int i,j,spamindex,temp;
	char buffer[500];
	cl->numfield=1;
	cl->maxCoord = 70;
	cl->q3mR=AllocMatR(s->Ato,7);
	cl->q3mI=AllocMatR(s->Ato,7);
	cl->q6mR=AllocMatR(s->Ato,13);
	cl->q6mI=AllocMatR(s->Ato,13);
	cl->r=AllocVecR(4);
	cl->C=AllocMatI(200*s->Ato,3);
	cl->oldclustSeb=AllocVecI(s->Ato+1);
	cl->nucleo=AllocMatI(s->Ato+1,3);		// Solid + liquid5 + liquid4
	cl->clustSeb=AllocVecI(s->Ato+1);
	cl->csizeSeb=AllocVecI(s->Ato+1);
	cl->q3=AllocVecR(s->Ato);
	cl->q6=AllocVecR(s->Ato);
	cl->neigh=AllocVecI(s->Ato);
	printf("cluster size, all allocated??\n");
	cl->ClTime=AllocMatI(s->Cl_limit+1,6);	// S + l5 + l4 + s_Max + l5_Max + l4_Max + nBond
	cl->coordHisto=AllocMatR(3,cl->maxCoord);
	cl->coordList=AllocMatI(s->Ato,cl->maxCoord); // list of particles, no. neigh and nerigh ids.
	cl->CryNeigh=AllocMatR(4000,3);
	cl->bondHisto=AllocMatR(10,3);
    cl->qqdistro3=AllocVecR(100);
    cl->qqdistro4=AllocVecR(100);
    cl->qqdistro5=AllocVecR(100);
    cl->qqdistro6=AllocVecR(100);
    cl->qq6distro3=AllocVecR(100);
    cl->qq6distro4=AllocVecR(100);
    cl->qq6distro5=AllocVecR(100);
    cl->qq6distro6=AllocVecR(100);
    int spam = (int)(0.33*s->OrPar_bin);
    cl->q3NCDistro1=AllocVecR(spam);
    cl->q3NCDistro2=AllocVecR(spam);
    cl->q6NCDistro1=AllocVecR(spam);
    cl->q6NCDistro2=AllocVecR(spam);
    cl->qqNeig=AllocVecR(10);
	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<7;j++)
		{
			cl->q3mR[i][j]=0.0;
			cl->q3mI[i][j]=0.0;
		}
		cl->q3[i]=0.0;
		for(j=0;j<13;j++)
		{
			cl->q6mR[i][j]=0.0;
			cl->q6mI[i][j]=0.0;
		}
		cl->q6[i]=0.0;
	}
	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<cl->maxCoord;j++)
		{
		cl->coordList[i][j]=0;
		if(i<3)
		{
		cl->coordHisto[i][j]=0.0;
		}
		}

	}
	for(i=0;i<4000;i++)
	{
		for(j=0;j<3;j++)
		{
		cl->CryNeigh[i][j]=0.0;
		}
	}
	for(i=0;i<10;i++)
	{
		for(j=0;j<3;j++)
		{
		cl->bondHisto[i][j]=0.0;
		}
	}
	for(i=0;i<100;i++)
    {
     cl->qqdistro3[i]=0.0;
     cl->qqdistro4[i]=0.0;
     cl->qqdistro5[i]=0.0;
     cl->qqdistro6[i]=0.0;
     cl->qq6distro3[i]=0.0;
     cl->qq6distro4[i]=0.0;
     cl->qq6distro5[i]=0.0;
     cl->qq6distro6[i]=0.0;
    }
    for(i=0;i<spam;i++)
    {
        cl->q3NCDistro1[i]=0.0;
        cl->q3NCDistro2[i]=0.0;
        cl->q6NCDistro1[i]=0.0;
        cl->q6NCDistro2[i]=0.0;
    }
    for(i=0;i<10;i++)
    {
        cl->qqNeig[i]=0.0;
    }
	cl->counter=0;
	cl->fileCounter=0;
}

void print_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	FILE *fp;	
	int t;
	char buffer[400];
	double deltaR, norma1,norma2,norma3;
	
	printf("entering print! cl size\n");
	sprintf(buffer,"%sClSize%g-%d.dat",s->output,s->Cl_cut,cl->fileCounter);
  	fp=fopen(buffer,"w");
  	printf("%s\n",buffer);
  	if(fp==NULL)
  	printf("null pointer problems cluster\n");
	for(t=0;t<s->Cl_limit;t++)
	{
		fprintf(fp,"%le	%d	%d	%d	%d	%d	%d\n",s->DeltaT*s->Step*t,cl->ClTime[t][0],cl->ClTime[t][1],cl->ClTime[t][2],cl->ClTime[t][3],cl->ClTime[t][4],cl->ClTime[t][5]);
	}
	fclose(fp);
	printf("done with cluster\n");
	
	// section to make a file and print the coord no histo 
	sprintf(buffer,"%sCoordNoClSizeSegg%g-%d.dat",s->output,s->Cl_r,cl->fileCounter);
	fp=fopen(buffer,"w");
	// code below prints a separate histogram of coordination numbers for LDL, HDL and crystal
	for(t=0;t<cl->maxCoord;t++) // for each coordination number
	{	
		cl->coordHisto[0][t] /= s->Cl_limit; // normalise against
		cl->coordHisto[1][t] /= s->Cl_limit; // normalise against
		cl->coordHisto[2][t] /= s->Cl_limit; // normalise against
		fprintf(fp,"%d %f %f %f\n",t,cl->coordHisto[0][t],cl->coordHisto[1][t],cl->coordHisto[2][t]);
	}
	fclose(fp);

	
	deltaR = 0.5*d->lBox[0]/cl->numfield;
	norma2=2.0*pow(cl->numfield,3.0)/(M_PI*cl->nucleo[0][1]*cl->nucleo[0][1]*(s->Cl_limit-1));
	norma3=2.0*pow(cl->numfield,3.0)/(M_PI*cl->nucleo[0][1]*cl->nucleo[0][2]*(s->Cl_limit-1));
	norma1=2.0*pow(cl->numfield,3.0)/(M_PI*cl->nucleo[0][1]*cl->nucleo[0][0]*(s->Cl_limit-1));
	//printf("zeros?? %f %f %f\n",norma1,norma2,norma3);
	//printf("%d %d\n",cl->nucleo[0][2],cl->numfield);
	//norma1 = 0.0;
	//norma2 = 0.0;
	sprintf(buffer,"%sLDL_Gr_%d.dat",s->output,cl->fileCounter);
	fp=fopen(buffer,"w");		
	for(t=0;t<cl->numfield;t++)
	{
		//norma1 = (pow((t+1),3)-pow(t,3))*pow(deltaR,3);
		//norma2 = (4/3)*M_PI*norma1*(s->Ato/(pow(d->lBox[0],3)));
		//cl->CryNeigh[t][0]/=s->Cl_limit*cl->nucleo[0][1]*norma2;
		//cl->CryNeigh[t][1]/=s->Cl_limit*cl->nucleo[0][2]*norma2;
		//printf("%f %f %f - printing??\n",cl->CryNeigh[t][1],cl->CryNeigh[t][2],cl->CryNeigh[t][0]);
		cl->CryNeigh[t][1]*=norma2/sqr(t+1);
		cl->CryNeigh[t][2]*=norma3/sqr(t+1);
		cl->CryNeigh[t][0]*=norma1/sqr(t+1);
		//printf("%f %f %f - what are we printing??\n",cl->CryNeigh[t][1],cl->CryNeigh[t][2],cl->CryNeigh[t][0]);
		//cl->CryNeigh[t][0]/=s->Cl_limit*cl->nucleo[0][1];
		//cl->CryNeigh[t][1]/=s->Cl_limit*cl->nucleo[0][2];
		fprintf(fp,"%le	%le %le %le\n",((double)t+0.5)*deltaR,cl->CryNeigh[t][0],cl->CryNeigh[t][1],cl->CryNeigh[t][2]);
	}
	fclose(fp);
    
    // print q3q3 histo 
    sprintf(buffer,"%sq3q3Histo%g.dat",s->output,s->Cl_r);
    fp=fopen(buffer,"w");
    double spamsum=0.0;
    for(t=0;t<100;t++)
    {
        spamsum=cl->qqdistro3[t]/s->Cl_limit + cl->qqdistro4[t]/s->Cl_limit + cl->qqdistro5[t]/s->Cl_limit + cl->qqdistro6[t]/s->Cl_limit;
        fprintf(fp,"%f %f %f %f %f %f\n",-1.0+2.0*t/100.0,cl->qqdistro3[t]/s->Cl_limit,cl->qqdistro4[t]/s->Cl_limit,cl->qqdistro5[t]/s->Cl_limit,cl->qqdistro6[t]/s->Cl_limit,spamsum);
    }
    fclose(fp);
    
    sprintf(buffer,"%sq6q6Histo%g.dat",s->output,s->Cl_r);
    fp=fopen(buffer,"w");
    spamsum=0.0;
    for(t=0;t<100;t++)
    {
        spamsum=cl->qq6distro3[t]/s->Cl_limit + cl->qq6distro4[t]/s->Cl_limit + cl->qq6distro5[t]/s->Cl_limit + cl->qq6distro6[t]/s->Cl_limit;
        fprintf(fp,"%f %f %f %f %f %f\n",-1.0+2.0*t/100.0,cl->qq6distro3[t]/s->Cl_limit,cl->qq6distro4[t]/s->Cl_limit,cl->qq6distro5[t]/s->Cl_limit,cl->qq6distro6[t]/s->Cl_limit,spamsum);
    }
    fclose(fp);
    
    sprintf(buffer,"%sq3NCDistro%g.dat",s->output,s->Cl_r);
    fp=fopen(buffer,"w");
    for(t=0;t<0.33*s->OrPar_bin;t++)
    {
        fprintf(fp,"%f %f %f\n",3.0*t/s->OrPar_bin,cl->q3NCDistro1[t]/s->Cl_limit,cl->q3NCDistro2[t]/s->Cl_limit); // q3 distributions of two kinds of particles
    }
    fclose(fp);
    
    sprintf(buffer,"%sq6NCDistro%g.dat",s->output,s->Cl_r);
    fp=fopen(buffer,"w");
    for(t=0;t<0.33*s->OrPar_bin;t++)
    {
        fprintf(fp,"%f %f %f\n",3.0*t/s->OrPar_bin,cl->q6NCDistro1[t]/s->Cl_limit,cl->q6NCDistro2[t]/s->Cl_limit); // q3 distributions of two kinds of particles
    }
    fclose(fp);
    
    sprintf(buffer,"%sqqNeig%g.dat",s->output,s->Cl_r);
    fp=fopen(buffer,"w");
    for(t=0;t<10;t++)
    {
        fprintf(fp,"%d %f\n",t,cl->qqNeig[t]/s->Cl_limit); // number of "bonded" neighbours 
    }
    fclose(fp);
	cl->fileCounter++;
	printf("done coord stats, exiting\n");
}

void Compute_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	int i,k,nBond,j,m,spamindex,temp;
	int sp1,spSolid,spLiq5,spLiq4,ipair;
	double q3q3R,q3q3I;
	double q6q6R,q6q6I;
	int change, cnum, cref, maxcsizeSeb, maxclustSeb;
	double cosTheta,phi;

	int NeiCount;
	double posR[4];
	
	//char buffer[500],buffer2[500],buffer3[500];	
	//FILE *fp, *fp2, *fp3;
	printf("entering cl size\n");
	
	
	
	// we need to open files to write the list of LDL, HDL and crystal particles and their positions
	// need to pay attention to appropriate indexing.
	// need to also make sure it's in format compatible with ovito
	// total 3 files.
	
	
	// we also need to open files to write the list of bonds between LDL, between HDl and between crystal particles.
	// need to pay mind to the fact that for each bond, we want to store the indices corresponding to the nucleo list, rather than in the actual full list.
	// total 3 files. 
	// NOTE: not needed if Ovito is used
	for(i=0;i<s->Ato;i++)
	{	
	//printf("sanity 1 cluster size %d\n",i);
		nBond=0;
		for(j=0;j<s->Ato;j++)
		{
			
			if(i!=j)
			{
				for(k=0;k<3;k++)
				{
					cl->r[k]=d->r[j][k]-d->r[i][k];
					cl->r[k]-=d->lBox[k]*rint(cl->r[k]/d->lBox[k]);
				}
				cl->r[3]=sqrt(sqr(cl->r[0])+sqr(cl->r[1])+sqr(cl->r[2]));
				if(cl->r[3]<s->Cl_r)			// Atomo dentro la prima shell
				{
					nBond++;
				
					cosTheta=cl->r[2]/cl->r[3];
					phi=atan2(cl->r[1],cl->r[0]);
					
					for(m=0;m<=3;m++)
					{
						cl->q3mR[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(m*phi);
						cl->q3mI[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(m*phi);
					}
					for(m=1;m<=3;m++)
					{
						cl->q3mR[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(-m*phi)*pow(-1,m);
						cl->q3mI[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(-m*phi)*pow(-1,m);
					}
					
					for(m=0;m<=6;m++)
					{
						cl->q6mR[i][m]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*cos(m*phi);
						cl->q6mI[i][m]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*sin(m*phi);
					}
					for(m=1;m<=6;m++)
					{
						cl->q6mR[i][m+6]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*cos(-m*phi)*pow(-1,m);
						cl->q6mI[i][m+6]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*sin(-m*phi)*pow(-1,m);
					}
				}
			}
		}
			
		for(m=0;m<7;m++)
		{
			cl->q3mR[i][m]/=nBond;
			cl->q3mI[i][m]/=nBond;
			cl->q3[i]+=cl->q3mR[i][m]*cl->q3mR[i][m]+cl->q3mI[i][m]*cl->q3mI[i][m];
		}
		cl->q3[i]=sqrt(cl->q3[i]*4.0*M_PI/7.0);
		for(m=0;m<13;m++)
		{
			cl->q6mR[i][m]/=nBond;
			cl->q6mI[i][m]/=nBond;
			cl->q6[i]+=cl->q6mR[i][m]*cl->q6mR[i][m]+cl->q6mI[i][m]*cl->q6mI[i][m];
		}
		cl->q6[i]=sqrt(cl->q6[i]*4.0*M_PI/13.0);
		//cl->q6[i]=sqrt(cl->q6[i]);
		//printf("%d %f %f\n",i,cl->q3[i],cl->q6[i]);
		cl->neigh[i]=nBond;								// Adesso ho il q3(i) ed il numero di bond
	}
	
	printf("Finiti q3 e bond\n");
	
	spLiq5=spLiq4=spSolid=0.0;
    int qqindex=0;
    int q3index=0;
    int qq6index=0;
    int q6index=0;
	for(i=0;i<s->Ato;i++)
	{				
		sp1=0;			
		for(j=0;j<s->Ato;j++)
		{
			if(i!=j)
			{
				for(k=0;k<3;k++)
				{
					cl->r[k]=d->r[j][k]-d->r[i][k];
					cl->r[k]-=d->lBox[k]*rint(cl->r[k]/d->lBox[k]);
				}
				cl->r[3]=sqrt(sqr(cl->r[0])+sqr(cl->r[1])+sqr(cl->r[2]));
				if(cl->r[3]<s->Cl_r)			// Atomo dentro la prima shell
				{			
					q3q3R=0.0;
					q3q3I=0.0;
					
					for(m=0;m<7;m++)
					{
						q3q3R+=(cl->q3mR[i][m]*cl->q3mR[j][m])+(cl->q3mI[i][m]*cl->q3mI[j][m]);
						q3q3I+=(cl->q3mI[i][m]*cl->q3mR[j][m])-(cl->q3mR[i][m]*cl->q3mI[j][m]);
					}
	
					if(q3q3R<-0.23)
					{
						sp1++;
					}
					qqindex= (int)((q3q3R+1.0)*50); // []-1.0,1.0] divided into 100 bins of 0.02
                    //printf("the qq index is %d for %f\n",qqindex,q3q3R);
                    if(qqindex<100 && qqindex>=0)
                    {
                        if(cl->neigh[i]==3)
                        cl->qqdistro3[qqindex]+=1.0;
                        
                        if(cl->neigh[i]==4)
                        cl->qqdistro4[qqindex]+=1.0;
                        
                        if(cl->neigh[i]==5)
                        cl->qqdistro5[qqindex]+=1.0;
                        
                        if(cl->neigh[i]==6)
                        cl->qqdistro6[qqindex]+=1.0;
                    }
                    
                    q6q6R=0.0;
					q6q6I=0.0;
					
					for(m=0;m<13;m++)
					{
						q6q6R+=(cl->q6mR[i][m]*cl->q6mR[j][m])+(cl->q6mI[i][m]*cl->q6mI[j][m]);
						q6q6I+=(cl->q6mI[i][m]*cl->q6mR[j][m])-(cl->q6mR[i][m]*cl->q6mI[j][m]);
					}
	                
	                //if(q6q6R>0.7)
					//{
					//	sp1++;
					//}
                    
                    qq6index= (int)((q6q6R+1.0)*50); // []-1.0,1.0] divided into 100 bins of 0.02
                    //printf("the qq index is %d for %f\n",qqindex,q3q3R);
                    if(qq6index<100 && qq6index>=0)
                    {
                        if(cl->neigh[i]==3)
                        cl->qq6distro3[qq6index]+=1.0;
                        
                        if(cl->neigh[i]==4)
                        cl->qq6distro4[qq6index]+=1.0;
                        
                        if(cl->neigh[i]==5)
                        cl->qq6distro5[qq6index]+=1.0;
                        
                        if(cl->neigh[i]==6)
                        cl->qq6distro6[qq6index]+=1.0;
                    }
                    
				}
			}
		}
		cl->qqNeig[sp1]+=1.0;
		if(sp1>=s->Cl_solid) // 3 or more bonds
		{	
			if(sp1<10)
			{
				cl->bondHisto[sp1][0]++;
			}
			spSolid++;
			cl->nucleo[spSolid][0]=i;
			// set current type of i here
		}
		else
		{
            // check the coordination number and q3 distro of those that are bonded to 2 or fewer
            q3index = (int) (cl->q3[i]*0.33*s->OrPar_bin);
            printf("q3 index is %d for %f\n",q3index,cl->q3[i]);
            if(cl->neigh[i]<=4)
            cl->q3NCDistro1[q3index]+=1.0;
            if(cl->neigh[i]>4)
            cl->q3NCDistro2[q3index]+=1.0;
            
            q6index = (int) (cl->q6[i]*0.33*s->OrPar_bin);
            printf("q6 index is %d for %f\n",q6index,cl->q6[i]);
            if(cl->neigh[i]<=4)
            cl->q6NCDistro1[q6index]+=1.0;
            if(cl->neigh[i]>4)
            cl->q6NCDistro2[q6index]+=1.0;

			if(cl->q3[i]>=0.6 && cl->neigh[i]==4)
			{
				if(sp1<10)
				{
					cl->bondHisto[sp1][2]++;
				}
				spLiq4++;
				cl->nucleo[spLiq4][2]=i; 
				// set current type of i here
			}
			else
			{
				if(sp1<10)
				{
					cl->bondHisto[sp1][1]++;
				}
				spLiq5++;
				cl->nucleo[spLiq5][1]=i; 
				// set current type of i here
			}
		}
		
		// insert code here for bond life time
		///////////////////////////////////////////////
		
		// first up, we have current type.
		// loop j goes from i+1 to s->Ato
		// next condition is if types are unchanged
			// if unchanged
			    // check type and if bonded ->t(i,j) increment or write+reset depending (3 wing "if" with one level nesting )
			// if changed
			   // check if bonded -> if yes, reset t(i,j) ; if no, do nothing
		// update type
				
	}
	cl->nucleo[0][0]=spSolid;
	cl->nucleo[0][1]=spLiq5;
	cl->nucleo[0][2]=spLiq4;
	printf("N solid, N HDl, N LDL - %d %d %d\n",spSolid,spLiq5,spLiq4);
	while(cl->numfield<20*sqrt(cl->nucleo[0][1]))
	{
		cl->numfield *=2;
	}
	cl->numfield/=2;
	printf("one step further\n");
	//cl->numfield = 200;
	//printf("%d numfield thing within\n",cl->numfield);

	
	if(cl->nucleo[0][0]+cl->nucleo[0][1]+cl->nucleo[0][2]!=s->Ato)
	{
		printf("Errore al tempo %d: solid:%d hdl:%d ldl:%d\n",cl->counter,cl->nucleo[0][0],cl->nucleo[0][1],cl->nucleo[0][2]);
		printf("La somma dovrebbe valere %d\n",s->Ato);
		exit(0);
	}
	//////////////// FILE WRITE PORTION
	// write all atoms here
	printf("maybe here?\n");
	// now print the first three lines of info
	printf("%d %d %d - bogey numbers\n",cl->nucleo[0][0]-1,cl->nucleo[0][1]-1,cl->nucleo[0][2]-1);
	
	// file write portion ended
	printf("finding clusters?\n");
	for(m=0;m<3;m++)
	{
		ipair=0; 
		for(i=1;i<=cl->nucleo[0][m];i++)
		{
			for(j=i+1;j<=cl->nucleo[0][m];j++)
			{
				for(k=0;k<3;k++)
				{
					cl->r[k]=d->r[cl->nucleo[j][m]][k]-d->r[cl->nucleo[i][m]][k]; 
					// need cl->nucleo[i][m]-1 because that number refers to index of ith atom
					cl->r[k]-=d->lBox[k]*rint(cl->r[k]/d->lBox[k]);
				}
				cl->r[3]=sqrt(sqr(cl->r[0])+sqr(cl->r[1])+sqr(cl->r[2]));
				
				if(cl->r[3]<s->Cl_cut)	// Le particelle solide sono define entro il SW cut-off
				{
					//printf("problems with this C ipair thingy %d %d %d %d\n?",ipair,i,j,m);
					ipair++;
					cl->C[ipair][1]=cl->nucleo[i][m];
					cl->C[ipair][2]=cl->nucleo[j][m];
				}
			}
		}
		
		for(i=1;i<=cl->nucleo[0][m];i++)
		{
			sp1=cl->nucleo[i][m];
			cl->clustSeb[sp1]=sp1;
		}
		
		change=1;
		while(change==1)
		{
			for(i=1;i<=cl->nucleo[0][m];i++)
			{
				sp1=cl->nucleo[i][m];
				cl->oldclustSeb[sp1]=cl->clustSeb[sp1];
			}
			
			for(i=1;i<=cl->nucleo[0][m];i++)
			{
				sp1=cl->nucleo[i][m];
				for(j=1;j<=ipair;j++)
				{
					if((cl->C[j][1]==sp1)||(cl->C[j][2]==sp1))
					{
						if(cl->C[j][1]==sp1)
						{
							spSolid=cl->C[j][2];
						}
						if(cl->C[j][2]==sp1)
						{
							spSolid=cl->C[j][1];
						}
						
						cl->clustSeb[sp1]=MIN(cl->clustSeb[sp1],cl->clustSeb[spSolid]);
						cl->clustSeb[sp1]=MIN(cl->clustSeb[sp1],cl->clustSeb[spSolid]);
					}
				}
			}
			
			change=0;
			for(i=1;i<=cl->nucleo[0][m];i++)
			{
				sp1=cl->nucleo[i][m];
				if(cl->oldclustSeb[sp1]!=cl->clustSeb[sp1])
				{
					change=1;
				}
			}
		}

		cnum=0;
		cref=-1;
		for(i=1;i<=s->Ato;i++)
		{
			cl->csizeSeb[i]=0;
		}
		
		for(i=1;i<=cl->nucleo[0][m];i++)
		{
			sp1=cl->nucleo[i][m];
			if(cl->clustSeb[sp1]>cref)
			{
                //printf("found a cluster\n");
				cnum++;
				cref=cl->clustSeb[sp1];
				cl->csizeSeb[cnum]=0;

				for(j=1;j<=cl->nucleo[0][m];++j)
				{
					spSolid=cl->nucleo[j][m];
                    
					if(cl->clustSeb[spSolid]==cref)
					{
						cl->clustSeb[spSolid]=cnum; 
						cl->csizeSeb[cnum]++; 
					}
				}
			}
		}
	
		maxcsizeSeb=0;
		maxclustSeb=0;
		for(i=1;i<=s->Ato;i++)
		{
            //printf("cluster index?? %d\n",cl->csizeSeb[i]);
			if(cl->csizeSeb[i]>maxcsizeSeb)
			{
				maxcsizeSeb=cl->csizeSeb[i]; 
				maxclustSeb=i;
			}
		}
		printf("%d %d index and max size\n",cl->counter,s->Cl_limit);
		printf("%d number of particles,\n",cl->nucleo[0][m]);
		cl->ClTime[cl->counter][m+3]=maxcsizeSeb;
		cl->ClTime[cl->counter][m]=cl->nucleo[0][m];
	}
	
	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<7;j++)
		{
			cl->q3mR[i][j]=0.0;
			cl->q3mI[i][j]=0.0;
		}
		for(j=0;j<13;j++)
		{
			cl->q6mR[i][j]=0.0;
			cl->q6mI[i][j]=0.0;
		}
		cl->q3[i]=0.0;
		cl->q6[i]=0.0;
	}
	// insert code to calc coordination number of all atoms here
	printf("old code cluster works??\n");
	// reset the coord/neighbour list
	for(i=0;i<s->Ato;i++)
	{	
		for(j=0;j<cl->maxCoord;j++)
		{
			cl->coordList[i][j]=0;
		}
	}
	// for each atom
	//printf("cluster coord entered\n");
	NeiCount =0;
	for(i=0;i<s->Ato;i++)
	{
		//printf("whussappening???? %d\n",j);
		// for each j other atoms
		NeiCount = 1;
		for(j=0;j<s->Ato;j++)
		{
			if(i!=j)
			{
				// if dist betw i and j < cutoff
				for(k=0;k<3;k++)
				{
					posR[k]=d->r[i][k]-d->r[j][k];
					posR[k]-=d->lBox[k]*rint(posR[k]/d->lBox[k]);
				}
				posR[3]=sqrt(sqr(posR[0])+sqr(posR[1])+sqr(posR[2]));
				if(posR[3]<s->cn_cutoff) // cn_cutoff is cut-off in anafile, read in structures.h
				{
					//printf("entered if for %d and %d\n",i,j);
					// store j as neigh of i, list[i][0] will be incremented for one more neighbour.
					cl->coordList[i][0]++; // one more neighbour for i
					cl->coordList[i][NeiCount]=j; // id of NeiCountth neighbour
					NeiCount++; // for the next neighbour
				
				}
			}
			
		}
		//printf("finished atom no %d %d\n",i,cl->coordList[i][0]);
		// now we know how many neighbours i has. We can update the histogram now
		// first find out which bin it will go in
		temp = cl->coordList[i][0]; // <-is an int
		// then increment the value at that bin index.

	}	
	printf("cluster coord no exited");
	// section to check coordination numbers of LDL, HDL and crystal
	for(k=0;k<3;k++) // for the three types
	{	
	//printf("next type happening %d\n",k);
	//printf("number of nuclei?? %d %d %d \n",cl->nucleo[0][0],cl->nucleo[0][1],cl->nucleo[0][2]);
		for(i=1;i<cl->nucleo[0][k];i++) // for each particle of kth type - REMEMBER that i should start from 1 bec 0 is total number
		{
			//printf("%d index problems for HDL??\n",i);
			spamindex = cl->nucleo[i][k]; //index of the atom
			temp = cl->coordList[spamindex][0]; // coordination number of that atom
			//printf("%d %d\n",spamindex,temp);
			if(temp<cl->maxCoord)
			{
				cl->coordHisto[k][temp]++; // one more atom with that coordination for this type
			}
			else
			printf("erroneous values %d has %f neighbours\n",spamindex,cl->coordHisto[k][temp]);

		}
	}
	
	
	cl->counter++;
	printf("%d\n",cl->counter);
	if(cl->counter==s->Cl_limit)
	{
		//Compute_Cl_correl(s,d,cl);
		printf("entering print\n");
		print_ClSize(s,d,cl);
	}
}

	
// insert code here for correlations
// we need a separate subroutine - make sure all args are there
// ALGO
// for each crystal particle
// loop over all HDL particles
// calc distance - bin it. (bin size - box length/100)
// make sure that the contribution from each HDL is divided by total number of crystal particles
// then loop over all LDL particles
// calc distance  - bin it.
// make sure that the contribution from each HDL is divided by total number of crystal particles 

void Compute_Cl_correl(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	int i,j,k,m,spi,spj,temp;
	double posR[4];
	cl->numfield=1;
 	while(cl->numfield<20*sqrt(cl->nucleo[0][1]))
	{
		cl->numfield *=2;
	}
	cl->numfield/=2;
	//cl->numfield = 200;
	printf("%d numfield thing\n",cl->numfield);
	for(i=1;i<cl->nucleo[0][1];i++)
	{	
		spi = cl->nucleo[i][1];
		// for each crystal particle
		for(m=0;m<3;m++) // for each particle of other types 
		{
			for(j=1;j<cl->nucleo[0][m];j++)// for every atom of that other type
			{
				spj = cl->nucleo[j][m];
				//find distance between the given crystal particle and that other particle
				if(spi!=spj)
				{
					for(k=0;k<3;k++)
					{
						posR[k]=d->r[spi][k]-d->r[spj][k];
						posR[k]-=d->lBox[k]*rint(posR[k]/d->lBox[k]);
					}
					posR[3]=sqrt(sqr(posR[0])+sqr(posR[1])+sqr(posR[2]))*2.0*cl->numfield/d->lBox[0];
					
					// bin it - add to histo with appropriate weight
					//temp = (int) 100.0*posR[3]/(1.414*d->lBox[0]); // find dist as frac of box length times nbins. Reverse transform while printing
					//cl->CryNeigh[temp][m-1] = cl->CryNeigh[temp][m-1] + 1.0/cl->nucleo[0][0]; // one more atom at that distance - 
					// contribution normalised against number of crystal atoms
					temp = rint(posR[3]);
					if(temp<=cl->numfield)
					{
						//printf("updating?? %d\n",temp);
						cl->CryNeigh[temp][m]++;
						// normalise in print	
					}
			
				}
				
			
			}
		}
	}	
}

