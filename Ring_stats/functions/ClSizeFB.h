//#include"YgCoordNo.h"
struct CLSIZE
{
	int counter;          
	int fileCounter;
	int *oldclustSeb,**nucleo,*clustSeb,*csizeSeb;
	int **C,*neigh;
	double **q3mR,**q3mI,*r,*q3, *q6mR, *q6mI;
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
	char buffer[150];
	cl->numfield=1;
	cl->maxCoord = 70;
	cl->q3mR=AllocMatR(s->Ato,7);
	cl->q3mI=AllocMatR(s->Ato,7);
	cl->q6mR=AllocVecR(13);
	cl->q6mI=AllocVecR(13);
	cl->r=AllocVecR(4);
	cl->C=AllocMatI(200*s->Ato,3);
	cl->oldclustSeb=AllocVecI(s->Ato+1);
	cl->nucleo=AllocMatI(s->Ato+1,3);		// Solid + liquid5 + liquid4
	cl->clustSeb=AllocVecI(s->Ato+1);
	cl->csizeSeb=AllocVecI(s->Ato+1);
	cl->q3=AllocVecR(s->Ato);
	cl->neigh=AllocVecI(s->Ato);
	printf("cluster size, all allocated??\n");
	cl->ClTime=AllocMatI(s->Cl_limit+1,6);	// S + l5 + l4 + s_Max + l5_Max + l4_Max + nBond
	cl->coordHisto=AllocMatR(3,cl->maxCoord);
	cl->coordList=AllocMatI(s->Ato,cl->maxCoord); // list of particles, no. neigh and nerigh ids.
	cl->CryNeigh=AllocMatR(4000,3);
	cl->bondHisto=AllocMatR(10,3);
	for(i=0;i<s->Ato;i++)
	{
		for(j=0;j<7;j++)
		{
			cl->q3mR[i][j]=0.0;
			cl->q3mI[i][j]=0.0;
		}
		cl->q3[i]=0.0;
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
	for(i=0;i<13;i++)
	{
		cl->q6mR[i]=0.0;
		cl->q6mI[i]=0.0;
	}
	cl->counter=0;
	cl->fileCounter=0;
}

void print_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	FILE *fp;	
	int t,spam;
	char buffer[150];
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
	sprintf(buffer,"%sBondHisto%g-%d.dat",s->output,s->Cl_r,cl->fileCounter);
	fp=fopen(buffer,"w");
	for(t=0;t<10;t++)
	{
		cl->bondHisto[t][0] /= s->Cl_limit;
		cl->bondHisto[t][1] /= s->Cl_limit;
		cl->bondHisto[t][2] /= s->Cl_limit;
		fprintf(fp,"%d %f %f %f\n",t,cl->bondHisto[t][0],cl->bondHisto[t][1],cl->bondHisto[t][2]);
	}

	
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

	sprintf(buffer,"%sAllAt_ovito_%d.dat",s->output,cl->fileCounter);
	fp = fopen(buffer,"w");
	fprintf(fp,"%d\n",s->Ato);
	fprintf(fp,"random text\n");
	
	for(t=1;t<=cl->nucleo[0][0];t++)
	{
		spam = cl->nucleo[t][0];
		fprintf(fp,"Si0 %f %f %f\n",d->r[spam][0],d->r[spam][1],d->r[spam][2]);
	}
	for(t=1;t<=cl->nucleo[0][1];t++)
	{
		spam = cl->nucleo[t][1];
		fprintf(fp,"Si1 %f %f %f\n",d->r[spam][0],d->r[spam][1],d->r[spam][2]);
	}
	for(t=1;t<=cl->nucleo[0][2];t++)
	{
		spam = cl->nucleo[t][2];
		fprintf(fp,"Si2 %f %f %f\n",d->r[spam][0],d->r[spam][1],d->r[spam][2]);
	}
	fclose(fp);
	
	sprintf(buffer,"%stCry_ovito_%d.dat",s->output,cl->fileCounter);
	fp = fopen(buffer,"w");
	fprintf(fp,"%d\n",cl->nucleo[0][0]);
	fprintf(fp,"random text\n");
	
	for(t=1;t<=cl->nucleo[0][0];t++)
	{
		spam = cl->nucleo[t][0];
		fprintf(fp,"Si0 %f %f %f\n",d->r[spam][0],d->r[spam][1],d->r[spam][2]);
	}
	fclose(fp);
	
	
	sprintf(buffer,"%st_ClCry_ovito_%d.dat",s->output,cl->fileCounter);
	fp = fopen(buffer,"w");
	fprintf(fp,"%d\n",cl->nucleo[0][0]);
	fprintf(fp,"random text\n");
	
	for(t=1;t<=cl->nucleo[0][0];t++)
	{
		spam = cl->nucleo[t][0];
		fprintf(fp,"Cl%d %f %f %f\n",cl->clustSeb[spam],d->r[spam][0],d->r[spam][1],d->r[spam][2]);
	}
	fclose(fp);
	
	cl->fileCounter++;
	printf("done coord stats, exiting\n");
}

void Compute_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	int i,k,nBond,j,m,spamindex,temp,mspam;
	int sp1,spSolid,spLiq5,spLiq4,ipair;
	double q3q3R,q3q3I;
	int change, cnum, cref, maxcsizeSeb, maxclustSeb;
	double cosTheta,phi;

	int NeiCount;
	double posR[4];
	
	char buffer[150],buffer2[150],buffer3[150],buffer4[150],buffer5[150];	
	FILE *fp, *fp2, *fp3, *fp4, *fp5;
	printf("entering cl size\n");
	sprintf(buffer,"%scry%g-%g.xyz",s->output,s->Cl_r,s->Cl_cut);
	sprintf(buffer2,"%sHDL%g-%g.xyz",s->output,s->Cl_r,s->Cl_cut);
	sprintf(buffer3,"%sLDL%g-%g.xyz",s->output,s->Cl_r,s->Cl_cut);
	printf("%scry%g-%g.xyz\n",s->output,s->Cl_r,s->Cl_cut);
	printf("%sHDL%g-%g.xyz\n",s->output,s->Cl_r,s->Cl_cut);
	printf("%sLDL%g-%g.xyz\n",s->output,s->Cl_r,s->Cl_cut);
	fp=fopen(buffer,"w");
	fp2=fopen(buffer2,"w");
	fp3=fopen(buffer3,"w");
	
	sprintf(buffer4,"%sq6mRhisto%g-%g.dat",s->output,s->Cl_r,s->Cl_cut);
	sprintf(buffer5,"%sq6mIhisto%g-%g.dat",s->output,s->Cl_r,s->Cl_cut);
	fp4=fopen(buffer4,"w");
	fp5=fopen(buffer5,"w");
	
	mspam = 0;
	
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
		for(j=0;j<13;j++)
		{
		cl->q6mR[j]=0.0;
		cl->q6mI[j]=0.0;
		}
		nBond=0;
		fprintf(fp4,"%d ",i);
		fprintf(fp5,"%d ",i);
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
					// make changes here to output the spherical harmonic thing for all particles.
					for(m=0;m<=6;m++)
					{
						cl->q6mR[m]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*cos(m*phi);
						cl->q6mI[m]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*sin(m*phi);
					}
					for(m=1;m<=6;m++)
					{
						cl->q6mR[m+6]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*cos(-m*phi)*pow(-1,m);
						cl->q6mI[m+6]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*sin(-m*phi)*pow(-1,m);
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
		cl->neigh[i]=nBond;								// Adesso ho il q3(i) ed il numero di bond
		
		
		/// for the q6 contribution thing
		for(m=0;m<13;m++)
		{	
			//cl->q6mR[m]/=nBond;
			//cl->q6mI[m]/=nBond;
			// now print to file
			fprintf(fp4," %f ",cl->q6mR[m]);
			fprintf(fp5," %f ",cl->q6mI[m]);
		}
		fprintf(fp4,"\n");
		fprintf(fp5,"\n");
		
	}
	fclose(fp4);
	fclose(fp5);
	printf("Finiti q3 e bond\n");
	
	spLiq5=spLiq4=spSolid=0.0;
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
				}
			}
		}
		
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
	fprintf(fp,"%d\n",cl->nucleo[0][0]-1); // N atoms cry	
	printf("really???\n");
	fprintf(fp,"random comment\n");
	fprintf(fp2,"%d\n",cl->nucleo[0][1]-1); // N atoms HDL
	fprintf(fp2,"random comment\n");
	fprintf(fp3,"%d\n",cl->nucleo[0][2]-1); // N atoms LDL
	fprintf(fp3,"random comment\n");

	printf("is it here?\n");
	// for each type separate loop because of file name shit
	for(i=1;i<cl->nucleo[0][0];i++) // for each atom of that type
	{	
		spamindex = cl->nucleo[i][0];
		for(k=0;k<3;k++)
		{
			cl->r[k]=d->r[spamindex][k];
			// need cl->nucleo[i][m]-1 because that number refers to index of ith atom
			cl->r[k]-=d->lBox[k]*rint(cl->r[k]/d->lBox[k]);
			cl->r[k]+=d->lBox[k]/2.0; // for ovito - make all positive
		}
		if(cl->r[0]>25.0 && cl->r[0]<43.0 && cl->r[1]>25.0 && cl->r[1]<43.0 && cl->r[2]>25.0 && cl->r[2]<43.0)
		{ // bounds on the if condition specific to 8000 particles
		fprintf(fp,"Si   %f   %f   %f\n",cl->r[0],cl->r[1],cl->r[2]);
		}
		
	}
	fclose(fp);
	for(i=1;i<cl->nucleo[0][1];i++) // for each atom of that type
	{	
		spamindex = cl->nucleo[i][1];
		for(k=0;k<3;k++)
		{
			cl->r[k]=d->r[spamindex][k]; 
			// need cl->nucleo[i][m]-1 because that number refers to index of ith atom
			cl->r[k]-=d->lBox[k]*rint(cl->r[k]/d->lBox[k]);
			cl->r[k]+=d->lBox[k]/2.0; // for ovito - make all positive
		}
		if(cl->r[0]>25.0 && cl->r[0]<43.0 && cl->r[1]>25.0 && cl->r[1]<43.0 && cl->r[2]>25.0 && cl->r[2]<43.0)
		{// bounds on the if condition specific to 8000 particles
		fprintf(fp2,"Si   %f   %f   %f\n",cl->r[0],cl->r[1],cl->r[2]);
		}
	}
	fclose(fp2);
	for(i=1;i<cl->nucleo[0][2];i++) // for each atom of that type
	{	
		spamindex = cl->nucleo[i][2];
		for(k=0;k<3;k++)
		{
			cl->r[k]=d->r[spamindex][k];
			// need cl->nucleo[i][m]-1 because that number refers to index of ith atom
			cl->r[k]-=d->lBox[k]*rint(cl->r[k]/d->lBox[k]);
			cl->r[k]+=d->lBox[k]/2.0; // for ovito - make all positive
		}
		if(cl->r[0]>25.0 && cl->r[0]<43.0 && cl->r[1]>25.0 && cl->r[1]<43.0 && cl->r[2]>25.0 && cl->r[2]<43.0)
		{// bounds on the if condition specific to 8000 particles
		fprintf(fp3,"Si   %f   %f   %f\n",cl->r[0],cl->r[1],cl->r[2]);
		}
		
	}
	fclose(fp3);
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
			//printf("%d how many of type %d\n",cl->nucleo[0][m],m);
			for(i=1;i<=cl->nucleo[0][m];i++)
			{
				sp1=cl->nucleo[i][m];
				//printf("%d %d\n",sp1,m);
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
						//printf("%d %d wtf man!\n",cl->clustSeb[cl->C[1][1]],cl->clustSeb[cl->C[1][2]]);
						if(cl->clustSeb[sp1]!=cl->clustSeb[spSolid])
						{
							mspam = 1;
							//printf("%d %d %d %d no really\n",sp1, spSolid,cl->clustSeb[sp1],cl->clustSeb[spSolid]);							
						}
						cl->clustSeb[sp1]=MIN(cl->clustSeb[sp1],cl->clustSeb[spSolid]);
						cl->clustSeb[spSolid]=MIN(cl->clustSeb[sp1],cl->clustSeb[spSolid]); // remember the change here
						if(mspam==1)
						//printf("%d %d %d %d no no really\n",sp1, spSolid,cl->clustSeb[sp1],cl->clustSeb[spSolid]);						
						
						mspam = 0;
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
			printf("%d %d\n",cl->clustSeb[sp1],m);
			if(cl->clustSeb[sp1]>cref)
			{
				cnum++;
				cref=cl->clustSeb[sp1];
				cl->csizeSeb[cnum]=0;

				for(j=1;j<=cl->nucleo[0][m];j++)
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
		printf("%d %d %d %d\n",cl->csizeSeb[1],cl->csizeSeb[2],cl->csizeSeb[3],m);
	
		maxcsizeSeb=0;
		maxclustSeb=0;
		for(i=1;i<=s->Ato;i++)
		{
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

