//#include"YgCoordNo.h"
struct CLSIZE
{
	int counter;          
	int fileCounter;
	int *oldclustSeb,**nucleo,*clustSeb,*csizeSeb;
	int **C,*neigh;
	double **q3mR,**q3mI,*r,*q3;
	int **ClTime;	// Max cluster size e numero di solid particles,
	double **coordHisto, **CryNeigh, **bondHisto;	
	//	and histogram of coord nos. double because we are normalising histograms, also crystal specific g(r), histogram of bonds
	int **coordList;
	int maxCoord;
	int numfield; // complicated shit for normalising a la g(r)
	double *clusterHisto, *Pcluster;
} clsize;

void Init_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	int i,j,spamindex,temp;
	char buffer[500];
	cl->numfield=1;
	cl->maxCoord = 70;
	cl->q3mR=AllocMatR(s->Ato,7);
	cl->q3mI=AllocMatR(s->Ato,7);
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
    cl->clusterHisto=AllocVecR(s->Ato+1);
    cl->Pcluster=AllocVecR(s->Ato+1);
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
		for(j=0;j<2;j++)
		{
		cl->CryNeigh[i][j]=0.0;
		}
	}
	for(i=0;i<10;i++)
	{
		for(j=0;j<2;j++)
		{
		cl->bondHisto[i][j]=0.0;
		}
	}
	for(i=0;i<s->Ato+1;i++)
    {
        cl->clusterHisto[i]=0.0;
        cl->Pcluster[i]=0.0;
    }
    for(i=0;i<s->Cl_limit+1;i++)
    {
        for(j=0;j<6;j++)
        {
            cl->ClTime[i][j]=0.0;
        }    
    }
	cl->counter=0;
	cl->fileCounter=0;
}

void print_ClSize(SIMDAT *s, DATASIM *d, CLSIZE *cl)
{
	FILE *fp;	
	int t;
	char buffer[500];
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
    fclose(fp);
    
    sprintf(buffer,"%sClusterHisto%g-%d.dat",s->output,s->Cl_r,cl->fileCounter);
    fp=fopen(buffer,"w");
    for(t=0;t<=s->Ato;t++)
    {
        if(cl->clusterHisto[t]>0.0)
        fprintf(fp,"%f %f\n",t+0.001,cl->clusterHisto[t]/cl->counter);
    }
    fclose(fp);
    
    sprintf(buffer,"%sP_Cluster%g-%d.dat",s->output,s->Cl_r,cl->fileCounter);
    fp=fopen(buffer,"w");
    for(t=0;t<=s->Ato;t++)
    {
        fprintf(fp,"%f %f\n",t+0.001,cl->Pcluster[t]/cl->counter);
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
	int change, cnum, cref, maxcsizeSeb, maxclustSeb;
	double cosTheta,phi;

	int NeiCount;
	double posR[4];
	
	

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
	}
	
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
	printf("N solid, N HDL, N LDL - %d %d %d\n",spSolid,spLiq5,spLiq4);
	


	
	if(cl->nucleo[0][0]+cl->nucleo[0][1]+cl->nucleo[0][2]!=s->Ato)
	{
		printf("Errore al tempo %d: solid:%d hdl:%d ldl:%d\n",cl->counter,cl->nucleo[0][0],cl->nucleo[0][1],cl->nucleo[0][2]);
		printf("La somma dovrebbe valere %d\n",s->Ato);
		exit(0);
	}
	
	printf("finding clusters?\n");
	for(m=0;m<1;m++)
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
        int spammer=0;
		while(change==1)
		{
            spammer+=1;
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
						cl->clustSeb[spSolid]=MIN(cl->clustSeb[sp1],cl->clustSeb[spSolid]);
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

		cnum=-1;
		cref=-1;
		for(i=0;i<=s->Ato;i++)
		{
			cl->csizeSeb[i]=0;
		}
//         for(i=0;i<s->Ato;i++)
//             printf("cluster index of %d is %d\n",i,cl->clustSeb[i]);
		for(i=1;i<=cl->nucleo[0][m];i++)
		{
			sp1=cl->nucleo[i][m];
//             if(m==0)
//             printf("cluster index re-check %d %d\n",cl->clustSeb[sp1],cref);
			if(cl->clustSeb[sp1]>cref)
			{
				cnum++;
//              printf("are we updating cref %d %d\n",cref,cl->clustSeb[sp1]);
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
        // update the full cluster distribution
        if(m==0)
		{	
           // printf("%d is the number of clusters \n",cnum+1);
		for(i=0;i<=cnum;i++) // take stats of all cluster sizes in the system - cnum clusters
		{
		cl->clusterHisto[cl->csizeSeb[i]]+=1.0; // how many clusters of size csizeSeb[index]
		cl->Pcluster[cl->csizeSeb[i]]+=1.0/s->Ato;
		//printf("clus size current %d %d %f\n",i,csizeSeb[i],box->currentHisto[csizeSeb[i]]);
		}
		cl->Pcluster[0]+=1.0*(s->Ato-cl->nucleo[0][0])/(1.0*s->Ato);
        }
		maxcsizeSeb=0;
		maxclustSeb=0;
		for(i=0;i<=s->Ato;i++)
		{
			if(cl->csizeSeb[i]>maxcsizeSeb)
			{
				maxcsizeSeb=cl->csizeSeb[i]; 
				maxclustSeb=i;
			}
		}
		//printf("%d %d index and max size\n",cl->counter,s->Cl_limit);
		//printf("%d number of particles,\n",cl->nucleo[0][m]);
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

	

