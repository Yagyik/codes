long int lround(double x);

void largestClus(struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,k,nbi,nBond,j,m;
	int sp1,sp2,ipair;
	int *oldclustSeb,**nucleo,*clustSeb,*csizeSeb;
	int **C,*neigh;
	double cosTheta, phi, *r;
	double **q3mR,**q3mI,*q3;
	double q3q3R,q3q3I;
	int change, cnum, cref, maxcsizeSeb, maxclustSeb;
	int spSolid,spLiq5,spLiq4;

	q3mR=calloc(p->nAto+1,sizeof(double *));
	q3mI=calloc(p->nAto+1,sizeof(double *));
	C=calloc(20*p->nAto,sizeof(int *));
	q3=calloc(p->nAto+1,sizeof(double));
	oldclustSeb=calloc(p->nAto+1,sizeof(int));
	
	nucleo=calloc(p->nAto+1,sizeof(int *));
	for(i=0;i<p->nAto+1;i++)
	{
		nucleo[i]=calloc(3,sizeof(int));		// solid-liq5-liq4
	}
	
	clustSeb=calloc(p->nAto+1,sizeof(int));
	csizeSeb=calloc(p->nAto+1,sizeof(int));
	neigh=calloc(p->nAto+1,sizeof(int));
	r=calloc(4,sizeof(double));
	for(i=0;i<=p->nAto;i++)
	{
		q3mR[i]=calloc(7,sizeof(double));
		q3mI[i]=calloc(7,sizeof(double));
	}
	for(i=0;i<20*p->nAto;i++)
	{
		C[i]=calloc(3,sizeof(int));
	}
	
	for(i=0;i<=p->nAto;i++)
	{
		q3[i]=0.0;
		for(k=0;k<7;k++)
		{
			q3mR[i][k]=0.0;
			q3mI[i][k]=0.0;
		}
	}

	for(i=1;i<=p->nAto;i++)
	{
		nBond=0;
		for(j=1;j<=box->NEIGLIST[i][0];j++)
		{
			r[0]=box->X[box->NEIGLIST[i][j]]-box->X[i];
			r[0]-=box->lx*lround(r[0]/box->lx);

			r[1]=box->Y[box->NEIGLIST[i][j]]-box->Y[i];
			r[1]-=box->ly*lround(r[1]/box->ly);

			r[2]=box->Z[box->NEIGLIST[i][j]]-box->Z[i];
			r[2]-=box->lz*lround(r[2]/box->lz);		

			r[3]=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	
			if(r[3]<p->q3Cut)
			{
				nBond++;

				cosTheta=r[2]/r[3];
				phi=atan2(r[1],r[0]);

				for(m=0;m<=3;m++)
				{
					q3mR[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(m*phi);
					q3mI[i][m]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(m*phi);
				}
				for(m=1;m<=3;m++)
				{
					q3mR[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*cos(-m*phi)*pow(-1,m);
					q3mI[i][m+3]+=gsl_sf_legendre_sphPlm(3,m,cosTheta)*sin(-m*phi)*pow(-1,m);
				}
			}
		}
		for(m=0;m<7;m++)
		{
			q3mR[i][m]/=nBond;
			q3mI[i][m]/=nBond;
			q3[i]+=q3mR[i][m]*q3mR[i][m]+q3mI[i][m]*q3mI[i][m];
		}
		q3[i]=sqrt(q3[i]*4.0*M_PI/7.0);		// Now I have q3 for each particle
		neigh[i]=nBond;							// Now I have number of bond for each particle
	}
	
	spLiq5=spLiq4=spSolid=0.0;
	for(i=1;i<=p->nAto;i++)
	{
		sp1=0;
		for(j=1;j<=box->NEIGLIST[i][0];j++)
		{
			r[0]=box->X[box->NEIGLIST[i][j]]-box->X[i];
			r[0]-=box->lx*lround(r[0]/box->lx);

			r[1]=box->Y[box->NEIGLIST[i][j]]-box->Y[i];
			r[1]-=box->ly*lround(r[1]/box->ly);

			r[2]=box->Z[box->NEIGLIST[i][j]]-box->Z[i];
			r[2]-=box->lz*lround(r[2]/box->lz);		

			r[3]=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	
			if(r[3]<p->q3Cut)
			{
				q3q3R=0.0;
				q3q3I=0.0;
				for(m=0;m<7;m++)
				{
					q3q3R+=(q3mR[i][m]*q3mR[box->NEIGLIST[i][j]][m])+(q3mI[i][m]*q3mI[box->NEIGLIST[i][j]][m]);
					q3q3I+=(q3mI[i][m]*q3mR[box->NEIGLIST[i][j]][m])-(q3mR[i][m]*q3mI[box->NEIGLIST[i][j]][m]);
				}
				if(q3q3R<-0.23)
				{
					sp1++;
				}
			}
		}		
	
		if(sp1>=p->clSize)
		{
			spSolid++;
			nucleo[spSolid][0]=i;			// Crystal-like
		}
		else
		{
			if(q3[i]>0.6 && neigh[i]==4)
			{
				spLiq4++;
				nucleo[spLiq4][2]=i;			// LDL-like
			}
			else
			{
				spLiq5++;
				nucleo[spLiq5][1]=i;			// HDL-like
			}
		}
	}
	nucleo[0][0]=spSolid;
	nucleo[0][1]=spLiq5;
	nucleo[0][2]=spLiq4;
	box->nSp=spSolid;
	box->nl5=spLiq5;
	box->nl4=spLiq4;
	
	if(spSolid+spLiq5+spLiq4!=p->nAto)
	{
		printf("Cluster Error! Solid:%d hdl:%d ldl:%d\n",spSolid,spLiq5,spLiq4);
		exit(-1);
	}
	
	for(m=0;m<3;m++)
	{
		ipair=0; 
		for(i=1;i<=nucleo[0][m];i++)
		{
			for(j=i+1;j<=nucleo[0][m];j++)
			{
				r[0]=box->X[nucleo[j][m]]-box->X[nucleo[i][m]];
				r[0]-=box->lx*lround(r[0]/box->lx);

				r[1]=box->Y[nucleo[j][m]]-box->Y[nucleo[i][m]];
				r[1]-=box->ly*lround(r[1]/box->ly);
	
				r[2]=box->Z[nucleo[j][m]]-box->Z[nucleo[i][m]];
				r[2]-=box->lz*lround(r[2]/box->lz);		

				r[3]=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
	
				if(r[3]<=p->clCut)
				{
					ipair++;
					C[ipair][1]=nucleo[i][m];
					C[ipair][2]=nucleo[j][m];
				}
			}
		}

		for(i=1;i<=nucleo[0][m];i++)
		{
			sp1=nucleo[i][m];
			clustSeb[sp1]=sp1;
		}
   
		change=1;
		while(change==1)
		{
			for(i=1;i<=nucleo[0][m];i++)
			{
				sp1=nucleo[i][m];
				oldclustSeb[sp1]=clustSeb[sp1];
			}

			for(i=1;i<=nucleo[0][m];i++)
			{
				sp1=nucleo[i][m];
				for(j=1;j<=ipair;j++)
				{
					if((C[j][1]==sp1)||(C[j][2]==sp1))
					{	
						if(C[j][1]==sp1)
						{
							sp2=C[j][2];
						}
						if(C[j][2]==sp1)
						{
							sp2=C[j][1];
						}
						clustSeb[sp1]=MIN(clustSeb[sp1],clustSeb[sp2]);
						clustSeb[sp2]=MIN(clustSeb[sp1],clustSeb[sp2]);
					}
				}				
			}

			change=0;
			for(i=1;i<=nucleo[0][m];i++)
			{
				sp1=nucleo[i][m];
				if(oldclustSeb[sp1]!=clustSeb[sp1])
				{
					change=1;
				}
			}
		}

		cnum=-1;
		cref=-1;
		for(i=0;i<=p->nAto;i++)
		{
			csizeSeb[i]=0;
		}

		for(i=1;i<=nucleo[0][m];i++)
		{
			sp1=nucleo[i][m];
			if(clustSeb[sp1]>cref)
			{
				cnum++;
				cref=clustSeb[sp1];
				csizeSeb[cnum]=0;				
	
				for(j=1;j<=nucleo[0][m];++j)
				{
					sp2=nucleo[j][m];
					if(clustSeb[sp2]==cref)
					{
						clustSeb[sp2]=cnum; 
						csizeSeb[cnum]++; 
					}
				}
			}
		}
		if(m==0)
		{
            //printf("wtf %d\n",csizeSeb[1]);
            for(i=0;i<NspUppLim;i++) // first reset for this time step/function call
            box->currentHisto[i]=0.0;
            
            for(i=0;i<=cnum;i++) // take stats of all cluster sizes in the system - cnum clusters
            {
            box->currentHisto[csizeSeb[i]]+=1.0/(1.0*(p->nAto-spSolid)); // how many clusters of size csizeSeb[index]
            //printf("clus size current %d %d %f\n",i,csizeSeb[i],box->currentHisto[csizeSeb[i]]);
            }
            box->currentHisto[0] = 1.0*(p->nAto-spSolid)/(1.0*(p->nAto-spSolid)); // the weight of zero clusters is the remaining non-crystalline atoms
		}
		
		maxcsizeSeb=0;
		maxclustSeb=0;
        //printf("%d %d m and maxcsize reset\n",m,maxcsizeSeb);
		for(i=0;i<=p->nAto;i++)
		{
            //printf("%d %d\n",i,csizeSeb[i]);
			if(csizeSeb[i]>maxcsizeSeb)
			{
				maxcsizeSeb=csizeSeb[i]; 
				maxclustSeb=i;
			}
		}
		//printf("%d %d m and maxcsize found\n",m,maxcsizeSeb);
		if(m==0)
		{
			box->nSp_max=1.0*maxcsizeSeb;
		}
		else
		{
			if(m==1)
			{
                //printf("%f %d\n",box->nl5_max,maxcsizeSeb);
				box->nl5_max=1.0*maxcsizeSeb;
			}
			else
			{
				box->nl4_max=1.0*maxcsizeSeb;
			}
		}
	}
	
	for(i=0;i<=p->nAto;i++)
	{	
		free(q3mR[i]);
		free(q3mI[i]);
	}
	for(i=0;i<20*p->nAto;i++)
	{
		free(C[i]);
	}
	for(i=0;i<=p->nAto;i++)
	{
		free(nucleo[i]);
	}
	
	free(nucleo);
	free(neigh);
	free(clustSeb);
	free(csizeSeb);
	free(C);
	free(q3mR);
	free(q3mI);
	free(r);
	free(oldclustSeb);
	free(q3);
}
