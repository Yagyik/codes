long int lround(double x);

double q6eval(struct SYSTEM *box,struct PARAINPUT *p)
{
	int i,k,nBond,j,m;
	double cosTheta,phi,*r;
	double *q6mR,*q6mI,q6;

	q6mR=calloc(13,sizeof(double));
	q6mI=calloc(13,sizeof(double));
	r=calloc(4,sizeof(double));

	nBond=0;
	for(i=0;i<13;i++)
	{
		q6mR[i]=0.0;
		q6mI[i]=0.0;
	}

	for(i=1;i<=p->nAto;i++)
	{
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

				for(m=0;m<=6;m++)
				{
					q6mR[m]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*cos(m*phi);
					q6mI[m]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*sin(m*phi);
				}
				for(m=1;m<=6;m++)
				{
					q6mR[m+6]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*cos(-m*phi)*pow(-1,m);
					q6mI[m+6]+=gsl_sf_legendre_sphPlm(6,m,cosTheta)*sin(-m*phi)*pow(-1,m);
				}
			}
		}
	}

	q6=0.0;
	for(m=0;m<13;m++)
	{
		q6mR[m]/=nBond;
		q6mI[m]/=nBond;
		q6+=q6mR[m]*q6mR[m]+q6mI[m]*q6mI[m];
	}

	q6*=4.0*M_PI/13.0;

	free(q6mR);
	free(q6mI);
	free(r);
	
	return(sqrt(q6));
}
