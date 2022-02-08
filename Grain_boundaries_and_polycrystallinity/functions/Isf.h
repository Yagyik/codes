struct ISF
{
	double *isf;		// Final ISF
	double **pos0;		// Initial position
	int nWaveVec;		// Number of Wavevectors
	int fileCounter;
	int counter; 
} IsfK;

void Init_Isf(SIMDAT *s, DATASIM *d, ISF *isf)
{
	int t;

	isf->pos0=AllocMatR(s->Ato,3); 
	isf->isf=AllocVecR(s->Isf_limit);
	for(t=0;t<s->Isf_limit;t++)
	{
		isf->isf[t]=0.0;
	}

	isf->counter=0;
	isf->fileCounter=0;
}

void print_Isf(SIMDAT *s, DATASIM *d, ISF *isf)
{	
	FILE *fp;
	int t;
	char buffer[500];

	sprintf(buffer,"%sIsfK%g_%d.dat",s->output,s->IsfK_value,isf->fileCounter);
	fp=fopen(buffer,"w");
	
	for(t=1;t<s->Isf_limit;t++)
	{
		fprintf(fp,"%le	%le\n",s->DeltaT*s->Step*t,isf->isf[t]);
	}
	fclose(fp);

	isf->fileCounter++;
	isf->counter=0;
	
}

void Compute_Isf(SIMDAT *s, DATASIM *d, ISF *isf)
{
	int a,x;
	int qx,qy,qz,qthr;
	double q;

	if(isf->counter==0)
	{
		for(a=0;a<s->Ato;a++)
		{
		for(x=0;x<3;x++)
		{
			isf->pos0[a][x]=d->r[a][x];
		}
		}
	}

	isf->nWaveVec=0;												// Reset number of K
	qthr=rint(d->lBox[0]*(s->IsfK_value+0.05)/(2.0*M_PI));			// Tolerance of 0.05

	for(qx=-qthr;qx<=qthr;qx++)
	{
	for(qy=-qthr;qy<=qthr;qy++)
	{
	for(qz=-qthr;qz<=qthr;qz++)
	{
		q=2.0*M_PI/d->lBox[0]*sqrt(qx*qx+qy*qy+qz*qz); 				// Modulus of K
		if(q<(s->IsfK_value+0.05) && q>(s->IsfK_value-0.05))		
		{
			isf->nWaveVec++;										// Counter number of K
			
			for(a=0;a<s->Ato;a++)
			{
				isf->isf[isf->counter]+=cos(qx*2.0*M_PI/d->lBox[0]*(d->r[a][0]-isf->pos0[a][0])+qy*2.0*M_PI/d->lBox[1]*(d->r[a][1]-isf->pos0[a][1])+qz*2.0*M_PI/d->lBox[2]*(d->r[a][2]-isf->pos0[a][2]));		// q dot DeltaR
			}
		}
	}
	}
	}

	isf->isf[isf->counter]/=s->Ato*isf->nWaveVec;
	
	isf->counter++; 
  
	if(isf->counter==s->Isf_limit) 
	{
		print_Isf(s,d,isf);
	}
}

